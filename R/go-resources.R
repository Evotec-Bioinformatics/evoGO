## This file includes functions for retrieving the annotations from
## external resources and to transform them into the input structures,
## which are required by the evoGO constructor.



##' Retrieve the ontologyIndex object from the obolibrary.
##'
##' @return ontologyIndex object
##' @export
getOBO <- function() {
  url <- "http://purl.obolibrary.org/obo/go/go-basic.obo"
  filename <- tempfile()

  curl::curl_download(url, filename)

  go <- ontologyIndex::get_ontology(filename,
    extract_tags = "everything",
    propagate_relationships = c(
      "is_a", "part_of", "regulates",
      "positively_regulates", "negatively_regulates"
    )
  )
  go
}


##' Create an input data.frame for evoGO constructor function.
##'
##' The data.frame includes "child" - "parent" relationships between the GO ids.
##' @param go an ontologyIndex object returned by getOBO function (also see ontologyIndex package)
##' @param rootName the name of the root elements from which to build
##' the hierarchy.
##'
##' @return a data.frame to be used as graphTable argument in evoGO function.
##' @export
ontologyIndex2graph <- function(go, rootName = "cellular_component") {
  parents <- go$parents
  ids <- go$id
  hasParent <- vapply(go$parents, length, numeric(1)) > 0
  rootIndex <- which(go$name == rootName)
  if (length(rootIndex) == 0) {
    stop(sprintf("no such root element found: %s", rootName))
  }
  if (length(rootIndex) > 1) {
    stop(sprintf("more than 1 element found for the root: %s", rootName))
  }
  root <- go$id[rootIndex]
  with_root <- ontologyIndex::get_descendants(go, root)
  o <- which((go$id %in% with_root) & hasParent)

  res <- lapply(o, function(i) {
    data.frame(child = go$id[[i]], parent = go$parents[[i]])
  })
  res <- do.call(rbind, res)
  rownames(res) <- NULL
  res
}


##' Get gene-GO mapping from biomaRt
##'
##' @param biomart the name of the biomart database (default "ensembl").
##' @param dataset the name of the dataset (default "hsapiens_gene_ensembl"),
##'   use `biomaRt::listDatasets` to see the available ones.
##'
##' @return an object to use as a geneSets argument for evoGO constructor
##'   function.
##' @export
getGeneSets <- function(biomart = "ensembl", dataset = "hsapiens_gene_ensembl") {
  ensembl <- biomaRt::useEnsembl(biomart = biomart, dataset = dataset)

  # Split the query into chunks to avoid timeouts
  gene_ids <- biomaRt::getBM(
    attributes = c("ensembl_gene_id"),
    mart = ensembl
  )
  gene_ids <- unname(unlist(gene_ids))
  gene_ids <- split(gene_ids, ceiling(seq_along(gene_ids) / 1000))

  gene_annot <- lapply(1:length(gene_ids), function(i) {
    message("Retrieving GO to gene mapping... (", i, "/", length(gene_ids), ")")
    mart <- biomaRt::getBM(
      attributes = c("go_id", "ensembl_gene_id"),
      filters = "ensembl_gene_id",
      values = gene_ids[[i]],
      mart = ensembl
    )
  })
  gene_annot <- do.call(rbind, gene_annot)

  empty <- gene_annot$go_id == "" | gene_annot$ensembl_gene_id == ""
  gene_annot <- gene_annot[!empty, , drop = FALSE]
  split(gene_annot$ensembl_gene_id, gene_annot$go_id)
}



##' Retrieve and save the latest annotation for species of interest
##'
##' @param species species identifier (example: 'hsapiens' for human)
##' @param database database to be used for acquiring of genesets. Only 'ensembl'
##' is currently supported.
##' @param nCores maximum number of cores to be used
##' @param save logical value indicating whether the result should be saved
##' @param path path to a directory where the result will be saved. If NULL, the file
##' will be saved to the package 'extdata' directory
##' @param deletePrevious logical value indicating whether previous annotation for
##' that species should be removed when newer is saved to save disk space. Only
##' used in case results are saved to 'extdata' directory
##'
##' @return list of evoGO objects containing the latest annotation (one per GO domain).
##' @export
getGOAnnotation <- function(species, database = "ensembl", nCores = 1, save = TRUE,
                            path = NULL, deletePrevious = TRUE) {

  # Check inputs
  assertthat::assert_that(isValidString(species),
    msg = "'species' should be a character vector with length of one"
  )
  assertthat::assert_that(isValidString(database),
    msg = "'database' should be a character vector with length of one"
  )
  assertthat::assert_that(database == "ensembl",
    msg = "Only Ensembl database is currently supported"
  )
  assertthat::assert_that(is.null(path) | isValidString(path),
    msg = "'path' should be a character vector with length of one"
  )

  # Get/create save directory
  if (is.null(path)) {
    package_dir <- path.package("evoGO")
    path <- file.path(package_dir, "extdata")
    if (!dir.exists(path)) {
      dir.create(path = path, recursive = TRUE, showWarnings = FALSE)
    }
  } else {
    assertthat::assert_that(dir.exists(path), msg = "Provided path does not exist")
  }
  if (save) {
    assertthat::assert_that(file.access(path, mode = 2) == 0,
      msg = paste0("No permission to write to ", path, ". Please specify other path.")
    )
  }

  # Get GO data
  message("Looking for latest annotation...")
  go <- getOBO()
  go_annotation <- as.data.frame(go[c("id", "name", "def")])
  goVersion <- strsplit(attr(go, "version")[2], "\\/")[[1]][2]
  out_list <- list()
  go_domains <- c("cellular_component", "biological_process", "molecular_function")

  if (database == "ensembl") {
    # Get geneset data
    message("Waiting for response from Ensembl...")
    ensembl <- biomaRt::useEnsembl(biomart = database)
    ensembl_datasets <- biomaRt::searchDatasets(ensembl, paste0(species, "_gene_ensembl"))
    assertthat::assert_that(!is.null(ensembl_datasets),
      msg = "Ensembl dataset is not found. Make sure that 'species' value format is correct (example: 'hsapiens' for human)"
    )
    dataset <- ensembl_datasets[1, "dataset"]
    dataset_version <- ensembl_datasets[1, "version"]

    # Get file name for the latest annotation and load it if exists
    file_name <- paste0(
      "evogo__go", gsub(
        pattern = "-", replacement = "",
        x = goVersion
      ), "__ensembl_", species, "_", dataset_version, ".rds"
    )
    full_file_name <- file.path(path, file_name)
    if (file.exists(full_file_name)) {
      message("Previously saved up-to-date annotation found")
      return(invisible(readRDS(full_file_name)))
    }

    # Create evoGO object with latest annotation
    genesets <- getGeneSets(biomart = database, dataset = dataset)
    attr(genesets, "Annotation_version") <- dataset_version

    for (domain in go_domains) {
      message("Retrieving graph for '", domain, "'...")

      # Prepare evoGO object
      this_graph <- ontologyIndex2graph(go, rootName = domain)
      attr(this_graph, "GO_version") <- goVersion
      attr(this_graph, "GO_domain") <- domain

      out_list[[domain]] <- evoGO(this_graph, genesets, annotation = go_annotation, nCores = nCores)
    }

    message("New annotation retrieved")
    message("GO release: ", goVersion)
    message("Ensembl version: ", ensembl_datasets[1, "description"])

    # Save results
    if (save) {
      saveRDS(out_list, full_file_name)
      message("Results were successfully saved")

      # Delete previous annotation file to save space if needed (only look in 'extdata')
      if (deletePrevious & path == file.path(path.package("evoGO"), "extdata")) {
        all_annotation_files <- list.files(path = path, pattern = "^evogo__")
        old_annotation_files <- all_annotation_files[
          grepl("\\.rds$", all_annotation_files) &
            grepl(species, all_annotation_files) &
            all_annotation_files != file_name
        ]
        if (length(old_annotation_files) > 0) {
          file.remove(file.path(path, old_annotation_files))
          message(
            "Older annotation files were deleted: ",
            paste0(old_annotation_files, collapse = ", ", ".")
          )
        }
      }
    }
  }

  return(invisible(out_list))
}


##' Load previously retrieved GO annotation as evoGO object
##'
##' @param species species identifier (example: 'hsapiens' for human)
##' @param database database to be used for acquiring of genesets. Only 'ensembl' is
##'  currently supported
##' @param goRelease GO release date in yyyy-mm-dd or yyyymmdd format (e.g., "2022-01-31").
##' If not specified, latest saved version is loaded
##' @param annoVersion gene annotation version (e.g., "GRCh38.p13" for Ensembl). If not specified,
##' latest saved version is loaded
##' @param path directory where evoGO annotations are stored. If NULL, the default package
##' directory is used ('extdata')
##'
##' @return list of evoGO objects containing the latest annotation (one per GO domain).
##' @export
loadGOAnnotation <- function(species, database = "ensembl", goRelease = NULL,
                             annoVersion = NULL, path = NULL) {

  # Check inputs
  assertthat::assert_that(isValidString(species),
    msg = "Argument 'species' should be a character vector with length of one"
  )
  assertthat::assert_that(isValidString(database),
    msg = "Argument 'database' should be a character vector with length of one"
  )
  assertthat::assert_that(database == "ensembl",
    msg = "Only Ensembl database is currently supported"
  )
  assertthat::assert_that(is.null(goRelease) | isValidString(goRelease),
    msg = "Argument 'goRelease' should be a character vector with length of one"
  )
  assertthat::assert_that(is.null(annoVersion) | isValidString(annoVersion),
    msg = "Argument 'annoVersion' should be a character vector with length of one"
  )
  assertthat::assert_that(is.null(path) | isValidString(path),
    msg = "Argument 'path' should be a character vector with length of one"
  )

  # Check path
  if (is.null(path)) {
    package_dir <- path.package("evoGO")
    path <- file.path(package_dir, "extdata")
  }
  assertthat::assert_that(file.access(path, mode = 4) == 0,
    msg = paste0("Path ", path, " is not available")
  )

  # Get file name
  goRelease <- ifelse(is.null(goRelease), "", goRelease)
  annoVersion <- ifelse(is.null(annoVersion), "", annoVersion)
  goVersion <- gsub(pattern = "-", replacement = "", x = goRelease)
  if (goVersion != "" & annoVersion != "") {
    file_name <- paste0("evogo__go", goVersion, "__", database, "_", species, "_", annoVersion, ".rds")
  } else {
    # Assume newer versions of GO data may only have same or newer versions of annotation and vice versa
    all_files <- list.files(path)
    selected_files <- all_files[
      grepl("^evogo__", all_files) &
        grepl("\\.rds$", all_files) &
        grepl(goVersion, all_files) &
        grepl(species, all_files) &
        grepl(annoVersion, all_files)
    ]
    assertthat::assert_that(length(selected_files) > 0, msg = "No annotation files found")

    # Also assume Ensembl genome build or patch number will not reach 100 in this century
    file_name <- sort(selected_files, decreasing = T)[1]
  }
  full_file_name <- file.path(path, file_name)
  assertthat::assert_that(file.access(full_file_name, mode = 4) == 0,
    msg = paste0("File ", full_file_name, " is not available")
  )
  readRDS(full_file_name)
}


##' Download example annotation for Homo sapiens
##'
##' This function obtains an example annotation and saves it to the package 'extdata'
##' directory (a default directory for stored annotations). The purpose of the function is to
##' provide the means for quick testing of the package functions. For performing a proper GO
##' enrichment analysis please acquire an up-to-date annotation using \code{\link{getGOAnnotation}}.
##' 
##' @param path directory where example evoGO annotation will be stored. If NULL, the default
##' package directory is used ('extdata').
##'
##' @return list of evoGO objects containing the example annotation (one per GO domain).
##' @export
exampleGOAnnotation <- function(path = NULL) {
  
  # Check the path
  assertthat::assert_that(is.null(path) | isValidString(path),
    msg = "'path' should be a character vector with length of one"
  )
  if (is.null(path)) {
    package_dir <- path.package("evoGO")
    path <- file.path(package_dir, "extdata")
    if (!dir.exists(path)) {
      dir.create(path = path, recursive = TRUE, showWarnings = FALSE)
    }
  } else {
    assertthat::assert_that(dir.exists(path), msg = "Provided path does not exist")
  }
  assertthat::assert_that(file.access(path, mode = 2) == 0,
    msg = paste0("No permission to write to ", path, ". Please specify other path.")
  )
  
  # Download the annotation
  url <- "https://evogodata.s3.eu-central-1.amazonaws.com/evogo__go20220701__ensembl_hsapiens_GRCh38.p13.rds"
  filename <- file.path(path, basename(url))
  curl::curl_download(url, filename)
  readRDS(filename)
}


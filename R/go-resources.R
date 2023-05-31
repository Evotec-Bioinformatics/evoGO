## This file includes functions for retrieving the annotations from
## external resources and to transform them into the input structures,
## which are required by the evoGO constructor.



##' Retrieve the ontologyIndex object from the obolibrary
##'
##' @param goRelease GO release date in yyyy-mm-dd format (e.g., "2022-01-31").
##' If not specified, latest available version is loaded.
##'
##' @return ontologyIndex object
##' @export
getOBO <- function(goRelease = NULL) {
  if (!is.null(goRelease)) {
    goRelease <- validateGoRelease(goRelease)
  }

  url <- ifelse(is.null(goRelease), "http://purl.obolibrary.org/obo/go/go-basic.obo",
    paste0("http://release.geneontology.org/", goRelease, "/ontology/go-basic.obo")
  )
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

##' Create an input data.frame for evoGO constructor function
##'
##' The data.frame includes "child" - "parent" relationships between the GO ids.
##' @param go an ontologyIndex object returned by getOBO function (also see ontologyIndex package).
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
##' use `biomaRt::listDatasets` to see the available ones.
##' @param ensemblRelease gene annotation version (e.g., 100 for Ensembl). If not specified,
##' latest available version is loaded.
##'
##' @return an object to use as a geneSets argument for evoGO constructor
##'   function.
##' @export
getGeneSets <- function(biomart = "ensembl", dataset = "hsapiens_gene_ensembl", ensemblRelease = NULL) {
  assertthat::assert_that(isValidEnsemblRelease(ensemblRelease))

  ensembl <- biomaRt::useEnsembl(biomart = biomart, dataset = dataset, version = ensemblRelease)

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

##' Parse GO annotation file information
##'
##' Extract information from GO annotation file name.
##'
##' There are three possible GO annotation file naming formats:
##' Previous evoGO v0.1.0 standard annotation format
##' [1] evogo__goYYYYMMDD__ensembl_SP_ENSEMBLV.rds
##' Current evoGO release formats for standard and custom annotation files
##' [2] evogo__goYYYYMMDD__ensemblXX__speciesSP.rds
##' [3] evogo__goYYYYMMDD__customXYZ.rds
##' This function parses a file name into named vector
##' c(goRelease:'YYYY-MM-DD', ensemblRelease:XX, customName:'XYZ', species:'SP', fileName:'file')
##' by splitting the file name at the double underscores '__'.
##' Note that for [1], the species and Ensembl release information ('SP_ENSEMBLV')
##' does not follow current format convention [2], so this information is labeled
##' together as 'XYZ'.
##' Where information is missing, NAs are introduced.
##'
##' @param file_name name of the RDS file which stores the GO annotation
##' @param database database used for acquiring of genesets. Only 'ensembl' is
##'  currently supported.
##' @return named vector
##' @noRd
parseFileInfo <- function(file_name, database = "ensembl") {
  basename <- gsub(pattern = "\\.rds$", "", file_name)
  # split file name at the '__'
  name_split <- strsplit(basename, split = "__")[[1]]
  # all [1], [2], [3] naming formats will have the goRelease
  # find chunk starting with 'go' and remove characters up to 3rd character
  # 'goYYYYMMDD' -> 'YYYY-MM-DD'
  go <- substring(name_split[startsWith(name_split, "go")], 3)
  go <- as.character(as.Date(go, "%Y%m%d"))
  # [1] and [3] will split into 3 chunks
  if (length(name_split) == 3) {
    ensembl <- NA
    if (startsWith(name_split[3], "custom")) {
      # 'customXYZ' -> 'XYZ'
      custom_name <- substring(name_split[3], 7)
    } else {
      # 'ensembl_SP_ENSEMBLV' -> 'SP_ENSEMBLV''
      custom_name <- name_split[3]
    }
    sp <- NA
    # [2] will split into 4 chunks
  } else {
    # 'ensemblXX' -> XX
    ensembl <- name_split[startsWith(name_split, database)]
    ensembl <- ifelse(length(ensembl) == 0, NA, substring(ensembl, 8))
    custom_name <- NA
    # 'speciesSP' -> 'SP'
    sp <- name_split[startsWith(name_split, "species")]
    sp <- ifelse(length(ensembl) == 0, NA, substring(sp, 8))
  }
  info <- c(go, ensembl, custom_name, sp, basename)
  names(info) <- c("goRelease", "ensemblRelease", "customName", "species", "fileName")

  info
}

##' List previously stored GO annotation files
##'
##' List of available files can be filtered using filter arguments.
##'
##' @param species species identifier (example: 'hsapiens' for human).
##' @param goRelease GO release date in yyyy-mm-dd format (e.g., "2022-01-31").
##' @param ensemblRelease numeric gene annotation version (e.g., 100 for Ensembl).
##' @param customName a string that was used to save custom annotation file.
##' @param database database used for acquiring of genesets. Only 'ensembl' is
##'  currently supported.
##' @param path path to the directory where the annotation files were saved.
##' @param returnLatest logical value indicating whether to return information for
##' the file with the latest GO and Ensembl releases (default FALSE). If the list of
##' files contains both standard and custom annotation files, the standard annotation
##' file with the latest GO and Ensembl releases will be selected as custom annotation
##' files do not use the Ensembl database.
##'
##' @return a data.frame containing information about previously stored GO annotations.
##' @export
listGOAnnotations <- function(species = NULL, goRelease = NULL, ensemblRelease = NULL, customName = NULL,
                              database = "ensembl", path = NULL, returnLatest = FALSE) {
  # Check inputs
  assertthat::assert_that(is.null(path) | isValidString(path),
    msg = "Argument 'path' should be a character vector with length of one"
  )
  assertthat::assert_that(is.null(species) | isValidString(species),
    msg = "Argument 'species' should be a character vector with length of one"
  )
  assertthat::assert_that(isValidString(database),
    msg = "Argument 'database' should be a character vector with length of one"
  )
  assertthat::assert_that(isValidEnsemblRelease(ensemblRelease))
  assertthat::assert_that(database == "ensembl",
    msg = "Only Ensembl database is currently supported"
  )
  assertthat::assert_that(is.null(customName) | isValidString(customName),
    msg = "Argument 'customName' should be a character vector with length of one"
  )

  if (!is.null(goRelease)) {
    goRelease <- validateGoRelease(goRelease)
  }

  if (is.null(path)) {
    package_dir <- path.package("evoGO")
    path <- file.path(package_dir, "extdata")
  }
  assertthat::assert_that(file.access(path, mode = 4) == 0,
    msg = paste0("Path ", path, " is not available")
  )

  annotation_files <- list.files(path, pattern = c("^evogo__", "\\.rds$"))
  assertthat::assert_that(length(annotation_files) != 0,
    msg = "No previously saved annotations found"
  )

  annotation_file_info <- list()
  for (i in 1:length(annotation_files)) {
    file_info <- parseFileInfo(annotation_files[i])
    annotation_file_info[[i]] <- file_info
  }

  annotation_file_info <- as.data.frame(do.call(rbind, annotation_file_info))
  annotation_file_info$ensemblRelease <- as.integer(annotation_file_info$ensemblRelease)

  if (!is.null(goRelease)) {
    annotation_file_info <- annotation_file_info[annotation_file_info$goRelease %in% goRelease, ]
  }
  if (!is.null(ensemblRelease)) {
    annotation_file_info <- annotation_file_info[annotation_file_info$ensemblRelease %in% ensemblRelease, ]
  }
  if (!is.null(customName)) {
    annotation_file_info <- annotation_file_info[annotation_file_info$customName %in% customName, ]
  }
  if (!is.null(species)) {
    annotation_file_info <- annotation_file_info[annotation_file_info$species %in% species, ]
  }

  assertthat::assert_that(nrow(annotation_file_info) != 0,
    msg = "No matching annotation files found. Run listGOAnnotations() to view available files"
  )

  filter_annotation_file_info <- annotation_file_info[order(annotation_file_info$goRelease, decreasing = TRUE), ]

  # Select the file with latest releases of GO (and Ensembl)
  if (returnLatest) {
    # if filtered list of files doesn't contain files with Ensembl release information, return file with latest GO release
    if (all(is.na(filter_annotation_file_info$ensemblRelease))) {
      latest_go <- max(filter_annotation_file_info$goRelease)
      latest_go_files <- filter_annotation_file_info[filter_annotation_file_info$goRelease %in% latest_go, ]
      assertthat::assert_that(nrow(latest_go_files) == 1,
        msg = paste("Multiple annotation files with the same releases:", "",
          paste(capture.output(latest_go_files), collapse = "\n"), "",
          "Please provide a more explicit 'customName' argument",
          sep = "\n"
        )
      )
      return(latest_go_files)
      # if any files in the filtered list contain Ensembl release information, select file based on latest GO and Ensembl releases
    } else {
      filter_annotation_file_info <- filter_annotation_file_info[!is.na(filter_annotation_file_info$ensemblRelease), ]
      latest_go <- max(filter_annotation_file_info$goRelease)
      latest_go_files <- filter_annotation_file_info[filter_annotation_file_info$goRelease %in% latest_go, ]
      latest_ensemlb <- max(na.omit(latest_go_files$ensemblRelease))
      latest_file <- latest_go_files[latest_go_files$ensemblRelease %in% latest_ensemlb, ]
      assertthat::assert_that(nrow(latest_file) == 1,
        msg = paste("Multiple annotation files with the same releases:", "",
          paste(capture.output(latest_file), collapse = "\n"), "",
          "Please provide more filter arguments",
          sep = "\n"
        )
      )
      return(latest_file)
    }
  }

  filter_annotation_file_info
}

##' Save GO annotation in RDS file
##'
##' Function saves provided GO annotation to specified path using specified name.
##' Previously stored GO annotation files can be removed to save disk space.
##'
##' @param annotation GO annotation to be saved.
##' @param path path to a directory where the result will be saved.
##' @param file_name name given to the RDS file which stores the GO annotation
##' @param deletePrevious logical value indicating whether previous annotations for
##' that species should be removed when newer is saved to save disk space (default FALSE).
##' @noRd
saveGOAnnotationFile <- function(annotation, path, file_name, deletePrevious) {
  saveRDS(annotation, file.path(path, file_name))
  message("Results were successfully saved")

  # Delete previous annotation file to save space if needed (only look in 'extdata')
  if (deletePrevious & path == file.path(path.package("evoGO"), "extdata")) {
    all_annotation_files <- list.files(path = path, pattern = "^evogo__")
    old_annotation_files <- all_annotation_files[
      grepl("\\.rds$", all_annotation_files) &
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

##' Retrieve and save a GO annotation for species of interest
##'
##' @param species species identifier (example: 'hsapiens' for human).
##' @param goRelease GO release date in yyyy-mm-dd format (e.g., "2022-01-31").
##' If not specified, latest version is retrieved.
##' @param ensemblRelease gene annotation version (e.g., 100 for Ensembl). If not specified,
##' latest version is retrieved.
##' @param database database to be used for acquiring of genesets. Only 'ensembl'
##' is currently supported.
##' @param nCores maximum number of cores to be used.
##' @param save logical value indicating whether the result should be saved (default TRUE).
##' @param path path to a directory where the result will be saved. If NULL, the file
##' will be saved to the package 'extdata' directory.
##' @param deletePrevious logical value indicating whether previous annotations for
##' that species should be removed when newer is saved to save disk space (default FALSE).
##'
##' @return list of evoGO objects containing the annotation (one per GO domain).
##' @export
getGOAnnotation <- function(species = NULL, goRelease = NULL, ensemblRelease = NULL,
                            database = "ensembl", nCores = 1, save = TRUE, path = NULL,
                            deletePrevious = FALSE) {
  # Check inputs
  assertthat::assert_that(isValidString(database),
    msg = "'database' should be a character vector with length of one"
  )
  assertthat::assert_that(database == "ensembl",
    msg = "Only Ensembl database is currently supported"
  )
  assertthat::assert_that(isValidEnsemblRelease(ensemblRelease))
  assertthat::assert_that(isValidString(species),
    msg = "'species' should be a character vector with length of one"
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
  message("Looking for the annotation...")
  go <- getOBO(goRelease = goRelease)
  go_annotation <- as.data.frame(go[c("id", "name", "def")])
  goVersion <- strsplit(attr(go, "version")[2], "\\/")[[1]][2]
  out_list <- list()
  go_domains <- c("cellular_component", "biological_process", "molecular_function")

  # Get geneset data
  message("Waiting for response from Ensembl...")
  ensembl <- biomaRt::useEnsembl(biomart = database)
  ensembl_datasets <- biomaRt::searchDatasets(ensembl, paste0(species, "_gene_ensembl"))
  assertthat::assert_that(!is.null(ensembl_datasets),
    msg = "Ensembl dataset is not found. Make sure that 'species' value format is correct (example: 'hsapiens' for human)"
  )

  dataset <- ensembl_datasets[1, "dataset"]

  if (is.null(ensemblRelease)) {
    ensembl_archive <- biomaRt::listEnsemblArchives()
    database_version <- sub(".*? ", "", ensembl_archive$name[2])
  } else {
    database_version <- ensemblRelease
  }

  # check if file exists
  file_name <- paste0(
    "evogo__go", gsub(
      pattern = "-", replacement = "",
      x = goVersion
    ), "__ensembl", database_version, "__species", species, ".rds"
  )
  full_file_name <- file.path(path, file_name)

  if (file.exists(full_file_name)) {
    message(paste0("Previously saved ", file_name, " annotation found and will be loaded"))
    return(invisible(readRDS(full_file_name)))
  }

  # Get gene-GO mapping
  genesets <- getGeneSets(biomart = database, dataset = dataset, ensemblRelease = ensemblRelease)
  attr(genesets, "Annotation") <- paste0("Ensembl ", database_version)

  for (domain in go_domains) {
    message("Retrieving graph for '", domain, "'...")

    # Prepare evoGO object
    this_graph <- ontologyIndex2graph(go, rootName = domain)
    attr(this_graph, "GO_version") <- goVersion
    attr(this_graph, "GO_domain") <- domain

    out_list[[domain]] <- evoGO(this_graph, genesets, annotation = go_annotation, nCores = nCores)
  }

  message(
    "New annotation retrieved\n",
    "GO release: ", goVersion, "\n",
    "Ensembl version: ", ensembl_datasets[1, "description"], "\n",
    "Ensembl release: ", database_version
  )

  # Save results
  if (save) {
    saveGOAnnotationFile(out_list, path, file_name, deletePrevious)
  }

  return(invisible(out_list))
}


##' Construct and save a GO annotation using a cutom gene-GO mapping
##'
##' @param customAnnotation a named list with character vectors of gene groups per GO term,
##' e.g. list(`GO:01` = c("ENS01", "ENS02"), `GO:02` = c("ENS03")), where
##' the names are the GO term identifiers.
##' The list items may contain the genes belonging to all children terms or
##' include only the genes with exclusion of the children terms.
##' See output of getGeneSets() as format reference.
##' @param customName a string that is used to store custom annotation file.
##' @param goRelease GO release date in yyyy-mm-dd format (e.g., "2022-01-31").
##' If not specified, latest version is retrieved.
##' @param nCores maximum number of cores to be used.
##' @param save logical value indicating whether the result should be saved (default TRUE).
##' @param path path to a directory where the result will be saved. If NULL, the file
##' will be saved to the package 'extdata' directory.
##' @param deletePrevious logical value indicating whether previous annotations for
##' that species should be removed when newer is saved to save disk space (default FALSE).
##'
##' @return list of evoGO objects containing the annotation (one per GO domain).
##' @export
getCustomGOAnnotation <- function(customAnnotation = NULL, customName = NULL, goRelease = NULL,
                                  nCores = 1, save = TRUE, path = NULL, deletePrevious = FALSE) {
  # Check inputs
  assertthat::assert_that(!is.null(customAnnotation) & !is.null(customName),
    msg = "Both 'customAnnotation' and 'customName' must be provided in order to use custom annotation"
  )
  assertthat::assert_that(is.list(customAnnotation) && length(names(customAnnotation)) > 0,
    msg = paste0(
      "Argument 'customAnnotation' should be a named list where each element is named after a GO term ",
      "(eg. \"GO:0000002\") and the element is a character vector containing the names of the annotation terms"
    )
  )
  assertthat::assert_that(isValidString(customName) & nchar(customName) > 0 & !grepl("__", customName, fixed = TRUE),
    msg = "'customName' should be a non-empty character vector with length of one and should not contain \"__\""
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
  message("Looking for the annotation...")
  go <- getOBO(goRelease = goRelease)
  go_annotation <- as.data.frame(go[c("id", "name", "def")])
  goVersion <- strsplit(attr(go, "version")[2], "\\/")[[1]][2]
  out_list <- list()
  go_domains <- c("cellular_component", "biological_process", "molecular_function")

  genesets <- customAnnotation
  customName <- sub(" ", "_", customName)
  attr(genesets, "Annotation") <- customName

  # check if file exists
  file_name <- paste0(
    "evogo__go", gsub(
      pattern = "-", replacement = "",
      x = goVersion
    ), "__custom", customName, ".rds"
  )
  full_file_name <- file.path(path, file_name)

  assertthat::assert_that(!file.exists(full_file_name),
    msg = paste0(
      "Annotation file ", file_name,
      " already exists. Please use a different 'cutomName' to save the annotation ",
      "or use loadGOAnnotation() to load the existing annotation"
    )
  )

  for (domain in go_domains) {
    message("Retrieving graph for '", domain, "'...")

    # Prepare evoGO object
    this_graph <- ontologyIndex2graph(go, rootName = domain)
    attr(this_graph, "GO_version") <- goVersion
    attr(this_graph, "GO_domain") <- domain

    out_list[[domain]] <- evoGO(this_graph, genesets, annotation = go_annotation, nCores = nCores)
  }

  message(
    "New annotation retrieved\n",
    "GO release: ", goVersion, " \n",
    "Custom annotation: ", customName
  )

  # Save results
  if (save) {
    saveGOAnnotationFile(out_list, path, file_name, deletePrevious)
  }

  return(invisible(out_list))
}


##' Load previously saved GO annotation
##'
##' Either a specific GO annotation can be loaded by using filter arguments or an
##' annotation will be selected based on the most recent available GO and Ensembl releases.
##' @param species species identifier (example: 'hsapiens' for human).
##' @param database database used for acquiring of genesets. Only 'ensembl' is
##'  currently supported.
##' @param goRelease GO release date in yyyy-mm-dd format (e.g., "2022-01-31").
##' @param ensemblRelease gene annotation version (e.g., 100 for Ensembl).
##' @param customName a string that was used to store custom annotation file.
##' @param path directory where GO annotations are stored. If NULL, the default package
##' directory is used ('extdata').
##'
##' @return list of evoGO objects containing the latest annotation (one per GO domain).
##' @export
loadGOAnnotation <- function(species = NULL, database = "ensembl", goRelease = NULL,
                             ensemblRelease = NULL, customName = NULL, path = NULL) {
  # Check inputs
  assertthat::assert_that(!is.null(species) || !is.null(customName),
    msg = "Either 'species' or 'customName' argument has to be provided"
  )
  assertthat::assert_that(is.null(species) != is.null(customName),
    msg = "Only one of 'species' and 'customName' arguments can be provided"
  )

  # Check path
  if (is.null(path)) {
    package_dir <- path.package("evoGO")
    path <- file.path(package_dir, "extdata")
  }
  assertthat::assert_that(file.access(path, mode = 4) == 0,
    msg = paste0("Path ", path, " is not available")
  )

  # Find files matching filter arguments
  file_info <- listGOAnnotations(
    species = species, database = database, goRelease = goRelease,
    ensemblRelease = ensemblRelease, custom = customName, path = path
  )

  if (nrow(file_info) == 1) {
    file_name <- file_info$fileName
    # Retrieve latest file
  } else if (nrow(file_info) > 1) {
    message(
      "Multiple annotation file matches:\n\n",
      paste0(capture.output(file_info), collapse = "\n"),
      "\n\nFile with latest releases will be loaded"
    )

    file_name <- listGOAnnotations(
      species = species, database = database, goRelease = goRelease,
      ensemblRelease = ensemblRelease, custom = customName, path = path,
      returnLatest = TRUE
    )$fileName
  }

  message(paste0(file_name, ".rds", " will be loaded"))

  full_file_name <- file.path(path, paste0(file_name, ".rds"))
  assertthat::assert_that(file.access(full_file_name, mode = 4) == 0,
    msg = paste0("File ", full_file_name, " is not available")
  )
  readRDS(full_file_name)
}


##' Download example GO annotation for Homo sapiens
##'
##' This function obtains an example annotation and saves it to the package 'extdata'
##' directory (a default directory for stored annotations). The purpose of the function is to
##' provide the means for quick testing of the package functions. For performing a proper GO
##' enrichment analysis please acquire an up-to-date annotation using \code{\link{getGOAnnotation}}.
##'
##' @param path directory where example GO annotation will be stored. If NULL, the default
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
  url <- "https://evogodata.s3.eu-central-1.amazonaws.com/evogo__go20230611__ensembl110__specieshsapiens.rds"
  filename <- file.path(path, basename(url))
  curl::curl_download(url, filename)
  readRDS(filename)
}

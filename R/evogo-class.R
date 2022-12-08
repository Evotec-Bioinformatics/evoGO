
##' Create an input object for the gene set analysis.
##'
##' @param graphTable a data.frame defining relationship between the GO terms.
##'   The first column must be called "child" and the second is "parent".
##' @param geneSets a list with character vectors of gene groups per a GO term,
##'   e.g. list(`GO:01` = c("ENS01", "ENS02"), `GO:02` = c("ENS03")), where
##'   the names are the GO term identifiers.
##'   The list items may contain the genes belonging to all children terms or
##'   include only the genes with exclusion of the children terms.
##' @param annotation a data.frame with first column "id" with the term identifiers
##' and columns with the annotation information
##' @param minGenes minimal size of a term to keep in the object
##' @param nCores number of cores to use
##' @export
evoGO <- function(graphTable, geneSets, annotation = NULL, minGenes = 3, nCores = 1) {
  if (!is.null(annotation)) {
    if (!"id" %in% colnames(annotation)) {
      stop("The annotation data.frame must have a columns named 'id' with GO term identifiers.")
    }
  }
  g <- new.env(hash = TRUE)
  g_stat <- new.env(hash = TRUE)

  g_stat$nodes <- union(unique(graphTable[, 2]), unique(graphTable[, 1]))
  g_stat$leaves <- setdiff(unique(graphTable[, 1]), unique(graphTable[, 2]))
  g_stat$root <- setdiff(unique(graphTable[, 2]), unique(graphTable[, 1]))

  # Set all attributes
  for (node in g_stat$nodes) {
    g[[node]] <- list(parents = c(), p = 1, degs = c(), weights = c(), ancestors = c())
  }

  # Set parents
  parent_list <- tapply(graphTable[, "parent"], graphTable[, "child"], function(x) {
    x
  })
  for (go in names(parent_list)) {
    g[[go]][["parents"]] <- parent_list[[go]]
  }

  # Set attributed genes (not weights yet)
  for (go_id in names(geneSets)) {
    g[[go_id]]$weights <- as.character(na.omit(geneSets[[go_id]]))
  }

  # Lengths of gene vectors
  g_stat$wt_len <- lengths(sapply(g, "[[", "weights", simplify = FALSE))

  # Remove empty leaves (and update vector with leave names on each iteration)
  nodes_to_remove <- intersect(names(g_stat$wt_len[g_stat$wt_len == 0]), g_stat$leaves)
  g_stat$nodes <- setdiff(g_stat$nodes, nodes_to_remove)
  g_stat$wt_len <- g_stat$wt_len[g_stat$nodes]
  rm(list = nodes_to_remove, envir = g)

  # Get ancestors and distances to root
  ancestors <- parallel::mcmapply(function(go, g) {
    distance <- 0
    parents <- list()
    res <- g[[go]][["parents"]]
    while (length(res) != 0) {
      distance <- distance + 1
      parents[[distance]] <- res
      res <- unlist(sapply(res, function(x) g[[x]][["parents"]]))
    }
    list(parents = unique(unlist(parents)), distance = distance)
  }, g_stat$nodes, MoreArgs = list(g), SIMPLIFY = FALSE, mc.cores = nCores)

  # Set all ancestors
  for (n in g_stat$nodes) {
    g[[n]]$ancestors <- ancestors[[n]]$parents
  }

  # Set all distances
  g_stat$distances <- unlist(sapply(ancestors, "[", "distance"), use.names = FALSE)
  names(g_stat$distances) <- names(ancestors)
  g_stat$distances <- g_stat$distances[g_stat$nodes]

  # Propagate genes
  for (d in sort(unique(g_stat$distances), decreasing = TRUE)) {
    nodes <- names(g_stat$distances[g_stat$distances == d])
    for (go in nodes) {
      for (p in g[[go]]$parents) {
        if (exists(p, envir = g)) {
          g[[p]]$weights <- union(g[[go]]$weights, g[[p]]$weights)
        }
      }
    }
  }
  # Set weights to 1
  for (go in g_stat$nodes) {
    g[[go]]$weights <- setNames(rep(1, length(g[[go]]$weights)), g[[go]]$weights)
  }
  g_stat$wt_len <- lengths(sapply(g, "[[", "weights", simplify = FALSE))

  # Remove the rest of terms where wt_len < minGenes
  nodes_to_remove <- names(g_stat$wt_len[g_stat$wt_len < minGenes])
  g_stat$nodes <- setdiff(g_stat$nodes, nodes_to_remove)
  g_stat$wt_len <- g_stat$wt_len[g_stat$nodes]
  g_stat$distances <- g_stat$distances[g_stat$nodes]
  rm(list = nodes_to_remove, envir = g)

  # Update lists of ancestors
  for (go in g_stat$nodes) {
    g[[go]]$parents <- NULL
    g[[go]]$ancestors <- intersect(g[[go]]$ancestors, g_stat$nodes)
  }

  # Get total number of genes
  g_stat$total_genes <- length(unique(names(unlist(unname(sapply(g, "[[", "weights"))))))

  res_list <- list(g = g, g_stat = g_stat, annotation = annotation)
  attr(res_list, "GO_version") <- attr(graphTable, "GO_version")
  attr(res_list, "GO_domain") <- attr(graphTable, "GO_domain")
  attr(res_list, "Annotation_version") <- attr(geneSets, "Annotation_version")
  class(res_list) <- "evoGO"
  res_list
}


##' @export
print.evoGO <- function(x, ...) {
  cat("evoGO object\n")
  if (!is.null(attr(x, "GO_version"))) {
    cat("GO release: ", attr(x, "GO_version"), "\n")
    cat("GO Data Archive: https://doi.org/10.5281/zenodo.1205166\n")
  }
  if (!is.null(attr(x, "Annotation_version"))) {
    cat("Gene annotation:", attr(x, "Annotation_version"), "\n")
  }
  cat(sprintf("Total genes: %d\n", x$g_stat$total_genes))
  cat(sprintf("Total GO terms: %d\n", length(x$g_stat$nodes)))
  root <- x$g_stat$root
  cat("Root element:\n")
  if (!is.null(x$annotation)) {
    i <- which(x$annotation$id == root)
    if (length(i) == 1) {
      cat(sprintf(
        "|  %s\n",
        strtrim(paste(x$annotation[i, ], collapse = " "), 30)
      ))
    }
  }
}

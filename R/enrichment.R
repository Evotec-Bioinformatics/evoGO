
#' Helper function that removes genes from graph nodes based on the given set of
#' possible gene IDs (updates 'weights' node attribute)
#'
#' @param allIDs character vector of gene IDs
#' @param data object containing GO graph
#' @param minGenes minimum number of genes attributed to GO term for this term
#'   to be kept for the further statistical testing
#' @param nCores maximum number of processes
#'
#' @return updated data object containing GO graph
setUniverse <- function(allIDs, data, minGenes, nCores = 1) {
  g <- data$g
  g_stat <- data$g_stat

  weights <- parallel::mcmapply(function(wts, universe) {
    wts[intersect(names(wts), universe)]
  }, sapply(g, "[[", "weights", simplify = FALSE), MoreArgs = list(allIDs), mc.cores = nCores)

  for (i in 1:length(weights)) {
    g[[names(weights[i])]]$weights <- weights[[i]]
  }

  # Update values and remove terms which became empty
  g_stat$wt_len <- lengths(sapply(g, "[[", "weights", simplify = FALSE))
  nodes_to_remove <- names(g_stat$wt_len[g_stat$wt_len < minGenes])
  g_stat$nodes <- setdiff(g_stat$nodes, nodes_to_remove)
  g_stat$wt_len <- g_stat$wt_len[g_stat$nodes]
  g_stat$distances <- g_stat$distances[g_stat$nodes]
  rm(list = nodes_to_remove, envir = g)
  g_stat$total_genes <- length(unique(names(unlist(unname(sapply(g, "[[", "weights"))))))

  list(g = g, g_stat = g_stat)
}


#' Helper function that assigns differentially expressed genes to graph nodes
#' (updates 'degs' node attribute)
#'
#' @param diffExprIDs character vector of gene IDs
#' @param data object containing GO graph
#' @param nCores maximum number of processes
#'
#' @return updated data object containing GO graph
setDEGs <- function(diffExprIDs, data, nCores = 1) {
  g <- data$g
  g_stat <- data$g_stat

  degs <- parallel::mcmapply(function(wts, degs) {
    intersect(names(wts), degs)
  }, sapply(g, "[[", "weights", simplify = FALSE),
  MoreArgs = list(unique(diffExprIDs)), mc.cores = nCores
  )

  for (i in seq_along(degs)) {
    g[[names(degs[i])]]$degs <- degs[[i]]
  }

  g_stat$deg_len <- lengths(sapply(g, "[[", "degs", simplify = FALSE))
  g_stat$deg_len <- g_stat$deg_len[g_stat$nodes]
  list(g = g, g_stat = g_stat)
}


#' Helper function that calculates (using one-sided Fisher exact test) and
#' assigns p values to nodes (updates 'p' node attribute)
#'
#' @param data object containing GO graph
#' @param nCores maximum number of processes
#'
#' @return updated igraph object containing GO graph
setInitPVals <- function(data, nCores = 1) {
  g <- data$g
  g_stat <- data$g_stat

  n_degs <- length(unique(unlist(unname(sapply(g, "[[", "degs")))))

  p_vals <- parallel::mcmapply(function(term_degs, term_genes, n_degs, n_all_genes) {
    stats::fisher.test(matrix(c(
      term_degs,
      term_genes - term_degs,
      n_degs - term_degs,
      n_all_genes - n_degs - term_genes + term_degs
    ),
    nrow = 2
    ), alternative = "greater")$p
  },
  g_stat$deg_len,
  g_stat$wt_len,
  MoreArgs = list(n_degs, g_stat$total_genes), mc.cores = nCores
  )

  for (i in seq_along(p_vals)) {
    g[[names(p_vals[i])]]$p <- p_vals[i]
  }

  list(g = g, g_stat = g_stat)
}

#' Helper function that calculates coefficients for downweighting of genes of
#' the current node (term) in the ancestor nodes
#'
#' @param go character, GO term ID
#' @param goWeights numeric, named vector containing gene weights for the
#'   current term
#' @param goP numeric, p value of the current term
#' @param k numeric, additional coefficient to be applied to weights when
#'   downweighting genes
#' @param g environment containing all nodes
#'
#' @return list of coefficients for genes of the ancestor terms
getCoefs <- function(go, goWeights, goP, k = 1, g) {
  ancestors <- g[[go]]$ancestors
  ancestor_p <- unlist(sapply(ancestors, function(x) g[[x]][["p"]], USE.NAMES = FALSE))
  ancestor_p[ancestor_p == 0] <- .Machine$double.xmin

  # Keep only less significant ancestors
  # but not those with p = 1, as their down-weigting has no effect
  ancestors <- names(ancestor_p[ancestor_p > goP & ancestor_p < 1])
  if (length(ancestors) == 0) {
    return(NULL)
  }

  # Get sets of coefficients for ancestor nodes
  mapply(function(a, w, coef, genes) {
    w[names(w)] <- 1
    w[intersect(names(w), genes)] <- coef
    w
  },
  ancestors,
  sapply(ancestors, function(x) g[[x]][["weights"]], simplify = FALSE),
  goP / ancestor_p[ancestors] * k,
  MoreArgs = list(names(goWeights)),
  SIMPLIFY = FALSE
  )
}


#' Function that downweights genes of the given term in the less significant
#' ancestor terms and calculates resulting p values for the whole GO term graph
#'
#' @param data object containing GO graph
#' @param k numeric, additional coefficient to be applied to weights when
#'   downweighting genes
#' @param nCores integer, maximum number of processes spawned during
#'   calculations
#'
#' @return updated graph object
getGOscores <- function(data, k = 1, nCores = 1) {
  g <- data$g
  g_stat <- data$g_stat

  # Iteration by distance to root
  for (d in sort(unique(g_stat$distances), decreasing = TRUE)) {

    # message('Working: level ', d)
    nodes <- names(g_stat$distances[g_stat$distances == d])
    nodes <- nodes[sapply(nodes, function(x) g[[x]][["p"]], USE.NAMES = FALSE) < 1]

    if (length(nodes) != 0) {
      # Get lists of gene weight coefficients
      new_coefs <- mapply(getCoefs,
        nodes,
        sapply(nodes, function(x) g[[x]][["weights"]], simplify = FALSE),
        sapply(nodes, function(x) g[[x]][["p"]], USE.NAMES = FALSE),
        MoreArgs = list(k, g),
        SIMPLIFY = FALSE
      )

      # Remove cases when no ancestors need an update
      new_coefs <- new_coefs[lengths(new_coefs) != 0]

      if (length(new_coefs) != 0) {

        # Get final coefs by multiplying coef sets for the same nodes
        new_coefs <- unlist(unname(new_coefs), recursive = FALSE)
        nodes_to_update <- unique(names(new_coefs))

        new_coefs <- mapply(function(go, coefs) {
          Reduce("*", coefs[names(coefs) == go])
        },
        nodes_to_update,
        MoreArgs = list(new_coefs), SIMPLIFY = FALSE
        )


        # Update weights
        for (i in 1:length(new_coefs)) {
          g[[names(new_coefs[i])]]$weights <- new_coefs[[i]] * g[[names(new_coefs[i])]]$weights
        }

        # Update p values
        updated_degs <- sapply(nodes_to_update, function(x) g[[x]][["degs"]], simplify = FALSE)
        updated_weights <- sapply(nodes_to_update, function(x) g[[x]][["weights"]], simplify = FALSE)

        degs_n <- lengths(updated_degs)
        genes_n <- lengths(updated_weights)

        deg_wts <- round(mapply(function(wts, degs) {
          sum(wts[degs])
        }, updated_weights, updated_degs))
        term_wts <- round(sapply(updated_weights, sum))

        total_degs <- length(unique(unlist(unname(sapply(g, "[[", "degs")))))

        p_vals <- parallel::mcmapply(
          function(deg_wts, term_wts, degs_n, genes_n, total_degs, total_genes) {
            stats::fisher.test(matrix(c(
              deg_wts,
              term_wts - deg_wts,
              total_degs - degs_n,
              total_genes - total_degs - genes_n + degs_n
            ),
            nrow = 2
            ), alternative = "greater")$p
          },
          deg_wts, term_wts, degs_n, genes_n,
          MoreArgs = list(total_degs, g_stat$total_genes),
          mc.cores = nCores
        )

        for (i in 1:length(p_vals)) {
          g[[names(p_vals[i])]]$p <- p_vals[i]
        }
      }
    }
  }
  list(g = g, g_stat = g_stat)
}

#' Perform GO enrichment analysis
#'
#' @param goAnnotation evoGO object or list of evoGO objects
#' @param deGenes character vector of differentially expressed genes
#' @param universe character vector of all possible genes for this dataset
#' @param domain character, name of GO domain: 'BP' for biological process, "CC" for
#'   cellular component, "MF" for molecular function.
#' @param minGenes integer indicating the minimal number of genes attributed to
#'   a GO term for this term to be considered during the analysis. This value should
#'   not be lower then 3. Recommended value: 5.
#' @param k numeric, additional coefficient to be applied to weights when
#'   down-weighting genes
#' @param nCores integer, maximum number of cores used during calculations.
#'   If 0, maximum available number of cores is used but not more than 8.
#'
#' @importFrom assertthat assert_that
#'
#' @return data.frame containing the following columns:
#' \itemize{
#'   \item{id - GO term ID}
#'   \item{name - GO term name}
#'   \item{def - GO term definition}
#'   \item{fisher.pvalue - p value obtained using classic overrepresentation analysis 
#'   relying on Fisher's exact test}
#'   \item{evogo.pvalue - p value provided by evoGO algorithm}
#'   \item{annotated - number of genes attributed to a corresponding GO term}
#'   \item{significant - number of differentially expressed genes for a corresponding GO term}
#'}
#' @export
calcGOenrichment <- function(goAnnotation, deGenes = character(), universe = character(),
                             domain = "BP", minGenes = 5, k = 1, nCores = 0) {

  # Annotation check and selection
  go_domains <- c(bp = "biological_process", cc = "cellular_component", mf = "molecular_function")
  assertthat::assert_that(isValidString(domain) & tolower(domain) %in% c("bp", "cc", "mf"),
    msg = "Invalid GO domain. 'BP', 'MF' and 'CC' are supported"
  )
  domain <- go_domains[tolower(domain)]

  assertthat::assert_that(!is.null(goAnnotation), msg = "Annotation cannot be NULL")
  if (inherits(goAnnotation, "evoGO")) {
    if (!is.null(attr(goAnnotation, "GO_domain"))) {
      assertthat::assert_that(attr(goAnnotation, "GO_domain") == domain,
        msg = "Incorrect annotation is provided for selected domain"
      )
    } else {
      warning("GO domain of evoGO object could not be verified")
    }
    evogo <- goAnnotation
  } else {
    names(goAnnotation) <- sapply(goAnnotation, function(x) attr(x, "GO_domain"))
    evogo <- goAnnotation[[domain]]
    assertthat::assert_that(length(evogo) != 0 && inherits(evogo, "evoGO"),
      msg = "Incorrect annotation is provided for selected domain"
    )
  }

  # Other input checks
  assertthat::assert_that(is.numeric(minGenes) && minGenes %% 1 == 0 && minGenes >= 3,
    msg = "Argument 'minGenes' is expected to be an integer with minimum value of 3."
  )
  assertthat::assert_that(is.numeric(nCores) && nCores %% 1 == 0 && length(nCores) == 1 && nCores >= 0,
    msg = "Argument 'nCores' is expected to be a positive integer or 0"
  )
  assertthat::assert_that(is.numeric(k) && length(k) == 1 && k >= 0 & k <= 1,
    msg = "Argument 'k' is expected to be a single numeric value where 0 <= k <= 1"
  )
  assertthat::assert_that(is.character(deGenes),
    msg = "Argument 'deGenes' is expected to be a character vector"
  )
  assertthat::assert_that(is.character(universe),
    msg = "Argument 'universe' is expected to be a character vector"
  )

  if (length(deGenes) == 0) {
    warning("The length of DE genes vector is 0.")
  }

  # Use environments to speed up calculations
  data <- evogo
  data$g <- copyEnvironment(evogo$g)
  data$g_stat <- copyEnvironment(evogo$g_stat)

  # Get default number of cores
  if (nCores == 0) {
    nCores <- min(8, parallel::detectCores())
  }

  # Check universe
  message("Setting universe...")
  if (length(universe) == 0) {
    warning("'universe' was not provided. Using default gene set")
  } else {
    data <- setUniverse(universe, data, minGenes, nCores = nCores)
  }
  data <- setDEGs(deGenes, data, nCores = nCores)

  # Attribute DEGs for terms and calculate p values
  message("Calculating scores...")
  data <- setInitPVals(data, nCores = nCores)
  init_pvalue <- vapply(mget(data$g_stat$nodes, data$g), "[[", "p", FUN.VALUE = numeric(1))
  data <- getGOscores(data, k = k, nCores = nCores)

  # Prepare the output
  g <- data$g
  g_stat <- data$g_stat
  go_tab <- data.frame(
    id = g_stat$nodes,
    fisher.pvalue = init_pvalue,
    evogo.pvalue = vapply(mget(g_stat$nodes, g), "[[", "p", FUN.VALUE = numeric(1))
  )
  ## add annotation of the GO terms if provided
  annotationColumns <- c()
  if (!is.null(evogo$annotation)) {
    go_tab <- merge(go_tab, evogo$annotation, by.x = "id", by.y = "id", all.x = TRUE)
    annotationColumns <- setdiff(colnames(evogo$annotation), "id")
  }
  go_tab$id <- as.character(go_tab$id)
  rownames(go_tab) <- go_tab$id
  go_tab <- go_tab[g_stat$nodes, , drop = FALSE]
  go_tab$annotated <- g_stat$wt_len[g_stat$nodes]
  go_tab$significant <- g_stat$deg_len[g_stat$nodes]
  resultColumns <- c("id", annotationColumns, "fisher.pvalue", "evogo.pvalue", "annotated", "significant")
  go_tab <- go_tab[order(go_tab$evogo.pvalue, decreasing = FALSE), resultColumns]
  rownames(go_tab) <- NULL
  message("Done")
  class(go_tab) <- c("evoGOResult", class(go_tab))
  go_tab
}

##' @export
print.evoGOResult <- function(x, ...) {
  ## Trim the column with the gene lists to a fixed width
  for (i in which(vapply(x, is.character, logical(1)))) {
    x[[i]] <- trimToWidth(x[[i]], 35)
  }
  NextMethod()
}

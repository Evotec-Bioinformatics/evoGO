context("calcGOenrichment")

test_that("significant ontology is the top one", {

  ## a binary tree with node position encoded in the names
  graph <- data.frame(
    child = c(
      "L", "R", "LL", "LR", "RL", "RR", "RRR"
    ),
    parent = c(
      "root", "root", "L", "L", "R", "R", "RR"
    )
  )
  attr(graph, "GO_domain") <- "biological_process"
  go_genes <- lapply(
    setNames(nm = unique(graph$child)),
    function(x) paste0(x, "_gene_", 1:10)
  )
  evogo <- evoGO(graph, go_genes, nCores = 1)
  universe <- unlist(go_genes)
  de_genes <- grep("RRR", universe, value = TRUE)
  res <- calcGOenrichment(evogo, de_genes, universe = universe)
  expect_equal(res$id[1], "RRR")
  ## all RRR are significant
  expect_equal(res$annotated[1], 10)
  expect_equal(res$significant[1], 10)
  ## no significant among left
  expect_equal(res$annotated[res$id == "L"], 30)
  expect_equal(res$significant[res$id == "L"], 0)
  ## half of RR are significant
  expect_equal(res$annotated[res$id == "RR"], 20)
  expect_equal(res$significant[res$id == "RR"], 10)

  ## If a parent ontology is all significant, the child is the next significant
  de_genes <- grep("RR", universe, value = TRUE)
  res <- calcGOenrichment(evogo, de_genes, universe = universe, nCores = 1)
  expect_equal(res$id[1], "RR")
  expect_equal(res$id[2], "RRR")

  ## handle gracefully empty gene sets
  res <- calcGOenrichment(evogo, de_genes, universe = universe, nCores = 1, minGenes = 100)
  expect_equal(nrow(res), 0)
})

test_that("correct graph object is created", {
  ## Here we create a binary tree with the path from the root encoded in the node names
  ##
  ##                     root
  ##                L           R
  ##
  ##             LL    LR     RL  RR
  ##
  ##                                RRR
  graph <- data.frame(
    child = c(
      "L", "R", "LL", "LR", "RL", "RR", "RRR"
    ),
    parent = c(
      "root", "root", "L", "L", "R", "R", "RR"
    )
  )

  go_genes <- lapply(
    setNames(nm = unique(graph$child)),
    function(x) paste0(x, 1:2)
  )
  evogo <- evoGO(graph, go_genes, minGenes = 1, nCores = 1)

  expect_length(evogo$g[["L"]]$ancestors, 0)
  expect_equal(evogo$g[["LL"]]$ancestors, c("L"))
  expect_equal(evogo$g[["RR"]]$ancestors, c("R"))
  expect_setequal(evogo$g[["RRR"]]$ancestors, c("RR", "R"))
  expect_mapequal(evogo$g[["RRR"]]$weights, c(RRR1 = 1, RRR2 = 1))
  expect_mapequal(evogo$g[["RR"]]$weights, c(RRR1 = 1, RRR2 = 1, RR1 = 1, RR2 = 1))
})

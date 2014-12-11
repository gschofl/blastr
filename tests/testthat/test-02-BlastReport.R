context("blastReport/accessors")

hla.xml <- system.file("extdata", "HLA-B-UTR.blast.xml", package = "blastr")
x <- blastReport(hla.xml)

test_that("Sequenctial subsetting works for BlastReport", {
  expect_is(x[1], "IterationList")
  expect_is(x[[1]], "Iteration")
  expect_is(x[[1]][1], "HitList")
  expect_is(x[[1]][[1]], "Hit")
  expect_is(x[[1]][[1]][1], "HspList")
  expect_is(x[[1]][[1]][[1]], "Hsp")
})

test_that("Subsetting an IterationList works as expected", {
  expect_is(x[1], "IterationList")
  expect_equal(length(x[1]), 1)
  expect_equal(length(x[1:2]), 2)
  expect_equal(length(x[1:3]), 2)
  expect_equal(length(x[3]), 0)
  expect_is(getIteration(x), "IterationList")
  expect_equal(length(getIteration(x)), 2)
  expect_equal(length(getIteration(x, 1)), 1)  
  expect_equal(length(getIteration(x, 3)), 0)
  expect_is(getIteration(x, 1), "Iteration")
  expect_is(getIteration(x, 1, drop = FALSE), "IterationList")
})

test_that("'getQueryID' returns query_id", {
  expect_equal(getQueryID(x), c("Query_1", "Query_2"))
  expect_equal(getQueryID(x[2]), "Query_2")
})

test_that("'getQueryDef' returns query_def", {
  expect_equal(getQueryDef(x), c("human HLA-B-5UTR, 100 bp", "human HLA-B-3UTR, 100 bp"))
})



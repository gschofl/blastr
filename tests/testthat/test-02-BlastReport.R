context("blastReport/accessors")

hla.xml <- system.file("extdata", "HLA-B-UTR.blast.xml", package = "blastr")
x <- blastReport(hla.xml)

test_that("Sequenctial subsetting works for BlastReport", {
  expect_is(x[1], "QueryList")
  expect_is(x[[1]], "Query")
  expect_is(x[[1]][1], "HitList")
  expect_is(x[[1]][[1]], "Hit")
  expect_is(x[[1]][[1]][1], "HspList")
  expect_is(x[[1]][[1]][[1]], "Hsp")
})

test_that("Subsetting an QueryList works as expected", {
  expect_is(x[1], "QueryList")
  expect_equal(length(x[1]), 1)
  expect_equal(length(x[1:2]), 2)
  expect_equal(length(x[1:3]), 2)
  expect_equal(length(x[3]), 0)
  expect_is(getQuery(x), "QueryList")
  expect_equal(length(getQuery(x)), 2)
  expect_equal(length(getQuery(x, 1)), 1)  
  expect_equal(length(getQuery(x, 3)), 0)
  expect_is(getQuery(x, 1), "Query")
  expect_is(getQuery(x, 1, drop = FALSE), "QueryList")
})

test_that("'getQueryID' returns query_id", {
  expect_equal(getQueryID(x), c("Query_1", "Query_2"))
  expect_equal(getQueryID(x[2]), "Query_2")
})

test_that("'getQueryDef' returns query_def", {
  expect_equal(getQueryDef(x), c("human HLA-B-5UTR, 100 bp", "human HLA-B-3UTR, 100 bp"))
})



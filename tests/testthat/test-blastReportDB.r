context("blastReport/getters")

db <- blastReportDBConnect("../pmp.db")
load("../pmp.RData")
  

test_that("Sequenctial subsetting works for blastReport", {
  expect_is(bl[1], "IterationList")
  expect_is(bl[[1]], "Iteration")
  expect_is(bl[[1]][1], "HitList")
  expect_is(bl[[1]][[1]], "Hit")
  expect_is(bl[[1]][[1]][1], "HspList")
  expect_is(bl[[1]][[1]][[1]], "Hsp")
})

test_that("Subsetting works for blastReportDB", {
  expect_equal(db_count(db[1], "query"), 1)
  expect_equal(db_count(db[1:3], "query"), 3)
})

test_that("'getQueryID' returns query_id", {
  expect_equal(getQueryID(db), getQueryID(db, 1:3))
  expect_equal(getQueryID(db, 1), getQueryID(db[[1]]))
  expect_equal(getQueryID(bl), c("Query_1", "Query_2", "Query_3"))
  expect_equal(getQueryID(bl[1:2]), c("Query_1", "Query_2"))
})

test_that("'getQueryDef' returns query_def", {
  expect_equal(getQueryDef(db), getQueryDef(db, 1:3))
  expect_equal(getQueryDef(db, 1), getQueryDef(db[[1]]))
})

test_that("'getQueryLen' returns query_len", {
  expect_equal(getQueryLen(db), getQueryLen(db, 1:3))
  expect_equal(getQueryLen(db, 1), getQueryLen(db[[1]]))
  expect_equal(getQueryLen(bl), c(1560, 153, 168))
})

test_that("'getHitID' returns hit_id", {
  expect_identical(getHitID(db), getHitID(db, 1:3))
  expect_is(getHitID(db), "list")
  expect_equal(getHitID(db[[1]]), 1:10)
})

test_that("'getBitscore' returns bitscore", {
  expect_is(getBitscore(db), "list")
  expect_is(getBitscore(db, 1:3), "list")
  expect_equivalent(getBitscore(db, 2), rep(277.202, 10))
})

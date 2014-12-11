context("blastReport/getters")

db <- blastReportDBConnect("../pmp.db")
load("../pmp.RData")
  
test_that("Subsetting a BlastReportDB returns a BlastReportDB", {
  expect_is(db[1], "BlastReportDB")
  expect_is(db[1:2], "BlastReportDB")
})

test_that("Subsetting works for blastReportDB", {
  expect_equal(db_count(db[1], "query"), 1)
  expect_equal(db_count(db[1:3], "query"), 3)
})

test_that("'getQueryID' returns query_id", {
  expect_equal(getQueryID(db), getQueryID(db, 1:3))
  expect_equal(getQueryID(db, 1), getQueryID(db[1]))
})

test_that("'getQueryDef' returns query_def", {
  expect_equal(getQueryDef(db), getQueryDef(db, 1:3))
  expect_equal(getQueryDef(db, 1), getQueryDef(db[1]))
})

test_that("'getQueryLen' returns query_len", {
  expect_equal(getQueryLen(db), getQueryLen(db, 1:3))
  expect_equal(getQueryLen(db, 1), getQueryLen(db[1]))
})


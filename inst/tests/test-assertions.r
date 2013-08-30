context("Testing Assertions")

test_that("'is.empty' returns TRUE on empty vectors of length one", {
  expect_that(is.empty(NULL), is_true())
  expect_that(is.empty(""), is_true())
  expect_that(is.empty(numeric()), is_true())
  
  expect_that(is.empty(list(NULL)), is_false())
  expect_that(is.empty(c("","")), is_false())
})

test_that("'are_empty' returns a logical vector", {
  expect_that(are_empty(NULL), is_true())
  expect_that(are_empty(""), is_true())
  expect_that(are_empty(numeric()), is_true())
  
  expect_that(are_empty(list(NULL)), is_true())
  expect_that(are_empty(list(NULL, NULL)), equals(c(TRUE,TRUE)))
  expect_that(are_empty(c("","")), equals(c(TRUE,TRUE)))
  expect_that(are_empty(list(NA,NULL,"",1)), equals(c(FALSE,TRUE,TRUE,FALSE)))
})

test_that("'all_empty' returns TRUE on empty vectors of all lengths", {
  expect_that(all_empty(NULL), is_true())
  expect_that(all_empty(""), is_true())
  expect_that(all_empty(numeric()), is_true())
  
  expect_that(all_empty(list(NULL, NULL)), is_true())
  expect_that(all_empty(c("","")), is_true())
})

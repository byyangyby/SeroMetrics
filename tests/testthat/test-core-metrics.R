titer <- c(3, 5, 6, 4, 2)
year  <- c(1968, 1972, 1977, 1987, 1997)

test_that("gmt returns numeric scalar", {
  result <- gmt(titer)
  expect_type(result, "double")
  expect_length(result, 1)
})

test_that("auc returns numeric scalar", {
  result <- auc(titer, year, input.log.trans = TRUE)
  expect_type(result, "double")
  expect_length(result, 1)
})

test_that("aty returns numeric scalar", {
  result <- aty(titer, year, input.log.trans = TRUE)
  expect_type(result, "double")
  expect_length(result, 1)
})

test_that("kurt returns numeric scalar", {
  result <- kurt(titer, year, input.log.trans = TRUE)
  expect_type(result, "double")
  expect_length(result, 1)
})

test_that("skew returns numeric scalar", {
  result <- skew(titer, year, input.log.trans = TRUE)
  expect_type(result, "double")
  expect_length(result, 1)
})

test_that("gini returns value between 0 and 1", {
  result <- gini(titer, input.log.trans = TRUE)
  expect_gte(result, 0)
  expect_lte(result, 1)
})

test_that("prop returns value between 0 and 1", {
  result <- prop(titer, threshold = 4)
  expect_gte(result, 0)
  expect_lte(result, 1)
})

test_that("prot_prop returns value between 0 and 1", {
  result <- prot_prop(titer)
  expect_gte(result, 0)
  expect_lte(result, 1)
})

test_that("max_titer returns max of input when log-transformed", {
  result <- max_titer(titer)
  expect_equal(result, max(titer))
})

test_that("width returns non-negative value", {
  result <- width(titer, year)
  expect_gte(result, 0)
})

test_that("var_trans = TRUE and var_trans = 1 give identical gmt results", {
  r1 <- gmt(titer, input.log.trans = TRUE)
  r2 <- gmt(titer, input.log.trans = 1)
  expect_equal(r1, r2)
})

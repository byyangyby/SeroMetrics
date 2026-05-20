make_df <- function(use_char_id = FALSE) {
  id <- if (use_char_id) c("P1", "P1", "P1", "P2", "P2", "P2") else c(1, 1, 1, 2, 2, 2)
  data.frame(
    id             = id,
    isolation_year = c(1968, 1972, 1977, 1968, 1972, 1977),
    titer          = c(3, 5, 6, 4, 2, 3),
    stringsAsFactors = FALSE
  )
}

test_that("gmt_data returns one row per participant (numeric IDs)", {
  df     <- make_df()
  result <- gmt_data(df, part_col = "id", val_col = "titer")
  expect_equal(nrow(result), 2)
  expect_true("GMT" %in% names(result))
})

test_that("gmt_data works with character participant IDs", {
  df     <- make_df(use_char_id = TRUE)
  result <- gmt_data(df, part_col = "id", val_col = "titer")
  expect_equal(nrow(result), 2)
  expect_setequal(result$id, c("P1", "P2"))
})

test_that("aty_data returns one row per participant", {
  df     <- make_df()
  result <- aty_data(df, part_col = "id", weight_col = "isolation_year",
                    val_col = "titer")
  expect_equal(nrow(result), 2)
})

test_that("auc_data returns one row per participant", {
  df     <- make_df()
  result <- auc_data(df, part_col = "id", weight_col = "isolation_year",
                    val_col = "titer")
  expect_equal(nrow(result), 2)
})

test_that("kurtosis_data: mode = 'max' no longer errors (roung regression)", {
  df <- data.frame(
    id             = c(1, 1, 1, 1),
    isolation_year = c(1968, 1968, 1972, 1972),
    titer          = c(3, 5, 4, 6)
  )
  expect_no_error(
    kurtosis_data(df, part_col = "id", weight_col = "isolation_year",
                 val_col = "titer", mode = "max")
  )
})

test_that("skewness_data: mode = 'max' no longer errors (roung regression)", {
  df <- data.frame(
    id             = c(1, 1, 1, 1),
    isolation_year = c(1968, 1968, 1972, 1972),
    titer          = c(3, 5, 4, 6)
  )
  expect_no_error(
    skewness_data(df, part_col = "id", weight_col = "isolation_year",
                 val_col = "titer", mode = "max")
  )
})

test_that("summarise_landscape: unknown metric throws informative error", {
  df <- make_df()
  expect_error(
    summarise_landscape(df, part_col = "id", weight_col = "isolation_year",
                val_col = "titer", required_metrics = c("GMT", "aty")),
    regexp = "Unknown metric"
  )
})

test_that("summarise_landscape: does not bleed Overall into caller environment", {
  df <- make_df()
  rm_if_exists <- function(x) if (exists(x, envir = .GlobalEnv)) rm(list = x, envir = .GlobalEnv)
  rm_if_exists("Overall")
  summarise_landscape(df, part_col = "id", weight_col = "isolation_year",
              val_col = "titer", required_metrics = c("GMT", "AUC"))
  expect_false(exists("Overall", envir = .GlobalEnv))
})

test_that("summarise_landscape returns correct column count for two metrics", {
  df     <- make_df()
  result <- summarise_landscape(df, part_col = "id", weight_col = "isolation_year",
                        val_col = "titer", required_metrics = c("GMT", "AUC"))
  expect_equal(ncol(result), 3)  # id + GMT + AUC
})

test_that("uniq deduplicates rows", {
  df <- data.frame(a = c(1, 1, 2), b = c("x", "x", "y"),
                   stringsAsFactors = FALSE)
  expect_equal(nrow(uniq(df)), 2)
})

test_that("uniq subset deduplicates selected columns only", {
  df <- data.frame(a = c(1, 1, 2), b = c("x", "x", "y"),
                   stringsAsFactors = FALSE)
  expect_equal(nrow(uniq(df, columns = "a")), 2)
})

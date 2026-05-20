## Run this script once to build the package data objects from the raw xlsx files.
## Requires: readxl, dplyr, tidyr, stringr, usethis

library(readxl)
library(dplyr)
library(tidyr)
library(stringr)

transform <- function(df) {
  df |>
    pivot_longer(
      cols = 7:ncol(df),
      names_to = "strain",
      values_to = "titer"
    ) |>
    mutate(strain = case_when(
      str_detect(strain, "_MDCK$") ~ str_remove(strain, "_MDCK"),
      TRUE ~ strain
    )) |>
    mutate(strain = case_when(
      str_detect(strain, "_EGG$") ~ str_remove(strain, "_EGG"),
      TRUE ~ strain
    )) |>
    mutate(strain = case_when(
      str_detect(strain, "/02$") ~ str_replace(strain, "/02$", "/2002"),
      TRUE ~ strain
    )) |>
    mutate(titer = as.character(titer)) |>
    mutate(titer = case_when(
      str_detect(titer, "\\*")    ~ str_replace(titer, "\\*", "10"),
      TRUE ~ titer
    )) |>
    mutate(titer = case_when(
      str_detect(titer, "<10")    ~ str_replace(titer, "<10", "5"),
      TRUE ~ titer
    )) |>
    mutate(titer = case_when(
      str_detect(titer, ">=1280") ~ str_replace(titer, ">=1280", "9"),
      TRUE ~ titer
    )) |>
    mutate(titer = as.numeric(titer)) |>
    mutate(
      isolation_year = case_when(
        str_detect(str_sub(strain, -4), "^\\d{4}$") ~ as.numeric(str_sub(strain, -4)),
        TRUE ~ as.numeric(str_sub(strain, -2)) + 1900
      )
    ) |>
    filter(!is.na(titer))
}

Fonville <- read_excel("data-raw/Fonville.xlsx") |> transform()
Lessler  <- read_excel("data-raw/Lessler.xlsx")  |> transform()

usethis::use_data(Fonville, Lessler, overwrite = TRUE)

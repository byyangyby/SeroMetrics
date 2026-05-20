## Run this script once to build the package data objects from the raw xlsx files.
## Requires: readxl, dplyr, tidyr, stringr, usethis

library(readxl)
library(dplyr)
library(tidyr)
library(stringr)

## Fonville is in wide format: columns 1-6 are sample metadata and columns
## 7+ are one column per strain. Pivot to long, clean the titer cell
## conventions, and extract the isolation year from the strain name.
transform_fonville <- function(df) {
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

## Lessler is already in long format: each row is one (id, strain) pair with
## a single (log-transformed) titer in the `titers` column. The strain name
## sits in `neut.against`, formatted like "1 A/HK/1968(H3N2)" — a leading
## index, then the strain string with the isolation year embedded.
transform_lessler <- function(df) {
  df |>
    mutate(
      strain = str_remove(neut.against, "^\\d+\\s+"),
      isolation_year = as.numeric(str_extract(neut.against, "\\d{4}")),
      titer = titers,
      is.vac = na_if(is.vac, "NA")
    ) |>
    select(id, age, is.vac, shift.age, strain, isolation_year, titer) |>
    filter(!is.na(titer))
}

Fonville <- read_excel("data-raw/Fonville.xlsx") |> transform_fonville()
Lessler  <- read_excel("data-raw/Lessler.xlsx")  |> transform_lessler()

usethis::use_data(Fonville, Lessler, overwrite = TRUE)

transform <- function(df){
  transformed_df <- df %>%
    pivot_longer(
      cols = 7:ncol(df),
      names_to = "strain",
      values_to = "titer"
    ) %>%
    mutate(strain = case_when(
      str_detect(strain, "_MDCK$") ~ str_remove(strain, "_MDCK"),
      TRUE ~ strain
    ))%>%
    mutate(strain = case_when(
      str_detect(strain, "_EGG$") ~ str_remove(strain, "_EGG"),
      TRUE ~ strain))%>%
    mutate(strain = case_when(
      str_detect(strain, "/02$") ~ str_replace(strain, "/02$", "/2002"),
      TRUE ~ strain))%>%
    #processing NA values, not sure if it is correct
    mutate(titer = case_when(
      str_detect(titer, "\\*") ~ str_replace(titer, "\\*", "10"),
      TRUE ~ titer))%>%
    mutate(titer = case_when(
      str_detect(titer, "<10") ~ str_replace(titer, "<10", "5"),
      TRUE ~ titer))%>%
    mutate(titer = case_when(
      str_detect(titer, ">=1280") ~ str_replace(titer, ">=1280", "9"),
      TRUE ~ titer))%>%
    mutate(titer = as.numeric(titer))
  
  
  transformed_df <- transformed_df %>%
    mutate(
      isolation_year = case_when(
        str_detect(str_sub(strain, -4), "^\\d{4}$") ~ as.numeric(str_sub(strain, -4)),
        TRUE ~ as.numeric(str_sub(strain, -2)) + 1900
      )
    )
  
  return(transformed_df)
}
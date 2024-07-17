# SeroMetrics

## I. Introduction 
### Data Format

#### -  participant id: mandatory. This should be the id of the sampled subject.

#### -  value: mandatory. This should be the investigated value, usually being titers.

#### -  distance: mandatory. This should be the absolute distance among virus strains, usually being isolation year or AA distance.

#### -  strain: optional. This should be the virus strains. When distance is not provided, one may need to derive isolation year from the information in virus strain.

#### -  group: optional. This should be the group of some sample of a participant. Sometimes, one participant may take tests on the same strains multiple times, and the samples collected from each time forms a sample group.

#### -  other variables: optional. There might be age, sex, and other information regarding one sample. At most times, these variables will not affect the final result of the metrics.

### Data Simulation (Sample codes)

```r

    library(dplyr)

    set.seed(12986)

    id <- rep(1:30, each = 2)

    unique_birth_year <- sample(1937:2014, size = 30, replace = TRUE)

    birth_year <- rep(unique_birth_year, each = 2)

    sample_year <- rep(c(2020, 2024), times = 30)

    df <- data.frame(id = id, birth_year = birth_year, sample_year = sample_year)

    unique_isolation_years <- sample(1950:2020, size = 20, replace = TRUE)

    df <- df %>%
      slice(rep(1:n(), each = 20)) %>%
      mutate(isolation_year = rep(unique_isolation_years, times = n() / 20))

    generate_log_titer <- function(birth_year, isolation_year) {
      if(isolation_year<1960){
        if(birth_year > 2005){
          return(sample(-1:2,1))
        }else if(birth_year >=1995 & birth_year < 2005){
          return(sample(2:5,1))
        }else if(birth_year >=1955 & birth_year < 1995){
          return(sample(5:6,1))
        }else{
          return(sample(4:6,1))
        }
      }else if(isolation_year >= 1960 & isolation_year<=1970){
        if(birth_year > 2005){
          return(sample(-1:2,1))
        }else if(birth_year >=1995 & birth_year < 2005){
          return(sample(3:5,1))
        }else if(birth_year >=1955 & birth_year < 1995){
          return(sample(6:7,1))
        }else{
          return(sample(5:6,1))
        }
      }else if(isolation_year >= 1970 & isolation_year<=1995){
        if(birth_year > 2005){
          return(sample(-1:2,1))
        }else if(birth_year >=1995 & birth_year < 2005){
          return(sample(2:6,1))
        }else if(birth_year >=1955 & birth_year < 1995){
          return(sample(4:5,1))
        }else{
          return(sample(3:4,1))
        }
      }else if(isolation_year >= 1995 & isolation_year<=2005){
        if(birth_year > 2005){
          return(sample(0:4,1))
        }else if(birth_year >=1995 & birth_year < 2005){
          return(sample(6:7,1))
        }else if(birth_year >=1955 & birth_year < 1995){
          return(sample(4:7,1))
        }else{
          return(sample(3:5,1))
        }
      }else if(isolation_year >= 2005){
        if(birth_year > 2005){
          if(isolation_year > birth_year){
            return(sample(5:8,1))
          }else{
            return(sample(2:6,1))
          }
        }else if(birth_year >=1995 & birth_year < 2005){
          return(sample(7:8,1))
        }else if(birth_year >=1955 & birth_year < 1995){
          return(sample(2:5,1))
        }else{
          return(sample(0:4,1))
        }
      }
    }

    df <- df %>%
      rowwise() %>%
      mutate(log_titer = generate_log_titer(birth_year, isolation_year)) %>%
      ungroup()

    head(df,10)
```
```r
    Output:
    ## # A tibble: 10 × 5
    ##       id birth_year sample_year isolation_year log_titer
    ##    <int>      <int>       <dbl>          <int>     <int>
    ##  1     1       1952        2020           1980         3
    ##  2     1       1952        2020           1952         4
    ##  3     1       1952        2020           1974         3
    ##  4     1       1952        2020           1975         4
    ##  5     1       1952        2020           1992         3
    ##  6     1       1952        2020           2001         3
    ##  7     1       1952        2020           2016         4
    ##  8     1       1952        2020           1978         3
    ##  9     1       1952        2020           1960         5
    ## 10     1       1952        2020           1969         5
```
---
## II. Package Build-Up

### Step 1: Clear up the global environment.

### Step 2: Download the R folder and store it in a folder called "SeroMetrics".

### Step 3: Run the codes below.

``` r
# Download necessary packages if you have not installed them in your PC
# install.packages("roxygen2")
# install.packages("devtools")
# install.packages("usethis")

library(roxygen2)
library(devtools)
library(usethis)

# input your own path to the SeroMetrics folder
create_package("path/to/SeroMetrics")
roxygenise("path/to/SeroMetrics")
build("path/to/SeroMetrics")
install("path/to/SeroMetrics")

library(SeroMetrics)
```
---
## III. Examples of real datasets

#### To run the examples, you may download the data folder into the SeroMetrics folder on your PC.

### 1. Lessler

#### This is a basic example with a simple dataset, in which for each participant, there is only one testing result on each virus strain, i.e.”group\_col” can be NULL, and all the virus strains have a different isolation year, i.e. “mode” can be NULL.

##### The first method is to use the overallData function directly.
```r
    library(readxl)
    library(stringr)
    library(DescTools)
    library(moments)
    library(dplyr)

    #load data
    df <- read_excel("SeroMetrics/data/Lessler.xlsx")
    #data pre-processing
    df$isolation_year <- as.numeric(str_extract(df$neut.against, "\\d{4}"))
    df$titers <- floor(df$titers)
    #show the basic structure of the dataset
    head(df,10)
```

```r
    Output:
    ## # A tibble: 10 × 8
    ##      age is.vac shift.age titers neut.against         id   loc isolation_year
    ##    <dbl> <chr>      <dbl>  <dbl> <chr>             <dbl> <dbl>          <dbl>
    ##  1    75 NA            34      5 1 A/HK/1968(H3N2)     1     1           1968
    ##  2    35 NA            -6      1 1 A/HK/1968(H3N2)     2     1           1968
    ##  3    71 1             30      4 1 A/HK/1968(H3N2)     3     1           1968
    ##  4    65 1             24      3 1 A/HK/1968(H3N2)     4     1           1968
    ##  5    64 1             23      5 1 A/HK/1968(H3N2)     5     1           1968
    ##  6    33 5             -8      2 1 A/HK/1968(H3N2)     6     1           1968
    ##  7    60 1             19      5 1 A/HK/1968(H3N2)     7     1           1968
    ##  8    43 1              2      5 1 A/HK/1968(H3N2)     8     1           1968
    ##  9    17 1            -24      1 1 A/HK/1968(H3N2)     9     1           1968
    ## 10    36 1             -5      1 1 A/HK/1968(H3N2)    10     1           1968
```

###### You can choose the metrics as you want. However, there are some restrictions when inputting the arguments. For details, check the explanation of the function overallData.
```r
    overallResult <- overallData(df, part_col = "id", weight_col = "isolation_year", val_col = "titers", var_trans = 1, gmt_output_trans = 1, required_metrics = c("AUC","ATY","GMT","Prot_prop","Width0","Width2","Width3"), col_names = c("auc","aty","gmt","prot_prop","Width10","Width40","Width80"))
    head(overallResult, 10)
```
```r
    Output:
    ##    id   auc      aty      gmt  prot_prop   Width10   Width40 Width80
    ## 1   1  63.0 1988.754 3.555556 0.11111111 0.8916667 0.0000000  0.0000
    ## 2   2  93.5 1985.372 2.777778 0.19444444 0.8812500 0.2687500  0.1000
    ## 3   3  55.5 1987.554 2.555556 0.12777778 0.9083333 0.0250000  0.0000
    ## 4   4  19.0 1980.289 1.444444 0.02222222 0.1750000 0.0000000  0.0000
    ## 5   5  80.5 1984.481 2.777778 0.16666667 0.8750000 0.0875000  0.0000
    ## 6   6 115.5 1990.405 3.555556 0.35000000 0.8000000 0.6416667  0.6000
    ## 7   7  26.0 1983.577 1.777778 0.08333333 0.3312500 0.0437500  0.0000
    ## 8   8  85.0 1988.106 4.222222 0.17777778 0.8666667 0.1750000  0.0000
    ## 9   9  66.5 1996.665 2.777778 0.23333333 0.4750000 0.2916667  0.2375
    ## 10 10 125.5 1989.727 4.555556 0.49444444 0.9650000 0.5345833  0.3725
```

##### The second method is to calculate the metrics one by one, and join them manually.
```r

    library(readxl)
    library(stringr)
    library(DescTools)
    library(moments)
    library(dplyr)

    #load data
    df <- read_excel("SeroMetrics/data/Lessler.xlsx")
    #data pre-processing
    df$isolation_year <- as.numeric(str_extract(df$neut.against, "\\d{4}"))
    df$titers <- floor(df$titers)
    #show the basic structure of the dataset
    head(df,10)
```
```r
    Output:
    ## # A tibble: 10 × 8
    ##      age is.vac shift.age titers neut.against         id   loc isolation_year
    ##    <dbl> <chr>      <dbl>  <dbl> <chr>             <dbl> <dbl>          <dbl>
    ##  1    75 NA            34      5 1 A/HK/1968(H3N2)     1     1           1968
    ##  2    35 NA            -6      1 1 A/HK/1968(H3N2)     2     1           1968
    ##  3    71 1             30      4 1 A/HK/1968(H3N2)     3     1           1968
    ##  4    65 1             24      3 1 A/HK/1968(H3N2)     4     1           1968
    ##  5    64 1             23      5 1 A/HK/1968(H3N2)     5     1           1968
    ##  6    33 5             -8      2 1 A/HK/1968(H3N2)     6     1           1968
    ##  7    60 1             19      5 1 A/HK/1968(H3N2)     7     1           1968
    ##  8    43 1              2      5 1 A/HK/1968(H3N2)     8     1           1968
    ##  9    17 1            -24      1 1 A/HK/1968(H3N2)     9     1           1968
    ## 10    36 1             -5      1 1 A/HK/1968(H3N2)    10     1           1968
```
###### Choose the metrics as you want.
```r
    atyResult <- atyData(df,part_col = c("id"),weight_col = c("isolation_year"), val_col = c("titers"), var_trans = 1)
    aucResult <- aucData(df,part_col = c("id"),weight_col = c("isolation_year"), val_col = c("titers"), var_trans = 1)
    giniResult <- giniData(df,part_col = c("id"), val_col = c("titers"), var_trans = 1, gini_col = "Gini")
    gmtResult <- gmtData(df,part_col = c("id"), val_col = c("titers"), var_trans = 1,output_trans = 1)
    kurtosisResult <- kurtosisData(df,part_col = c("id"), weight_col = c("isolation_year"), val_col = c("titers"), var_trans = 1)
    max_titerResult <- max_titerData(df,part_col = c("id"), val_col = c("titers"), var_trans = 1,output_trans = 1)
    propResult2 <- propData(df,part_col = c("id"), val_col = c("titers"), var_trans = 1,threshold_trans = 1,threshold = 2, prop_col = "Prop2") #the proportion of titers above 2
    propResult3 <- propData(df,part_col = c("id"), val_col = c("titers"), var_trans = 1,threshold_trans = 1,threshold = 3, prop_col = "Prop3") #the proportion of titers above 3
    prot_propResult <- prot_propData(df,part_col = c("id"), val_col = c("titers")) #using default weight vector
    skewnessResult <- skewnessData(df,part_col = c("id"), weight_col = c("isolation_year"), val_col = c("titers"), var_trans = 1)
    widthResult <- widthData(df,part_col = c("id"), weight_col = c("isolation_year"), val_col = c("titers"), var_trans = 1, threshold_trans = 1,threshold = 2, width_col = "Width2") #the width of titers above 2

    overallResult <- atyResult %>%
      left_join(aucResult, by = "id") %>%
      left_join(giniResult, by = "id") %>%
      left_join(gmtResult, by = "id") %>%
      left_join(kurtosisResult, by = "id") %>%
      left_join(max_titerResult, by = "id") %>%
      left_join(propResult2, by = "id") %>%
      left_join(propResult3, by = "id") %>%
      left_join(prot_propResult, by = "id") %>%
      left_join(skewnessResult, by = "id") %>%
      left_join(widthResult, by = "id")

    head(overallResult,10)
```
```r
    Output:
    ##    id      ATY   AUC      Gini      GMT Kurtosis Max_titer     Prop2     Prop3  Prot_prop   Skewness    Width2
    ## 1   1 1988.754  63.0 0.1640625 3.555556 1.687270         5 1.0000000 0.8888889 0.11111111 -0.3850070 0.0000000
    ## 2   2 1985.372  93.5 0.2426471 2.777778 1.601629         5 0.7777778 0.4444444 0.19444444  0.3521357 0.2687500
    ## 3   3 1987.554  55.5 0.1718750 2.555556 1.545701         4 0.8888889 0.3333333 0.12777778 -0.4852686 0.0250000
    ## 4   4 1980.289  19.0 0.1477273 1.444444 2.134154         3 0.3333333 0.1111111 0.02222222  0.9432719 0.0000000
    ## 5   5 1984.481  80.5 0.1985294 2.777778 1.668338         5 0.8888889 0.4444444 0.16666667  0.2886074 0.0875000
    ## 6   6 1990.405 115.5 0.2073171 3.555556 2.516994         5 0.8888889 0.5555556 0.35000000 -0.7167264 0.6416667
    ## 7   7 1983.577  26.0 0.2200000 1.777778 1.191072         5 0.4444444 0.1111111 0.08333333  0.3490311 0.0437500
    ## 8   8 1988.106  85.0 0.1184211 4.222222 1.512408         5 1.0000000 0.8888889 0.17777778 -0.3804884 0.1750000
    ## 9   9 1996.665  66.5 0.2941176 2.777778 2.364882         6 0.6666667 0.3333333 0.23333333 -0.5062919 0.2916667
    ## 10 10 1989.727 125.5 0.2150000 4.555556 1.815939         7 0.8888889 0.7777778 0.49444444 -0.7370504 0.5345833
```
---
### 2. Fonville

#### This is a more complicated example with the dataset called “Fonville”, in which for each participant, there are multiple testing results, one for each sample year, i.e. “group\_col” is sample year, and there are different virus strains share the same isolation\_year, i.e. “mode” is necessary. Besides, the dataset has a particular structure with each virus strain forming one column, so data transformation should be done first.

##### The first method is to use the overallData function directly.
```r
    library(readxl)
    library(stringr)
    library(DescTools)
    library(moments)
    library(dplyr)
    library(tidyverse) #this is for data transformation

    #This transformation R-script specifically works for the Fonville dataset. If one's dataset has a structure different than the standard one, one may need to write the transformation codes by oneself, and use the package after transforming the data.
    source('SeroMetrics/data/transform.R')

    #load valid data
    df <- read_excel("SeroMetrics/data/Fonville.xlsx", n_max = 324)
    #data transformation
    transformed_df <- transform(df)

    head(transformed_df,10)
```
```r
    Output:
    ## # A tibble: 10 × 9
    ##    `Subject Number` `Sample Year` `Year of Birth` Sample  `PCR Results` `Row in Fonville Fig S15` strain      titer isolation_year
    ##               <dbl>         <dbl>           <dbl> <chr>   <chr>                             <dbl> <chr>       <dbl>          <dbl>
    ##  1                1          2007            1957 H3 PCR+ <NA>                                 63 BI/16190/68    10           1968
    ##  2                1          2007            1957 H3 PCR+ <NA>                                 63 BI/21793/72    10           1972
    ##  3                1          2007            1957 H3 PCR+ <NA>                                 63 BI/1761/76     10           1976
    ##  4                1          2007            1957 H3 PCR+ <NA>                                 63 BI/2271/76     10           1976
    ##  5                1          2007            1957 H3 PCR+ <NA>                                 63 NL/233/82      10           1982
    ##  6                1          2007            1957 H3 PCR+ <NA>                                 63 NL/620/89      10           1989
    ##  7                1          2007            1957 H3 PCR+ <NA>                                 63 NL/823/92      10           1992
    ##  8                1          2007            1957 H3 PCR+ <NA>                                 63 NL/179/93     160           1993
    ##  9                1          2007            1957 H3 PCR+ <NA>                                 63 JO/33/94       80           1994
    ## 10                1          2007            1957 H3 PCR+ <NA>                                 63 SD/9/93        10           1993
```

###### Choose the metrics as you want.
```r
    overallResult <- overallData(transformed_df, part_col = "Subject Number", weight_col = "isolation_year", val_col = "titer", group_col = "Sample Year", var_trans = 0, max_titer_output_trans = 1, mode = "mean", required_metrics = c("gini_coefficient","Kurtosis","Skewness","Prop0","Prop2","Prop3"), col_names = c("Gini","Kurt","Skew","prop0","prop2","prop3"))
    head(overallResult, 10)
```
```r
    Output:
    ##    Subject Number Sample Year      Gini     Kurt       Skew     prop0     prop2     prop3
    ## 1               1        2007 0.2250316 5.431954 -1.2607987 0.8771930 0.3333333 0.2631579
    ## 2               1        2008 0.2461859 1.960863  0.4710923 1.0000000 0.4385965 0.3333333
    ## 3               1        2009 0.2415966 1.878263  0.4097566 1.0000000 0.4210526 0.3157895
    ## 4               1        2010 0.2403342 4.798016 -1.1516109 0.7543860 0.3508772 0.2105263
    ## 5               1        2011 0.2407295 5.031745 -1.1495281 0.8070175 0.4210526 0.3157895
    ## 6               1        2012 0.1727029 5.417621 -1.1458572 0.9649123 0.2456140 0.1578947
    ## 7               2        2007 0.2097812 5.290075 -1.3655160 0.8947368 0.2631579 0.2105263
    ## 8               2        2008 0.2664053 1.904847 -0.1764279 1.0000000 0.5964912 0.4912281
    ## 9               2        2009 0.2621916 1.918048 -0.1420216 1.0000000 0.6491228 0.4736842
    ## 10              2        2010 0.2060224 5.711489 -1.5264905 0.9824561 0.3859649 0.2982456
```

##### The second method is to calculate the metrics one by one, and join them manually.
```r
    library(readxl)
    library(stringr)
    library(DescTools)
    library(moments)
    library(dplyr)
    library(tidyverse) #this is for data transformation

    #This transformation R-script specifically works for the Fonville dataset. If one's dataset has a structure different than the standard one, one may need to write the transformation codes by oneself, and use the package after transforming the data.
    source('SeroMetrics/data/transform.R')

    #load valid data
    df <- read_excel("SeroMetrics/data/Fonville.xlsx", n_max = 324)
    #data transformation
    transformed_df <- transform(df)

    head(transformed_df,10)
```
```r
    Output:
    ## # A tibble: 10 × 9
    ##    `Subject Number` `Sample Year` `Year of Birth` Sample  `PCR Results` `Row in Fonville Fig S15` strain      titer isolation_year
    ##               <dbl>         <dbl>           <dbl> <chr>   <chr>                             <dbl> <chr>       <dbl>          <dbl>
    ##  1                1          2007            1957 H3 PCR+ <NA>                                 63 BI/16190/68    10           1968
    ##  2                1          2007            1957 H3 PCR+ <NA>                                 63 BI/21793/72    10           1972
    ##  3                1          2007            1957 H3 PCR+ <NA>                                 63 BI/1761/76     10           1976
    ##  4                1          2007            1957 H3 PCR+ <NA>                                 63 BI/2271/76     10           1976
    ##  5                1          2007            1957 H3 PCR+ <NA>                                 63 NL/233/82      10           1982
    ##  6                1          2007            1957 H3 PCR+ <NA>                                 63 NL/620/89      10           1989
    ##  7                1          2007            1957 H3 PCR+ <NA>                                 63 NL/823/92      10           1992
    ##  8                1          2007            1957 H3 PCR+ <NA>                                 63 NL/179/93     160           1993
    ##  9                1          2007            1957 H3 PCR+ <NA>                                 63 JO/33/94       80           1994
    ## 10                1          2007            1957 H3 PCR+ <NA>                                 63 SD/9/93        10           1993
```
###### Choose the metrics as you want.
```r
    atyResult <- atyData(transformed_df,part_col = c("Subject Number"),weight_col = c("isolation_year"), val_col = c("titer"), group_col = c("Sample Year"),var_trans = 0, aty_col = "aty")
    aucResult <- aucData(transformed_df,part_col = c("Subject Number"),weight_col = c("isolation_year"), val_col = c("titer"), group_col = c("Sample Year"),var_trans = 0,mode="mean", auc_col = "auc")
    giniResult <- giniData(transformed_df,part_col = c("Subject Number"), val_col = c("titer"), group_col = c("Sample Year"),var_trans = 0)
    gmtResult <- gmtData(transformed_df,part_col = c("Subject Number"), val_col = c("titer"), group_col = c("Sample Year"),var_trans = 0,output_trans = 1)
    kurtosisResult <- kurtosisData(transformed_df,part_col = c("Subject Number"),weight_col = c("isolation_year"), val_col = c("titer"), group_col = c("Sample Year"),var_trans = 0)
    max_titerResult <-max_titerData(transformed_df,part_col = c("Subject Number"), val_col = c("titer"),group_col = c("Sample Year"),var_trans = 0,output_trans = 1)
    propResult <- propData(transformed_df,part_col = c("Subject Number"), val_col = c("titer"), group_col = c("Sample Year"),var_trans = 0, threshold_trans = 1,threshold = 2) #the proportion of titers above 2
    prot_propResult <- prot_propData(transformed_df,part_col = c("Subject Number"),val_col = c("titer"),group_col = c("Sample Year"), min.y = 0, weight = c(0,0.3,0.5,0.8,0.9))
    skewnessResult <- skewnessData(transformed_df,part_col = c("Subject Number"),weight_col = c("isolation_year"), val_col = c("titer"), group_col = c("Sample Year"),var_trans = 0)
    widthResult2 <- widthData(transformed_df,part_col = c("Subject Number"),weight_col = c("isolation_year"), val_col = c("titer"), group_col = c("Sample Year"),var_trans = 0,threshold_trans = 1,threshold = 2,mode="mean",width_col = "Width2") #the width of titers above 2
    widthResult3 <- widthData(transformed_df,part_col = c("Subject Number"),weight_col = c("isolation_year"), val_col = c("titer"), group_col = c("Sample Year"),var_trans = 0,threshold_trans = 1,threshold = 3,mode="mean",width_col = "Width3") #the width of titers above 3

    overallResult <- atyResult %>%
      left_join(aucResult, by = c("Subject Number","Sample Year")) %>%
      left_join(giniResult, by = c("Subject Number","Sample Year")) %>%
      left_join(gmtResult, by = c("Subject Number","Sample Year")) %>%
      left_join(kurtosisResult, by = c("Subject Number","Sample Year")) %>%
      left_join(max_titerResult, by = c("Subject Number","Sample Year")) %>%
      left_join(propResult, by = c("Subject Number","Sample Year")) %>%
      left_join(prot_propResult, by = c("Subject Number","Sample Year")) %>%
      left_join(skewnessResult, by = c("Subject Number","Sample Year")) %>%
      left_join(widthResult2, by = c("Subject Number","Sample Year"))%>%
      left_join(widthResult3, by = c("Subject Number","Sample Year"))

    head(overallResult,10)
```
```r
    Output:
    ##    Subject Number Sample Year      aty      auc gini_coefficient       GMT Kurtosis Max_titer Proportion Prot_prop   Skewness
    ## 1               1        2007 1992.256 77.07798        0.2250316 0.9649123 5.822214         5  0.3333333       0.9 -0.9739634
    ## 2               1        2008 1999.937 47.78503        0.2461859 1.6140351 1.646947         5  0.4385965       0.9  0.2071604
    ## 3               1        2009 1999.986 46.80581        0.2415966 1.5789474 1.609897         5  0.4210526       0.9  0.1802292
    ## 4               1        2010 1992.265 74.21888        0.2403342 0.8245614 5.732402         4  0.3508772       0.9 -1.0324670
    ## 5               1        2011 1992.826 81.71422        0.2407295 1.1228070 5.458734         5  0.4210526       0.9 -0.9726327
    ## 6               1        2012 1992.221 73.83704        0.1727029 0.7368421 5.197581         4  0.2456140       0.9 -0.9195227
    ## 7               2        2007 1992.647 76.59212        0.2097812 0.8947368 5.898057         4  0.2631579       0.9 -1.1269794
    ## 8               2        2008 2002.108 66.79414        0.2664053 2.3157895 1.892472         6  0.5964912       0.9 -0.3756182
    ## 9               2        2009 2001.715 67.87382        0.2621916 2.3684211 1.812409         6  0.6491228       0.9 -0.2995806
    ## 10              2        2010 1994.853 88.58177        0.2060224 1.4736842 6.350863         5  0.3859649       0.9 -1.4286519
    ##       Width2        Width3
    ## 1  0.2294611  1.032853e-01
    ## 2  0.1387480 -5.287760e-15
    ## 3  0.1293311 -5.287760e-15
    ## 4  0.1650965  4.217115e-02
    ## 5  0.2339210  9.654034e-02
    ## 6  0.2107012  3.102572e-03
    ## 7  0.1870788  6.899400e-02
    ## 8  0.3015671  9.849466e-02
    ## 9  0.2786792  1.152329e-01
    ## 10 0.2842462  1.059375e-01
```

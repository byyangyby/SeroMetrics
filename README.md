# SeroMetrics

`SeroMetrics` is an R package designed to quantify an individual’s
*antigenic landscape*, which comprises hemagglutination inhibition (HAI)
titers against a panel of influenza viruses isolated at various times or
in different antigenic spaces. Visualization functions are currently
being planned and developed, and will be released soon.

Summary metrics and application examples are adapted from the following
publications:

-   Yang B, Lessler J, Zhu H, Jiang CQ, Read JM, Hay JA, Kwok KO, Shen
    R, Guan Y, Riley S, Cummings DA. [Life course exposures continually
    shape antibody profiles and risk of seroconversion to
    influenza](https://journals.plos.org/plospathogens/article?id=10.1371/journal.ppat.1008635).
    PLoS pathogens. 2020 Jul 23;16(7):e1008635.

-   Yang B, Leung NHL, Fox A, Sullivan SG, Ho F, Barr IG, Cheng SMS,
    Wong SS, Levine MZ, Iuliano AD, Peiris M, Thompson MG, Cowling BJ.
    Evaluating Antibody Breadth Following Various Enhanced Influenza
    Vaccines in Older Adults. *In preparation*.

------------------------------------------------------------------------

## 1. Installation

------------------------------------------------------------------------

## 2. Data structure

<table>
<colgroup>
<col style="width: 18%" />
<col style="width: 18%" />
<col style="width: 63%" />
</colgroup>
<thead>
<tr class="header">
<th style="text-align: left;">Variable</th>
<th style="text-align: left;">Type</th>
<th style="text-align: left;">Explanation</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td style="text-align: left;"><em>id</em></td>
<td style="text-align: left;">mandatory</td>
<td style="text-align: left;">This should be the unique id of the antibody profile, e.g.,
each sampled participant or participant-sample-year.</td>
</tr>
<tr class="even">
<td style="text-align: left;"><em>value</em></td>
<td style="text-align: left;">mandatory</td>
<td style="text-align: left;">This should be the investigated value,
usually being titers.</td>
</tr>
<tr class="odd">
<td style="text-align: left;"><em>distance</em></td>
<td style="text-align: left;">mandatory</td>
<td style="text-align: left;">This should be the absolute distance among
virus strains, usually being isolation year of the or amino acid (AA)
distance between the tested and referenced virus.</td>
</tr>
<tr class="even">
<td style="text-align: left;"><em>strain</em></td>
<td style="text-align: left;">optional</td>
<td style="text-align: left;">This should be the virus strains. When
distance is not provided, one may need to derive isolation year from the
information in virus strain.</td>
</tr>
<tr class="odd">
<td style="text-align: left;"><em>group</em></td>
<td style="text-align: left;">optional</td>
<td style="text-align: left;">This should be the group of some sample of
a participant. Sometimes, one participant may take tests on the same
strains multiple times, and the samples collected from each time forms a
sample group.</td>
</tr>
<tr class="even">
<td style="text-align: left;"><em>other variables</em></td>
<td style="text-align: left;">optional</td>
<td style="text-align: left;">There might be age, sex, and other
information regarding the participant. At most times, these variables
will not affect the final result of the metrics.</td>
</tr>
</tbody>
</table>

### 2.1. Simulated data

The following shows a simulated antibody landscape from the following
participant:

-   `id` in the study is 25

-   `birth_year` born in the year of 1952

-   `sample_year` antiserum was collected in the year of 2020

-   `isolation_year` the isolation year of the tested virus; there are
    40 rows in the sample profile, meaning that we tested the antiserum
    against 40 tested viruses

-   `log_titer` log-transformed antibody titer

<!-- -->

    head(profile, 5)

    ## # A tibble: 5 × 5
    ##      id birth_year sample_year isolation_year log_titer
    ##   <int>      <int>       <dbl>          <int>     <int>
    ## 1    25       1985        2020           1980         4
    ## 2    25       1985        2020           1952         6
    ## 3    25       1985        2020           1974         4
    ## 4    25       1985        2020           1975         5
    ## 5    25       1985        2020           1992         5

### 2.2. Transform distance variable

The following shows examples to compute the `distance` variable using
different references, which will be further used for metric calculation.
In this example, the antibody landscape is arranged by isolation time of
the tested virus.

If the antibody landscape is arranged in the antigenic space, the
`distance` variable should be computed as the difference in AA
substitutions between tested virus and the reference virus.
distance\_null will be unable to calculate in such situation.

    profile = profile %>%
      mutate(
        distance_1 = isolation_year - birth_year, # align to the year of birth
        distance_2 = isolation_year - sample_year, # align to the antiserum collection year
        distance_null = isolation_year # without reference; raw value for isolation year
      )

    profile %>%
      select(-id) %>%
      head()

    ## # A tibble: 6 × 7
    ##   birth_year sample_year isolation_year log_titer distance_1 distance_2 distance_null
    ##        <int>       <dbl>          <int>     <int>      <int>      <dbl>         <int>
    ## 1       1985        2020           1980         4         -5        -40          1980
    ## 2       1985        2020           1952         6        -33        -68          1952
    ## 3       1985        2020           1974         4        -11        -46          1974
    ## 4       1985        2020           1975         5        -10        -45          1975
    ## 5       1985        2020           1992         5          7        -28          1992
    ## 6       1985        2020           2001         7         16        -19          2001

------------------------------------------------------------------------

## 3. Example of applications

### 3.1. Lessler *et al*. 2012 *PLoS Pathogens*

This is an example using the data set published by [Lessler *et al.*,
2012 *PLoS
Pathogens*](https://journals.plos.org/plospathogens/article?id=10.1371/journal.ppat.1002802).

In this data set, cross-sectional antisera were collected, indicating
that for each participant, only one titer was available for each tested
virus. In this situation, the `group_col` argument can be set to *NULL*,
and since all virus strains have different isolation years, the
`mode` argument can be set to *NULL*.

        #show the basic structure of the dataset
        head(lessler,10)

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

#### 3.1.1. Calculate multiple metrics at once

Using `overallData()` function to simultaneously calculate multiple
metrics specified in the `required_metrics` argument. The output is a
matrix with each individual’s antibody profile, in this case for each
participant.

    overallResult <- overallData(data = lessler, 
                                 part_col = "id", # unique id of the participant
                                 weight_col = "isolation_year", # specify the column name for weights.
                                 val_col = "titers", # specify the column name for the titer value
                                 var_trans = 1, 
                                 gmt_output_trans = 1, 
                                 required_metrics = c("AUC",
                                                      "ATY",
                                                      "GMT",
                                                      "Prot_prop",
                                                      "Width0",
                                                      "Width2",
                                                      "Width3"), # metrics to be computed
                                 col_names = c("auc",
                                               "aty",
                                               "gmt",
                                               "prot_prop",
                                               "Width10",
                                               "Width40",
                                               "Width80") # column names for the calculated metrics in the output matrix

                                 )

    head(overallResult, 10)

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

#### 3.1.2. Calculate individual metric

You can also calculate a single metric using the corresponding function
in the package. This generates a two-column matrix with a unique ID and
values for the chosen metric, which can be manually integrated into an
overall matrix.

        # average titer year
        atyResult <- atyData(
          data = lessler,
          part_col = c("id"),
          weight_col = c("isolation_year"), 
          val_col = c("titers"), 
          var_trans = 1
        )
        
        # area under the curve
        aucResult <- aucData(
          data = lessler,
          part_col = c("id"),
          weight_col = c("isolation_year"), 
          val_col = c("titers"), 
          var_trans = 1
        )

        # Gini coefficient   
        giniResult <- giniData(
          data = lessler,
          part_col = c("id"), 
          val_col = c("titers"), 
          var_trans = 1
        )
        
         # Geometric mean titer   
        gmtResult <- gmtData(
          data = lessler,
          part_col = c("id"), 
          val_col = c("titers"), 
          var_trans = 1,
          output_trans = 1
        )
        
        # Kurtosis
        kurtosisResult <- kurtosisData(
          data = lessler,
          part_col = c("id"), 
          weight_col = c("isolation_year"), 
          val_col = c("titers"), 
          var_trans = 1
        )
        
        # Maximum titer
        max_titerResult <- max_titerData(
          data = lessler,
          part_col = c("id"), 
          val_col = c("titers"), 
          var_trans = 1,
          output_trans = 1
        )
        
        # The proportion of titers above 40
        propResult2 <- propData(
          data = lessler,
          part_col = c("id"), 
          val_col = c("titers"), 
          var_trans = 1,
          threshold_trans = 1,
          threshold = 2, 
          prop_col = "Prop2"
        )
        
        # The proportion of titers above 80
        propResult3 <-  propData(
          data = lessler,
          part_col = c("id"), 
          val_col = c("titers"), 
          var_trans = 1,
          threshold_trans = 1,
          threshold = 3, 
          prop_col = "Prop3"
        )
        
        # Using default weight vector
        prot_propResult <- prot_propData(
          data = lessler,
          part_col = c("id"), 
          val_col = c("titers")
        )
        
        # Skewness
        skewnessResult <- skewnessData(
          data = lessler,
          part_col = c("id"), 
          weight_col = c("isolation_year"), 
          val_col = c("titers"), 
          var_trans = 1
        )
        
         # The width of titers above 40
        widthResult <- widthData(
          data = lessler,
          part_col = c("id"), 
          weight_col = c("isolation_year"), 
          val_col = c("titers"), 
          var_trans = 1, 
          threshold_trans = 1,
          threshold = 2, 
          width_col = "Width2"
        )

        # Combine calculated metrics 
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

    ##    id      ATY   AUC gini_coefficient      GMT Kurtosis Max_titer     Prop2     Prop3  Prot_prop
    ## 1   1 1988.754  63.0        0.1640625 3.555556 1.687270         5 1.0000000 0.8888889 0.11111111
    ## 2   2 1985.372  93.5        0.2426471 2.777778 1.601629         5 0.7777778 0.4444444 0.19444444
    ## 3   3 1987.554  55.5        0.1718750 2.555556 1.545701         4 0.8888889 0.3333333 0.12777778
    ## 4   4 1980.289  19.0        0.1477273 1.444444 2.134154         3 0.3333333 0.1111111 0.02222222
    ## 5   5 1984.481  80.5        0.1985294 2.777778 1.668338         5 0.8888889 0.4444444 0.16666667
    ## 6   6 1990.405 115.5        0.2073171 3.555556 2.516994         5 0.8888889 0.5555556 0.35000000
    ## 7   7 1983.577  26.0        0.2200000 1.777778 1.191072         5 0.4444444 0.1111111 0.08333333
    ## 8   8 1988.106  85.0        0.1184211 4.222222 1.512408         5 1.0000000 0.8888889 0.17777778
    ## 9   9 1996.665  66.5        0.2941176 2.777778 2.364882         6 0.6666667 0.3333333 0.23333333
    ## 10 10 1989.727 125.5        0.2150000 4.555556 1.815939         7 0.8888889 0.7777778 0.49444444
    ##      Skewness    Width2
    ## 1  -0.3850070 0.0000000
    ## 2   0.3521357 0.2687500
    ## 3  -0.4852686 0.0250000
    ## 4   0.9432719 0.0000000
    ## 5   0.2886074 0.0875000
    ## 6  -0.7167264 0.6416667
    ## 7   0.3490311 0.0437500
    ## 8  -0.3804884 0.1750000
    ## 9  -0.5062919 0.2916667
    ## 10 -0.7370504 0.5345833

------------------------------------------------------------------------

### 3.2. Fonville *et al*. 2014 *Science*

This is an example using the data set published by [*Fonville* et al.,
2014 *Science*](https://www.science.org/doi/10.1126/science.1256427).

In this dataset, longitudinal antisera were collected from the same
population over six years. For each participant, multiple titer values
are available for each tested virus, corresponding to antiserum
collected during each visit.

In this case, the `group_col` is the sample year, and different viruses
share the same `isolation_year`, making the `mode` necessary.

Additionally, the raw data has a specific structure where each tested
virus forms one column. Therefore, data transformation (from a wide
table to a long table) should be done first.

        # data transformation
        transformed_df <- transform(fonville)

    ## Warning: There was 1 warning in `mutate()`.
    ## ℹ In argument: `isolation_year = case_when(...)`.
    ## Caused by warning:
    ## ! NAs introduced by coercion

        # show the basic structure of the dataset
        head(transformed_df,10)

    ## # A tibble: 10 × 9
    ##    `Subject Number` `Sample Year` `Year of Birth` Sample  `PCR Results` Row in Fonville Fig S1…¹ strain
    ##               <dbl>         <dbl>           <dbl> <chr>   <chr>                            <dbl> <chr> 
    ##  1                1          2007            1957 H3 PCR+ <NA>                                63 BI/16…
    ##  2                1          2007            1957 H3 PCR+ <NA>                                63 BI/21…
    ##  3                1          2007            1957 H3 PCR+ <NA>                                63 BI/17…
    ##  4                1          2007            1957 H3 PCR+ <NA>                                63 BI/22…
    ##  5                1          2007            1957 H3 PCR+ <NA>                                63 NL/23…
    ##  6                1          2007            1957 H3 PCR+ <NA>                                63 NL/62…
    ##  7                1          2007            1957 H3 PCR+ <NA>                                63 NL/82…
    ##  8                1          2007            1957 H3 PCR+ <NA>                                63 NL/17…
    ##  9                1          2007            1957 H3 PCR+ <NA>                                63 JO/33…
    ## 10                1          2007            1957 H3 PCR+ <NA>                                63 SD/9/…
    ## # ℹ abbreviated name: ¹​`Row in Fonville Fig S15`
    ## # ℹ 2 more variables: titer <dbl>, isolation_year <dbl>

#### 3.2.1. Calculate multiple metrics at once

Using `overallData()` function to simultaneously calculate multiple
metrics specified in the `required_metrics` argument. The output is a
matrix with each individual’s antibody profile, in this case for each
participant.

     overallResultFonville <- overallData(
       transformed_df, 
       part_col = "Subject Number", 
       weight_col = "isolation_year", 
       val_col = "titer", 
       group_col = "Sample Year", 
       var_trans = 0, 
       max_titer_output_trans = 1, 
       mode = "mean", 
       required_metrics = c("gini_coefficient","Kurtosis","Skewness","Prop0","Prop2","Prop3"),
       col_names = c("Gini","Kurt","Skew","prop0","prop2","prop3")
     )
     
        head(overallResultFonville, 10)

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

#### 3.2.2. Calculate individual metric

Similarly, you can also calculate a single metric using the
corresponding function in the package.

        atyResult <- atyData(
          transformed_df,
          part_col = c("Subject Number"),
          weight_col = c("isolation_year"), 
          val_col = c("titer"), 
          group_col = c("Sample Year"),
          var_trans = 0, 
          aty_col = "aty"
        )

        aucResult <- aucData(
          transformed_df,
          part_col = c("Subject Number"),
          weight_col = c("isolation_year"), 
          val_col = c("titer"), 
          group_col = c("Sample Year"),
          var_trans = 0,
          mode="mean", 
          auc_col = "auc"
        )
        
        giniResult <- giniData(
          transformed_df,
          part_col = c("Subject Number"), 
          val_col = c("titer"), 
          group_col = c("Sample Year"),
          var_trans = 0
        )
        
        gmtResult <- gmtData(
          transformed_df,
          part_col = c("Subject Number"), 
          val_col = c("titer"), 
          group_col = c("Sample Year"),
          var_trans = 0,
          output_trans = 1
        )
        
        kurtosisResult <- kurtosisData(
          transformed_df,
          part_col = c("Subject Number"),
          weight_col = c("isolation_year"), 
          val_col = c("titer"), 
          group_col = c("Sample Year"),
          var_trans = 0
        )
        
        max_titerResult <-max_titerData(
          transformed_df,
          part_col = c("Subject Number"), 
          val_col = c("titer"),
          group_col = c("Sample Year"),
          var_trans = 0,
          output_trans = 1
        )
        
        #the proportion of titers above 40
        propResult <- propData(
          transformed_df,
          part_col = c("Subject Number"),
          val_col = c("titer"), 
          group_col = c("Sample Year"),
          var_trans = 0,
          threshold_trans = 1,
          threshold = 2
        ) 
        
        prot_propResult <- prot_propData(
          transformed_df,
          part_col = c("Subject Number"),
          val_col = c("titer"),
          group_col = c("Sample Year"), 
          min.y = 0, 
          weight = c(0,0.3,0.5,0.8,0.9)
        )
        
        skewnessResult <- skewnessData(
          transformed_df,
          part_col = c("Subject Number"),
          weight_col = c("isolation_year"), 
          val_col = c("titer"), 
          group_col = c("Sample Year"),
          var_trans = 0
        )
        
        #the width of titers above 40
        widthResult2 <- widthData(
          transformed_df,
          part_col = c("Subject Number"),
          weight_col = c("isolation_year"), 
          val_col = c("titer"), 
          group_col = c("Sample Year"),
          var_trans = 0,
          threshold_trans = 1,
          threshold = 2,
          mode="mean",
          width_col = "Width2"
        ) 
        
        #the width of titers above 80
        widthResult3 <- widthData(
          transformed_df,
          part_col = c("Subject Number"),
          weight_col = c("isolation_year"), 
          val_col = c("titer"), 
          group_col = c("Sample Year"),
          var_trans = 0,
          threshold_trans = 1,
          threshold = 3,
          mode="mean",
          width_col = "Width3"
        )
        
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

    ##    Subject Number Sample Year      aty      auc gini_coefficient       GMT Kurtosis Max_titer
    ## 1               1        2007 1992.256 77.07798        0.2250316 0.9649123 5.822214         5
    ## 2               1        2008 1999.937 47.78503        0.2461859 1.6140351 1.646947         5
    ## 3               1        2009 1999.986 46.80581        0.2415966 1.5789474 1.609897         5
    ## 4               1        2010 1992.265 74.21888        0.2403342 0.8245614 5.732402         4
    ## 5               1        2011 1992.826 81.71422        0.2407295 1.1228070 5.458734         5
    ## 6               1        2012 1992.221 73.83704        0.1727029 0.7368421 5.197581         4
    ## 7               2        2007 1992.647 76.59212        0.2097812 0.8947368 5.898057         4
    ## 8               2        2008 2002.108 66.79414        0.2664053 2.3157895 1.892472         6
    ## 9               2        2009 2001.715 67.87382        0.2621916 2.3684211 1.812409         6
    ## 10              2        2010 1994.853 88.58177        0.2060224 1.4736842 6.350863         5
    ##    Proportion Prot_prop   Skewness    Width2        Width3
    ## 1   0.3333333       0.9 -0.9739634 0.2294611  1.032853e-01
    ## 2   0.4385965       0.9  0.2071604 0.1387480 -5.287760e-15
    ## 3   0.4210526       0.9  0.1802292 0.1293311 -5.287760e-15
    ## 4   0.3508772       0.9 -1.0324670 0.1650965  4.217115e-02
    ## 5   0.4210526       0.9 -0.9726327 0.2339210  9.654034e-02
    ## 6   0.2456140       0.9 -0.9195227 0.2107012  3.102572e-03
    ## 7   0.2631579       0.9 -1.1269794 0.1870788  6.899400e-02
    ## 8   0.5964912       0.9 -0.3756182 0.3015671  9.849466e-02
    ## 9   0.6491228       0.9 -0.2995806 0.2786792  1.152329e-01
    ## 10  0.3859649       0.9 -1.4286519 0.2842462  1.059375e-01

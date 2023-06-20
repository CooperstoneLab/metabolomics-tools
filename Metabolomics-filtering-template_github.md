Filtering for Untargeted Metabolomics Data
================
The Cooperstone lab
Lots of times

- [Load libraries](#load-libraries)
- [Read in Data](#read-in-data)
- [RT filter if necessary](#rt-filter-if-necessary)
- [Cleaning up data](#cleaning-up-data)
  - [Create mz_rt](#create-mz_rt)
  - [Clean up file names](#clean-up-file-names)
  - [Remove unwanted samples](#remove-unwanted-samples)
- [Start filtering](#start-filtering)
  - [CV function](#cv-function)
  - [Counting QCs](#counting-qcs)
  - [Filter on QC CV](#filter-on-qc-cv)
  - [Merge back the rest of the data](#merge-back-the-rest-of-the-data)
  - [Process blanks](#process-blanks)
  - [Save your file](#save-your-file)

Below is a series of code put together by the Cooperstone lab (Emma
Bilbrey, Jenna Miller, JL Hartman, Daniel Quiroz Moreno, and Jessica
Cooperstone) to be used as a template for filtering untargeted
metabolomics data after deconvolution with software like MZmine. This
approach should work for many datasets, though may need to be slightly
modified depending on what you are trying to do.

## Load libraries

``` r
library(tidyverse)
library(janitor) # if you want to clean_names()
```

Once you get deconvoluted data from MZmine or similar programs, you need
to wrangle your data in such a way that you can conduct your analysis on
it.

# Read in Data

First we want to read in our raw data. The code here is to read in data
directly from MZmine, so if you are using a different program for
deconvolution, you might need to make some light adjustments. The sample
data here was collected by JL Hartman.

``` r
metabdata <- read_csv(file = "./data/20210831_BNB_neg_fullscan_final_4528feat.csv",
                      col_names = TRUE) # has headers
```

    ## Rows: 4528 Columns: 120
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: ","
    ## dbl (120): row m/z, row retention time, BNBpilot_neg_1002_1_013.mzML Peak he...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

``` r
# replaces zeroes with NAs
metabdata <- replace(metabdata, metabdata == 0, NA)

head(metabdata)
```

    ## # A tibble: 6 × 120
    ##   `row m/z` `row retention time` BNBpilot_neg_1002_1_01…¹ BNBpilot_neg_1005_1_…²
    ##       <dbl>                <dbl>                    <dbl>                  <dbl>
    ## 1      215.                0.682                 7679077                7447086.
    ## 2      132.                0.626                 6205564.               4927318.
    ## 3      146.                0.638                 4960052.               2120450.
    ## 4      223.                0.681                 4479680                2886648.
    ## 5      277.                0.649                 4078063                4136307 
    ## 6      191.                0.711                 2988186.               3989021.
    ## # ℹ abbreviated names: ¹​`BNBpilot_neg_1002_1_013.mzML Peak height`,
    ## #   ²​`BNBpilot_neg_1005_1_016.mzML Peak height`
    ## # ℹ 116 more variables: `BNBpilot_neg_1028_7_042.mzML Peak height` <dbl>,
    ## #   `BNBpilot_neg_1006_1_017.mzML Peak height` <dbl>,
    ## #   `BNBpilot_neg_1011_3_023.mzML Peak height` <dbl>,
    ## #   `BNBpilot_neg_1012_3_024.mzML Peak height` <dbl>,
    ## #   `BNBpilot_neg_1008_3_020.mzML Peak height` <dbl>, …

Note there is no metadata included in this file. Just m/z, retention
time, and a column for each sample, where values are peak heights. We
are using peak height instead of peak area because it is less dependent
on bad peak shape which you get sometimes with metabolomics.

# RT filter if necessary

You might have deconvoluted more data than you plan to use in your
analysis. For example, you may want to exclude the first bit and last
bit of your run, since you do not expect to have good reproducibility in
those areas.

Here, we are filtering to only include features that elute between
0.5-7.5 min of this 10 min run.

``` r
metabdata_RTfilt <- metabdata %>%
  filter(between(`row retention time`, 0.5, 7.5))

# did it work?
range(metabdata_RTfilt$`row retention time`)
```

    ## [1] 0.5086431 7.4997881

# Cleaning up data

## Create mz_rt

This creates a unique identifier for each feature using its
mass-to-charge ratio (m/z) and retention time (RT).

``` r
MZ_RT <- metabdata_RTfilt %>%
  mutate(mz = round(metabdata_RTfilt$`row m/z`, digits = 4), # Decrease number of decimals for m/z & rt
         rt = round(metabdata_RTfilt$`row retention time`, digits = 3),
         .before=1,
         .keep="unused") %>%
  unite(mz_rt, c(mz, rt), remove=TRUE) #%>% # Combine m/z & rt with _ in between

# remove old RT and MZ columns
MZ_RT <- MZ_RT %>%
  select(-`row m/z`, -`row retention time`) 
```

## Clean up file names

We are using `gsub()` to replace strings (i.e. characters) in our sample
names. Here is some useful info about
[gsub](http://www.endmemo.com/r/gsub.php), and two different tutorials
[here](https://www.youtube.com/watch?v=4ZokHoF99DY) and
[here](https://www.youtube.com/watch?v=r4Sh7H6wzPA). You will likely
need to change this code to suit your purposes.

``` r
# remove stuff from the end of file names, ".mzML Peak height"
newcolumnnames <- gsub(".mzML.*","", colnames(MZ_RT))
colnames(MZ_RT) <- newcolumnnames

# remove stuff from the front of file names "BNB_pilot_neg"
newcolumnnames <- gsub("^BNBpilot_neg_","", colnames(MZ_RT))
colnames(MZ_RT) <- newcolumnnames
```

## Remove unwanted samples

There is a solvent blank sample present in this dataset, which we will
not be using, so we will remove it now.

``` r
MZ_RT <- MZ_RT %>%
  select(-SB3againagain_0_3)
```

# Start filtering

### CV function

Since base R does not have a function to calculate coefficient of
variance, let’s write one.

``` r
cv <- function(x){
        (sd(x)/mean(x))
}
```

## Counting QCs

Subset QCs and filter features to keep only those that are present in
100% of QCs. You could change this parameter based on your data.

``` r
# check dimensions of current df
dim(MZ_RT)
```

    ## [1] 4525  118

``` r
MZ_RT_QCs <- MZ_RT %>%
  select(mz_rt, contains("QC")) %>% # select QCs
  filter(rowSums(is.na(.)) <= 1) # remove rows that have 1 or more NAs
```

``` r
# check dimensions of QCs filtered df
dim(MZ_RT_QCs)
```

    ## [1] 4319   18

``` r
# how many features got removed with this filtering?
nrow(MZ_RT) - nrow(MZ_RT_QCs)
```

    ## [1] 206

## Filter on QC CV

Here we are removing features that have a CV of more than 30% in the
QCs. The rationale is that if a feature cannot be reproducibly measured
in samples that are all the same, it should not be included in our
analysis.

``` r
# calculate CV row-wise (1 means row-wise)
QC_CV <- apply(MZ_RT_QCs[, 2:ncol(MZ_RT_QCs)], 1, cv)

# bind the CV vector back to the QC df
MZ_RT_QCs_CV <- cbind(MZ_RT_QCs, QC_CV)

# filter for keeping features with QC_CV <= 0.30 (or 30%)
MZ_RT_QCs_CVfilt <- MZ_RT_QCs_CV %>%
  filter(QC_CV <= 0.30)
```

How many features did I remove with this CV filtering?

``` r
nrow(MZ_RT_QCs) - nrow(MZ_RT_QCs_CVfilt)
```

    ## [1] 308

## Merge back the rest of the data

MZ_RT_QCs_CVfilt only contains the QCs, We want to keep only the rows
that are present in this df, and then merge back all of the other
samples present in MZ_RT. We will do this by creating a vector that has
the mz_rt features we want to keep, and then using `filter()` and `%in%`
to keep only features that are a part of this list.

``` r
dim(MZ_RT_QCs_CVfilt)
```

    ## [1] 4011   19

``` r
dim(MZ_RT)
```

    ## [1] 4525  118

``` r
# make a character vector of the mz_rt features we want to keep
# i.e., the ones that passed our previous filtering steps
features_to_keep <- as.character(MZ_RT_QCs_CVfilt$mz_rt)

MZ_RT_filt <- MZ_RT %>%
  filter(mz_rt %in% features_to_keep)

dim(MZ_RT_filt)
```

    ## [1] 4011  118

You should have the same number of features in MZ_RT_QCs_CVfilt as you
do in your new filtered df MZ_RT_filt.

``` r
all.equal(nrow(MZ_RT_QCs_CVfilt), nrow(MZ_RT_filt))
```

    ## [1] TRUE

## Process blanks

We want to remove features that are present in our process blanks as
they are not coming from compounds present in our samples. In this
dataset, the sample (there is only one, typically you would have at
least 3 process blanks to average) representing this process blank (a
sample that includes all the extraction materials, minus the sample,
here the tomato was replaced by mass with water) has “PB” in the sample
name.

``` r
# grab the name of the column/sample that is the process blank
grep("PB", colnames(MZ_RT_filt), value = TRUE)
```

    ## [1] "PB1again_0_6"

Calculate the average value across the QCs, then remove features that
are not at least 10x higher in the QCs than in the process blank. To do
this we will use
[`apply()`](https://www.rdocumentation.org/packages/base/versions/3.6.2/topics/apply).

`apply(X, MARGIN, FUN,...)` where X is your df, MARGIN is 1 for
row-wise, and 2 for col-wise, and FUN is your function

``` r
# pull avg peak height across QCs
avg_height_QC <- apply(MZ_RT_QCs_CVfilt[, 2:ncol(MZ_RT_QCs_CVfilt)], 1, mean)

# bind back to rest of data
MZ_RT_filt_QC_avg <- cbind(MZ_RT_filt, avg_height_QC)

# check dimensions
dim(MZ_RT_filt_QC_avg)
```

    ## [1] 4011  119

Pull the name of your process blank, and make a new column that
indicates how many fold higher your peak height is in your average QC vs
your process blank.

``` r
# pull name of process blank 
grep("PB", colnames(MZ_RT_filt), value = TRUE)
```

    ## [1] "PB1again_0_6"

``` r
# make a new column that has a value of how many fold higher peak height is
# in QCs as compared to PB
# here there is only  one PB, but would be better to have > 1 but this is ok
# then you can avg your PBs together and do the same thing
MZ_RT_filt_PB <- MZ_RT_filt_QC_avg %>% 
  mutate(fold_higher_in_QC = avg_height_QC/PB1again_0_6)
```

We want to keep features that are at least 10x higher in QCs than
process blanks, and we also want to keep NAs, because an NA indicates
that a feature absent in the process blanks (i.e., you get an NA because
you’re trying to divide by zero).

``` r
# keep features that are present at least 10x higher in QCs vs PB
# or, keep NAs because those are absent in blank
MZ_RT_filt_PBremoved <- MZ_RT_filt_PB %>%
  filter(fold_higher_in_QC > 10 | is.na(fold_higher_in_QC)) %>%
  select(-fold_higher_in_QC, -avg_height_QC, -PB1again_0_6)
  
dim(MZ_RT_filt_PBremoved)
```

    ## [1] 3948  117

## Save your file

Now you have a list of features present in your samples after filtering
for CV in QCs, and removing all the extraneous columns we added to help
us do this, along with removing any process blanks.

``` r
write_csv(MZ_RT_filt_PBremoved,
          "data/BNB_neg_filtered_template.csv")
```
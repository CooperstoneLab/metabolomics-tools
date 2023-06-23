Filtering for Untargeted Metabolomics Data
================
The Cooperstone lab

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
library(janitor) # if you want to clean_names() or get_dupes()
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
         rt = round(metabdata_RTfilt$`row retention time`, digits = 4),
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

Do we have any duplicate features?

``` r
get_dupes(MZ_RT_filt_PBremoved, mz_rt)
```

    ##             mz_rt dupe_count 1002_1_013 1005_1_016 1028_7_042 1006_1_017
    ## 1 239.0776_0.6677          2   18398.26  40836.883  25596.115   34077.79
    ## 2 239.0776_0.6677          2   18398.26  40836.883  25596.115   34077.79
    ## 3 270.9855_0.7435          2   31931.09  35832.715  31064.273   27352.09
    ## 4 270.9855_0.7435          2   31931.09  35832.715  31064.273   27352.09
    ## 5  677.1727_4.303          2    1029.13   4552.599   9868.488   18662.72
    ## 6  677.1727_4.303          2    1029.13   4552.599   9868.488   18662.72
    ##   1011_3_023 1012_3_024 1008_3_020 1003_1_014 1004_1_015 1010_3_022 1026_7_040
    ## 1   24453.42   23670.80   30323.73   60174.62   30499.17   34480.30   26423.76
    ## 2   24453.42   23670.80   30323.73   60174.62   30499.17   34480.30   26423.76
    ## 3   26064.35   25829.42   53138.13   44746.21   31656.12   50254.85   27008.75
    ## 4   26064.35   25829.42   53138.13   44746.21   31656.12   50254.85   27008.75
    ## 5   19399.51   11100.51   11389.76   11889.22   32046.98   14976.43   10299.67
    ## 6   19399.51   11100.51   11389.76   11889.22   32046.98   14976.43   10299.67
    ##   1031_9_046 1025_7_039 1007_1_018 1020_5_033 1016_5_029 1014_3_026 1019_5_032
    ## 1  26000.303   26258.80   32385.83   27665.63  27169.596   30924.73  32761.664
    ## 2  26000.303   26258.80   32385.83   27665.63  27169.596   30924.73  32761.664
    ## 3  28121.031   33599.97   35116.06   28601.47  26155.508   28363.31  37751.633
    ## 4  28121.031   33599.97   35116.06   28601.47  26155.508   28363.31  37751.633
    ## 5   3429.637   11847.22   36846.87   17553.60   8155.577   16408.99   6993.417
    ## 6   3429.637   11847.22   36846.87   17553.60   8155.577   16408.99   6993.417
    ##   1013_3_025 1029_9_044 1023_7_037 1018_5_031 1022_7_036 1021_5_034 1032_9_047
    ## 1 40810.3750  16785.041   22237.80   25636.65   31415.17   30469.16  22598.445
    ## 2 40810.3750  16785.041   22237.80   25636.65   31415.17   30469.16  22598.445
    ## 3 41156.9062  33820.180   40601.19   25459.74   35669.09   28086.34  27039.873
    ## 4 41156.9062  33820.180   40601.19   25459.74   35669.09   28086.34  27039.873
    ## 5   674.3444   9710.306   13211.57   15605.22   18984.73   30304.05   6306.936
    ## 6   674.3444   9710.306   13211.57   15605.22   18984.73   36257.68   6306.936
    ##   1009_3_021 1039_11_054 1035_9_049 1046_13_062 1030_9_045 1040_11_055
    ## 1  41088.746   25952.139  28047.129    22712.89  31025.541    21571.04
    ## 2  41088.746   25952.139  28047.129    22712.89  31025.541    21571.04
    ## 3   4337.607   31361.414  39152.219    33796.79  25816.801    27523.80
    ## 4   4337.607   31361.414  39152.219    33796.79  25816.801    27523.80
    ## 5  16678.168    3359.229   6854.434    11038.56   7249.560    10789.50
    ## 6  16678.168    3359.229   6854.434    11038.56   7814.693    10789.50
    ##   1033_31_132 1038_11_053 1041_11_056 1036_9_050 1034_9_048 1052_15_069
    ## 1    50582.26   13343.054   18304.398   20093.10   31052.56   30123.566
    ## 2    50582.26   13343.054   18304.398   20093.10   31052.56   30123.566
    ## 3    43392.93   23035.375   24987.023   31239.41   28920.15   31446.117
    ## 4    43392.93   23035.375   24987.023   31239.41   28920.15   31446.117
    ## 5    33364.76    9680.562    5536.901   11632.77   14842.04    6412.228
    ## 6    33364.76    9680.562    5536.901   11632.77   14842.04    6860.068
    ##   1050_13_066 1048_13_064 1037_11_052 1053_15_070 1054_15_071 1049_13_065
    ## 1    28652.63   71729.469    22287.97   31470.201   25497.828    29934.35
    ## 2    28652.63   71729.469    22287.97   31470.201   25497.828    29934.35
    ## 3    31721.76   35101.863    28286.07   40586.426   38719.844    36615.30
    ## 4    31721.76   35101.863    28286.07   40586.426   38719.844    36615.30
    ## 5    12644.80    4791.601    11417.55    7280.908    7410.219    22088.54
    ## 6    12644.80    4791.601    11417.55    7280.908    7410.219    22088.54
    ##   1051_15_068 1043_11_058 1045_13_061 1047_13_063 1061_17_079 1059_17_077
    ## 1    27290.59    29320.68    36906.17    38293.96   41689.906    32799.84
    ## 2    27290.59    29320.68    36906.17    38293.96   41689.906    32799.84
    ## 3    46237.28    39160.27    20905.79    35087.88   50573.789    37116.84
    ## 4    46237.28    39160.27    20905.79    35087.88   50573.789    37116.84
    ## 5    28865.32     6468.92    18683.06    14712.74    9504.914    17715.73
    ## 6    28865.32     6468.92    18683.06    14712.74    9504.914    17715.73
    ##   1068_19_087 1015_5_028 1066_19_085 1057_15_074 1044_13_060 1074_21_094
    ## 1    33285.69   44214.22    26814.53    31772.42   74542.961    22350.04
    ## 2    33285.69   44214.22    26814.53    31772.42   74542.961    22350.04
    ## 3    35486.48   14713.38    30550.87    20521.46   41070.328    44747.61
    ## 4    35486.48   14713.38    30550.87    20521.46   41070.328    44747.61
    ## 5    12941.47         NA    13736.32    30811.87    1263.334    15031.01
    ## 6    12941.47         NA    13736.32    30811.87    1263.334    15031.01
    ##   1056_15_073 1063_17_081 1069_19_088 1070_19_089 1071_19_090 1073_21_093
    ## 1    32365.32    55035.33    24627.73   31706.551    30007.19    31194.64
    ## 2    32365.32    55035.33    24627.73   31706.551    30007.19    31194.64
    ## 3    31104.39    26444.88    35722.28   39341.266    34431.98    35499.24
    ## 4    31104.39    26444.88    35722.28   39341.266    34431.98    35499.24
    ## 5    22873.80    15065.03    17330.65    9536.247    13046.53    12894.24
    ## 6    22873.80    15065.03    17330.65    9536.247    13046.53    12894.24
    ##   1065_19_084 1077_21_097 1060_17_078 1055_15_072 1075_21_095 1081_23_102
    ## 1    37810.50   17524.598  44080.1406    40505.07    40108.16   22763.008
    ## 2    37810.50   17524.598  44080.1406    40505.07    40108.16   22763.008
    ## 3    35451.66   31593.078  42624.1016    39685.69    32341.06   34898.891
    ## 4    35451.66   31593.078  42624.1016    39685.69    32341.06   34898.891
    ## 5    37387.27    7165.357    768.6859    33270.22    22451.71    9066.538
    ## 6    37387.27    7165.357    768.6859    33270.22    22451.71    9066.538
    ##   1089_25_111 1079_23_100 1058_17_076 1101_29_125 1088_25_110 1083_23_104
    ## 1   25207.078    34094.61    51771.14    24695.00   20457.891    38854.05
    ## 2   25207.078    34094.61    51771.14    24695.00   20457.891    38854.05
    ## 3   34448.879    36356.45    14175.11    41887.39   33767.055    18140.17
    ## 4   34448.879    36356.45    14175.11    41887.39   33767.055    18140.17
    ## 5    1428.239    19680.36          NA    11796.10    5596.355    15186.31
    ## 6    1428.239    19680.36          NA    11796.10    5596.355    15186.31
    ##   1082_23_103 1067_19_086 1090_25_112 1084_23_105 1102_29_126 1091_25_113
    ## 1    30509.62   30878.744    33287.71    57551.87    22418.34    25103.63
    ## 2    30509.62   30878.744    33287.71    57551.87    22418.34    25103.63
    ## 3    26506.68    8079.695    29109.50    35573.06    52447.38    25609.38
    ## 4    26506.68    8079.695    29109.50    35573.06    52447.38    25609.38
    ## 5    29676.77    5860.513    15811.58    31383.47     5321.65    12534.20
    ## 6    29676.77    5860.513    15811.58    31383.47     5321.65    12534.20
    ##   1087_25_109 NC28173_7_041 1072_21_092 1094_27_117 1095_27_118 1092_25_114
    ## 1    33995.05     26381.053    35625.85    24713.68    37406.63    33851.84
    ## 2    33995.05     26381.053    35625.85    24713.68    37406.63    33851.84
    ## 3    39055.54     27348.859    37164.93    30007.59    38345.87    45905.78
    ## 4    39055.54     27348.859    37164.93    30007.59    38345.87    45905.78
    ## 5    12029.20      1765.003          NA    25372.91    13026.90    32147.29
    ## 6    12029.20      1765.003          NA    25372.91    13026.90    30455.05
    ##   OH8243_7_038 1099_27_122 1078_21_098 1097_27_120 1096_27_119 1086_25_108
    ## 1    25054.664    31101.32    40683.86    26057.97    32389.38  81419.9375
    ## 2    25054.664    31101.32    40683.86    26057.97    32389.38  81419.9375
    ## 3    31343.670    31371.69    43608.37    55376.84    33137.48  44680.5430
    ## 4    31343.670    31371.69    43608.37    55376.84    33137.48  44680.5430
    ## 5     2930.516    12671.18          NA    21073.95    13812.59    652.7218
    ## 6     2930.516    12671.18          NA    21073.95    13812.59    652.7218
    ##   1104_29_128 1093_27_116 QC4_4_027 1105_29_129 QC6_6_035 QC10_10_051
    ## 1    38956.54   28699.266  36469.80    42033.24  32471.32    31703.25
    ## 2    38956.54   28699.266  36469.80    42033.24  32471.32    31703.25
    ## 3    24576.62    7526.663  22872.97    58895.51  21509.74    22965.00
    ## 4    24576.62    7526.663  22872.97    58895.51  21509.74    22965.00
    ## 5    15926.47          NA  15498.30    12539.64  13342.73    13693.94
    ## 6    15926.47          NA  15498.30    12539.64  13342.73    13693.94
    ##   1106_29_130 QC16_16_075 1098_27_121 QC12_12_059 QC8_8_043 QC2_2_019
    ## 1    37241.99    33429.03  49566.6992    32092.48  32396.12  36491.33
    ## 2    37241.99    33429.03  49566.6992    32092.48  32396.12  36491.33
    ## 3    50551.69    20884.95  39650.4727    21365.26  23788.25  23923.99
    ## 4    50551.69    20884.95  39650.4727    21365.26  23788.25  23923.99
    ## 5    18593.43    13461.82    668.6479    12833.02  13190.19  14376.58
    ## 6    18593.43    13461.82    668.6479    12833.02  13190.19  14376.58
    ##   1108_31_134 QC14_14_067 QC20_20_091 QC18_18_083 QC28_28_123 QC24_24_107
    ## 1    46188.59    32283.70    31491.38    31135.02    33043.91    33227.86
    ## 2    46188.59    32283.70    31491.38    31135.02    33043.91    33227.86
    ## 3    50958.50    23258.71    22374.13    18739.85    24087.93    25921.85
    ## 4    50958.50    23258.71    22374.13    18739.85    24087.93    25921.85
    ## 5    14798.74    14302.13    14951.52    12926.49    12153.80    13067.70
    ## 6    14798.74    14302.13    14951.52    12926.49    12153.80    13067.70
    ##   1100_29_124 QC26_26_115 QC0_0_011 QC30_30_131 TAINAN_31_135 QC22_22_099
    ## 1    54715.12    32171.46  36954.71    33248.15      43252.49    32531.70
    ## 2    54715.12    32171.46  36954.71    33248.15      43252.49    32531.70
    ## 3    18300.66    26181.96  24362.04    26241.43      28649.83    22154.52
    ## 4    18300.66    26181.96  24362.04    26241.43      28649.83    22154.52
    ## 5          NA    12592.44  14645.66    12130.53      21955.80    12486.18
    ## 6          NA    12592.44  14645.66    13500.00      21955.80    12486.18
    ##   QC32_32_136
    ## 1    33799.56
    ## 2    33799.56
    ## 3    25507.86
    ## 4    25507.86
    ## 5    15458.10
    ## 6    15458.10

Yes, let’s remove them.

``` r
MZ_RT_filt_PBremoved_nodupes <- MZ_RT_filt_PBremoved %>%
  rownames_to_column(var = "row_number")

# find the row numbers with the duplicates
get_dupes(MZ_RT_filt_PBremoved_nodupes, mz_rt)
```

    ##             mz_rt dupe_count row_number 1002_1_013 1005_1_016 1028_7_042
    ## 1 239.0776_0.6677          2       1457   18398.26  40836.883  25596.115
    ## 2 239.0776_0.6677          2       3785   18398.26  40836.883  25596.115
    ## 3 270.9855_0.7435          2        307   31931.09  35832.715  31064.273
    ## 4 270.9855_0.7435          2       3459   31931.09  35832.715  31064.273
    ## 5  677.1727_4.303          2       1403    1029.13   4552.599   9868.488
    ## 6  677.1727_4.303          2       3887    1029.13   4552.599   9868.488
    ##   1006_1_017 1011_3_023 1012_3_024 1008_3_020 1003_1_014 1004_1_015 1010_3_022
    ## 1   34077.79   24453.42   23670.80   30323.73   60174.62   30499.17   34480.30
    ## 2   34077.79   24453.42   23670.80   30323.73   60174.62   30499.17   34480.30
    ## 3   27352.09   26064.35   25829.42   53138.13   44746.21   31656.12   50254.85
    ## 4   27352.09   26064.35   25829.42   53138.13   44746.21   31656.12   50254.85
    ## 5   18662.72   19399.51   11100.51   11389.76   11889.22   32046.98   14976.43
    ## 6   18662.72   19399.51   11100.51   11389.76   11889.22   32046.98   14976.43
    ##   1026_7_040 1031_9_046 1025_7_039 1007_1_018 1020_5_033 1016_5_029 1014_3_026
    ## 1   26423.76  26000.303   26258.80   32385.83   27665.63  27169.596   30924.73
    ## 2   26423.76  26000.303   26258.80   32385.83   27665.63  27169.596   30924.73
    ## 3   27008.75  28121.031   33599.97   35116.06   28601.47  26155.508   28363.31
    ## 4   27008.75  28121.031   33599.97   35116.06   28601.47  26155.508   28363.31
    ## 5   10299.67   3429.637   11847.22   36846.87   17553.60   8155.577   16408.99
    ## 6   10299.67   3429.637   11847.22   36846.87   17553.60   8155.577   16408.99
    ##   1019_5_032 1013_3_025 1029_9_044 1023_7_037 1018_5_031 1022_7_036 1021_5_034
    ## 1  32761.664 40810.3750  16785.041   22237.80   25636.65   31415.17   30469.16
    ## 2  32761.664 40810.3750  16785.041   22237.80   25636.65   31415.17   30469.16
    ## 3  37751.633 41156.9062  33820.180   40601.19   25459.74   35669.09   28086.34
    ## 4  37751.633 41156.9062  33820.180   40601.19   25459.74   35669.09   28086.34
    ## 5   6993.417   674.3444   9710.306   13211.57   15605.22   18984.73   30304.05
    ## 6   6993.417   674.3444   9710.306   13211.57   15605.22   18984.73   36257.68
    ##   1032_9_047 1009_3_021 1039_11_054 1035_9_049 1046_13_062 1030_9_045
    ## 1  22598.445  41088.746   25952.139  28047.129    22712.89  31025.541
    ## 2  22598.445  41088.746   25952.139  28047.129    22712.89  31025.541
    ## 3  27039.873   4337.607   31361.414  39152.219    33796.79  25816.801
    ## 4  27039.873   4337.607   31361.414  39152.219    33796.79  25816.801
    ## 5   6306.936  16678.168    3359.229   6854.434    11038.56   7249.560
    ## 6   6306.936  16678.168    3359.229   6854.434    11038.56   7814.693
    ##   1040_11_055 1033_31_132 1038_11_053 1041_11_056 1036_9_050 1034_9_048
    ## 1    21571.04    50582.26   13343.054   18304.398   20093.10   31052.56
    ## 2    21571.04    50582.26   13343.054   18304.398   20093.10   31052.56
    ## 3    27523.80    43392.93   23035.375   24987.023   31239.41   28920.15
    ## 4    27523.80    43392.93   23035.375   24987.023   31239.41   28920.15
    ## 5    10789.50    33364.76    9680.562    5536.901   11632.77   14842.04
    ## 6    10789.50    33364.76    9680.562    5536.901   11632.77   14842.04
    ##   1052_15_069 1050_13_066 1048_13_064 1037_11_052 1053_15_070 1054_15_071
    ## 1   30123.566    28652.63   71729.469    22287.97   31470.201   25497.828
    ## 2   30123.566    28652.63   71729.469    22287.97   31470.201   25497.828
    ## 3   31446.117    31721.76   35101.863    28286.07   40586.426   38719.844
    ## 4   31446.117    31721.76   35101.863    28286.07   40586.426   38719.844
    ## 5    6412.228    12644.80    4791.601    11417.55    7280.908    7410.219
    ## 6    6860.068    12644.80    4791.601    11417.55    7280.908    7410.219
    ##   1049_13_065 1051_15_068 1043_11_058 1045_13_061 1047_13_063 1061_17_079
    ## 1    29934.35    27290.59    29320.68    36906.17    38293.96   41689.906
    ## 2    29934.35    27290.59    29320.68    36906.17    38293.96   41689.906
    ## 3    36615.30    46237.28    39160.27    20905.79    35087.88   50573.789
    ## 4    36615.30    46237.28    39160.27    20905.79    35087.88   50573.789
    ## 5    22088.54    28865.32     6468.92    18683.06    14712.74    9504.914
    ## 6    22088.54    28865.32     6468.92    18683.06    14712.74    9504.914
    ##   1059_17_077 1068_19_087 1015_5_028 1066_19_085 1057_15_074 1044_13_060
    ## 1    32799.84    33285.69   44214.22    26814.53    31772.42   74542.961
    ## 2    32799.84    33285.69   44214.22    26814.53    31772.42   74542.961
    ## 3    37116.84    35486.48   14713.38    30550.87    20521.46   41070.328
    ## 4    37116.84    35486.48   14713.38    30550.87    20521.46   41070.328
    ## 5    17715.73    12941.47         NA    13736.32    30811.87    1263.334
    ## 6    17715.73    12941.47         NA    13736.32    30811.87    1263.334
    ##   1074_21_094 1056_15_073 1063_17_081 1069_19_088 1070_19_089 1071_19_090
    ## 1    22350.04    32365.32    55035.33    24627.73   31706.551    30007.19
    ## 2    22350.04    32365.32    55035.33    24627.73   31706.551    30007.19
    ## 3    44747.61    31104.39    26444.88    35722.28   39341.266    34431.98
    ## 4    44747.61    31104.39    26444.88    35722.28   39341.266    34431.98
    ## 5    15031.01    22873.80    15065.03    17330.65    9536.247    13046.53
    ## 6    15031.01    22873.80    15065.03    17330.65    9536.247    13046.53
    ##   1073_21_093 1065_19_084 1077_21_097 1060_17_078 1055_15_072 1075_21_095
    ## 1    31194.64    37810.50   17524.598  44080.1406    40505.07    40108.16
    ## 2    31194.64    37810.50   17524.598  44080.1406    40505.07    40108.16
    ## 3    35499.24    35451.66   31593.078  42624.1016    39685.69    32341.06
    ## 4    35499.24    35451.66   31593.078  42624.1016    39685.69    32341.06
    ## 5    12894.24    37387.27    7165.357    768.6859    33270.22    22451.71
    ## 6    12894.24    37387.27    7165.357    768.6859    33270.22    22451.71
    ##   1081_23_102 1089_25_111 1079_23_100 1058_17_076 1101_29_125 1088_25_110
    ## 1   22763.008   25207.078    34094.61    51771.14    24695.00   20457.891
    ## 2   22763.008   25207.078    34094.61    51771.14    24695.00   20457.891
    ## 3   34898.891   34448.879    36356.45    14175.11    41887.39   33767.055
    ## 4   34898.891   34448.879    36356.45    14175.11    41887.39   33767.055
    ## 5    9066.538    1428.239    19680.36          NA    11796.10    5596.355
    ## 6    9066.538    1428.239    19680.36          NA    11796.10    5596.355
    ##   1083_23_104 1082_23_103 1067_19_086 1090_25_112 1084_23_105 1102_29_126
    ## 1    38854.05    30509.62   30878.744    33287.71    57551.87    22418.34
    ## 2    38854.05    30509.62   30878.744    33287.71    57551.87    22418.34
    ## 3    18140.17    26506.68    8079.695    29109.50    35573.06    52447.38
    ## 4    18140.17    26506.68    8079.695    29109.50    35573.06    52447.38
    ## 5    15186.31    29676.77    5860.513    15811.58    31383.47     5321.65
    ## 6    15186.31    29676.77    5860.513    15811.58    31383.47     5321.65
    ##   1091_25_113 1087_25_109 NC28173_7_041 1072_21_092 1094_27_117 1095_27_118
    ## 1    25103.63    33995.05     26381.053    35625.85    24713.68    37406.63
    ## 2    25103.63    33995.05     26381.053    35625.85    24713.68    37406.63
    ## 3    25609.38    39055.54     27348.859    37164.93    30007.59    38345.87
    ## 4    25609.38    39055.54     27348.859    37164.93    30007.59    38345.87
    ## 5    12534.20    12029.20      1765.003          NA    25372.91    13026.90
    ## 6    12534.20    12029.20      1765.003          NA    25372.91    13026.90
    ##   1092_25_114 OH8243_7_038 1099_27_122 1078_21_098 1097_27_120 1096_27_119
    ## 1    33851.84    25054.664    31101.32    40683.86    26057.97    32389.38
    ## 2    33851.84    25054.664    31101.32    40683.86    26057.97    32389.38
    ## 3    45905.78    31343.670    31371.69    43608.37    55376.84    33137.48
    ## 4    45905.78    31343.670    31371.69    43608.37    55376.84    33137.48
    ## 5    32147.29     2930.516    12671.18          NA    21073.95    13812.59
    ## 6    30455.05     2930.516    12671.18          NA    21073.95    13812.59
    ##   1086_25_108 1104_29_128 1093_27_116 QC4_4_027 1105_29_129 QC6_6_035
    ## 1  81419.9375    38956.54   28699.266  36469.80    42033.24  32471.32
    ## 2  81419.9375    38956.54   28699.266  36469.80    42033.24  32471.32
    ## 3  44680.5430    24576.62    7526.663  22872.97    58895.51  21509.74
    ## 4  44680.5430    24576.62    7526.663  22872.97    58895.51  21509.74
    ## 5    652.7218    15926.47          NA  15498.30    12539.64  13342.73
    ## 6    652.7218    15926.47          NA  15498.30    12539.64  13342.73
    ##   QC10_10_051 1106_29_130 QC16_16_075 1098_27_121 QC12_12_059 QC8_8_043
    ## 1    31703.25    37241.99    33429.03  49566.6992    32092.48  32396.12
    ## 2    31703.25    37241.99    33429.03  49566.6992    32092.48  32396.12
    ## 3    22965.00    50551.69    20884.95  39650.4727    21365.26  23788.25
    ## 4    22965.00    50551.69    20884.95  39650.4727    21365.26  23788.25
    ## 5    13693.94    18593.43    13461.82    668.6479    12833.02  13190.19
    ## 6    13693.94    18593.43    13461.82    668.6479    12833.02  13190.19
    ##   QC2_2_019 1108_31_134 QC14_14_067 QC20_20_091 QC18_18_083 QC28_28_123
    ## 1  36491.33    46188.59    32283.70    31491.38    31135.02    33043.91
    ## 2  36491.33    46188.59    32283.70    31491.38    31135.02    33043.91
    ## 3  23923.99    50958.50    23258.71    22374.13    18739.85    24087.93
    ## 4  23923.99    50958.50    23258.71    22374.13    18739.85    24087.93
    ## 5  14376.58    14798.74    14302.13    14951.52    12926.49    12153.80
    ## 6  14376.58    14798.74    14302.13    14951.52    12926.49    12153.80
    ##   QC24_24_107 1100_29_124 QC26_26_115 QC0_0_011 QC30_30_131 TAINAN_31_135
    ## 1    33227.86    54715.12    32171.46  36954.71    33248.15      43252.49
    ## 2    33227.86    54715.12    32171.46  36954.71    33248.15      43252.49
    ## 3    25921.85    18300.66    26181.96  24362.04    26241.43      28649.83
    ## 4    25921.85    18300.66    26181.96  24362.04    26241.43      28649.83
    ## 5    13067.70          NA    12592.44  14645.66    12130.53      21955.80
    ## 6    13067.70          NA    12592.44  14645.66    13500.00      21955.80
    ##   QC22_22_099 QC32_32_136
    ## 1    32531.70    33799.56
    ## 2    32531.70    33799.56
    ## 3    22154.52    25507.86
    ## 4    22154.52    25507.86
    ## 5    12486.18    15458.10
    ## 6    12486.18    15458.10

``` r
MZ_RT_filt_PBremoved_nodupes <- MZ_RT_filt_PBremoved_nodupes %>%
  filter(!row_number %in% c(3785, 3459, 3887)) # remove one of the dupes for each feature

dim(MZ_RT_filt_PBremoved_nodupes)  
```

    ## [1] 3945  118

``` r
# get rid of row_number
MZ_RT_filt_PBremoved_nodupes <- MZ_RT_filt_PBremoved_nodupes %>%
  select(-row_number)
```

## Save your file

Now you have a list of features present in your samples after filtering
for CV in QCs, and removing all the extraneous columns we added to help
us do this, along with removing any process blanks.

``` r
write_csv(MZ_RT_filt_PBremoved_nodupes,
          "data/BNB_neg_filtered_template.csv")
```

# Adherence to Canada’s Food Guide 2019 and nutrient intakes
Didier Brassard

- [Introduction](#introduction)
- [1. Requirements and file structure](#requirements-and-file-structure)
  - [1.1 Data and software](#data-and-software)
  - [1.2 SAS macros](#sas-macros)
  - [1.3 Structure](#structure)
- [Reference](#reference)
- [Session Info](#session-info)

# Introduction

This repository presents the analysis code for the study *Relationship
between adherence to the 2019 Canada’s Food Guide recommendations on
healthy food choices and nutrient intakes in older adults* published in
the Journal of Nutrition (Brassard & Chevalier, 2023).

The general objective was to assess the relationship between adherence
to CFG recommendations on healthy food choices, measured using the
HEFI-2019 (Brassard et al, 2022a; 2022b), and intake of key nutrients in
adults 65 years and older from the Canadian Community Health Survey
(CCHS) 2015 - Nutrition.

# 1. Requirements and file structure

## 1.1 Data and software

Complete raw data file from the Canadian Community Health Survey are
available publicly [upon request to Statistics
Canada](https://www150.statcan.gc.ca/n1/en/catalogue/82M0024X2018001).

The main analyses were executed in SAS (v9.4; maintenance release
9.04.01M7P08052020) in a Windows 10 64-bits environment. The manuscript
file was successfully executed using Quarto and R (version 4.2.2) on
macOS Big Sur 10.16 (64-bit; see complete session information details
below).

## 1.2 SAS macros

SAS macros made by the National Cancer Institute (NCI) were used to
perform measurement error correction (Zhang et al., 2011). The macros
are [available on the NCI
website](https://prevention.cancer.gov/research-groups/biometry/measurement-error-impact/software-measurement-error/several-regularly-consumed-or-0),
but are included in the present repository in the `/Macros/` folder.

The HEFI-2019 scoring algorithm SAS macro was used and it is also
[available online](https://github.com/didierbrassard/hefi2019) and as
well as included in the `/Macros/` folder.

## 1.3 Structure

Key steps of the main analysis are presented in separate `.sas` files.
Each file would need to be executed in sequential order for successful
execution.

Summary statistics output of each step are provided. 

### Folders

...

# Reference

Brassard D, Chevalier S. Relationship between adherence to the 2019
Canada’s Food Guide recommendations on healthy food choices and nutrient
intakes in older adults. J Nutr 2023. doi:
[10.1016/j.tjnut.2023.07.005](http://doi.org/10.1016/j.tjnut.2023.07.005).

> **Pre-print version of this manuscript**: Brassard D, Chevalier S.
> Relationship between adherence to the 2019 Canada’s Food Guide
> recommendations on healthy food choices and nutrient intakes in older
> adults. medRxiv 2023:2023.02.13.23285868. doi:
> [10.1101/2023.02.13.23285868](https://doi.org/10.1101/2023.02.13.23285868).

Brassard D, Elvidge Munene LA, St-Pierre S, et al. *Development of the
Healthy Eating Food Index (HEFI)-2019 measuring adherence to Canada’s
Food Guide 2019 recommendations on healthy food choices*. Appl Physiol
Nutr Metab 2022a;47:595-610. doi:
[10.1139/apnm-2021-0415](https://doi.org/10.1139/apnm-2021-0415).

Brassard D, Elvidge Munene LA, St-Pierre S, et al. *Evaluation of the
Healthy Eating Food Index (HEFI)-2019 measuring adherence to Canada’s
Food Guide 2019 recommendations on healthy food choices*. Appl Physiol
Nutr Metab 2022b;47(5):582-94. doi:
[10.1139/apnm-2021-0416](https://doi.org/10.1139/apnm-2021-0416).

Zhang S, Midthune D, Guenther PM, et al. *A New Multivariate Measurement
Error Model with Zero-Inflated Dietary Data, and Its Application to
Dietary Assessment*. Ann Appl Stat 2011;5(2B):1456-87. doi:
[10.1214/10-AOAS446](https://doi.org/10.1214/10-AOAS446).

# Session Info

<details>
<summary>
Expand for details
</summary>

    [1] "2023-07-28 15:32:46 EDT"

    R version 4.2.2 (2022-10-31)
    Platform: x86_64-apple-darwin17.0 (64-bit)
    Running under: macOS Big Sur ... 10.16

    Matrix products: default
    BLAS:   /Library/Frameworks/R.framework/Versions/4.2/Resources/lib/libRblas.0.dylib
    LAPACK: /Library/Frameworks/R.framework/Versions/4.2/Resources/lib/libRlapack.dylib

    locale:
    [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8

    attached base packages:
    [1] stats     graphics  grDevices utils     datasets  methods   base     

    loaded via a namespace (and not attached):
     [1] compiler_4.2.2    fastmap_1.1.1     cli_3.6.1         tools_4.2.2      
     [5] htmltools_0.5.5   rstudioapi_0.15.0 yaml_2.3.7        rmarkdown_2.23   
     [9] knitr_1.43        jsonlite_1.8.7    xfun_0.39         digest_0.6.33    
    [13] rlang_1.1.1       evaluate_0.21    

</details>

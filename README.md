# Adherence to Canada’s Food Guide 2019 and nutrient intakes
Didier Brassard

- [Introduction](#introduction)
  - [Quick links](#quick-links)
- [1. Requirements and file structure](#requirements-and-file-structure)
  - [1.1 Data and software](#data-and-software)
  - [1.2 SAS macros and R functions](#sas-macros-and-r-functions)
    - [National Cancer Institute (NCI)
      method](#national-cancer-institute-nci-method)
    - [HEFI-2019 macro](#hefi-2019-macro)
    - [In-house SAS macro](#in-house-sas-macro)
    - [In-house R functions](#in-house-r-functions)
  - [1.3 Structure](#structure)
    - [Folders](#folders)
- [2. Description of analyses and
  codes](#description-of-analyses-and-codes)
- [Reference](#reference)
- [Session Info](#session-info)

# Introduction

This repository presents the analysis code for the study *Relationship
between adherence to the 2019 Canada’s Food Guide recommendations on
healthy food choices and nutrient intakes in older adults* published in
the Journal of Nutrition ([Brassard & Chevalier,
2023](https://authors.elsevier.com/a/1hW0W_WUPSMQ4)).

The general objective was to assess the relationship between adherence
to CFG recommendations on healthy food choices, measured using the
HEFI-2019 (Brassard et al, 2022a; 2022b), and intake of key nutrients in
adults 65 years and older from the Canadian Community Health Survey
(CCHS) 2015 - Nutrition.

## Quick links

- [**Journal of Nutrition**: Free PDF of the article (until Sept. 19,
  2023)](https://authors.elsevier.com/a/1hW0W_WUPSMQ4)
- [**medRxiv**: Pre-print version of the
  study](https://doi.org/10.1101/2023.02.13.23285868)
- [**GitHub**: Main article code, text, figures and
  tables](https://didierbrassard.github.io/hefi2019_nutrient/9.1-HEFI19-NUT_ARTICLE.html)
- [**GitHub**: Supplemental material code, text, figures and
  tables](https://didierbrassard.github.io/hefi2019_nutrient/9.2-HEFI19-NUT_OSM.html)
- [**OSF**: Raw and processed data used for this
  study](https://osf.io/6na42/)

# 1. Requirements and file structure

## 1.1 Data and software

All data used for this project are [available on OSF in a
repository](https://osf.io/6na42/). To run analyses, all data from the
OSF repository should first be downloaded and placed in the same folder
as the project’s codes. More precisely, data in the `/data/` and in the
`/NCI/` folder from OSF should be put in the same folder as the codes
`1.0-Data_preparation.sas`, `2.1-NCI_multi_pro_mcmc.sas`, etc.

The main analyses were executed in SAS studio (v3.81). The manuscript
and supplemental file were successfully executed using Quarto and R
(version 4.2.2) on macOS Big Sur 10.16 (64-bit; see complete session
information details below). Analyses were conducted partly in SAS since
the National Cancer Institute (NCI) method is only available in SAS
codes at the moment.

A copy of all raw data used in this study can also be obtained
elsewhere:

- Free sugars estimate of foods ([Rana et
  al. 2021](https://www.mdpi.com/2072-6643/13/5/1471));
- HEFI-2019 categories and grams equivalent to one reference amount
  ([Health Canada,
  2023](https://open.canada.ca/data/en/dataset/29892c85-2ff5-484c-873c-f494ffba6e1b))
- Canadian Community Health Survey 2015 - Nutrition ([Statistics Canada,
  2018](https://www150.statcan.gc.ca/n1/en/catalogue/82M0024X2018001))

## 1.2 SAS macros and R functions

### National Cancer Institute (NCI) method

SAS macros made by the NCI were used to perform measurement error
correction (Zhang et al., 2011). The macros are [available on the NCI
website](https://prevention.cancer.gov/research-groups/biometry/measurement-error-impact/software-measurement-error/several-regularly-consumed-or-0),
but are included in the present repository in the `/Macros/` folder.

- `boxcox_survey.macro.v1.2.sas`
- `std_cov_boxcox24hr_conday_minamt_macro_v2.0.sas`
- `multivar_mcmc_macro_v2.1.sas`
- `multivar_distrib_macro_v2.1.sas`
- `percentiles_survey.macro.v1.1.sas`

Please note that small modifications were made vs. the original macros
for efficiency purpose, i.e., a relative path is used in
`multivar_mcmc_macro_v2.1.sas` to save trace plots and the default
cut-point probability was modified in
`percentiles_survey.macro.v1.1.sas` to $Pr(X< x)$ instead of
$Pr(X\le x)$.

### HEFI-2019 macro

The HEFI-2019 scoring algorithm SAS macro `hefi2019.scoring.macro.sas`
was used ([available in a GitHub
repository](https://github.com/didierbrassard/hefi2019)) and is included
in the `/Macros/` folder.

### In-house SAS macro

The `boot_auxiliary.sas` file includes a suite of macros used to observe
and analyze bootstrap data. These macros are used in
7.0-Bootstrap_variance.sas to calculate bootstrap variance.

### In-house R functions

Details are provided within the R script for the two functions made for
this project: `statDistrib.R` and `post_mcmc_auxiliary.R`.
`statDistrib.R` performs the same operations as
`percentiles_survey.macro.v1.1.sas` except cut-point probabilities.
`post_mcmc_auxiliary.R` load data generated by the
`multivar_distrib_macro_v2.1.sas` macros (MCMC output) and calculates
energy-adjusted correlations for HEFI-2019 dietary constituents.

## 1.3 Structure

Key steps of the main analysis are performed in separate `.sas` files.
Each file would need to be executed in sequential order for successful
execution, as indicated by the numerical prefix.

The flowchart below illustrates the relationship between codes, data and
folders in more details.

``` mermaid
%%{init: {'theme': 'neutral' } }%%
flowchart TB
  data1[(Data/Raw)]
  data2[(Data/Processed)]
  data3[(NCI/)]
  code1[1-Data_preparation.sas]
  code2[2, 3, ... 6-NCI_multi_*.sas]
  code3[7.0-Bootstrap_variance.sas]
  qmd1(<b>10-HEFI19NUT_ARTICLE.qmd)
  qmd2(<b>20-HEFI19NUT_OSM.qmd)
  out("<b>Text, figures, tables")
  data1 --> code1 --> data2
  data2 --> code2 --> data3
  data3 --> code3 
  code3 --> data3
  data1 --> code2
  macros1["Macros/*.sas"]
  macros1 --> code2
  data2 --> qmd1
  data3 --> qmd1
  data4[(Data/Results)]
  data4 --> qmd2
  qmd1 --> data4
  data4 --> qmd1
  data3 --> qmd2
  qmd1 --> out
  qmd2 --> out
  
```

Summary statistics output of each step are provided in the OSF
repository. Thus, the *10_HEFI19NUT_ARTICLE.qmd* Quarto document can be
used to generate the manuscript including tables and figures directly,
without having to run all analyses beforehand. The
*20_HEFI19NUT_OSM.qmd* Quarto document generates the online supplemental
material file.

### Folders

Available on GitHub (here):

- **Macros/:** SAS macros (i.e., scripts) needed for main analysis based
  on `.sas` files 1 to 7. R functions for repetitive steps are also
  included in this folder;
- **docs/**: text output of the quarto documents.

[Available on OSF](https://osf.io/6na42/):

- **Data/:** raw datafile (unedited), processed data and additional
  results generated by the Quarto documents;
- **NCI/:** results data generated by the NCI methods (see `.sas` files
  1 to 7);
- **Manuscript/:** manuscript tables and figures generated by the Quarto
  document.

# 2. Description of analyses and codes

The flowchart below presents a generic overview of the main analyses in
this study. Complete details about analyses are provided in the article.

``` mermaid
%%{init: {'theme': 'neutral' } }%%
flowchart TB
  A["1) Measurement error correction<br>(<I>NCI multivariate method</I>)"]
  B["2) Continuous relationship<br>(<I>Linear regression</I>)"]
  C["3) Nutrient intake adequacy<br>(<I>Logistic regression</I>)"]
  D["4) Bootstrap variance estimation<br>(<I>bootstrap replicate weights</I>)"]
  A-->B
  B-->C
  C-->D
  D-."500<br>repetitions".->A
```

<div id="zevrmshmfn" style="padding-left:0px;padding-right:0px;padding-top:10px;padding-bottom:10px;overflow-x:auto;overflow-y:auto;width:auto;height:auto;">
  &#10;  <table class="gt_table" data-quarto-disable-processing="false" data-quarto-bootstrap="false">
  <thead>
    <tr class="gt_heading">
      <td colspan="3" class="gt_heading gt_title gt_font_normal gt_bottom_border" style>Brief description of SAS and QMD codes</td>
    </tr>
    &#10;    <tr class="gt_col_headings">
      <th class="gt_col_heading gt_columns_bottom_border gt_left" rowspan="1" colspan="1" scope="col" id="Name">Name</th>
      <th class="gt_col_heading gt_columns_bottom_border gt_left" rowspan="1" colspan="1" scope="col" id="Description">Description</th>
      <th class="gt_col_heading gt_columns_bottom_border gt_center" rowspan="1" colspan="1" scope="col" id="Link">Link</th>
    </tr>
  </thead>
  <tbody class="gt_table_body">
    <tr><td headers="Name" class="gt_row gt_left">1.0-Data_preparation.sas</td>
<td headers="Description" class="gt_row gt_left"><div class='gt_from_md'><p>Based on CCHS 2015 - Nutrition PUMF data, calculate intakes of HEFI-2019 dietary constituents per respondent and recall. Prepare an input data for the NCI MCMC multivariate method. Look at proportion of respondents with zero intake</p>
</div></td>
<td headers="Link" class="gt_row gt_center"><a href="1.0-Data_preparation.sas">Open code</a></td></tr>
    <tr><td headers="Name" class="gt_row gt_left">2.1-NCI_multi_pro_mcmc.sas</td>
<td headers="Description" class="gt_row gt_left"><div class='gt_from_md'><p>Apply the multivariate measurement error correction model to HEFI-2019 dietary constituents and <strong>protein intakes</strong></p>
</div></td>
<td headers="Link" class="gt_row gt_center"><a href="2.1-NCI_multi_pro_mcmc.sas">Open code</a></td></tr>
    <tr><td headers="Name" class="gt_row gt_left">2.2-NCI_multi_pro_res.sas</td>
<td headers="Description" class="gt_row gt_left"><div class='gt_from_md'><p>Estimate usual intakes distribution, relationship between variables and prevalence of inadequacy among pseudo-individual (example code of original sample only)</p>
</div></td>
<td headers="Link" class="gt_row gt_center"><a href="2.2-NCI_multi_pro_res.sas">Open code</a></td></tr>
    <tr><td headers="Name" class="gt_row gt_left">3.1-NCI_multi_ca_vit_d_mcmc.sas</td>
<td headers="Description" class="gt_row gt_left"><div class='gt_from_md'><p>Apply the multivariate measurement error correction model to HEFI-2019 dietary constituents, <strong>calcium and vitamin D intakes</strong></p>
</div></td>
<td headers="Link" class="gt_row gt_center"><a href="3.1-NCI_multi_ca_vit_d_mcmc.sas">Open code</a></td></tr>
    <tr><td headers="Name" class="gt_row gt_left">3.2-NCI_multi_ca_vit_d_res.sas</td>
<td headers="Description" class="gt_row gt_left"><div class='gt_from_md'><p>Estimate usual intakes distribution, relationship between variables and prevalence of inadequacy among pseudo-individual (example code of original sample only)</p>
</div></td>
<td headers="Link" class="gt_row gt_center"><a href="3.2-NCI_multi_ca_vit_d_res.sas">Open code</a></td></tr>
    <tr><td headers="Name" class="gt_row gt_left">4.1-NCI_multi_miscA_mcmc.sas</td>
<td headers="Description" class="gt_row gt_left"><div class='gt_from_md'><p>Apply the multivariate measurement error correction model to HEFI-2019 dietary constituents, <strong>iron, zinc, vitamin B6 and B12</strong></p>
</div></td>
<td headers="Link" class="gt_row gt_center"><a href="4.1-NCI_multi_miscA_mcmc.sas">Open code</a></td></tr>
    <tr><td headers="Name" class="gt_row gt_left">4.2-NCI_multi_miscA_res.sas</td>
<td headers="Description" class="gt_row gt_left"><div class='gt_from_md'><p>Estimate usual intakes distribution, relationship between variables and prevalence of inadequacy among pseudo-individual (example code of original sample only)</p>
</div></td>
<td headers="Link" class="gt_row gt_center"><a href="4.2-NCI_multi_miscA_res.sas">Open code</a></td></tr>
    <tr><td headers="Name" class="gt_row gt_left">5.1-NCI_multi_miscB_mcmc.sas</td>
<td headers="Description" class="gt_row gt_left"><div class='gt_from_md'><p>Apply the multivariate measurement error correction model to HEFI-2019 dietary constituents, <strong>folate, magnesium, fibers and potassium</strong></p>
</div></td>
<td headers="Link" class="gt_row gt_center"><a href="5.1-NCI_multi_miscB_mcmc.sas">Open code</a></td></tr>
    <tr><td headers="Name" class="gt_row gt_left">5.2-NCI_multi_miscB_res.sas</td>
<td headers="Description" class="gt_row gt_left"><div class='gt_from_md'><p>Estimate usual intakes distribution, relationship between variables and prevalence of inadequacy among pseudo-individual (example code of original sample only)</p>
</div></td>
<td headers="Link" class="gt_row gt_center"><a href="5.2-NCI_multi_miscB_res.sas">Open code</a></td></tr>
    <tr><td headers="Name" class="gt_row gt_left">6.1-NCI_multi_vit_a_mcmc.sas</td>
<td headers="Description" class="gt_row gt_left"><div class='gt_from_md'><p>Apply the multivariate measurement error correction model to HEFI-2019 dietary constituents and <strong>vitamin A</strong></p>
</div></td>
<td headers="Link" class="gt_row gt_center"><a href="6.1-NCI_multi_vit_a_mcmc.sas">Open code</a></td></tr>
    <tr><td headers="Name" class="gt_row gt_left">6.2-NCI_multi_vit_a_res.sas</td>
<td headers="Description" class="gt_row gt_left"><div class='gt_from_md'><p>Estimate usual intakes distribution, relationship between variables and prevalence of inadequacy among pseudo-individual (example code of original sample only)</p>
</div></td>
<td headers="Link" class="gt_row gt_center"><a href="6.2-NCI_multi_vit_a_res.sas">Open code</a></td></tr>
    <tr><td headers="Name" class="gt_row gt_left">7.0-Bootstrap_variance.sas</td>
<td headers="Description" class="gt_row gt_left"><div class='gt_from_md'><p>Parametric bootstrap variance estimation (i.e., 95%CI) based on the previous analyses repeated in 500 bootstrap replicate weights. The repetition is not shown, but statistical estimates obtained with each bootstrap sample are available within the <em>NCI/</em> folders  on OSF</p>
</div></td>
<td headers="Link" class="gt_row gt_center"><a href="7.0-Bootstrap_variance.sas">Open code</a></td></tr>
    <tr><td headers="Name" class="gt_row gt_left">9.1-HEFI19-NUT_ARTICLE.qmd</td>
<td headers="Description" class="gt_row gt_left"><div class='gt_from_md'><p>Quarto document used to generate the text, figures and tables of the article</p>
</div></td>
<td headers="Link" class="gt_row gt_center"><a href="9.1-HEFI19-NUT_ARTICLE.qmd">Open code</a></td></tr>
    <tr><td headers="Name" class="gt_row gt_left">9.2-HEFI19-NUT_OSM.qmd</td>
<td headers="Description" class="gt_row gt_left"><div class='gt_from_md'><p>Quarto document used to generate the online supplemental material</p>
</div></td>
<td headers="Link" class="gt_row gt_center"><a href="9.2-HEFI19-NUT_OSM.qmd">Open code</a></td></tr>
  </tbody>
  &#10;  <tfoot class="gt_footnotes">
    <tr>
      <td class="gt_footnote" colspan="3"> CCHS-2015, Canadian Community Health Survey 2015; HEFI-2019, Healthy Eating Food Index 2019; MCMC, Markov Chain Monte Carlo; NCI, National Cancer Institute; PUMF, Public Use Microdata File</td>
    </tr>
  </tfoot>
</table>
</div>

Of note, due to model complexity and sample size, the main analysis is
computationally intensive. For example, a **single run** of codes
4.1-NCI_multi_miscA_mcmc.sas and 4.2-NCI_multi_miscA_res.sas on a
dedicated SAS studio server took approximately 97 minutes. To obtain
correct variance due to the multistep modelling approach and survey
data, this analysis had to be repeated **500** times with each bootstrap
weight replicate.

# Reference

Brassard D, Chevalier S. *Relationship between adherence to the 2019
Canada’s Food Guide recommendations on healthy food choices and nutrient
intakes in older adults*. J Nutr 2023. doi:
[10.1016/j.tjnut.2023.07.005](https://authors.elsevier.com/a/1hW0W_WUPSMQ4).

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

Rana H, Mallet M-C, Gonzalez A, et al. Free Sugars Consumption in
Canada. Nutrients 2021;13(5):1471. doi:
[10.3390/nu13051471](https://doi.org/10.3390/nu13051471)

Zhang S, Midthune D, Guenther PM, et al. *A New Multivariate Measurement
Error Model with Zero-Inflated Dietary Data, and Its Application to
Dietary Assessment*. Ann Appl Stat 2011;5(2B):1456-87. doi:
[10.1214/10-AOAS446](https://doi.org/10.1214/10-AOAS446).

# Session Info

<details>
<summary>
Expand for details
</summary>

    [1] "2023-08-01 12:16:28 EDT"

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

    other attached packages:
    [1] gt_0.9.0    purrr_1.0.1 tidyr_1.3.0 dplyr_1.1.2

    loaded via a namespace (and not attached):
     [1] rstudioapi_0.15.0 xml2_1.3.5        knitr_1.43        magrittr_2.0.3   
     [5] tidyselect_1.2.0  here_1.0.1        R6_2.5.1          rlang_1.1.1      
     [9] fastmap_1.1.1     fansi_1.0.4       tools_4.2.2       xfun_0.39        
    [13] utf8_1.2.3        cli_3.6.1         withr_2.5.0       commonmark_1.9.0 
    [17] htmltools_0.5.5   yaml_2.3.7        rprojroot_2.0.3   digest_0.6.33    
    [21] tibble_3.2.1      lifecycle_1.0.3   sass_0.4.7        vctrs_0.6.3      
    [25] glue_1.6.2        evaluate_0.21     rmarkdown_2.23    compiler_4.2.2   
    [29] pillar_1.9.0      generics_0.1.3    jsonlite_1.8.7    markdown_1.7     
    [33] pkgconfig_2.0.3  

</details>

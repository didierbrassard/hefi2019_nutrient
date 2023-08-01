# Adherence to Canada’s Food Guide 2019 and nutrient intakes
Didier Brassard

- [Introduction](#introduction)
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
`percentiles_survey.macro.v1.1.sas` to $Pr(X<x)$ instead of
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
  data3 <--> code3 
  data1 --> code2
  macros1["Macros/*.sas"]
  macros1 --> code2
  data2 --> qmd1
  data3 --> qmd1
  data4[(Data/Results)]
  data4 --> qmd2
  qmd1 <--> data4
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
- **Manuscript/:** manuscript text, tables and figures generated by the
  Quarto document;
- **Docs/**: text output of the quarto documents.

[Available on OSF](https://osf.io/6na42/):

- **NCI/:** results data generated by the NCI methods (see `.sas` files
  1 to 7);
- **Data/:** raw datafile (unedited), processed data and additional
  results generated by the Quarto documents.

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

list(Name = c("1.0-Data_preparation.sas", "2.1-NCI_multi_pro_mcmc.sas", "2.2-NCI_multi_pro_res.sas", "3.1-NCI_multi_ca_vit_d_mcmc.sas", "3.2-NCI_multi_ca_vit_d_res.sas", "4.1-NCI_multi_miscA_mcmc.sas", "4.2-NCI_multi_miscA_res.sas", "5.1-NCI_multi_miscB_mcmc.sas", "5.2-NCI_multi_miscB_res.sas", "6.1-NCI_multi_vit_a_mcmc.sas", "6.2-NCI_multi_vit_a_res.sas", "7.0-Bootstrap_variance.sas", "10-HEFI19-NUT_ARTICLE.qmd", "20-HEFI19-NUT_OSM.qmd"), Description = c("Based on CCHS 2015 - Nutrition PUMF data, calculate intakes of HEFI-2019 dietary constituents per respondent and recall. Prepare an input data for the NCI MCMC mutlivaraite method. Look at proportion of respondents with zero intake", 
"Apply the multivariate measurement error correction model to HEFI-2019 dietary constituents and **protein intakes**", "Estimate usual intakes distribution, relationship between variables and prevalence of inadequacy among pseudo-individual (example code of original sample only)", "Apply the multivariate measurement error correction model to HEFI-2019 dietary constituents, **calcium and vitamin D intakes**", "Estimate usual intakes distribution, relationship between variables and prevalence of inadequacy among pseudo-individual (example code of original sample only)", 
"Apply the multivariate measurement error correction model to HEFI-2019 dietary constituents, **iron, zinc, vitamin B6 and B12**", "Estimate usual intakes distribution, relationship between variables and prevalence of inadequacy among pseudo-individual (example code of original sample only)", "Apply the multivariate measurement error correction model to HEFI-2019 dietary constituents, **folate, magnesium, fibers and potassium**", "Estimate usual intakes distribution, relationship between variables and prevalence of inadequacy among pseudo-individual (example code of original sample only)", 
"Apply the multivariate measurement error correction model to HEFI-2019 dietary constituents and **vitamin A**", "Estimate usual intakes distribution, relationship between variables and prevalence of inadequacy among pseudo-individual (example code of original sample only)", "Parametric bootstrap variance estimation (i.e., 95%CI) based on the previous analyses repeated in 500 bootstrap replicate weights. The repetition is not shown, but statistical estimates obtained with each bootstrap sample are available within the *NCI/* folders  on OSF", 
"Quarto document used to generate the text, figures and tables of the article", "Quarto document used to generate the online supplemental material"), Link = list("[Open code](1.0-Data_preparation.sas)", "[Open code](2.1-NCI_multi_pro_mcmc.sas)", "[Open code](2.2-NCI_multi_pro_res.sas)", "[Open code](3.1-NCI_multi_ca_vit_d_mcmc.sas)", "[Open code](3.2-NCI_multi_ca_vit_d_res.sas)", "[Open code](4.1-NCI_multi_miscA_mcmc.sas)", "[Open code](4.2-NCI_multi_miscA_res.sas)", "[Open code](5.1-NCI_multi_miscB_mcmc.sas)", 
    "[Open code](5.2-NCI_multi_miscB_res.sas)", "[Open code](6.1-NCI_multi_vit_a_mcmc.sas)", "[Open code](6.2-NCI_multi_vit_a_res.sas)", "[Open code](7.0-Bootstrap_variance.sas)", "[Open code](10-HEFI19-NUT_ARTICLE.qmd)", "[Open code](20-HEFI19-NUT_OSM.qmd)")) list(var = c("Name", "Description", "Link"), type = c("default", "default", "default"), column_label = list("Name", "Description", "Link"), column_align = c("left", "left", "center"), column_width = list(NULL, NULL, NULL), hidden_px = list(NULL, NULL, NULL)) list(rownum_i = 1:14, row_id = c(NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA), group_id = c(NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA), group_label = list(NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL), indent = c(NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA), built_group_label = c(NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA)) character(0) list(title = "Brief description of SAS and QMD codes", subtitle = NULL, preheader = NULL) list(vars = list(), spanner_label = list(), spanner_id = character(0), spanner_level = integer(0), gather = logical(0), built = character(0)) list(label = NULL) list(locname = "none", grpname = NA, colname = NA, locnum = 0, rownum = NA, colnum = NA, footnotes = list("CCHS-2015, Canadian Community Health Survey 2015; HEFI-2019, Healthy Eating Food Index 2019; MCMC, Markov Chain Monte Carlo; NCI, National Cancer Institute; PUMF, Public Use Microdata File"), placement = "auto") list() list(list(func = list(html = function (x) 
{
    md_to_html(x, md_engine = md_engine)
}, latex = function (x) 
{
    markdown_to_latex(x, md_engine = md_engine)
}, rtf = function (x) 
{
    markdown_to_rtf(x)
}, word = function (x) 
{
    markdown_to_xml(x)
}, default = function (x) 
{
    sub("\n$", "", vapply(x, FUN.VALUE = character(1), USE.NAMES = FALSE, commonmark::markdown_text))
}), cols = "Description", rows = 1:14, compat = NULL)) list() list(locname = character(0), grpname = character(0), colname = character(0), locnum = numeric(0), rownum = integer(0), colnum = integer(0), styles = list()) list() list(parameter = c("table_id", "table_caption", "table_width", "table_layout", "table_margin_left", "table_margin_right", "table_background_color", "table_additional_css", "table_font_names", "table_font_size", "table_font_weight", "table_font_style", "table_font_color", "table_font_color_light", "table_border_top_include", "table_border_top_style", "table_border_top_width", "table_border_top_color", "table_border_right_style", "table_border_right_width", "table_border_right_color", "table_border_bottom_include", 
"table_border_bottom_style", "table_border_bottom_width", "table_border_bottom_color", "table_border_left_style", "table_border_left_width", "table_border_left_color", "heading_background_color", "heading_align", "heading_title_font_size", "heading_title_font_weight", "heading_subtitle_font_size", "heading_subtitle_font_weight", "heading_padding", "heading_padding_horizontal", "heading_border_bottom_style", "heading_border_bottom_width", "heading_border_bottom_color", "heading_border_lr_style", "heading_border_lr_width", 
"heading_border_lr_color", "column_labels_background_color", "column_labels_font_size", "column_labels_font_weight", "column_labels_text_transform", "column_labels_padding", "column_labels_padding_horizontal", "column_labels_vlines_style", "column_labels_vlines_width", "column_labels_vlines_color", "column_labels_border_top_style", "column_labels_border_top_width", "column_labels_border_top_color", "column_labels_border_bottom_style", "column_labels_border_bottom_width", "column_labels_border_bottom_color", 
"column_labels_border_lr_style", "column_labels_border_lr_width", "column_labels_border_lr_color", "column_labels_hidden", "row_group_background_color", "row_group_font_size", "row_group_font_weight", "row_group_text_transform", "row_group_padding", "row_group_padding_horizontal", "row_group_border_top_style", "row_group_border_top_width", "row_group_border_top_color", "row_group_border_right_style", "row_group_border_right_width", "row_group_border_right_color", "row_group_border_bottom_style", 
"row_group_border_bottom_width", "row_group_border_bottom_color", "row_group_border_left_style", "row_group_border_left_width", "row_group_border_left_color", "row_group_default_label", "row_group_as_column", "table_body_hlines_style", "table_body_hlines_width", "table_body_hlines_color", "table_body_vlines_style", "table_body_vlines_width", "table_body_vlines_color", "table_body_border_top_style", "table_body_border_top_width", "table_body_border_top_color", "table_body_border_bottom_style", "table_body_border_bottom_width", 
"table_body_border_bottom_color", "data_row_padding", "data_row_padding_horizontal", "stub_background_color", "stub_font_size", "stub_font_weight", "stub_text_transform", "stub_border_style", "stub_border_width", "stub_border_color", "stub_indent_length", "stub_row_group_background_color", "stub_row_group_font_size", "stub_row_group_font_weight", "stub_row_group_text_transform", "stub_row_group_border_style", "stub_row_group_border_width", "stub_row_group_border_color", "summary_row_padding", "summary_row_padding_horizontal", 
"summary_row_background_color", "summary_row_text_transform", "summary_row_border_style", "summary_row_border_width", "summary_row_border_color", "grand_summary_row_padding", "grand_summary_row_padding_horizontal", "grand_summary_row_background_color", "grand_summary_row_text_transform", "grand_summary_row_border_style", "grand_summary_row_border_width", "grand_summary_row_border_color", "footnotes_font_size", "footnotes_padding", "footnotes_padding_horizontal", "footnotes_background_color", "footnotes_margin", 
"footnotes_border_bottom_style", "footnotes_border_bottom_width", "footnotes_border_bottom_color", "footnotes_border_lr_style", "footnotes_border_lr_width", "footnotes_border_lr_color", "footnotes_marks", "footnotes_spec_ref", "footnotes_spec_ftr", "footnotes_multiline", "footnotes_sep", "source_notes_padding", "source_notes_padding_horizontal", "source_notes_background_color", "source_notes_font_size", "source_notes_border_bottom_style", "source_notes_border_bottom_width", "source_notes_border_bottom_color", 
"source_notes_border_lr_style", "source_notes_border_lr_width", "source_notes_border_lr_color", "source_notes_multiline", "source_notes_sep", "row_striping_background_color", "row_striping_include_stub", "row_striping_include_table_body", "container_width", "container_height", "container_padding_x", "container_padding_y", "container_overflow_x", "container_overflow_y", "ihtml_active", "ihtml_use_pagination", "ihtml_use_pagination_info", "ihtml_use_sorting", "ihtml_use_search", "ihtml_use_filters", 
"ihtml_use_resizers", "ihtml_use_highlight", "ihtml_use_compact_mode", "ihtml_use_text_wrapping", "ihtml_use_page_size_select", "ihtml_page_size_default", "ihtml_page_size_values", "ihtml_pagination_type", "page_orientation", "page_numbering", "page_header_use_tbl_headings", "page_footer_use_tbl_notes", "page_width", "page_height", "page_margin_left", "page_margin_right", "page_margin_top", "page_margin_bottom", "page_header_height", "page_footer_height", "quarto_disable_processing", "quarto_use_bootstrap"
), scss = c(FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, 
TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, 
FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE), category = c("table", "table", "table", "table", "table", "table", "table", "table", "table", "table", "table", "table", "table", "table", "table", "table", "table", "table", "table", "table", "table", "table", "table", "table", "table", "table", "table", "table", "heading", "heading", "heading", "heading", "heading", 
"heading", "heading", "heading", "heading", "heading", "heading", "heading", "heading", "heading", "column_labels", "column_labels", "column_labels", "column_labels", "column_labels", "column_labels", "table_body", "table_body", "table_body", "column_labels", "column_labels", "column_labels", "column_labels", "column_labels", "column_labels", "column_labels", "column_labels", "column_labels", "column_labels", "row_group", "row_group", "row_group", "row_group", "row_group", "row_group", "row_group", 
"row_group", "row_group", "row_group", "row_group", "row_group", "row_group", "row_group", "row_group", "row_group", "row_group", "row_group", "row_group", "row_group", "table_body", "table_body", "table_body", "table_body", "table_body", "table_body", "table_body", "table_body", "table_body", "table_body", "table_body", "table_body", "data_row", "data_row", "stub", "stub", "stub", "stub", "stub", "stub", "stub", "stub", "stub", "stub", "stub", "stub", "stub", "stub", "stub", "summary_row", "summary_row", 
"summary_row", "summary_row", "summary_row", "summary_row", "summary_row", "grand_summary_row", "grand_summary_row", "grand_summary_row", "grand_summary_row", "grand_summary_row", "grand_summary_row", "grand_summary_row", "footnotes", "footnotes", "footnotes", "footnotes", "footnotes", "footnotes", "footnotes", "footnotes", "footnotes", "footnotes", "footnotes", "footnotes", "footnotes", "footnotes", "footnotes", "footnotes", "source_notes", "source_notes", "source_notes", "source_notes", "source_notes", 
"source_notes", "source_notes", "source_notes", "source_notes", "source_notes", "source_notes", "source_notes", "row", "row", "row", "container", "container", "container", "container", "container", "container", "interactive", "interactive", "interactive", "interactive", "interactive", "interactive", "interactive", "interactive", "interactive", "interactive", "interactive", "interactive", "interactive", "interactive", "page", "page", "page", "page", "page", "page", "page", "page", "page", "page", 
"page", "page", "quarto", "quarto"), type = c("value", "value", "px", "value", "px", "px", "value", "values", "values", "px", "value", "value", "value", "value", "logical", "value", "px", "value", "value", "px", "value", "logical", "value", "px", "value", "value", "px", "value", "value", "value", "px", "value", "px", "value", "px", "px", "value", "px", "value", "value", "px", "value", "value", "px", "value", "value", "px", "px", "value", "px", "value", "value", "px", "value", "value", "px", "value", 
"value", "px", "value", "logical", "value", "px", "value", "value", "px", "px", "value", "px", "value", "value", "px", "value", "value", "px", "value", "value", "px", "value", "value", "logical", "value", "px", "value", "value", "px", "value", "value", "px", "value", "value", "px", "value", "px", "px", "value", "px", "value", "value", "value", "px", "value", "px", "value", "px", "value", "value", "value", "px", "value", "px", "px", "value", "value", "value", "px", "value", "px", "px", "value", "value", 
"value", "px", "value", "px", "px", "px", "value", "px", "value", "px", "value", "value", "px", "value", "values", "values", "values", "logical", "value", "px", "px", "value", "px", "value", "px", "value", "value", "px", "value", "logical", "value", "value", "logical", "logical", "px", "px", "px", "px", "overflow", "overflow", "logical", "logical", "logical", "logical", "logical", "logical", "logical", "logical", "logical", "logical", "logical", "values", "values", "value", "value", "logical", "logical", 
"logical", "value", "value", "value", "value", "value", "value", "value", "value", "logical", "logical"), value = list(NA, NA, "auto", "fixed", "auto", "auto", "#FFFFFF", character(0), c("system-ui", "Segoe UI", "Roboto", "Helvetica", "Arial", "sans-serif", "Apple Color Emoji", "Segoe UI Emoji", "Segoe UI Symbol", "Noto Color Emoji"), "16px", "normal", "normal", "#333333", "#FFFFFF", TRUE, "solid", "2px", "#A8A8A8", "none", "2px", "#D3D3D3", TRUE, "solid", "2px", "#A8A8A8", "none", "2px", "#D3D3D3", 
    NA, "center", "125%", "initial", "85%", "initial", "4px", "5px", "solid", "2px", "#D3D3D3", "none", "1px", "#D3D3D3", NA, "100%", "normal", "inherit", "5px", "5px", "none", "1px", "#D3D3D3", "solid", "2px", "#D3D3D3", "solid", "2px", "#D3D3D3", "none", "1px", "#D3D3D3", FALSE, NA, "100%", "initial", "inherit", "8px", "5px", "solid", "2px", "#D3D3D3", "none", "1px", "#D3D3D3", "solid", "2px", "#D3D3D3", "none", "1px", "#D3D3D3", NA, FALSE, "solid", "1px", "#D3D3D3", "none", "1px", "#D3D3D3", "solid", 
    "2px", "#D3D3D3", "solid", "2px", "#D3D3D3", "8px", "5px", NA, "100%", "initial", "inherit", "solid", "2px", "#D3D3D3", "5px", NA, "100%", "initial", "inherit", "solid", "2px", "#D3D3D3", "8px", "5px", NA, "inherit", "solid", "2px", "#D3D3D3", "8px", "5px", NA, "inherit", "double", "6px", "#D3D3D3", "90%", "4px", "5px", NA, "0px", "none", "2px", "#D3D3D3", "none", "2px", "#D3D3D3", "numbers", "^i", "^i", TRUE, " ", "4px", "5px", NA, "90%", "none", "2px", "#D3D3D3", "none", "2px", "#D3D3D3", TRUE, 
    " ", "rgba(128,128,128,0.05)", FALSE, FALSE, "auto", "auto", "0px", "10px", "auto", "auto", FALSE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, 10, c(10, 25, 50, 100), "numbers", "portrait", FALSE, FALSE, FALSE, "8.5in", "11.0in", "1.0in", "1.0in", "1.0in", "1.0in", "0.5in", "0.5in", FALSE, FALSE)) list() list(locale = NULL) FALSE

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

    [1] "2023-08-01 12:04:43 EDT"

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
    [13] utf8_1.2.3        cli_3.6.1         withr_2.5.0       htmltools_0.5.5  
    [17] yaml_2.3.7        rprojroot_2.0.3   digest_0.6.33     tibble_3.2.1     
    [21] lifecycle_1.0.3   vctrs_0.6.3       glue_1.6.2        evaluate_0.21    
    [25] rmarkdown_2.23    compiler_4.2.2    pillar_1.9.0      generics_0.1.3   
    [29] jsonlite_1.8.7    pkgconfig_2.0.3  

</details>

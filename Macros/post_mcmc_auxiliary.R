# *********************************************************************** #
#                                                                         #
#          Auxiliary functions to load and prepare MCMC output            #
#                                                                         #
# *********************************************************************** #

# ********************************************** #
# Function 1: <read_mcmc>                        #
# ********************************************** #

#' Read MCMC Data and Extract Relevant Variables
#'
#' The \code{read_mcmc} function is designed to read NCI Markov Chain Monte Carlo (MCMC) data and extract relevant variables.
#'
#' @param nut_vars A character vector specifying the names of the nutrient variables that need to be extracted from the MCMC data. These variables represent the nutrients of interest.
#' @param nut_label A character string indicating the label or identifier for the nutrient being analyzed. It is used to construct file paths for reading the MCMC data.
#' @param replicnum An optional numeric parameter representing the replication number. It defaults to 0. If \code{replicnum} is 0, the function reads MCMC data from the 'Model' folder. Otherwise, it reads MCMC data from the 'Bootlib' folder.
#'
#' @return A data frame containing the relevant variables extracted from the MCMC data.
#'
#' @details The \code{nut_vars} parameter is a character vector specifying the names of the nutrient variables to be extracted. The order of variables in this vector must be consistent with the MCMC models.
#'
#' @importFrom haven read_sas
#' @importFrom dplyr rename_with mutate select all_of
#'
#' @export

read_mcmc <- function(nut_vars, nut_label, replicnum=0){

  # assign folder | value of replicnum
  if(replicnum==0) {
    outlib <- 'Model'
  } else {
    outlib <- 'Bootlib'
  }
  
  # message
  if(replicnum==0){
    message("read_mcmc: ",file.path(path, 'NCI', paste0('MCMC_', nut_label), outlib , paste0('mc_t_distrib_out', replicnum, '.sas7bdat')))
  } 
  
  # create <hefi_vars> vector
  hefi_vars <- 
    c('wg', 'pfpb', 'otherbevs', 'milk', 'rg', 'vf', 'otherfoods',
      'pfab', 'water', 'mufa', 'pufa', 'sfa', 'freesugars', 'sodium', 'energy')
  
  # note: the actual variable names in <nut_vars> do not matter
  # !!AS LONG AS!! the order is consistent with MCMC models (see log files)
  
  # combine <mc_t> variables
  all_vars <- 
    c(hefi_vars, nut_vars)
  
  # load SAS file, rename variables, calculate weight and keep relevant variables
  if(replicnum==0) tictoc::tic(msg="Step 1 - reading MCMC") # log time for replicnum 0 only
  mcmc <- 
    haven::read_sas(file.path(path, 'NCI', paste0('MCMC_', nut_label), outlib , paste0('mc_t_distrib_out', replicnum, '.sas7bdat'))) |>
    dplyr::rename_with(~ all_vars, matches("^mc_t")) |>
    dplyr::mutate(
      # divide mcmc weight
      weight_nw_sumw_div = weight_nw_sumw/500
    ) |>
    # for efficiency purpose, select only essential variables
    dplyr::select('weight_nw_sumw_div', 'wg', 'pfpb', 'otherbevs', 'milk', 'rg',
                  'vf', 'otherfoods', 'pfab', 'water', 'energy', all_of(nut_vars))
  if(replicnum==0) tictoc::toc()
  return(mcmc)
}

# ********************************************** #
# Function 2.1: <calculate_residuals>            #
# ********************************************** #

#' Calculate Energy-Adjusted Residuals
#'
#' The \code{calculate_residuals} function calculates (weighted) energy-adjusted via the residual method values for a given variable using a linear regression model with residuals as the output.
#'
#' @param indata A data frame containing the data for regression analysis.
#' @param x A character string representing the name of the variable to adjust.
#' @param energy A character string indicating the name of the energy variable used in the regression model. It defaults to "energy".
#'
#' @return A tibble containing the energy-adjusted residuals for the specified variable (\code{x}).
#'
#' @details The function fits a linear regression model with \code{x} as the dependent variable and \code{energy} as the independent variable. It adjusts the values of \code{x} based on the regression model residuals, providing energy-adjusted values.
#'
#' @importFrom stats lm residuals setNames
#' @importFrom tibble as_tibble
#'
#' @family Energy-Adjusted Residuals
#' @export

calculate_residuals <- function(indata, x, energy="energy") {
  
  # note: the weight variable 'weight_nw_sumw_div' is hardcoded
  
  if (!x %in% names(indata)) {
    stop("Dependent variable not found in the data frame.")
  }
  model <- lm(paste(x, "~", energy), weight = weight_nw_sumw_div, data = indata)
  residuals <-  tibble::as_tibble(residuals(model))
  residuals <- stats::setNames(residuals, paste0(x,"_res"))
  return(residuals)
}

# ********************************************** #
# Function 2.2: <corr_mcmc_adj>                  #
# ********************************************** #

#' Calculate Energy-Adjusted, Weighted Correlation based on MCMC Data
#'
#' The \code{corr_mcmc_adj} function performs a two-step calculation to estimate correlations using MCMC data. The function calculates the energy-adjusted, weighted correlation using the \code{cov.wt()} function. It takes the following parameters:
#'
#' @param mcmc_data A data frame containing the MCMC data.
#' @param nut_vars A character vector specifying the names of the nutrient variables of interest.
#' @param replicnum An optional numeric parameter representing the replication number. It defaults to 0.
#'
#' @details The function performs the following steps:
#' 1. Identifies a vector of HEFI (Healthy Eating Food Index) foods and adds the "_res" suffix to create a vector of all HEFI food residuals and nutrient variables.
#' 2. Loops through the HEFI foods to calculate their residuals and merge them with the input data.
#' 3. Calculates the energy-adjusted, weighted correlation using the \code{cov.wt()} function from the \pkg{stats} package.
#'
#' @return A tibble containing the energy-adjusted, weighted correlation matrix for the specified nutrient variables.
#'
#' @seealso
#' \code{\link[stats]{cov.wt}} from the \pkg{stats} package for calculating energy-adjusted, weighted correlation.
#' \code{\link[purrr]{map}} from the \pkg{purrr} package for applying functions to elements of a list.
#' \code{\link[purrr]{list_cbind}} from the \pkg{purrr} package for combining data frames element-wise.
#' \code{\link[dplyr]{cbind}} from the \pkg{dplyr} package for column-wise binding of data frames.
#'
#' @importFrom stats cov.wt
#' @importFrom purrr map list_cbind
#' @importFrom dplyr cbind
#'
#' @family MCMC Correlation
#' @export

corr_mcmc_adj <- function(mcmc_data, nut_vars, replicnum=0){
  
  # note: the weight variable 'weight_nw_sumw_div' is hardcoded
  # note: the adjustment variable 'energy' is hardcoded
  # note: the variable for HEFI <hefi_foods> are hardcoded
  
  # 1) Identification
  # A) Vector of HEFI foods
  hefi_foods <- c('vf', 'wg', 'pfab', 'pfpb', 'milk', 'water', 'otherbevs', 'rg', 'otherfoods')
  # B) Add the <_res> suffix for easy id
  hefi_foods_res <- paste0(hefi_foods,"_res")
  # C) Vector of all hefi foods residuals AND nutrient variables
  all_vars <- c(hefi_foods_res, nut_vars)
  
  ## 2) loop through to get all residuals of <hefi_foods>, merge both input data and calculated residuals
  if(replicnum==0) tictoc::tic(msg="Step 2 - adj. and correlations") # log time for replicnum 0 only
  mcmc_data <-
    mcmc_data %>% 
    purrr::map(.x     = hefi_foods, 
               .f     = calculate_residuals,
               indata = .) |>
    purrr::list_cbind() |>
    cbind(mcmc_data)
  
  ## 3) Calculate energy-adjusted, weighted correlation using <cov.wt>
  weighted_corr <- 
    stats::cov.wt(x   = mcmc_data[,all_vars],
                  wt  = mcmc_data$weight_nw_sumw_div,
                  cor = TRUE)
  corr_matrix <- 
    data.frame(replicate = replicnum,
               index     = seq.int(nrow(weighted_corr$cor)),
               name      = all_vars,
               weighted_corr$cor) |>
    tibble::as_tibble()
  if(replicnum==0) tictoc::toc()
  return(corr_matrix)
}

# ********************************************** #
# WRAPPER of 1, 2.1, 2.2: <read_n_corr>          #
# ********************************************** #

#' Read MCMC Data, Calculate Correlations, and Tidy the Output
#'
#' The \code{read_n_corr} function reads NCI Markov Chain Monte Carlo (MCMC) data, calculates correlations using the \code{corr_mcmc_adj} function, and optionally tidies the output by removing autocorrelation.
#'
#' @param nut_vars A character vector specifying the names of the nutrient variables of interest.
#' @param nut_label A character string indicating the label or identifier for the nutrient being analyzed. It is used to construct file paths for reading the MCMC data.
#' @param replicnum An optional numeric parameter representing the replication number. It defaults to 0.
#' @param tidy A logical parameter indicating whether to tidy the output by removing autocorrelation. It defaults to \code{TRUE}.
#'
#' @return A data frame containing the correlations for the specified nutrient variables.
#'
#' @details The function performs the following steps:
#' 1. Reads MCMC data using the \code{read_mcmc} function for the specified nutrient variables.
#' 2. Calculates correlations using the \code{corr_mcmc_adj} function from the read MCMC data and the nutrient variables.
#' 3. Optionally tidies the output by removing autocorrelation if \code{tidy} is \code{TRUE}.
#'
#' @seealso
#' \code{\link{read_mcmc}} for reading MCMC data.
#' \code{\link{corr_mcmc_adj}} for calculating energy-adjusted, weighted correlation using MCMC data.
#' \code{\link[dplyr]{select}} from the \pkg{dplyr} package for selecting specific variables from a data frame.
#' \code{\link[dplyr]{slice}} from the \pkg{dplyr} package for selecting specific rows from a data frame.
#'
#' @importFrom dplyr select slice
#'
#' @export


read_n_corr <- function(nut_vars, nut_label, replicnum=0, tidy=TRUE){
  output <-  
    # A) ead mcmc data 
    read_mcmc(nut_vars  = nut_vars,
              nut_label = nut_label) |>
    # B) calculate correlations
    corr_mcmc_adj(nut_vars)
  
  # C) keep only most relevant?
  if(tidy==TRUE){
    output <-
      output |>
      # keep id variables, nutrient columns
      select(1:3, all_of(nut_vars)) |>
      # remove autocorrelation 
      slice(1:(nrow(output)-length(nut_vars)))
  }
  # end of code
  return(output) 
}

#' Example Parent Trial data file used in PRBpower
#'
#' This dataset is contrived and used only to demonstrate how to use of PRBpower.
#'
#' \describe{
#'   \item{nrx}{(num) treatment code 0=standard treatment, 1=experimental(targeted treatment)}
#'   \item{survtime}{(num) Time at risk of the study event}
#'   \item{survstat}{(num) an indicator for whether the inidividual evented; 0=censored, 1=evented }
#' }
#'
#' @source These data were generated using random number generators.
#'
#' @docType data
#' @keywords datasets
#' @name parentTrial
#' @usage data(parentTrial)
#' @format ## `parentTrial`
#' A data frame with 500 rows (one row per study participant) and 5 columns.
NULL

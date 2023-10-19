
#' Pre-execution check of the user specified program parameters.
#'
#' @param ppos Proportion of subjects expected to be M+ at time=0.
#' @param HR An array of length 2 containing the biomarker hazard ratio for the
#' standard treatment group and the experimental (targeted) treatment groups.
#' @param nsim The number of simulated trials to  be generated. Set to 0 if no simulations.
#' @param alpha The critical p-value for rejecting the null hypotheses that the ratio
#' of the biomarker (or treatment) hazard ratios = 1.
#' @param dfp A data frame with a row for each individual in the parent trial.
#' This data frame must include the variables: survtime, survstat, and nrx.
#'   survtime (dbl) is the time at risk of the event of interest.
#'   Survstat (dbl) an indicator for whether the subject evented (=1) or is censored (=0).
#'   nrx (dbl) assigned treatment: (0=standard treatment, 1=experimental (targeted) treatment).
#'
#' @return NULL.
#'
#' @export
#'
#' @examples
#' data(parentTrial)
#' # Recodes to make data consistent with program requirements.
#' parentTrial$nrx  <- ifelse(parentTrial$nrx == 1, 0, 1)        # recode $rx (1,2) -> (0,1)
#' parentTrial$survstat <- ifelse((parentTrial$survstat == 1), 1, 0) # censor ind. 1=evt; 0=cens
#' paramsCheck(ppos=0.50, HR=c(1,2), nsim=1000, alpha=0.05, dfp=parentTrial)
#'
paramsCheck <- function (ppos, HR, nsim, alpha, dfp) {
  if(missing(ppos)) stop("Specify ppos: the expected proportion of biomarker positive")
  if(ppos < 0.001 | ppos > 0.999) stop("invalid value ppos: proportion of biomarker positive")
  if(missing(alpha)) stop("Specify alpha: type I error")
  if(alpha <= 0 | alpha >= 1) stop("invalid value alpha")
  if(sum(is.na(HR)) > 0) stop("Specify the biomarker HR for each of the 2 treatment groups")
  if(length(HR) != 2) stop("Specify the TWO biomarker HR; one for each treatment")
  if(missing(nsim)) stop("Specify nsim: The number of simulated trials to be performed")
  if(!exists("dfp")) stop("Specify dfp, the parent trial dataframe")
  # if(!is.data.frame(dfp)) stop("dfp should be defined as a dataframe")
  vlist <- colnames(dfp)
  if(!("survtime" %in% vlist)) stop("Parent trial dataframe must have the variable: survtime ")
  if(!("survstat" %in% vlist)) stop("Parent trial dataframe must have the variable: survstat ")
  if(!("nrx" %in% vlist)) stop("Parent trial dataframe must have the variable: nrx ")
  if (sum(ifelse(dfp$nrx %in% c(0,1), 0, 1)) > 0) stop("valid values for nrx are only (0,1)")
  if (sum(ifelse(dfp$survstat %in% c(0,1), 0, 1)) > 0)
    stop("valid values for survstat are only (0,1)")
} # end of function paramsCheck

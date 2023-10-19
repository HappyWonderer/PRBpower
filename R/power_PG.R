
#' Calculates the Peterson-George Power of a Prospective-Retrospective Biomarker Study
#'
#' @param EE The expected number of events returned from eventsTable_exp().
#' @param HR An array of length 2 specifying the hypothesized biomarker hazard ratio
#' for each Treatment Group
#' @param alpha The critical p-value boundary for rejecting the null hypothesis which
#' is that the ratio of the biomarker hazard ratios (M+:M-) for individuals in the experimental
#' treatment group and the standard treatment group equals 1.
#'
#' @return (num) The calculated power
#'
#' @export
#'
#' @importFrom stats pnorm
#' @importFrom stats qnorm
#'
#' @examples
#' library(survival)
#' library(dplyr)
#' library(stats)
#' # read data from the parent clinical trial
#' data(parentTrial)
#' dfp <- parentTrial %>% select(nrx, survtime, survstat) %>%
#' arrange(., nrx,survtime)
#' dfp$nrx  <- ifelse(dfp$nrx==1, 0, 1)          # recode df$rx (1,2) -> (0,1)
#' dfp$survstat <- ifelse((dfp$survstat == 1), 1, 0) # recode censor indicator. 1=event; 0=censor
#' dfp$rx   <- factor(dfp$nrx, levels=c(0,1), labels=c("std", "Exp"),ordered=TRUE)
#' kmdef <- Surv(dfp$survtime, dfp$survstat == 1)
#' kmfit <- survfit(kmdef ~ dfp$rx )
#' ppos  <- 0.30    # proportion of subjects expected to be M+ at the beginning of the study.
#' HR    <- c(1,2)  # biomarker hazard ratios in the M- and M+ subgroups.
#' dfm <- dissect_km(ppos=ppos, HR=HR, kmfit=kmfit)
#' dfm <- ppos_calc(dfm, ppos=ppos)
#' dfx <- join_parent(dfp=dfp, dfm=dfm)
#' evt.exp <- eventsTable_exp(dfx)
#' alpha <- 0.05
#' power_PG(EE=evt.exp, HR=HR, alpha=alpha)
#'
power_PG <- function (EE, HR, alpha) {
  var <- sum(1/EE[1:2, 2:3])
  se  <- sqrt(var)
  z   <- qnorm(1-alpha/2)  # two-sided test
  pow <- 1 - pnorm(-log(HR[2]/HR[1])/se + z) + pnorm(-log(HR[2]/HR[1])/se - z)
  return (pow)
}

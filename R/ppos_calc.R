#' Computes the probability of being biomarker-positive (M+) at each event/censor
#' time in the parent study.
#'
#' @param xx This object is returned from the function dissect_km().
#' @param ppos The proportion of subjects expected to be M+ at time=0 (baseline)
#'
#' @return A modified version of the input data frame, xx, that includes an estimate
#' of being M+ at each event/censor time in the parent study.
#'
#' @export
#'
#' @import dplyr
#'
#' @examples
#' library(survival)
#' library(dplyr)
#' # read data from the parent clinical trial
#' data(parentTrial)
#' dfp <- parentTrial %>% select(nrx, survtime, survstat) %>%
#' arrange(., nrx,survtime)
#' dfp$nrx  <- ifelse(dfp$nrx==1, 0, 1)          # recode df$rx (1,2) -> (0,1)
#' dfp$survstat <- ifelse((dfp$survstat == 1), 1, 0) # recode censor indicator. 1=event; 0=censor
#' dfp$rx   <- factor(dfp$nrx, levels=c(0,1), labels=c("std", "Exp"),ordered=TRUE)
#' kmdef <- survival::Surv(dfp$survtime, dfp$survstat == 1)
#' kmfit <- survival::survfit(kmdef ~ dfp$rx )
#' ppos  <- 0.30    # proportion of subjects expected to be M+ at the beginning of the study.
#' HR    <- c(1,2)  # biomarker hazard ratios in the M- and M+ subgroups.
#' dfm <- dissect_km(ppos=ppos, HR=HR, kmfit=kmfit)
#' ppos_calc(xx=dfm, ppos=ppos)
#'
ppos_calc <- function (xx, ppos) {
  xx <- arrange(xx, .data$nrx, .data$survival) %>%
    mutate (xx, lag_survival = dplyr::lag(.data$survival),
            lag_surv_pos = dplyr::lag(.data$surv_pos),
            survstat = NA,
            ppos = NA,
            BM   = 0)
  xx[1,"lag_survival"] <- 1  # define missing lagged value
  xx[1,"lag_surv_pos"] <- 1  # define missing lagged value
  # Need to recover the censoring indicator which survival::survfit() does not retain.
  # if xx$survival = its lagged value then this indicates these are censored individual(s).
  xx$survstat <- ifelse(xx$survival == xx$lag_survival, 0, 1)
  pp1 <- ppos * (xx$lag_surv_pos-xx$surv_pos)/(xx$lag_survival-xx$survival)
  pp0 <- ppos * (xx$surv_pos/xx$survival)
  xx$ppos <- ifelse(xx$survstat == 0, pp0, pp1)
  return(xx)
}  # end of function ppos_calc ==============================================

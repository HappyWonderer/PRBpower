
#' Join Parent Data to the Probability of Biomarker-positive Data
#'
#' Use to join Parent data frame to the Probability of Biomarker-positive data frame based
#' on: time at risk of the study event, mmor indicator, and treatment group.
#'
#' @param dfp Data frame containing the data from the parent trial.
#' @param dfm Data frame containing the probability of being biomarker-positive.
#' The probability of being marker positive is computed and returned by the
#' function ppos_calc().
#'
#' @details The probability of being biomarker-positive is a function of a person's
#' study treatment, time at risk of the event, and censor status. The probabilities of
#' being biomarker-positive is computed by the function ppos_calc() for all unique
#' combinations of these three variables in the Parent trial.  The join_parent()
#' function combines the parent-trial data to the probability of being biomarker-positive
#' based on matching each patient's trial data to the appropriate probability of
#' being biomarker-positive based on the three variables mentioned above.
#'
#' @return  A joined data frame.
#'
#' @export
#'
#' @import dplyr
#' @import survival
#'
#' @examples
#' library(survival)
#' library(dplyr)
#' # read data from the parent clinical trial
#' data(parentTrial)
#' dfp <- parentTrial %>% select(nrx, survtime, survstat) %>%
#'   arrange( .data$nrx, .data$survtime)
#' dfp$nrx  <- ifelse(dfp$nrx==1, 0, 1)          # recode df$rx (1,2) -> (0,1)
#' dfp$survstat <- ifelse((dfp$survstat == 1), 1, 0) # recode censor indicator. 1=event; 0=censor
#' dfp$rx   <- factor(dfp$nrx, levels=c(0,1), labels=c("std", "Exp"),ordered=TRUE)
#' kmdef <- Surv(dfp$survtime, dfp$survstat == 1)
#' kmfit <- survfit(kmdef ~ dfp$rx )
#' ppos  <- 0.30    # proportion of subjects expected to be M+ at the beginning of the study.
#' HR    <- c(1,2)  # biomarker hazard ratios in the M- and M+ subgroups.
#' dfm <- dissect_km(ppos=ppos, HR=HR, kmfit=kmfit)
#' dfm <- ppos_calc(dfm, ppos=ppos)
#' join_parent(dfp=dfp, dfm=dfm)
#'
#'
join_parent <- function (dfp=dfp, dfm=dfm) {
  dfx <- dfp %>%
    left_join(., dfm, by=c('survtime'='time', 'nrx'='nrx', "survstat"="survstat"))
  # define labels for summaries
  dfx$rx   <- factor(dfx$nrx,  levels=c(0,1), labels=c("Std", "Exp"), ordered=TRUE)
  dfx$BM   <- factor(dfx$BM,   levels=c(0,1), labels=c("Neg", "Pos"), ordered=TRUE)
  dfx$cens <- factor(dfx$survstat, levels=c(0,1), labels=c("Censor", "Died"), ordered=TRUE)
  return (dfx)
}

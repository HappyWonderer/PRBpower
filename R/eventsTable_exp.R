
#' Computes the Expected Number of Events by Treatment Group and Biomarker Status
#'
#' @param dfx The data frame returned by join_parent()
#'
#' @return The expected number of events by treatment group and biomarker status
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
#' dfp <- parentTrial %>% select(.data$nrx, .data$survtime, .data$survstat) %>%
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
#' dfx <- join_parent(dfp=dfp, dfm=dfm)
#' eventsTable_exp(dfx)
#'

eventsTable_exp <- function (dfx) {
  eventsTable.exp <- dfx %>%
    filter( .data$cens=="Died") %>%
    group_by( .data$rx) %>%
    summarize(
              BM.neg = sum(1 - .data$ppos),
              BM.pos = sum(.data$ppos),
              total = n())
  return(eventsTable.exp)
}

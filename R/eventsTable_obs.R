
#' Computes the Number of Evented and Censored Individuals by Treatment Group
#' in the Parent Trial
#'
#' @param dfx The data frame returned by join_parent()
#'
#' @return The number of evented and censored individuals by treatment group
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
#' eventsTable_obs(dfx)
#'

eventsTable_obs <- function (dfx) {
  eventsTable.obs <- dfx %>%
    mutate ( nevt = ifelse (.data$cens == "Died", 1, 0),
             ncen = ifelse (.data$cens == "Censor", 1, 0)) %>%
    group_by( .data$rx) %>%
    summarize(
              rx.cen = sum(.data$ncen),
              rx.evt = sum(.data$nevt),
              total = n())
  colnames(eventsTable.obs) <- c("Rx", "Censored", "Evented", "Total")
  return(eventsTable.obs)
}

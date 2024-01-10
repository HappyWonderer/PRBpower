#' Computes the probability of being biomarker-positive (M+) at each event/censor
#' time in the parent study.
#'
#' @param dfm This object is returned from the function dissect_km().
#' @param ppos The proportion of subjects expected to be M+ at time=0 (baseline)
#'
#' @return A modified version of the input data frame, dfm, that includes an estimate
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
#' ppos_calc(dfm=dfm, ppos=ppos)
#'
ppos_calc <- function (dfm, ppos) {
  dfm <- arrange(dfm, .data$nrx, .data$survival) %>%
    mutate (dfm,
            lag_survival = dplyr::lag(.data$survival),
            lag_surv_pos = dplyr::lag(.data$surv_pos),
            lag_surv_neg = dplyr::lag(.data$surv_neg),
            survstat = NA,
            ppos = NA,
            BM   = 0)
  dfm[1,"lag_survival"] <- 1  # define missing lagged value
  dfm[1,"lag_surv_pos"] <- 1  # define missing lagged value
  dfm[1,"lag_surv_neg"] <- 1  # define missing lagged value
  # Need to recover the censoring indicator which survival::survfit() does not retain.
  # if dfm$survival = its lagged value then this indicates this a censored individual(s).
  dfm$survstat <- ifelse(dfm$survival == dfm$lag_survival, 0, 1)

  lik_pos <-   ppos   * dfm$surv_pos * (dfm$lag_surv_pos - dfm$surv_pos)/dfm$lag_surv_pos
  lik_neg <- (1-ppos) * dfm$surv_neg * (dfm$lag_surv_neg - dfm$surv_neg)/dfm$lag_surv_neg

  # pp_u is prob of marker-pos among uncensored at T=t; pp_c is among censored
  pp_u <- lik_pos/(lik_pos + lik_neg)
  pp_c <- ppos * dfm$surv_pos / (ppos * dfm$surv_pos + (1-ppos) * dfm$surv_neg)
  dfm$ppos <- ifelse(dfm$survstat == 0, pp_c, pp_u)
  return(dfm)
}  # end of function ppos_calc ==============================================

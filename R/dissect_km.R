#' Estimates the Survival For Biomarker positive (M+) and -negative (M-) subjects
#'
#' @param ppos The proportion of subjects assumed to be biomarker positive at time=0.
#' @param HR An array with length=2 containing the biomarker hazard ratios under the
#' alternative hypothesis for each treatment group in the parent trial.
#' @param kmfit The output from function survival::survfit comparing the treatment
#' groups in the parent trial.
#'
#' @return A data frame containing K-M survival estimates for M+ and M- subjects
#' in each treatment group in the parent trial.
#'
#' @export
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
#' kmdef <- Surv(dfp$survtime, dfp$survstat == 1)
#' kmfit <- survfit(kmdef ~ dfp$rx )
#' dissect_km(ppos=0.30, HR=c(1,2), kmfit)
#'

dissect_km <- function (ppos, HR, kmfit) {
  ns <- kmfit$strata # get the number of distinct event times in each trtmt grp.
  if (length(ns) != 2) stop ("There should be 2 and only 2 treatment groups")

  HR <- 1/HR  # This Newton-Raphson procedure works best with HR^(-1)
  xx   <- matrix(NA, sum(ns), 5)   # pre-allocate result matrix
  xx[ , 1]   <-  rep(c(0,1), ns)   # repeat 1 and 2 ns[1:2] times each
  xx[ , 2:3] <-  cbind(unlist(kmfit$time), unlist(kmfit$surv))

  for (nrx in 0:1) {
    nt_start <- ifelse ((nrx == 0), 1, ns[1]+1)
    nt_end   <- ifelse ((nrx == 0), ns[1], sum(ns))
    nt_len   <- nt_end - nt_start + 1

    # print(cbind(nt_start, nt_end, nt_len))

    # loop over all time points in K-M
    for (i in (nt_start:nt_end)) {
      sbar = xx[i, 3]   # Kaplan-Meier estimates of S(t) for the treatment group j.
      if (sbar <= 0) { # if S(t)=0 no need to dissect; if S(t) < 0 ** error **
        xx[i, 3:5] = c(0, 0, 0) # if s(t)=0 then s_plus=0 s_minus=0
      } else {         # if sbar > 0
        s = sbar       # initial guess of the root
        inc = 1        # initialize N-R increment

        # begin Newton-Raphson iteration

        while (abs(inc) > 0.0005) {
          fx  = ppos*s + (1-ppos)*s**HR[nrx+1] - sbar  # Function whose roots are to be determined;
          fx_d= ppos + HR[nrx+1]*(1-ppos)*s**(HR[nrx+1]-1)  # 1st derivative of the function wrt s1;

          inc = fx/fx_d
          inc = ifelse (inc < s, inc, inc/2)  # if inc causes out of range vaule then use inc/2
          s=s - inc  # revised estimate of the root
           # print (cbind(sbar, fx, fx_d, inc, s ))
        }  # end of while loop

        xx[i,4] = s
        xx[i,5] = s^HR[nrx+1]
      } # end of if (sbar > 0)
    } # end of for loop over i distinct event times
  } # end loop over nrx treatment groups
  colnames(xx) <- c("nrx", "time", "survival", "surv_pos", "surv_neg")
  xx <- as.data.frame(xx)   # convert matrix x to a data frame
  return (xx)
}  # end of function dissect_km =============================================

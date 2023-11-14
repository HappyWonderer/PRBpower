
#' Compute relative hazards for the groups:
#' (rx.std, BM.neg), (rx.std, BM.pos), (rx.exp, BM.neg), (rx.exp, BM.pos)
#'
#' @param df.sim The is the data frame return by sim_gen()
#'
#' @return A matrix containing the relative hazards for each treatment group and marker status.
#' The row and columns are labeled treatment and biomarker, respectively.
#'
#' @details This function computes the mean of the coefficients over all simulated trials.
#' The relative hazards of each treatment x biomarker status is computed by summing the
#' appropriate coefficients in the proportional hazards model and exponentiating the sum.
#'
#' @export
#'
#' @examples
#' library(survival)
#' library(dplyr)
#' # read data from the parent clinical trial
#' data(parentTrial)
#' dfp <- parentTrial %>% select(nrx, survtime, survstat) %>%
#'  arrange(., nrx,survtime)
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
#' nsim <- 100
#' df.sim <- sim_gen(dfx=dfx, nsim=nsim)
#' calc_RH(df.sim)
#'
#'
calc_RH <- function(df.sim) {
  nsim <- length(df.sim)
  RH <- matrix(NA, nsim, 9)
  xx   <- matrix(c(0, 0, 0,  1, 0, 0,  1, 0, 0,
                   0, 1, 0,  1, 1, 1,  1, 0, 1,
                   0, 1, 0,  0, 1, 1,  0, 0, 1), 9, 3, byrow=TRUE)
  for (i in 1:nsim) {
    RH[i, ] <- exp(xx %*% rbind(df.sim$coef.nrx[i], df.sim$coef.BM[i], df.sim$coef.nrx.BM[i]))
  }

  RHm <- colMeans(RH)
  RHm   <- matrix(RHm, 3, 3, byrow=FALSE)
  colnames(RHm) <- c("Std", "Exp", "Exp:Std"); rownames(RHm) <- c("Neg", "Pos", "Pos:Neg")
  return (RHm)
}

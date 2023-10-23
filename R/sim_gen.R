
#' Analyze the Parent Clincal Trial Data with Simulated Biomarker Status.
#'
#' @param dfx Data Frame returned from join_parent().
#' @param nsim Number of simulated trials to be generated and analyzed.
#'
#' @return Data Frame summarizing the results from a proportional hazards model of each simulated trial.
#'
#' @details The probability that a trial subject is biomarker positive is a function of her time
#' at risk of experiencing the study event, whether she is censored or evented, and her study treatment.
#' The probability of being marker positive for all unique values of these three variables is
#' calculated in ppos_calc(). The appropriate probability is then linked to each person in the parent
#' trial with join_parent(). This function generates nsim trials in which each person's biomarker status
#' is randomly generated proportional to the expected probability of being biomarker positive.
#'
#' Each simulated trial is then analyzed with a proportional hazards (PH) model which includes covariates for
#' patient's assigned treatment, simulated biomarker status, and an interaction term for these two covariates.
#'
#' This function returns a data frame summarizing the results of each simulated trial.
#'
#' @export
#'
#' @import survival
#' @importFrom utils txtProgressBar
#' @importFrom utils setTxtProgressBar
#' @importFrom stats runif
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
#' ppos  <- 0.30    # proportion of subjects expected to be M+ at the beginning of the study.
#' HR    <- c(1,2)  # biomarker hazard ratios in the M- and M+ subgroups.
#' dfm <- dissect_km(ppos=ppos, HR=HR, kmfit=kmfit)
#' dfm <- ppos_calc(dfm, ppos=ppos)
#' dfx <- join_parent(dfp=dfp, dfm=dfm)
#' nsim <- 100
#' sim_gen(dfx=dfx, nsim=nsim)
#'
sim_gen <- function (dfx, nsim) {
  nsiz <- nrow(dfx)
  rez <- matrix(NA, nsim, 9)
  colnames(rez) <- c("coef.nrx", "coef.BM", "coef.nrx.BM",
                     "se.nrx", "se.BM", "se.nrx.BM",
                     "p.nrx", "p.BM", "p.nrx.BM")
  pbar <- ifelse (interactive(), TRUE, FALSE)
  if (pbar) pb <- txtProgressBar(style=3, min=1, max=nsim)
  print ("starting simulations")
  for (isim in 1:nsim) {
    dfr <- dfx %>%
      mutate (BM = ifelse(runif(nsiz, min=0, max=1) < dfx$ppos, 1, 0))
    res <- summary(coxph(Surv(survtime, cens=="Died") ~ nrx + BM + BM*nrx, data = dfr))
    rez[isim, 1:3] <- res$coefficients[1:3]
    rez[isim, 4:6] <- res$coefficients [7:9]
    rez[isim, 7:9] <- res$coefficients [13:15]
    if (pbar) setTxtProgressBar(pb,isim)
  }  # end of for loop isim
  print ("Simulations complete")
  rez <- as.data.frame(rez)
  return(rez)
} # end of function sim_gen

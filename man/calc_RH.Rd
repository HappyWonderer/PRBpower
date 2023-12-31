% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/calc_RH.R
\name{calc_RH}
\alias{calc_RH}
\title{Compute relative hazards for the groups:
(rx.std, BM.neg), (rx.std, BM.pos), (rx.exp, BM.neg), (rx.exp, BM.pos)}
\usage{
calc_RH(df.sim)
}
\arguments{
\item{df.sim}{The data frame return by sim_gen()}
}
\value{
A matrix containing the relative hazards for each treatment group and marker status.
The row and columns are labeled treatment and biomarker, respectively.
}
\description{
Compute relative hazards for the groups:
(rx.std, BM.neg), (rx.std, BM.pos), (rx.exp, BM.neg), (rx.exp, BM.pos)
}
\details{
This function computes the mean of the coefficients over all simulated trials.
The relative hazards of each treatment x biomarker status is computed by combining the
appropriate coefficients and then exponentiating.
}
\examples{
library(survival)
library(dplyr)
# read data from the parent clinical trial
data(parentTrial)
dfp <- parentTrial \%>\% select(nrx, survtime, survstat) \%>\%
 arrange(., nrx,survtime)
dfp$nrx  <- ifelse(dfp$nrx==1, 0, 1)          # recode df$rx (1,2) -> (0,1)
dfp$survstat <- ifelse((dfp$survstat == 1), 1, 0) # recode censor indicator. 1=event; 0=censor
dfp$rx   <- factor(dfp$nrx, levels=c(0,1), labels=c("std", "Exp"),ordered=TRUE)
kmdef <- Surv(dfp$survtime, dfp$survstat == 1)
kmfit <- survfit(kmdef ~ dfp$rx )
ppos  <- 0.30    # proportion of subjects expected to be M+ at the beginning of the study.
HR    <- c(1,2)  # biomarker hazard ratios in the M- and M+ subgroups.
dfm <- dissect_km(ppos=ppos, HR=HR, kmfit=kmfit)
dfm <- ppos_calc(dfm, ppos=ppos)
dfx <- join_parent(dfp=dfp, dfm=dfm)
nsim <- 100
df.sim <- sim_gen(dfx=dfx, nsim=nsim)
calc_RH(df.sim)


}

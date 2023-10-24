#'
#' Commonly Created Plots for Prospective-Retrospective Biomarker Study Planning.
#'
#' @param dfx The data frame returned by join_parent().
#' @param plotType 1= Survival by treatment group.  2= Survival of the standard
#' treatment group and the expected survival of biomarker-negative subgroups
#' within the standard treatment group. 3= Survival of experimental (Targeted)
#' treatment group and expected survival of the biomarker-positive and
#' biomarker-negative subgroups within this treatment group.
#' @param ppos Expected proportion of the study population that are biomarker-positive
#' at time=0.
#' @param HR an array of length=2 with the biomarker hazard ratios for the standard
#' treatment group and the experimental treatment group.
#'
#' @return A ggplot2 figure object.
#'
#' @export
#'
#' @import ggplot2
#' @importFrom dplyr filter
#' @importFrom dplyr mutate
#' @importFrom dplyr select

#'
#' @examples
#' library(survival)
#' library(dplyr)
#' # read data from the parent clinical trial
#' data(parentTrial)
#' dfp <- parentTrial %>% select( .data$nrx, .data$survtime, .data$survstat) %>%
#'   arrange(.data$nrx, .data$survtime)
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
#' plot_PRB(dfx=dfx, plotType=3, ppos=ppos, HR=HR)
#'
plot_PRB <- function (dfx, plotType, ppos, HR ) {

  # Helper sub function  ==============================
  xpose <- function (dfx, selx_nrx) {     # date transpose to tydy-ize the dataframe

    dfx <- dfx %>%                       # select one treatment group
      filter(.data$nrx == selx_nrx)

    t1 <- dfx %>%
      mutate( gtype = "Entire Treatment Group") %>%
      select (.data$nrx, .data$survtime, .data$survival, .data$gtype)

    t2 <- dfx %>%
      mutate (survival = .data$surv_neg,
              gtype = "Biomarker[-]") %>%
      select (.data$nrx, .data$survtime, .data$survival, .data$gtype)

    t3 <- dfx %>%
      mutate (survival = .data$surv_pos,
              gtype = "Biomarker[+]") %>%
      select (.data$nrx, .data$survtime, .data$survival, .data$gtype)

    return(bind_rows(t1, t2, t3))
  } # end function xpose ==========================

  if (plotType==1) {

    dfx$rx <- factor(dfx$nrx, levels=c(0,1), labels=c("Std", "Exp"), ordered=TRUE)
    plt <- ggplot (data=dfx) +
      geom_step(mapping=aes(x=.data$survtime, y=.data$survival, group=.data$rx, color=.data$rx)) +
      scale_color_manual(values = c("blue","darkred"))  +
      labs(color = "Treatment Group") +
      xlab("Time since study enrollment") +
      ylab("Proportion Event-Free") +
      ggtitle("Survival by Treatment Group") +
      theme(legend.position = c(0.8, 0.8))   +
      theme(panel.background = element_rect(fill='transparent'),
            plot.background = element_rect(fill='transparent', color=NA),
            axis.line = element_line(color = "black"))
    return(plt)

  } else if (plotType %in% c(2,3)) {
    # survival of biomarker[+] and biomarker[-] Standard treatment group
    warnMsg=""
    if (plotType == 2) {
      textMsg=paste("p+ = ", ppos, "\nBiomarker HR = ", HR[1])
      if (HR[1]==1) warnMsg="When HR=1 these curves are superimposed"
      gtitle <- "Expected Survival for Standard Treatment Group\nand Biomarker Subgroups"
      df.xpose <- xpose(dfx=dfx, selx_nrx=0)
      # This factor statement is needed to control the order of the legend and mapping.
      df.xpose$gtype = factor(df.xpose$gtype, levels=c("Entire Treatment Group",
                                                       "Biomarker[-]", 'Biomarker[+]'), ordered=TRUE)
    }  else  {    # if plotType != 2
      textMsg=paste("p+ = ", ppos, "\nBiomarker HR = ", HR[2])
      if (HR[2]==1) warnMsg="When HR=1 these curves are superimposed"
      gtitle <- "Expected Survival for Experimental Treatment Group\nand Biomarker Subgroups"
      df.xpose <- xpose(dfx=dfx, selx_nrx=1)
      df.xpose$gtype = factor(df.xpose$gtype, levels=c("Entire Treatment Group",
                                                       "Biomarker[-]", 'Biomarker[+]'), ordered=TRUE)
    }  # end of if ... else (plotType = ...)

    # note: lty and color need to be mapped to the same variable in order to control
    # the linetype and color simultaneously with the scale_*_manual statements.  Moreover
    # the scale_*_manual must specify the same legend title in order to produce one
    # combined legend.
    plt <- ggplot (data=df.xpose) +
      geom_step(mapping=aes(x=.data$survtime, y=.data$survival,
                            lty=.data$gtype, color=.data$gtype), linewidth=1) +
      xlab("Time since study enrollment") +
      ylab("Proportion Event-Free") +
      ggtitle(gtitle) +
      scale_color_manual("Group/subgroup", values=c("gray70", "darkblue", "darkblue")) +
      scale_linetype_manual("Group/subgroup", values=c("solid", "dashed", "dotdash")) +
      theme(legend.position = c(0.8, 0.8))   +
      geom_text(aes(x=105, y=0.6, label=textMsg), size=3.5, fontface="plain", hjust=0) +
      geom_text(aes(x=0, y=0.05, label=warnMsg), size=3.0, fontface="plain", hjust=0, alpha=0.50) +
      theme(panel.background = element_rect(fill='transparent'),
            plot.background = element_rect(fill='transparent', color=NA),
            legend.key = element_rect(fill = "transparent", colour = "transparent"),
            legend.key.width = unit(1,"cm"),
            axis.line = element_line(color = "black"))

    return(plt)
  }  # end of if (plotType %in% c(2,3)) ..else if
}  # end of function plot_PRB


<!-- README.md is generated from README.Rmd. Please edit that file -->

# PRBpower

## Nonparametric approaches to Calculating Power for a Prospective-Retrospective Predictive Biomarker Study with a Time to an Event Endpoint

### Introduction

Modern phase 3 clinical trials often collect and archive biologic
specimens from the study participants. These specimens may be earmarked
for identifying individuals who are more likely to benefit from
receiving the experimental treatment vs those for whom standard
treatment can stll be recommended. often the study objectives that
depend on the collected specimens may not be incorporated into the
original study design, because the specific hypothesis regarding which
biomarker or biomarkers need to be evaluated is not known or the
appropriate laboratory methods for assessing the biomarker(s) may not
yet have been determined.

The clinical trial may be completed before the hypotheses concerning the
biomarkers mature. In this case, the clinical data has been collected
and recorded. These data are fixed and not random variables. Only each
individual’s biomarker status is considered random. Since the clinical
data was collected in a prospective fashion, but the analysis of the
biomarker is performed in a retrospective fashion it has been suggested
that these studies be referred to as prospective-retrospective studies
(Simon et al,, JNCI 2009).

Since the archived specimens are limited and considered precious, the
study proposals seeking access to the specimens are frequently
rigorously reviewed for scientific merit and feasibility. The purpose of
this library is to provide the tools needed to compute the statistical
power for evaluating the predictive nature of a biomarker.

The procedures can be used to compute power as proposed by Peterson and
George (Cont. Clin. Tri. 1993) for prospective trials. Also, Power can
be computed using Monte-Carlo simulation.

### User Supplied Input Required

The user must specify the expected prevalence of biomarker-positive
individuals in the study population, the hypothesized biomarker hazard
ratios (biomarker-positive:biomarker-negative) for those who were
assigned to received the standard treatment and for those assigned to
the experimental treatment in the parent trial, and the allocated type I
error for the biomarker study.

Also required is a data file with one record per study participant with
at least 3 variables: 1) assigned study treatment, 2) time at risk of
the study event, and 3) An indicator for whether the individual
experienced the event (or is right-censored).

## Installation

Install the development version of PRBpower from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("HappyWonderer/PRBpower")
```

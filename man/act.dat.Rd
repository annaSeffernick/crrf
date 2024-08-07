% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/act.dat.R
\docType{data}
\name{act.dat}
\alias{act.dat}
\title{Pediatric Adrenocortical Tumor (ACT) data}
\format{
\subsection{\code{act.dat}}{

A data frame with 48 rows and 27 columns:
\describe{
\item{Onco_Case, Sentrix_ID}{Patient IDs}
\item{Gender}{Biological sex}
\item{Age_years}{Age at diagnosis in years}
\item{Clinical_Presentation}{Tumor type at initial diagnosis}
\item{Stage}{Cancer stage}
\item{Wieneke_Classification}{Classification of disease}
\item{Initial_Treatment}{Initial treatment for ACT}
\item{Relapse_YN}{Indicator of relapse}
\item{Mortality}{Vital status}
\item{Deceased_From_Disease}{Indicator that death was due to ACT}
\item{Follow_Up_Days}{Length of follow-up in days}
\item{Tumor_Weight_Grams}{Weight of tumor in grams}
\item{Tumor_Size_cm}{Size of tumor in centimeters.}
\item{Ki67_Percent}{Percentage of Ki-67-positive tumor cells.}
\item{P53_IHC}{P53 Immunohistochemistry}
\item{P53_IHC_Percent}{P53 Immunohistochemistry Percentage}
\item{BC_IHC}{BC Immunohistochemistry}
\item{BC_IHC_Percent}{BC Immunohistochemistry Percentage}
\item{TP53_Mutation_Category}{Mutational status of TP53}
\item{CTNNB1_Mutation_Category}{Mutational status of CTNNB1}
\item{ATRX_Mutation_Category}{Mutational status of ATRX}
\item{X11p15_status}{11p15 status, LOH means loss of heterozygosity}
\item{Methylation_Group}{Methylation-based classification developed in Clay et al. 2019 JCO Precision Oncology}
\item{RFS}{Relapse-free survival: possibly censored time until relapse or death}
\item{RFS_status}{Event type indicator for RFS: Alive, Deceased, Relapse}
}
}
}
\source{
\url{https://doi.org/10.1200/PO.19.00163}
}
\usage{
act.dat
}
\description{
Dataset from the supplementary materials of Clay et al. 2019 JCO Precision Oncology.
}
\keyword{datasets}

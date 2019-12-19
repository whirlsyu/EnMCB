#' ExpressionSet of demo dataset.
#'
#' A ExpressionSet containing the 10020 CpGs beta value of 455 samples.
#'
#' @format ExpressionSet:
#' \describe{
#'   \item{assays}{10020 features, 455 samples}
#'   \item{colData}{Data descriptions.}
#'   ...
#' }
#' @usage data(demo_set)
"demo_set"


#' Survival data of demo dataset.
#'
#' A Surv containing survival value of 455 samples.
#'
#' @format Surv data:
#' @usage data(demo_survival_data)
"demo_survival_data"

#' MCB information.
#'
#' A dataset containing the number and other attributes of 94
#' MCBs.
#'
#' @format A data frame with 94 rows and 8 variables:
#' \describe{
#'   \item{MCB_no}{MCB code}
#'   \item{start}{start point of this MCB in the chromosome.}
#'   ...
#' }
#' @usage data(demo_MCBinformation)
"demo_MCBinformation"


#' CpGs information on Illumina Infinium Human Methylation 450K.
#'
#' A dataset containing the number and other attributes of all CpG information in
#' Illumina Infinium Human Methylation 450K.
#'
#' @format A data frame with 94 rows and 8 variables:
#' \describe{
#'   \item{Composite Element REF}{CpG code}
#'   \item{Chromosome}{Chromosome code for CpG.}
#'   \item{Start}{Start point of CpGs in the chromosome.}
#'   ...
#' }
#' @usage data(Illumina_Infinium_Human_Methylation_450K)
"Illumina_Infinium_Human_Methylation_450K"





#' Expression matrix of demo dataset.
#' @description 
#' A Expression matrix containing the 10020 CpGs beta value of 455 samples in TCGA lung Adenocarcinoma dataset.
#' This will call form create_demo() function.
#' @format ExpressionSet:
#' \describe{
#'   \item{rownames}{rownames of 10020 CpG features}
#'   \item{colnames}{colnames of 455 samples}
#'   \item{realdata}{Real data matrix for demo.}
#' }
#' @usage data(demo_data)
"demo_data"


#' Survival data of demo dataset.
#' @description
#' A Surv containing survival value of 455 samples in TCGA lung Adenocarcinoma dataset.
#'
#' @format Surv data created by Surv() function in survival package.
#' This data have two unnamed arguments, they will match time and event.
#' @usage data(demo_survival_data)
"demo_survival_data"

#' MCB information.
#' @description
#' A dataset containing the number and other attributes of 94 MCBs; 
#' This results was created by the identification function IdentifyMCB.
#' This data used for metricMCB function. 
#'
#' @format A data frame with 94 rows and 8 variables:
#' \describe{
#'   \item{MCB_no}{MCB code}
#'   \item{start}{Start point of this MCB in the chromosome.}
#'   \item{end}{End point of this MCB in the chromosome.}
#'   \item{CpGs}{All the CpGs probe names in the MCB.}
#'   \item{location}{Start, end point and the chromosome number of this MCB.}
#'   \item{chromosomes}{the chromosome number of this MCB.}
#'   \item{length}{the length of bps of this MCB in the chromosome.}
#'   \item{CpGs_num}{number of CpG probes of this MCB.}
#' }
#' @usage data(demo_MCBinformation)
"demo_MCBinformation"






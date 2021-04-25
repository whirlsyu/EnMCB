#' @title univariate and multivariate survival analysis using coxph
#' @author Xin Yu
#'
#' @param dataframe Clinic data and covariates ready to be tested. Rows are variables and columns are samples.
#' @param y_surv Survival function contain survival data, usually are obtained form Surv() function in survival package.
#' @param digits Integer indicating the number of decimal places.
#' @param asnumeric indicator that the data will be (True) / not (False) transformed into numeric. Default is true.
#' @export
#'
#' @examples
#' data(demo_survival_data)
#' data('demo_data',package = "EnMCB")
#' demo_set<-demo_data$realdata
#' res<-univ_coxph(demo_set,demo_survival_data)
#' @return Object of class \code{matrix} with results.
#' 
univ_coxph <- function(dataframe,y_surv,digits=4,asnumeric=TRUE) {

  #
  rownames_dataframe<-rownames(dataframe)
  res_coxph<-NULL
  for (i in seq_len(nrow(dataframe))) {
    if (asnumeric) {
      covariates <- as.numeric(dataframe[i,])
    }else{
      covariates <- unlist(dataframe[i,])
    }

    univ_models<-survival::coxph(y_surv ~ covariates)
    univ_models<-summary(univ_models)
    p.value<-signif(univ_models$wald["pvalue"], digits)
    wald.test<-signif(univ_models$wald["test"], digits)
    beta<-signif(univ_models$coef[1], digits);#coeficient beta
    HR <-signif(univ_models$coef[2], digits);#exp(beta)
    HR.confint.lower <- signif(univ_models$conf.int[,"lower .95"],digits)
    HR.confint.upper <- signif(univ_models$conf.int[,"upper .95"],digits)
    HR <- paste0(HR, " (",
                 HR.confint.lower, "-", HR.confint.upper, ")")
    res<-c(i,beta, HR, wald.test, p.value)
    names(res)<-c("label","beta", "HR (95% CI for HR)", "wald.test",
                  "p.value")
    res_coxph<-rbind(res_coxph,res)
  }
  colnames(res_coxph)<-c("label","beta", "HR (95% CI for HR)", "wald.test",
                         "p.value")
  res_coxph[,'label']<-rownames_dataframe
  # Source:Benjamini, Y., and Hochberg, Y. (1995).
  # Controlling the false discovery rate: a practical and powerful approach
  # to multiple testing. Journal of the Royal Statistical Society Series B,
  # 57, 289–300. http://www.jstor.org/stable/2346101.
  BH_adjust_p=stats::p.adjust(as.numeric(res_coxph[,"p.value"]),method = "BH")
  # final results
  res_coxph<-cbind(res_coxph,BH_adjust_p)
  res_coxph
}

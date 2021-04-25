#' @title multivariate survival analysis using coxph
#' @author Xin Yu
#'
#' @param dataframe Clinic data and covariates ready to be tested. Note that Rows are samples and columns are variables.
#' @param y_surv Survival function contain survival data, usually are obtained form Surv() function in survival package.
#' @param digits Integer indicating the number of decimal places.
#' @param asnumeric indicator that the data will be (True) / not (False) transformed into numeric. Default is true.
#' @export
#'
#' @examples
#' data(demo_survival_data)
#' data('demo_data',package = "EnMCB")
#' demo_set<-demo_data$realdata
#' res<-multi_coxph(t(demo_set),demo_survival_data)
#' @return Object of class \code{matrix} with results.
#' 
multi_coxph <- function(dataframe,y_surv,digits=4,asnumeric=TRUE) {
  dataframe <- as.data.frame(dataframe)
  if (asnumeric) {
    for (i in seq(ncol(dataframe))) {
      if (!is.numeric(dataframe[1,i]))  dataframe[,i] <-as.numeric(as.factor(dataframe[,i] ))
    }
  }
  #print(head(dataframe))
  covariates <- paste(colnames(dataframe),collapse = ' + ')
  formula_str<- paste('y_surv',covariates,sep = ' ~ ')
  #cat(formula_str)
  multi_models<-survival::coxph(as.formula(formula_str),data = dataframe)
  multi_models %>% finalfit::fit2df()
}
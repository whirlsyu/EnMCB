#' @title Trainging stacking ensemble model for Methylation Correlation Block
#'
#' @description Method for training a stacking ensemble model for Methylation Correlation Block.
#' @export
#' @param single_res Methylation Correlation Block information returned by the IndentifyMCB function.
#' @param training_set methylation matrix used for training the model in the analysis.
#' @param Surv_training Survival function contain the survival information for training.
#' @param testing_set methylation matrix used for testing the model in the analysis.
#' @param Surv_testing Survival function contain the survival information for testing.
#' @author Xin Yu
#' @keywords methylation ensemble stacking
#' @usage ensemble_model(single_res,training_set,Surv_training,testing_set,Surv_testing)
#' @return Object of class \code{list} with elements (XXX repesents the model you choose):
#'  \tabular{ll}{
#'    \code{cox} \tab Model object for the cox model at first level. \cr
#'    \code{svm} \tab Model object for the svm model at first level. \cr
#'    \code{lasso} \tab Model object for the lasso model at first level. \cr
#'    \code{stacking} \tab Model object for the stacking model. \cr
#'  }
#' @references
#'  Xin Yu et al. 2019 Predicting disease progression in lung adenocarcinoma patients based on methylation correlated blocks using ensemble machine learning classifiers (under review)
#'
ensemble_model <- function(single_res, training_set, Surv_training, testing_set=NULL, Surv_testing=NULL) {
  if (dim(single_res)[1]>dim(single_res)[2]) {
    single_res<-t(as.matrix(single_res))
  }
  if (nrow(single_res)!=1){
    stop(errorCondition("ERROR: Results information (single_res) should only contain one single Methylation Correlation Block."))
  }
  all_related_CpGs<-strsplit(paste(single_res[,"CpGs"],collapse = " ")," ")[[1]]
  related_training<-training_set[rownames(training_set) %in% all_related_CpGs,]
  if (!is.null(testing_set)) {
    related_testing<-testing_set[rownames(testing_set) %in% all_related_CpGs,]
  }else{
    related_testing<-NULL
  }
  cox <- EnMCB::metricMCB(MCBset = single_res,
                          training_set = related_training,
                          Surv = Surv_training,
                          testing_set = related_testing,
                          Surv.new = Surv_testing,
                          Method = "cox",
                          silent = T)
  svm<- EnMCB::metricMCB(MCBset = single_res,
                         training_set = related_training,
                         Surv = Surv_training,
                         testing_set = related_testing,
                         Surv.new = Surv_testing,
                         Method = "svm",
                         silent = T)
  lasso<- EnMCB::metricMCB(MCBset = single_res,
                           training_set = related_training,
                           Surv = Surv_training,
                           testing_set = related_testing,
                           Surv.new = Surv_testing,
                           Method = "lasso",
                           silent = T)
  data<-rbind(cox$MCB_cox_matrix_training,
              svm$MCB_svm_matrix_training,
              lasso$MCB_lasso_matrix_training
  )
  rownames(data)<-c('cox','svm','lasso')
  data<-t(data)
  data_f<-as.data.frame(data)

  try(univ_models<-rms::cph(formula = Surv_training ~. ,data=data_f) )
  if (is.null(univ_models)) {
    stop(errorCondition("Ensemble model can't be created, please check your data..."))
  }else{
    res<-list(cox=cox$best_cox_model,
              svm$best_svm_model,
              lasso$best_lasso_model,
              univ_models)
    names(res)<-c("cox","svm","lasso","stacking")
    return(res)
  }
}

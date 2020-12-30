#' @title fitting function using stacking ensemble model for Methylation Correlation Block
#'
#' @description predict is a generic function for predictions from the results of stacking ensemble model fitting functions.
#' The function invokes particular methods which is the ensemble model described in the reference.
#' @param ensemble_model ensemble model which built by ensemble_model() function
#' @param prediction_data A vector, matrix, list, or data frame containing the predictions (input).
#' @param mutiple_results Boolean vector, True for including the single model results.
#' @references
#' Xin Yu et al. 2019 Predicting disease progression in lung adenocarcinoma patients based on methylation correlated blocks using ensemble machine learning classifiers (under review)
#' @export
#' @return Object of numeric class \code{double}
#' @examples 
#' library(survival)
#' #import datasets
#' data(demo_survival_data)
#' datamatrix<-create_demo()
#' data(demo_MCBinformation)
#' #select MCB with at least 3 CpGs.
#' demo_MCBinformation<-demo_MCBinformation[demo_MCBinformation[,"CpGs_num"]>2,]
#' trainingset<-colnames(datamatrix) %in% sample(colnames(datamatrix),0.6*length(colnames(datamatrix)))
#' testingset<-!trainingset
#' #select one MCB
#' select_single_one=1
#' em<-ensemble_model(t(demo_MCBinformation[select_single_one,]),
#'     training_set=datamatrix[,trainingset],
#'     Surv_training=demo_survival_data[trainingset])
#'
#' em_prediction_results<-ensemble_prediction(ensemble_model = em,
#' prediction_data = datamatrix[,testingset])
#'
ensemble_prediction <- function(ensemble_model,prediction_data, mutiple_results = FALSE) {
  if (mutiple_results) {
    return(ensemble_prediction.m(ensemble_model,prediction_data))
  }
  prediction_data<-prediction_data[rownames(prediction_data) %in% ensemble_model$cox$cox_model$CpGs,]
  if (nrow(prediction_data)!=length(rownames(prediction_data))) {
    stop("ERROR: The predition data and the model have wrong dimensions!")
  }
  svm<-stats::predict(ensemble_model$svm$svm_model, data.frame(t(prediction_data)))$predicted
  cox<-stats::predict(ensemble_model$cox$cox_model, data.frame(t(prediction_data)))
  enet<-stats::predict(ensemble_model$enet$enet_model,t(prediction_data),s=ensemble_model$enet$`corrected_lambda(min)`)
  coxboost<-stats::predict(ensemble_model$coxboost$coxboost_model, t(prediction_data))[,1]
  data<-rbind(cox,
              svm,
              t(enet),
              coxboost
  )
  rownames(data)<-c('cox','svm','enet','coxboost')
  data<-t(data)
  data_f<-as.data.frame(data)
  stats::predict(ensemble_model$stacking, data_f)
}

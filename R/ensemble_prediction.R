#' @title fitting function using stacking ensemble model for Methylation Correlation Block
#'
#' @description predict is a generic function for predictions from the results of stacking ensemble model fitting functions.
#' The function invokes particular methods which is the ensemble model described in the reference.
#' @param ensemble_model ensemble model which built by ensemble_model() function
#' @param predition_data A vector, matrix, list, or data frame containing the predictions (input).
#' @references
#' Xin Yu et al. 2019 Predicting disease progression in lung adenocarcinoma patients based on methylation correlated blocks using ensemble machine learning classifiers (under review)
#' @export
#' @return Object of numeric class \code{double}
#' @examples 
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
#'     Surv_training=datamatrix[trainingset])
#'
#' em_prediction_results<-ensemble_prediction(ensemble_model = em,
#' predition_data = datamatrix[,testingset])
#'
ensemble_prediction <- function(ensemble_model,predition_data) {
  predition_data<-predition_data[rownames(predition_data) %in% names(coef(ensemble_model$cox$cox_model)),]
  if (nrow(predition_data)!=length(rownames(predition_data))) {
    stop("ERROR: The predition data and the model have wrong dimensions!")
  }
  svm<-stats::predict(ensemble_model$svm$svm_model, data.frame(t(predition_data)))$predicted
  cox<-stats::predict(ensemble_model$cox$cox_model, data.frame(t(predition_data)))
  lasso<-stats::predict(ensemble_model$lasso$`lasso model`,t(predition_data),s=ensemble_model$lasso$`corrected lambda(min)`)
  data<-rbind(cox,
              svm,
              t(lasso)
  )
  rownames(data)<-c('cox','svm','lasso')
  data<-t(data)
  data_f<-as.data.frame(data)
  stats::predict(ensemble_model$stacking, data_f)
}

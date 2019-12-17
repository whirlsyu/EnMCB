#' @title Fast calculation of AUC for ROC using parallel strategy
#'
#' @description This function is used to create time-dependent ROC curve from censored survival data using the Kaplan-Meier (KM) or Nearest Neighbor Estimation (NNE) method of Heagerty, Lumley and Pepe, 2000
#'
#'
#' @param test_matrix Test matrix used in the analysis. Colmuns are samples, rows are markers.
#' @param y_surv Survival information created by Surv function in survival package.
#' @param predict_time Time point of the ROC curve.
#' @param roc_method Method for fitting joint distribution of (marker,t), either of KM or NNE, the default method is NNE.
#'
#' @author Xin Yu
#' @return
#' This will retrun a numeric vector contains AUC results for each row in test_matrix.
#' @export
#'
#'
fast_roc_calculation<- function(test_matrix,y_surv,predict_time=5,roc_method="NNE") {
  requireNamespace(survivalROC)
  requireNamespace(parallel)
  no_cores <- parallel::detectCores()
  cl <- parallel::makeCluster(getOption("cl.cores", no_cores));
  survival_calculation<-function(x,y_surv,predict_time){
    requireNamespace(survivalROC)
    survivalROC::survivalROC(Stime=y_surv[,1],
                             status=y_surv[,2],
                             marker =x,
                             lambda = NULL,
                             predict.time = predict_time,method = "NNE",span =0.25*length(y_surv[,1])^(-0.20))$AUC}
  res <- parallel::parApply(cl, test_matrix,1, survival_calculation,y_surv,predict_time)
  parallel::stopCluster(cl)
  res
}

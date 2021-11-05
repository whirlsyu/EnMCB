#' @title Trainging stacking ensemble model for Methylation Correlation Block
#'
#' @description Method for training a stacking ensemble model for Methylation Correlation Block.
#' @export
#' @param single_res Methylation Correlation Block information returned by the IndentifyMCB function.
#' @param training_set methylation matrix used for training the model in the analysis.
#' @param Surv_training Survival function contain the survival information for training.
#' @param testing_set methylation matrix used for testing the model in the analysis.
#' @param Surv_testing Survival function contain the survival information for testing.
#' @param ensemble_type Secondary model use for ensemble, one of Cox and C-index.
#' @author Xin Yu
#' @keywords methylation ensemble stacking
#' @usage ensemble_model(single_res,training_set,Surv_training,testing_set,
#' Surv_testing,ensemble_type)
#' @return Object of class \code{list} with elements (XXX repesents the model you choose):
#'  \tabular{ll}{
#'    \code{cox} \tab Model object for the cox model at first level. \cr
#'    \code{svm} \tab Model object for the svm model at first level. \cr
#'    \code{enet} \tab Model object for the enet model at first level. \cr
#'    \code{mboost} \tab Model object for the mboost model at first level. \cr
#'    \code{stacking} \tab Model object for the stacking model. \cr
#'  }
#' @references
#'  Xin Yu et al. 2019 Predicting disease progression in lung adenocarcinoma patients based on methylation correlated blocks using ensemble machine learning classifiers (under review)
#' @examples
#' #import datasets
#' library(survival)
#' data(demo_survival_data)
#' datamatrix<-create_demo()
#' data(demo_MCBinformation)
#' #select MCB with at least 3 CpGs.
#' demo_MCBinformation<-demo_MCBinformation[demo_MCBinformation[,"CpGs_num"]>2,]
#' trainingset<-colnames(datamatrix) %in% sample(colnames(datamatrix),0.6*length(colnames(datamatrix)))
#' select_single_one=1
#' em<-ensemble_model(t(demo_MCBinformation[select_single_one,]),
#'     training_set=datamatrix[,trainingset],
#'     Surv_training=demo_survival_data[trainingset])
#'
#'
ensemble_model <- function(single_res, 
                           training_set, 
                           Surv_training, 
                           testing_set=NULL, 
                           Surv_testing=NULL,
                           ensemble_type = "Standard Cox") {
  if (dim(single_res)[1]>dim(single_res)[2]) {
    single_res<-t(as.matrix(single_res))
  }
  if (nrow(single_res)!=1){
    stop(errorCondition("ERROR: Results information (single_res) should only contain one single Methylation Correlation Block."))
  }
  rz=!(is.na(Surv_training)|Surv_training[,1]==0)
  Surv_training<-Surv_training[rz]
  all_related_CpGs<-strsplit(paste(single_res[,"CpGs"],collapse = " ")," ")[[1]]
  related_training<-training_set[rownames(training_set) %in% all_related_CpGs,rz]
  if (!is.null(testing_set)) {
    related_testing<-testing_set[rownames(testing_set) %in% all_related_CpGs,]
  }else{
    related_testing<-NULL
  }
  cox <- tryCatch(EnMCB::metricMCB(MCBset = single_res,
                          training_set = related_training,
                          Surv = Surv_training,
                          testing_set = related_testing,
                          Surv.new = Surv_testing,
                          Method = "cox",
                          silent = TRUE),error = function(e){NULL})
  svm<- tryCatch(EnMCB::metricMCB(MCBset = single_res,
                         training_set = related_training,
                         Surv = Surv_training,
                         testing_set = related_testing,
                         Surv.new = Surv_testing,
                         Method = "svm",
                         silent = TRUE),error = function(e){NULL})
  enet<- tryCatch(EnMCB::metricMCB(MCBset = single_res,
                           training_set = related_training,
                           Surv = Surv_training,
                           testing_set = related_testing,
                           Surv.new = Surv_testing,
                           Method = "enet",
                           silent = TRUE),error = function(e){NULL})
  mboost<- tryCatch(EnMCB::metricMCB(MCBset = single_res,
                                   training_set = related_training,
                                   Surv = Surv_training,
                                   testing_set = related_testing,
                                   Surv.new = Surv_testing,
                                   Method = "mboost",
                                   silent = TRUE),error = function(e){NULL})
  data<-rbind(cox$MCB_cox_matrix_training,
              svm$MCB_svm_matrix_training,
              enet$MCB_enet_matrix_training,
              mboost$MCB_mboost_matrix_training
  )
  rownames(data)<-c('cox','svm','enet','mboost')
  data<-t(data)
  data_f<-as.data.frame(data)
  if (ensemble_type == "Standard Cox"){
    univ_models<-tryCatch(rms::cph(formula = Surv_training ~ cox + svm + enet + mboost ,data=data_f),error=function(e){NULL} )
  }else if (ensemble_type == "C-index"){
    #This will let the prediction function use the lambda.1st as the default lambda, which maximum the C-index.
    #It returns the linear model predictions.
    univ_models<-tryCatch(glmnet::cv.glmnet(x=as.matrix(data_f[!is.na(Surv_training),]),
                                            Surv_training[!is.na(Surv_training)],
                                            family='cox',type.measure = 'C'),
                          error=function(e){NULL} )
  }else{
    warning("unsupported ensemble type! must be one of Cox and C-index")
    stop(0)
  }
  if (is.null(univ_models)) {
    stop(errorCondition("Ensemble model can't be created, please check your data..."))
  }else{
    res<-list(cox$best_cox_model,
              svm$best_svm_model,
              enet$best_enet_model,
              mboost$best_mboost_model,
              univ_models)
    names(res)<-c("cox","svm","enet","mboost","stacking")
    return(res)
  }
}

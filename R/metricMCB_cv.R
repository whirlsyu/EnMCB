#' @title Calculation of model AUC for Methylation Correlation Blocks using cross validation
#'
#' @description To enable quantitative analysis of the methylation patterns
#' within individual Methylation Correlation Blocks across many samples, a single metric to
#' define the methylated pattern of multiple CpG sites within each block.
#' Compound scores which calculated all CpGs within individual Methylation Correlation Blocks by SVM model
#' were used as the compound methylation values of Methylation Correlation Blocks.
#' @usage metricMCB.cv(MCBset,data_set,Surv,nfold,
#' Method,predict_time,alpha,n_mstop,n_nu,theta,silent)
#' @export
#' @param MCBset Methylation Correlation Block information returned by the IndentifyMCB function.
#' @param data_set methylation matrix used for training the model in the analysis.
#' @param Surv Survival function contain the survival information for training.
#' @param nfold fold used in the cross validation precedure.
#' @param Method model used to calculate the compound values for multiple Methylation correlation blocks. Options include "svm", "cox", "coxboost", and "enet". The default option is SVM method.
#' @param predict_time time point of the ROC curve used in the AUC calculations, default is 5 years.
#' @param alpha The elasticnet mixing parameter, with 0 ≤ alpha ≤ 1. alpha=1 is the lasso penalty, and alpha=0 the ridge penalty. It works only when "enet" Method is selected.
#' @param n_mstop an integer giving the number of initial boosting iterations. If mstop = 0, the offset model is returned. It works only when "coxboost" Method is selected.
#' @param n_nu a double (between 0 and 1) defining the step size or shrinkage parameter in coxboost model. It works only when "coxboost" Method is selected.
#' @param theta penalty used in the penalized coxph model, which is theta/2 time sum of squared coefficients. default is 1. It works only when "cox" Method is selected.
#' @param silent Ture indicates that processing information and progress bar will be shown.
#' @author Xin Yu
#' @keywords Methylation Correlation
#' @examples
#' #import datasets
#' data(demo_survival_data)
#' datamatrix<-create_demo()
#' data(demo_MCBinformation)
#' #select MCB with at least 3 CpGs.
#' demo_MCBinformation<-demo_MCBinformation[demo_MCBinformation[,"CpGs_num"]>2,]
#'
#' trainingset<-colnames(datamatrix) %in% sample(colnames(datamatrix),0.6*length(colnames(datamatrix)))
#' testingset<-!trainingset
#' #create the results using Cox regression. 
#' mcb_cox_res<-metricMCB.cv(MCBset = demo_MCBinformation,
#'                data_set = datamatrix,
#'                Surv = demo_survival_data,
#'                Method = "cox")
#'
#' @return Object of class \code{list} with elements (XXX will be replaced with the model name you choose):
#'  \tabular{ll}{
#'    \code{MCB_matrix} \tab Prediction results of model. \cr
#'    \code{auc_results} \tab AUC results for each model. \cr
#'  }
#' @references
#' Xin Yu et al. 2019 Predicting disease progression in lung adenocarcinoma patients based on methylation correlated blocks using ensemble machine learning classifiers (under review)
#' @importFrom mboost glmboost predict.glmboost CoxPH boost_control
#' 
metricMCB.cv<-function(
  MCBset,
  data_set,
  Surv,
  nfold=10,
  Method=c("svm","cox","enet","coxboost")[1],
  predict_time = 5,
  alpha = 0.5,
  n_mstop = 500,
  n_nu = 0.1,
  theta = 1,
  silent=FALSE
){
  requireNamespace("stats")
  if (!silent) {
    cat("Start anaylsis, this may take a while...\n")
    show_bar=nrow(MCBset)>1
  }else{
    show_bar=FALSE
  }
  if (show_bar) {
    bar<-utils::txtProgressBar(min = 1,max = nrow(MCBset),char = "#",style = 3)
  }
  if (is.null(Surv)) {
    stop(paste("You must have a survival function to train the data."))
  }
  if (is.integer0(grep("MCB_no|CpGs",colnames(MCBset)))){
    stop(paste("Methylation Correlation Block information in your result must have columns of MCB_no and CpGs. Please check your results."))
  }
  na_or_zero_data<-(is.na(Surv)|Surv[,1]==0)
  if (sum(na_or_zero_data)>0){
    data_set<-data_set[,!na_or_zero_data]
    Surv<-Surv[!na_or_zero_data]
    warning("survival data contains NAs or zero survival times, NAs or data with zero survival times are remove automaticly.")
  }
  # constuction of MCB Method matrix for SVM
  MCB_matrix<-matrix(0,nrow = nrow(MCBset),ncol = ncol(data_set))
  colnames(MCB_matrix)<-colnames(data_set)
  rownames(MCB_matrix)<-as.numeric(MCBset[,'MCB_no'])
  if (!Method %in% c("svm","cox","enet","coxboost")){
    stop(paste("Method:",Method,"is not supported, see hlep files for the details.",collapse = " "))
  }
  sp<-sample(1:ncol(data_set),replace = F)
  order_sp<-order(sp)
  data_set<-data_set[,sp]
  folds <- cut(seq(1,ncol(data_set)),breaks=nfold,labels=FALSE)
  #if it has a independent test set create the test_set res set
  FunctionResults<-NULL
  best_auc<-0
  best_model<-NULL
  mcb_model_res<-NULL
  for (mcb in seq_len(nrow(MCBset))) {
    #if (nrow(MCBset)>1){}
    if (show_bar&!silent) {
      utils::setTxtProgressBar(bar, mcb)
    }
    write_MCB<-rep(NA,5)
    #save the mcb number
    write_MCB[1]<-as.numeric(MCBset[mcb,'MCB_no'])
    write_MCB[2]<-MCBset[mcb,'CpGs']
    # build temp variable for saving the results.
    # MCB number
    # aquire information for CpG sites in MCB
    CpGs<-strsplit(MCBset[mcb,'CpGs']," ")[[1]]
    #cat(CpGs)
    for (i in seq(unique(folds))) {
      rz<- which(folds==i,arr.ind=TRUE)
      data_used_for_training<-data.frame(t(data_set[CpGs,-rz]))
      data_used_for_testing <-data.frame(t(data_set[CpGs,rz]))
      # train a svm model
      times = Surv[-rz]
      if (Method=="svm") {
        model<-tryCatch(survivalsvm::survivalsvm(times ~ .,
                                                 data_used_for_training,
                                                 gamma.mu = 0.1,
                                                 type = "regression"),
                        error = function(e){warning(paste('SVR can not be built, error occurs:', e));return(NULL)})
      }else if(Method=="cox"){
        data_used_for_training = data.frame(allvars = as.ridgemat(data_used_for_training))
        if (length(CpGs)<20){
          model<-tryCatch(survival::coxph(times ~ allvars, 
                                data_used_for_training),
                error = function(e){return(NULL)})
        }else{
          model<-tryCatch(ridge_model(times, data_used_for_training, theta),
                          error = function(e){return(NULL)})
        }
        if (!is.null(model)){
          model$CpGs <- CpGs
          model <- as.mcb.coxph.penal(model)
        }
      }else if(Method=="enet"){
        model<-tryCatch( glmnet::cv.glmnet(data_used_for_training,
                                           times,
                                           #cox model in enet was used, note that here penalty were used.
                                           family="cox",
                                           alpha=alpha,
                                           # The elasticnet mixing parameter, with 0≤α≤ 1. The penalty is defined as
                                           # (1-alpha)/2||beta||_2^2+alpha||beta||_1
                                           # alpha=1 is the lasso penalty, and alpha=0 the ridge penalty.
                                           # type.measure = "AUC"
                                           type.measure= "deviance"
                                           # It uses AUC as the criterion for 10-fold cross-validation.
                                           #foldid = 10
        ),error = function(e){return(NULL)})
        correctional_value=1
        while ( sum(stats::coef(model, s = model$lambda.min-0.001*(correctional_value-1))>0)<1 &
                (model$lambda.min-0.001*(correctional_value-1))>0 ) {
          correctional_value=correctional_value*1.25
        }
        lambda_min_corrected<-model$lambda.min-0.001*(correctional_value-1)
      }else if (Method=="coxboost"){
        model <- tryCatch(mboost::glmboost(y=times,x=data_used_for_training,family=mboost::CoxPH(),
                                                    control=mboost::boost_control(mstop=n_mstop,nu=n_nu)),error = NULL)
      }
      if (!is.null(model)) {
        #predictions
        if (Method=="svm") MCB_matrix[mcb,rz]<-stats::predict(model,data_used_for_testing)$predicted
        if (Method=="cox") MCB_matrix[mcb,rz]<-stats::predict(model,data_used_for_testing)
        if (Method=="enet"){
          #if you use lambda.1se instead, the penalty of enet would be larger, leading that most of covariates were removed form the final model.
          MCB_matrix[mcb,rz]<-stats::predict(model,data_used_for_testing,s=lambda_min_corrected)
        }
        if (Method=="coxboost")MCB_matrix[mcb,rz] <- stats::predict(model, data_used_for_testing)[,1]
      }else{
        MCB_matrix[mcb,rz]<-NA
      }
    }
    MCB_matrix[mcb,]<-MCB_matrix[mcb,order_sp]
    if (sum(is.na(MCB_matrix[mcb,])) == 0){
      AUC_value<-survivalROC::survivalROC(Stime = Surv[,1],
                                          status = Surv[,2],
                                          marker = MCB_matrix[mcb,],
                                          predict.time = predict_time,
                                          method = "NNE",
                                          span =0.25*length(Surv)^(-0.20))$AUC
      cindex<-survival::survConcordance(Surv ~ MCB_matrix[mcb,])
      write_MCB[4]<-cindex$concordance
      write_MCB[5]<-cindex$std.err
    }else{
      write_MCB[3:5]<-NA
    }
    mcb_model_res<-rbind(mcb_model_res,write_MCB)
  }
  colnames(mcb_model_res)<-c("MCB_no","CpGs","auc","C-index","C-index_SE")
  FunctionResults$MCB_matrix<-MCB_matrix
  FunctionResults$auc_results<-mcb_model_res
  return(FunctionResults)
}
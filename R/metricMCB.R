#' @title Calculation of the metric matrix for Methylation Correlation Block
#'
#' @description To enable quantitative analysis of the methylation patterns
#' within individual Methylation Correlation Blocks across many samples, a single metric to
#' define the methylated pattern of multiple CpG sites within each block.
#' Compound scores which calculated all CpGs within individual Methylation Correlation Blocks by linear, SVM or elastic-net model
#' Predict values were used as the compound methylation values of Methylation Correlation Blocks.
#' @usage metricMCB(MCBset,training_set,Surv,testing_set,
#' Surv.new,Method,predict_time,ci,silent,alpha,n_mstop,n_nu,theta)
#' @export
#' @param training_set methylation matrix used for training the model in the analysis.
#' @param testing_set methylation matrix used in the analysis. This can be missing then training set itself will be used as testing set.
#' @param MCBset Methylation Correlation Block information returned by the IndentifyMCB function.
#' @param Surv Survival function contain the survival information for training.
#' @param Surv.new Survival function contain the survival information for testing.
#' @param Method model used to calculate the compound values for multiple Methylation correlation blocks. Options include "svm" "cox" "coxboost" and "enet". The default option is SVM method.
#' @param predict_time time point of the ROC curve used in the AUC calculations, default is 5 years.
#' @param ci if True, the confidence intervals for AUC under area under the receiver operating characteristic curve will be calculated. This will be time consuming. default is False.
#' @param silent True indicates that processing information and progress bar will be shown.
#' @param alpha The elasticnet mixing parameter, with 0 ≤ alpha ≤ 1. alpha=1 is the lasso penalty, and alpha=0 the ridge penalty. It works only when "enet" Method is selected.
#' @param n_mstop an integer giving the number of initial boosting iterations. If mstop = 0, the offset model is returned. It works only when "coxboost" Method is selected.
#' @param n_nu a double (between 0 and 1) defining the step size or shrinkage parameter in coxboost model. It works only when "coxboost" Method is selected.
#' @param theta penalty used in the penalized coxph model, which is theta/2 time sum of squared coefficients. default is 1. It works only when "cox" Method is selected.
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
#' mcb_cox_res<-metricMCB(MCBset = demo_MCBinformation,
#'                training_set = datamatrix[,trainingset],
#'                Surv = demo_survival_data[trainingset],
#'                testing_set = datamatrix[,testingset],
#'                Surv.new = demo_survival_data[testingset],
#'                Method = "cox"
#'                )
#'
#' @return Object of class \code{list} with elements (XXX will be replaced with the model name you choose):
#'  \tabular{ll}{
#'    \code{MCB_XXX_matrix_training} \tab Prediction results of model for training set. \cr
#'    \code{MCB_XXX_matrix_test_set} \tab Prediction results of model for test set. \cr
#'    \code{XXX_auc_results} \tab AUC results for each model. \cr
#'    \code{best_XXX_model} \tab Model object for the model with best AUC. \cr
#'    \code{maximum_auc} \tab Maximum AUC for the whole generated models. \cr
#'  }
#' @references
#' Xin Yu et al. 2019 Predicting disease progression in lung adenocarcinoma patients based on methylation correlated blocks using ensemble machine learning classifiers (under review)
#' @importFrom mboost glmboost predict.glmboost CoxPH boost_control
metricMCB<-function(
  MCBset,
  training_set,
  Surv,
  testing_set=NULL,
  Surv.new=NULL,
  Method=c("svm","cox","enet","coxboost")[1],
  predict_time = 5,
  ci=FALSE,
  silent=FALSE,
  alpha = 0.5,
  n_mstop = 500,
  n_nu = 0.1,
  theta = 1
  ){
  requireNamespace("stats")
  # load private functions
  bs_ci <- function(data, indices, predict.time = predict_time) { 
    d <- data[indices,] # allows boot to select sample 
    surv.res = survivalROC(Stime = d$survival, 
                           status = d$survival_status, 
                           marker = d$marker, 
                           predict.time = predict.time, method = "NNE",
                           span = 0.25*NROW(d)^(-0.20))
    return(surv.res$AUC)
  }
  # end load
  
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
  #private functions
  create_doc<-function(ci){
    if (ci){
      write_MCB<-rep(NA,5)
      names(write_MCB)<-c('MCB_no','AUC_train','95_CI_train','AUC_test','95_CI_test')
      return(write_MCB)
    }else{
      write_MCB<-rep(NA,3)
      names(write_MCB)<-c('MCB_no','AUC_train','AUC_test')
      return(write_MCB)
    }
  }
  if (!Method %in% c("svm","cox","enet","coxboost")){
    stop(paste("Method:",Method,"is not supported, see hlep files for the details.",collapse = " "))
  }else if (Method=="svm") {
    # constuction of MCB Method matrix for SVM
    MCB_svm_matrix_training<-matrix(0,nrow = nrow(MCBset),ncol = ncol(training_set))
    colnames(MCB_svm_matrix_training)<-colnames(training_set)
    rownames(MCB_svm_matrix_training)<-as.numeric(MCBset[,'MCB_no'])
    #if it has a independent test set create the test_set res set
    if (!is.null(testing_set)) {
      MCB_svm_matrix_test_set<-matrix(0,nrow = nrow(MCBset),ncol = ncol(testing_set))
      colnames(MCB_svm_matrix_test_set)<-colnames(testing_set)
      rownames(MCB_svm_matrix_test_set)<-as.numeric(MCBset[,'MCB_no'])
    }else{
      MCB_svm_matrix_test_set<-NULL
    }
    FunctionResults<-NULL
    rz=!(is.na(Surv)|Surv[,1]==0)
    times=Surv[rz]
    best_auc<-0
    best_model<-NULL
    mcb_SVM_res<-NULL
    for (mcb in seq_len(nrow(MCBset))) {
      #if (nrow(MCBset)>1){}
      if (show_bar&!silent) {
        utils::setTxtProgressBar(bar, mcb)
      }
      write_MCB<-create_doc(ci)
      #save the mcb number
      write_MCB['MCB_no']<-as.numeric(MCBset[mcb,'MCB_no'])
      # build temp variable for saving the results.
      # MCB number
      # aquire information for CpG sites in MCB
      CpGs<-strsplit(MCBset[mcb,'CpGs']," ")[[1]]
      data_used_for_training<-data.frame(t(training_set[CpGs,rz]))
      # train a svm model
      svm_model <- tryCatch(survivalsvm::survivalsvm(times ~ ., data_used_for_training, gamma.mu = 0.1,type = "regression"),error = NULL)
      #predictions
      if (!is.null(svm_model)) {
        MCB_svm_matrix_training[mcb,]<-stats::predict(svm_model, data.frame(t(training_set[CpGs,])))$predicted
        auc_and_ci = calculate_auc_ci(survival = times,marker = MCB_svm_matrix_training[mcb,rz],predict_time,ci)
        write_MCB['AUC_train']<-auc_and_ci$AUC
        if (ci) write_MCB['95_CI_train']<-auc_and_ci$CI95
        #if it has a independent test set
        if (!is.null(testing_set)){
          MCB_svm_matrix_test_set[mcb,]<-stats::predict(svm_model, data.frame(t(testing_set[CpGs,])))$predicted
          auc_and_ci = calculate_auc_ci(Surv.new,marker = MCB_svm_matrix_test_set[mcb,],predict_time,ci)
          write_MCB['AUC_test']<-auc_and_ci$AUC
          if (ci) write_MCB['95_CI_test']<-auc_and_ci$CI95
          if ((write_MCB['AUC_train']+write_MCB['AUC_test'])>best_auc){
            best_auc<-write_MCB['AUC_train']+write_MCB['AUC_test']
            best_model<-list(mcb,svm_model)
          }
          #if it does not have a independent test set
        }else{
          write_MCB<-write_MCB[1:3]
          if (write_MCB['AUC_train']>best_auc){
            best_auc<-write_MCB['AUC_train']
            best_model<-list(mcb,svm_model)
          }
        }
      }else{
        stop("This svmr model can not be built.")
      }
      mcb_SVM_res<-rbind(mcb_SVM_res,write_MCB)
    }
    colnames(mcb_SVM_res)<-c("MCB_no","training_set_auc","test_set_auc")
    names(best_model)<-c("MCB_no","svm_model")
    FunctionResults$MCB_svm_matrix_training<-MCB_svm_matrix_training
    FunctionResults$MCB_svm_matrix_test_set<-MCB_svm_matrix_test_set
    FunctionResults$svm_auc_results<-mcb_SVM_res
    FunctionResults$maximum_auc<-best_auc
    FunctionResults$best_svm_model<-best_model
  }else if(Method=="cox"){
    # constuction of MCB Method matrix for cox
    MCB_cox_matrix_training<-matrix(0,nrow = nrow(MCBset),ncol = ncol(training_set))
    colnames(MCB_cox_matrix_training)<-colnames(training_set)
    rownames(MCB_cox_matrix_training)<-as.numeric(MCBset[,'MCB_no'])
    #if it has a independent test set create the test_set res set
    if (!is.null(testing_set)) {
      MCB_cox_matrix_test_set<-matrix(0,nrow = nrow(MCBset),ncol = ncol(testing_set))
      colnames(MCB_cox_matrix_test_set)<-colnames(testing_set)
      rownames(MCB_cox_matrix_test_set)<-as.numeric(MCBset[,'MCB_no'])
    }else{
      MCB_cox_matrix_test_set<-NULL
    }
    FunctionResults<-NULL
    rz=!(is.na(Surv)|Surv[,1]==0)
    times=Surv[rz]
    best_auc<-0
    best_model<-NULL
    mcb_cox_res<-NULL
    for (mcb in seq_len(nrow(MCBset))) {
      if (show_bar&!silent){utils::setTxtProgressBar(bar, mcb)}
      write_MCB<-create_doc(ci)
      #save the mcb number
      write_MCB['MCB_no']<-as.numeric(MCBset[mcb,'MCB_no'])
      # build temp variable for saving the results.
      # MCB number
      # aquire information for CpG sites in MCB
      CpGs<-strsplit(MCBset[mcb,'CpGs']," ")[[1]]
      data_used_for_training<-data.frame(allvars = as.ridgemat(t(training_set[CpGs,rz])))
      # train a ridge cox model
      if (length(CpGs)<=2){
        ridge_models<-tryCatch(survival::coxph(times ~ allvars,data=data_used_for_training),error = NULL)
      }else{
        ridge_models<-tryCatch(ridge_model(times,data_used_for_training,theta = theta),error = NULL)        
      }
      #predictions
      if (!is.null(ridge_models)) {
        ridge_models$CpGs <- CpGs
        ridge_models <- as.mcb.coxph.penal(ridge_models)
        MCB_cox_matrix_training[mcb,]<-predict.mcb.coxph.penal(ridge_models, data.frame(t(training_set[CpGs,])))
        auc_and_ci = calculate_auc_ci(survival = times,marker = MCB_cox_matrix_training[mcb,rz],predict_time,ci)
        write_MCB['AUC_train']<-auc_and_ci$AUC
        if (ci) write_MCB['95_CI_train']<-auc_and_ci$CI95
        #if it has a independent test set
        if (!is.null(testing_set)){
          MCB_cox_matrix_test_set[mcb,]<-predict.mcb.coxph.penal(ridge_models, data.frame(t(testing_set[CpGs,])))
          auc_and_ci = calculate_auc_ci(Surv.new,marker = MCB_cox_matrix_test_set[mcb,],predict_time,ci)
          write_MCB['AUC_test']<-auc_and_ci$AUC
          if (ci) write_MCB['95_CI_test']<-auc_and_ci$CI95
          if ((write_MCB['AUC_train']+write_MCB['AUC_test'])>best_auc){
            best_auc<-write_MCB['AUC_train']+write_MCB['AUC_test']
            best_model<-list(mcb,ridge_models)
          }
          #if it does not have a independent test set
        }else{
          write_MCB<-write_MCB[1:3]
          if (write_MCB['AUC_train']>best_auc){
            best_auc<-write_MCB['AUC_train']
            best_model<-list(mcb,ridge_models)
          }
        }
      }else{
        stop("This coxph model can not be built.")
      }
      mcb_cox_res<-rbind(mcb_cox_res,write_MCB)
    }
    colnames(mcb_cox_res)<-c("MCB_no","training_set_auc","test_set_auc")
    names(best_model)<-c("MCB_no","cox_model")
    FunctionResults$MCB_cox_matrix_training<-MCB_cox_matrix_training
    FunctionResults$MCB_cox_matrix_test_set<-MCB_cox_matrix_test_set
    FunctionResults$cox_auc_results<-mcb_cox_res
    FunctionResults$maximum_auc<-best_auc
    FunctionResults$best_cox_model<-best_model
  }else if (Method=="enet") {
    # constuction of MCB Method matrix for enet
    MCB_enet_matrix_training<-matrix(0,nrow = nrow(MCBset),ncol = ncol(training_set))
    colnames(MCB_enet_matrix_training)<-colnames(training_set)
    rownames(MCB_enet_matrix_training)<-as.numeric(MCBset[,'MCB_no'])
    #if it has a independent test set create the test_set res set
    if (!is.null(testing_set)) {
      MCB_enet_matrix_test_set<-matrix(0,nrow = nrow(MCBset),ncol = ncol(testing_set))
      colnames(MCB_enet_matrix_test_set)<-colnames(testing_set)
      rownames(MCB_enet_matrix_test_set)<-as.numeric(MCBset[,'MCB_no'])
    }else{
      MCB_enet_matrix_test_set<-NULL
    }
    FunctionResults<-NULL
    rz=!(is.na(Surv)|Surv[,1]==0)
    times=Surv[rz]
    best_auc<-0
    best_model<-NULL
    mcb_enet_res<-NULL
    for (mcb in seq_len(nrow(MCBset))) {
      if (show_bar&!silent){utils::setTxtProgressBar(bar, mcb)}
      write_MCB<-create_doc(ci)
      #save the mcb number
      write_MCB['MCB_no']<-as.numeric(MCBset[mcb,'MCB_no'])
      # build temp variable for saving the results.
      # MCB number
      # aquire information for CpG sites in MCB
      CpGs<-strsplit(MCBset[mcb,'CpGs']," ")[[1]]
      data_used_for_training<-t(training_set[CpGs,rz])
      # train a enet model
      enet_model <- tryCatch(glmnet::cv.glmnet(data_used_for_training,
                                                         times,
                                                         #cox model in enet was used, note that here cox and enet penalty were used.
                                                         family="cox",
                                                         #alpha = 0.5
                                                         alpha=alpha,
                                                         # The elasticnet mixing parameter, with 0≤α≤ 1. The penalty is defined as
                                                         # (1-alpha)/2||beta||_2^2+alpha||beta||_1
                                                         # alpha=1 is the lasso penalty, and alpha=0 the ridge penalty.
                                                         # type.measure = "AUC"
                                                         type.measure= "deviance"
                                                         # It uses AUC as the criterion for 10-fold cross-validation.
                                                         #foldid = 10
      ),error = NULL)
      #predictions
      if (!is.null(enet_model)) {
        correctional_value=1
        while ( sum(stats::coef(enet_model, s = enet_model$lambda.min-0.001*(correctional_value-1))>0)<1 &
               (enet_model$lambda.min-0.001*(correctional_value-1))>0 ) {
          correctional_value=correctional_value*1.25
        }
        lambda_min_corrected<-enet_model$lambda.min-0.001*(correctional_value-1)
        #if you use lambda.1se instead, the penalty of enet would be larger, leading that most of covariates were removed form the final model.
        MCB_enet_matrix_training[mcb,]<-stats::predict(enet_model,t(training_set[CpGs,]),s=lambda_min_corrected)
        auc_and_ci = calculate_auc_ci(survival = times,marker = MCB_enet_matrix_training[mcb,rz],predict_time,ci)
        write_MCB['AUC_train']<-auc_and_ci$AUC
        if (ci) write_MCB['95_CI_train']<-auc_and_ci$CI95
         #if it has a independent test set
        if (!is.null(testing_set)){
          # lambda.min was used.
          MCB_enet_matrix_test_set[mcb,]<-stats::predict(enet_model, t(testing_set[CpGs,]),s=lambda_min_corrected)
          auc_and_ci = calculate_auc_ci(Surv.new,marker = MCB_enet_matrix_test_set[mcb,],predict_time,ci)
          write_MCB['AUC_test']<-auc_and_ci$AUC
          if (ci) write_MCB['95_CI_test']<-auc_and_ci$CI95
          if ((write_MCB['AUC_train']+write_MCB['AUC_test'])>best_auc){
            best_auc<-write_MCB['AUC_train']+write_MCB['AUC_test']
            best_model<-list(mcb,enet_model,lambda_min_corrected)
          }
          #if it does not have a independent test set
        }else{
          write_MCB<-write_MCB[1:3]
          if (write_MCB['AUC_train']>best_auc){
            best_auc<-write_MCB['AUC_train']
            best_model<-list(mcb,enet_model,lambda_min_corrected)
          }
        }
      }else{
        stop("This enet model can not be built.")
      }
      mcb_enet_res<-rbind(mcb_enet_res,write_MCB)
    }
    colnames(mcb_enet_res)<-c("MCB_no","training_set_auc","test_set_auc")
    names(best_model)<-c("MCB_no","enet_model","corrected_lambda(min)")
    FunctionResults$MCB_enet_matrix_training<-MCB_enet_matrix_training
    FunctionResults$MCB_enet_matrix_test_set<-MCB_enet_matrix_test_set
    FunctionResults$enet_auc_results<-mcb_enet_res
    FunctionResults$maximum_auc<-best_auc
    FunctionResults$best_enet_model<-best_model
  }else if (Method=="coxboost") {
    # constuction of MCB Method matrix for CoxBoost
    MCB_coxboost_matrix_training<-matrix(0,nrow = nrow(MCBset),ncol = ncol(training_set))
    colnames(MCB_coxboost_matrix_training)<-colnames(training_set)
    rownames(MCB_coxboost_matrix_training)<-as.numeric(MCBset[,'MCB_no'])
    #if it has a independent test set create the test_set res set
    if (!is.null(testing_set)) {
      MCB_coxboost_matrix_test_set<-matrix(0,nrow = nrow(MCBset),ncol = ncol(testing_set))
      colnames(MCB_coxboost_matrix_test_set)<-colnames(testing_set)
      rownames(MCB_coxboost_matrix_test_set)<-as.numeric(MCBset[,'MCB_no'])
    }else{
      MCB_coxboost_matrix_test_set<-NULL
    }
    FunctionResults<-NULL
    rz=!(is.na(Surv)|Surv[,1]==0)
    times=Surv[rz]
    best_auc<-0
    best_model<-NULL
    mcb_coxboost_res<-NULL
    for (mcb in seq_len(nrow(MCBset))) {
      #if (nrow(MCBset)>1){}
      if (show_bar&!silent) {
        utils::setTxtProgressBar(bar, mcb)
      }
      write_MCB<-create_doc(ci)
      #save the mcb number
      write_MCB['MCB_no']<-as.numeric(MCBset[mcb,'MCB_no'])
      # build temp variable for saving the results.
      # MCB number
      # aquire information for CpG sites in MCB
      CpGs<-strsplit(MCBset[mcb,'CpGs']," ")[[1]]
      data_used_for_training<-t(training_set[CpGs,rz])
      # train a coxboost model
      coxboost_model <- tryCatch(mboost::glmboost(y=times,x=data_used_for_training,family=mboost::CoxPH(),
                                          control=mboost::boost_control(mstop=n_mstop,nu=n_nu)),error = NULL)
      #predictions
      if (!is.null(coxboost_model)) {
        MCB_coxboost_matrix_training[mcb,]<-stats::predict(coxboost_model, t(training_set[CpGs,]))[,1]
        auc_and_ci = calculate_auc_ci(survival = times,marker = MCB_coxboost_matrix_training[mcb,rz],predict_time,ci)
        write_MCB['AUC_train']<-auc_and_ci$AUC
        if (ci) write_MCB['95_CI_train']<-auc_and_ci$CI95
        #if it has a independent test set
        if (!is.null(testing_set)){
          MCB_coxboost_matrix_test_set[mcb,]<-stats::predict(coxboost_model, t(testing_set[CpGs,]))[,1]
          auc_and_ci = calculate_auc_ci(Surv.new,marker = MCB_coxboost_matrix_test_set[mcb,],predict_time,ci)
          write_MCB['AUC_test']<-auc_and_ci$AUC
          if (ci) write_MCB['95_CI_test']<-auc_and_ci$CI95
          if ((write_MCB['AUC_train']+write_MCB['AUC_test'])>best_auc){
            best_auc<-write_MCB['AUC_train']+write_MCB['AUC_test']
            best_model<-list(mcb,coxboost_model)
          }
          #if it does not have a independent test set
        }else{
          write_MCB<-write_MCB[1:3]
          if (write_MCB['AUC_train']>best_auc){
            best_auc<-write_MCB['AUC_train']
            best_model<-list(mcb,coxboost_model)
          }
        }
      }else{
        stop("This CoxBoost model can not be built.")
      }
      mcb_coxboost_res<-rbind(mcb_coxboost_res,write_MCB)
    }
    colnames(mcb_coxboost_res)<-c("MCB_no","training_set_auc","test_set_auc")
    names(best_model)<-c("MCB_no","coxboost_model")
    FunctionResults$MCB_coxboost_matrix_training<-MCB_coxboost_matrix_training
    FunctionResults$MCB_coxboost_matrix_test_set<-MCB_coxboost_matrix_test_set
    FunctionResults$coxboost_auc_results<-mcb_coxboost_res
    FunctionResults$maximum_auc<-best_auc
    FunctionResults$best_coxboost_model<-best_model
  }
  return(FunctionResults)
}



#' @title Calculation of the metric matrix for Methylation Correlation Block
#'
#' @description To enable quantitative analysis of the methylation patterns
#' within individual Methylation Correlation Blocks across many samples, a single metric to
#' define the methylated pattern of multiple CpG sites within each block.
#' Compound scores which calculated all CpGs within individual Methylation Correlation Blocks by SVM model
#' were used as the compound methylation values of Methylation Correlation Blocks.
#' @usage metricMCB(MCBset,training_set,Surv,testing_set,Surv.new,Method,silent)
#' @export
#' @param training_set methylation matrix used for training the model in the analysis.
#' @param testing_set methylation matrix used in the analysis. This can be missing then training set itself will be used as testing set.
#' @param MCBset Methylation Correlation Block information returned by the IndentifyMCB function.
#' @param Surv Survival function contain the survival information for training.
#' @param Surv.new Survival function contain the survival information for testing.
#' @param Method model used to calculate the compound values for multiple Methylation correlation blocks. Options include "svm" "cox" and "lasso". The default option is SVM method.
#' @param silent Ture indicates that processing information and progress bar will be shown.
#' @author Xin Yu
#' @keywords Methylation Correlation
#' @examples
#' #Not run! Remove # to run the example.
#'  require(EnMCB)
#' #import exprs function
#' data(demo_survival_data)
#' datamatrix<-create_demo()
#' data(demo_MCBinformation)
#' #select MCB with at least 3 CpGs.
#' demo_MCBinformation<-demo_MCBinformation[demo_MCBinformation[,"CpGs_num"]>2,]
#'
#' trainingset<-colnames(datamatrix) %in% sample(colnames(datamatrix),0.6*length(colnames(datamatrix)))
#' testingset<-!trainingset
#' #create the results using Cox regression. Remove # to run the example.
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
#'
metricMCB<-function(
  MCBset,
  training_set,
  Surv,
  testing_set=NULL,
  Surv.new=NULL,
  Method=c("svm","cox","lasso")[1],
  silent=F
  ){
  requireNamespace("stats")
  if (!silent) {
    cat("Start anaylsis, this may take a while...\n")
    show_bar=nrow(MCBset)>1
  }else{
    show_bar=F
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
  if (!Method %in% c("svm","cox","lasso")){
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
    for (mcb in 1:nrow(MCBset)) {
      #if (nrow(MCBset)>1){}
      if (show_bar&!silent) {
        utils::setTxtProgressBar(bar, mcb)
      }
      write_MCB<-c(NA,NA,NA)
      #save the mcb number
      write_MCB[1]<-as.numeric(MCBset[mcb,'MCB_no'])
      # build temp variable for saving the results.
      # MCB number
      # aquire information for CpG sites in MCB
      CpGs<-strsplit(MCBset[mcb,'CpGs']," ")[[1]]
      data_used_for_training<-data.frame(t(training_set[CpGs,rz]))
      # train a svm model
      svm_model<-NULL
      try(svm_model <- survivalsvm::survivalsvm(times ~ ., data_used_for_training, gamma.mu = 0.1,type = "regression"),silent = T)
      #predictions
      if (!is.null(svm_model)) {
        MCB_svm_matrix_training[mcb,]<-stats::predict(svm_model, data.frame(t(training_set[CpGs,])))$predicted
        write_MCB[2]<-survivalROC::survivalROC(Stime = times[,1],status = times[,2],marker = MCB_svm_matrix_training[mcb,rz],predict.time = 5,method = "NNE",span =0.25*length(times)^(-0.20)  )$AUC
        #if it has a independent test set
        if (!is.null(testing_set)){
          MCB_svm_matrix_test_set[mcb,]<-stats::predict(svm_model, data.frame(t(testing_set[CpGs,])))$predicted
          write_MCB[3]<-survivalROC::survivalROC(Stime = Surv.new[,1],status = Surv.new[,2],marker = MCB_svm_matrix_test_set[mcb,],predict.time = 5,method = "NNE",span =0.25*length(Surv.new)^(-0.20) )$AUC
          if (abs(write_MCB[2]+write_MCB[3]-1)>best_auc){
            best_auc<-abs(write_MCB[2]+write_MCB[3]-1)
            best_model<-list(mcb,svm_model)
          }
          #if it does not have a independent test set
        }else{
          if (abs(write_MCB[2]-0.5)>best_auc){
            best_auc<-abs(write_MCB[2]-0.5)+0.5
            best_model<-list(mcb,svm_model)
          }
        }
      }
      mcb_SVM_res<-rbind(mcb_SVM_res,write_MCB)
    }
    cat("\n")
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
    for (mcb in 1:nrow(MCBset)) {
      if (show_bar&!silent){utils::setTxtProgressBar(bar, mcb)}
      write_MCB<-c(NA,NA,NA)
      #save the mcb number
      write_MCB[1]<-as.numeric(MCBset[mcb,'MCB_no'])
      # build temp variable for saving the results.
      # MCB number
      # aquire information for CpG sites in MCB
      CpGs<-strsplit(MCBset[mcb,'CpGs']," ")[[1]]
      data_used_for_training<-data.frame(t(training_set[CpGs,rz]))
      # train a cox model
      univ_models<-NULL
      try(univ_models<-survival::coxph(times ~.,data=data_used_for_training),silent = T)
      #predictions
      if (!is.null(univ_models)) {
        MCB_cox_matrix_training[mcb,]<-stats::predict(univ_models, data.frame(t(training_set[CpGs,])))
        write_MCB[2]<-survivalROC::survivalROC(Stime = times[,1],status = times[,2],marker = MCB_cox_matrix_training[mcb,rz],predict.time = 5,method = "NNE",span = 0.25*length(times)^(-0.20))$AUC
        #if it has a independent test set
        if (!is.null(testing_set)){
          MCB_cox_matrix_test_set[mcb,]<-stats::predict(univ_models, data.frame(t(testing_set[CpGs,])))
          write_MCB[3]<-survivalROC::survivalROC(Stime = Surv.new[,1],status = Surv.new[,2],marker = MCB_cox_matrix_test_set[mcb,],predict.time = 5,method = "NNE",span =0.25*length(Surv.new)^(-0.20))$AUC
          if (abs(write_MCB[2]+write_MCB[3]-1)>best_auc){
            best_auc<-abs(write_MCB[2]+write_MCB[3]-1)
            best_model<-list(mcb,univ_models)
          }
          #if it does not have a independent test set
        }else{
          if (abs(write_MCB[2]-0.5)>best_auc){
            best_auc<-abs(write_MCB[2]-0.5)+0.5
            best_model<-list(mcb,univ_models)
          }
        }
      }
      mcb_cox_res<-rbind(mcb_cox_res,write_MCB)
    }
    cat("\n")
    colnames(mcb_cox_res)<-c("MCB_no","training_set_auc","test_set_auc")
    names(best_model)<-c("MCB_no","cox_model")
    FunctionResults$MCB_cox_matrix_training<-MCB_cox_matrix_training
    FunctionResults$MCB_cox_matrix_test_set<-MCB_cox_matrix_test_set
    FunctionResults$cox_auc_results<-mcb_cox_res
    FunctionResults$maximum_auc<-best_auc
    FunctionResults$best_cox_model<-best_model
  }else if (Method=="lasso") {
    # constuction of MCB Method matrix for lasso
    MCB_lasso_matrix_training<-matrix(0,nrow = nrow(MCBset),ncol = ncol(training_set))
    colnames(MCB_lasso_matrix_training)<-colnames(training_set)
    rownames(MCB_lasso_matrix_training)<-as.numeric(MCBset[,'MCB_no'])
    #if it has a independent test set create the test_set res set
    if (!is.null(testing_set)) {
      MCB_lasso_matrix_test_set<-matrix(0,nrow = nrow(MCBset),ncol = ncol(testing_set))
      colnames(MCB_lasso_matrix_test_set)<-colnames(testing_set)
      rownames(MCB_lasso_matrix_test_set)<-as.numeric(MCBset[,'MCB_no'])
    }else{
      MCB_lasso_matrix_test_set<-NULL
    }
    FunctionResults<-NULL
    rz=!(is.na(Surv)|Surv[,1]==0)
    times=Surv[rz]
    best_auc<-0
    best_model<-NULL
    mcb_lasso_res<-NULL
    for (mcb in 1:nrow(MCBset)) {
      if (show_bar&!silent){utils::setTxtProgressBar(bar, mcb)}
      write_MCB<-c(NA,NA,NA)
      #save the mcb number
      write_MCB[1]<-as.numeric(MCBset[mcb,'MCB_no'])
      # build temp variable for saving the results.
      # MCB number
      # aquire information for CpG sites in MCB
      CpGs<-strsplit(MCBset[mcb,'CpGs']," ")[[1]]
      data_used_for_training<-t(training_set[CpGs,rz])
      # train a lasso model
      lasso_model<-NULL
      try(lasso_model <- glmnet::cv.glmnet(data_used_for_training,
                                                         times,
                                                         #cox model in lasso was used, note that here cox and lasso penalty were used.
                                                         family="cox",
                                                         alpha=0.5,
                                                         # The elasticnet mixing parameter, with 0≤α≤ 1. The penalty is defined as
                                                         # (1-alpha)/2||beta||_2^2+alpha||beta||_1
                                                         # alpha=1 is the lasso penalty, and alpha=0 the ridge penalty.
                                                         # type.measure = "AUC"
                                                         type.measure= "deviance"
                                                         # It uses AUC as the criterion for 10-fold cross-validation.
                                                         #foldid = 10
      ),silent = T)
      #predictions
      if (!is.null(lasso_model)) {
        correctional_value=1
        while ( sum(stats::coef(lasso_model, s = lasso_model$lambda.min-0.001*(correctional_value-1))>0)<1 &
               (lasso_model$lambda.min-0.001*(correctional_value-1))>0 ) {
          correctional_value=correctional_value*1.25
        }
        lambda_min_corrected<-lasso_model$lambda.min-0.001*(correctional_value-1)
        #if you use lambda.1se instead, the penalty of lasso would be larger, leading that most of covariates were removed form the final model.
        MCB_lasso_matrix_training[mcb,]<-stats::predict(lasso_model,t(training_set[CpGs,]),s=lambda_min_corrected)
        write_MCB[2]<-survivalROC::survivalROC(Stime = times[,1],status = times[,2],marker = MCB_lasso_matrix_training[mcb,rz],predict.time = 5,method = "NNE",span = 0.25*length(times)^(-0.20) )$AUC
        #if it has a independent test set
        if (!is.null(testing_set)){
          # lambda.min was used.
          MCB_lasso_matrix_test_set[mcb,]<-stats::predict(lasso_model, t(testing_set[CpGs,]),s=lambda_min_corrected)
          write_MCB[3]<-survivalROC::survivalROC(Stime = Surv.new[,1],status = Surv.new[,2],marker = MCB_lasso_matrix_test_set[mcb,],predict.time = 5,method = "NNE",span = 0.25*length(Surv.new)^(-0.20))$AUC
          if (abs(write_MCB[2]+write_MCB[3]-1)>best_auc){
            best_auc<-abs(write_MCB[2]+write_MCB[3]-1)
            best_model<-list(mcb,lasso_model,lambda_min_corrected)
          }
          #if it does not have a independent test set
        }else{
          if (abs(write_MCB[2]-0.5)>best_auc){
            best_auc<-abs(write_MCB[2]-0.5)+0.5
            best_model<-list(mcb,lasso_model,lambda_min_corrected)
          }
        }
      }
      mcb_lasso_res<-rbind(mcb_lasso_res,write_MCB)
    }
    cat("\n")
    colnames(mcb_lasso_res)<-c("MCB_no","training_set_auc","test_set_auc")
    names(best_model)<-c("MCB_no","lasso model","corrected lambda(min)")
    FunctionResults$MCB_lasso_matrix_training<-MCB_lasso_matrix_training
    FunctionResults$MCB_lasso_matrix_test_set<-MCB_lasso_matrix_test_set
    FunctionResults$lasso_auc_results<-mcb_lasso_res
    FunctionResults$maximum_auc<-best_auc
    FunctionResults$best_lasso_model<-best_model
  }
  return(FunctionResults)
}

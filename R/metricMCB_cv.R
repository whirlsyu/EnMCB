#' @title Calculation of model AUC for Methylation Correlation Blocks using cross validation
#'
#' @description To enable quantitative analysis of the methylation patterns
#' within individual Methylation Correlation Blocks across many samples, a single metric to
#' define the methylated pattern of multiple CpG sites within each block.
#' Compound scores which calculated all CpGs within individual Methylation Correlation Blocks by SVM model
#' were used as the compound methylation values of Methylation Correlation Blocks.
#' @usage metricMCB(MCBset,training_set,Surv,testing_set,Surv.new,Method,silent)
#' @export
#' @param MCBset Methylation Correlation Block information returned by the IndentifyMCB function.
#' @param data_set methylation matrix used for training the model in the analysis.
#' @param Surv Survival function contain the survival information for training.
#' @param nfold fold used in the cross validation precedure.
#' @param Method model used to calculate the compound values for multiple Methylation correlation blocks. Options include "svm" "cox" and "lasso". The default option is SVM method.
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
#' mcb_cox_res<-metricMCB(MCBset = demo_MCBinformation,
#'                data_set = datamatrix,
#'                Surv = demo_survival_data,
#'                Method = "cox")
#'
#' @return Object of class \code{list} with elements (XXX will be replaced with the model name you choose):
#'  \tabular{ll}{
#'    \code{MCB_xxx_matrix} \tab Prediction results of model
#'    \code{XXX_auc_results} \tab AUC results for each model. \cr
#'    \code{best_XXX_model} \tab Model object for the model with best AUC. \cr
#'    \code{maximum_auc} \tab Maximum AUC for the whole generated models. \cr
#'  }
#' @references
#' Xin Yu et al. 2019 Predicting disease progression in lung adenocarcinoma patients based on methylation correlated blocks using ensemble machine learning classifiers (under review)
#'
metricMCB.cv<-function(
  MCBset,
  data_set,
  Surv,
  nfold=10,
  Method=c("svm","cox","lasso")[1],
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
  # constuction of MCB Method matrix for SVM
  MCB_matrix<-matrix(0,nrow = nrow(MCBset),ncol = ncol(data_set))
  colnames(MCB_matrix)<-colnames(data_set)
  rownames(MCB_matrix)<-as.numeric(MCBset[,'MCB_no'])
  if (sum(is.na(Surv)|Surv[,1]==0)>0){
    stop("survival data contain NA, please remove NA")
  }
  
  if (!Method %in% c("svm","cox","lasso")){
    stop(paste("Method:",Method,"is not supported, see hlep files for the details.",collapse = " "))
  }
    sp<-sample(1:ncol(data_set),replace = F)
    order_sp<-order(sp)
    data_set<-data_set[,sp]
    folds <- cut(seq(1,ncol(data_set)),breaks=ncol(data_set)/nfold,labels=FALSE)
    #if it has a independent test set create the test_set res set
    FunctionResults<-NULL
    best_auc<-0
    best_model<-NULL
    mcb_SVM_res<-NULL
    for (mcb in seq_len(nrow(MCBset))) {
      #if (nrow(MCBset)>1){}
      if (show_bar&!silent) {
        utils::setTxtProgressBar(bar, mcb)
      }
      write_MCB<-c(NA,NA,NA)
      #save the mcb number
      write_MCB[1]<-as.numeric(MCBset[mcb,'MCB_no'])
      write_MCB[2]<-MCBset[mcb,'CpGs']
      # build temp variable for saving the results.
      # MCB number
      # aquire information for CpG sites in MCB
      CpGs<-strsplit(MCBset[mcb,'CpGs']," ")[[1]]
      for (i in seq(unique(folds))) {
        rz<- which(folds==i,arr.ind=TRUE)
        data_used_for_training<-data.frame(t(data_set[CpGs,-rz]))
        data_used_for_testing <-data.frame(t(data_set[CpGs,rz]))
        # train a svm model
        times = Surv
          if (Method=="svm") {
          model<-tryCatch(survivalsvm::survivalsvm(times ~ ., 
                                                     data_used_for_training, 
                                                     gamma.mu = 0.1,
                                                     type = "regression"),
                            error = function(e){return(NULL)})
          }else if(Method=="cox"){
            model<-tryCatch(survival::coxph(times ~ ., 
                                            data_used_for_training),
                            error = function(e){return(NULL)})
          }else if(Method=="lasso"){
            model<-tryCatch( glmnet::cv.glmnet(data_used_for_training,
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
            ),error = function(e){return(NULL)})
            correctional_value=1
            while ( sum(stats::coef(lasso_model, s = lasso_model$lambda.min-0.001*(correctional_value-1))>0)<1 &
                    (lasso_model$lambda.min-0.001*(correctional_value-1))>0 ) {
              correctional_value=correctional_value*1.25
            }
            lambda_min_corrected<-lasso_model$lambda.min-0.001*(correctional_value-1)
          }
        
        if (!is.null(model)) {
          #predictions
          MCB_matrix[mcb,rz]<-stats::predict(model,data_used_for_testing)$predicted
          if (Method=="lasso"){
            #if you use lambda.1se instead, the penalty of lasso would be larger, leading that most of covariates were removed form the final model.
            MCB_matrix[mcb,rz]<-stats::predict(model,data_used_for_testing,s=lambda_min_corrected)
          }
        }else{
          stop('Model could not be built. check your dataset.')
        }
      }
      MCB_matrix[mcb,]<-MCB_matrix[mcb,order_sp]
      write_MCB[3]<-survivalROC::survivalROC(Stime = Surv[,1],
                                             status = Surv[,2],
                                             marker = MCB_svm_matrix_training[mcb,],
                                             predict.time = 5,
                                             method = "NNE",
                                             span =0.25*length(times)^(-0.20))$AUC
      mcb_SVM_res<-rbind(mcb_SVM_res,write_MCB)
      
    }
    colnames(mcb_SVM_res)<-c("MCB_no","CpGs","auc")
    FunctionResults$MCB_matrix<-MCB_matrix
    FunctionResults$auc_results<-mcb_SVM_res
  return(FunctionResults)
}

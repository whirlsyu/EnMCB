# load private functions
get450kAnno = function() {
  
  ca = BiocFileCache::BiocFileCache()
  
  q = BiocFileCache::bfcquery(ca, "450kanno.ilmn12.hg19")
  
  if (length(q$rpath)>0) return(readRDS(rev(q$rpath)[1])) # if multiple, use last
  
  if (!requireNamespace("minfi")) stop("install the minfi package to use this function")
  
  anno = minfi::getAnnotation("IlluminaHumanMethylation450kanno.ilmn12.hg19")
  
  tf = tempfile()
  
  saveRDS(anno, tf)
  
  BiocFileCache::bfcadd(ca, rname="IlluminaHumanMethylation450kanno.ilmn12.hg19", fpath=tf,
                        
                        action="move")  # for future
  
  anno
  
}

print_as_data <- function(variables,file) {
  sink(file)
  print(variables)
  sink()
}

ridge_model <- function(times,data_used_for_training,theta = 1){
  survival::coxph(times ~ survival::ridge(allvars,theta = theta),data=data_used_for_training)
}

#' @title ridge matrix
#' @description as.matrix attempts to turn its argument 
#' @param x data vector
#' @usage as.ridgemat(x)
#' @export
as.ridgemat <- function(x) {
  class(x) <- c("ridgemat", class(x))
  x
}

as.mcb.coxph.penal <- function(x){
  class(x) <- c("mcb.coxph.penal", class(x))
  x
}

removeobjectclass<- function(x){
  class(x) <- setdiff(class(x),"mcb.coxph.penal")
  x
}

#' @title data frame ridge matrix
#' @description data frame ridge matrix
#' @param x data vector
#' @param ... other parameters pass to as.data.frame.model.matrix()
#' @method as.data.frame ridgemat
#' @export
as.data.frame.ridgemat <- function(x, ...) {
  as.data.frame.model.matrix(as.matrix(x), ...)
}


#'@title predict coxph penal using MCB
#'@author Xin Yu
#'@description Compute fitted values and regression terms for a model fitted by coxph
#'@param object the results of a coxph fit.
#'@param newdata Optional new data at which to do predictions. If absent predictions are for the data frame used in the original fit. 
#'When coxph has been called with a formula argument created in another context, i.e., coxph has been called within another function and the formula was passed as an argument to that function, 
#'there can be problems finding the data set. See the note below.
#'@param ... other parameters pass to predict.coxph
#'@method predict mcb.coxph.penal
#'@return prediction values of regression.
#'@export
predict.mcb.coxph.penal <- function(object,newdata,...){
  internames = colnames(newdata) %in% object$CpGs
  if (sum(internames) == length(object$CpGs)){
    newdata <- data.frame(allvars = as.ridgemat(newdata[,internames]))
    object = removeobjectclass(object)
    return(predict(object,newdata,...))
  }else{
    stop('newdata do not contain sufficient CpGs in the model')
  }
}


test_logrank_p_in_KM<-function(mcb_matrix,y_surv){
  p.val_all<-NULL
  for (mcb in rownames(mcb_matrix)){
    Test<-mcb_matrix[mcb,]
    group_sur<-factor(1*(Test>median(Test)),
                      levels=0:1,
                      labels=c("High risk","Low risk"))
    diff_tcga<-survival::survdiff(y_surv ~ group_sur)
    p.val <- 1 - stats::pchisq(diff_tcga$chisq, length(diff_tcga$n) - 1)
    p.val_all<-rbind(p.val_all,p.val)
  }
  names(p.val_all)<-rownames(mcb_matrix)
  p.val_all
}


give_me_the_sen_spc<-function(pred,print=TRUE) {
  tpr<-(pred@tp[[1]]/(pred@tp[[1]]+pred@fn[[1]])) #sensitivity
  fpr<-(pred@fp[[1]]/(pred@tn[[1]]+pred@fp[[1]])) #specificity
  best_one<-which.min(sqrt( (1-tpr)^2+ fpr^2 ))
  if (print == TRUE)  cat(paste(" sensitivity:",tpr[best_one],"\n","specificity:",1-fpr[best_one],"\n"))
  pred@cutoffs[[1]][best_one]#specificity
}

random_matrix <- function(x) {
  x_random=matrix(x[sample(seq_along(x),length(x))],ncol = ncol(x),nrow = nrow(x))
  rownames(x_random)<-rownames(x)
  colnames(x_random)<-colnames(x)
  x_random
}


#'@title create demo matrix
#'@author Xin Yu
#'@description Demo matrix for methylation matrix.
#'@param model Two options, 'all' or 'short' for creating full dataset or very brief demo.
#'@return This function will generate a demo data.
#'@export
#'@examples 
#'demo_set<-create_demo()
create_demo<-function(model=c('all','short')[1]){
  utils::data(demo_data)
  rmax<-matrix(data = rnorm(length(demo_data$rownames)*length(demo_data$colnames), mean=0.5, sd=0.08),
           nrow = length(demo_data$rownames),
           ncol=length(demo_data$colnames),dimnames = list(demo_data$rownames,demo_data$colnames)
           )
  rmax[rownames(demo_data$realdata),]<-demo_data$realdata
  if (model=='short') {
    rmax<-rmax[ceiling(nrow(rmax)*0.8):nrow(rmax),]
  }
  rmax
}



as.numeric_matrix<-function(met_matrixf){
  met_matrixff<-matrix(as.numeric(met_matrixf),
                       nrow = NROW(met_matrixf),
                       ncol = NCOL(met_matrixf))
  rownames(met_matrixff)<-rownames(met_matrixf)
  colnames(met_matrixff)<-colnames(met_matrixf)
  met_matrixff
}

makeEset<-function(met_matrix,cli_data){
  all_col<-NULL
  for (i in seq_len(ncol(met_matrix))) {
    if (length(grep(pattern=tolower(colnames(met_matrix)[i]),cli_data))>0) {
      loca<-grep(pattern=tolower(colnames(met_matrix)[i]),cli_data)[1]
      col<-loca%/%nrow(cli_data)+1
      row<-loca%%nrow(cli_data)
      all_col<-c(all_col,col)
    }
  }
  cli_data_f<-cli_data[,all_col]
  colnames(cli_data_f)<-colnames(met_matrix)
  eSet<-make_eset(met_matrix,data.frame(t(cli_data_f)))
  eSet
}


multiple_time_ROC <- function(Test,y_surv,genesel) {
  requireNamespace("prognosticROC") #version 0.7
  sroclong_all<-NULL
  for (n_time in c(3,5,7)) {
    ROC_res= survivalROC::survivalROC(Stime=y_surv[,1],
                         status=y_surv[,2],
                         marker =as.numeric(Test),
                         lambda = NULL,window = "asymmetric",
                         predict.time = n_time,method = "NNE",span =0.25*length(y_surv[,1])^(-0.20))#
    # AUC at 7 years: 0·670 (0·578–0·762)
    sroclong_all<-ROCdata_save(sroclong_all,ROC_res,mark = paste("AUC at",n_time,"years:",round(ROC_res$AUC,2),collapse = " "))
  }
  gg<-ggplot(sroclong_all, aes(x = FPF, y = TPF, label = c, color = groups))+
    coord_cartesian(ylim=c(0, 1.05),xlim=c(0, 1.05),expand=FALSE)+
    geom_roc(labels = FALSE, stat = "identity") + style_roc(theme = theme_light())+
    theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(),axis.line = element_line(colour = "black"))+
    theme(legend.position=c(0.8,0.35), legend.justification=c(1,1),legend.title = element_blank())+
    geom_abline(slope=1, colour="black")
  ggsave(filename = paste("Time ROC of ",genesel,".jpeg",sep=""),plot = gg,device ="jpeg" ,
         path = getwd(),dpi = 300,units = "in",width = 5, height = 4.5,
         limitsize=FALSE)
}

draw_pheatmap_matrix<-function(eSet_matrix,group,savename="draw_pheatmap_matrix",scale = FALSE){
  requireNamespace(pheatmap)
  requireNamespace(ggplot2)
  group<-as.character(group)
  df<-as.numeric_matrix(eSet_matrix)
  df<-df[,order(group,decreasing = TRUE)]
  group<-group[order(group,decreasing = TRUE)]
  df[is.na(df)]<-0
  aka2 = data.frame(ID = group)
  rownames(aka2)<-colnames(df)
  aka3 = list(ID = c('TRUE' = "blue", 'FALSE'="red"))

  df<-df[apply(df, 1,stats::sd)>0,]
  if (scale == TRUE) {
    df=scale(df)
  }
   gg<-pheatmap(df,
           annotation_col = aka2,
           annotation_colors = aka3[1],
           annotation_legend = FALSE,
           show_colnames = FALSE, show_rownames = FALSE, cluster_rows = TRUE,
           cluster_cols = FALSE, legend = TRUE,
           clustering_distance_rows = "euclidean", border_color = FALSE)
   ggplot2::ggsave(filename = paste(savename,".jpeg",sep=""),plot = gg,device ="jpeg" ,
          path = getwd(),dpi = 300,units = "in",width = 10, height = 5,
          limitsize=FALSE)
}


give_me_the_lasso_gene<-function(cvfit,s=NULL,plot_figure=TRUE){
  if (plot_figure ==TRUE) {
    plot(cvfit)
    # cat("min: ",cvfit$lambda.min,"\n")
    # cat("1se: ",cvfit$lambda.1se,"\n")
  }
  if (is.null(s)) {
    s=cvfit$lambda.min
  }
  Coefficients <- coef(cvfit, s = s) #0.20
  Active.Index <- which(as.logical(Coefficients != 0))
  Active.Coefficients <- Coefficients[Active.Index]
  #print(Active.Coefficients)
  gene_lasso<-Coefficients@Dimnames[[1]][Active.Index]
  if (plot_figure ==TRUE) {
    plot(cvfit$glmnet.fit,xvar="lambda")
  }
  list(gene_lasso=gene_lasso,Active.Coefficients=Active.Coefficients)
}


plot_ROC<-function(perf){
  sroclong<-data.frame(TPF = perf@y.values[[1]], FPF =perf@x.values[[1]],
                       c = perf@alpha.values[[1]] #,
                       #group = rep("M stage", length(Mayo5[["FP"]]))
  )
  ggplot(sroclong, aes(x = FPF, y = TPF, label = c
                       #, color = group
  )) +
  geom_roc(labels = FALSE, stat = "identity") + style_roc(theme = theme_grey)
}

ROC_mutiple_clinical<-function(test_frame,y_surv,genesel="title",ntime=5){
    requireNamespace("prognosticROC") #version 0.7
    sroclong_all<-NULL
    for (n in seq_len(ncol(test_frame))) {
      ROC_res= survivalROC::survivalROC(Stime=y_surv[,1],
                                        status=y_surv[,2],
                                        marker =as.numeric(test_frame[,n]),
                                        lambda = NULL,window = "asymmetric",
                                        predict.time = ntime,method = "NNE",span =0.25*length(y_surv[,1])^(-0.20))#
      sroclong_all<-ROCdata_save(sroclong_all,ROC_res,mark = paste(ntime,"year AUC at",colnames(test_frame)[n],round(ROC_res$AUC,2),collapse = " "))
    }
    gg<-ggplot2::ggplot(sroclong_all, aes(x = FPF, y = TPF, label = c, color = groups))+
      coord_cartesian(ylim=c(0, 1.05),xlim=c(0, 1.05),expand=FALSE)+
      geom_roc(labels = FALSE, stat = "identity") + style_roc(theme = theme_light())+
      theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank(),axis.line = element_line(colour = "black"))+
      theme(legend.position=c(0.8,0.35), legend.justification=c(1,1),legend.title = element_blank())+
      geom_abline(slope=1, colour="black")
    ggplot2::ggsave(filename = paste("Time ROC of ",genesel,".jpeg",sep=""),plot = gg,device ="jpeg" ,
           path = getwd(),dpi = 300,units = "in",width = 5, height = 4.5,
           limitsize=FALSE)
  }

ROC_threshold <- function(predict, response) {
  perf <-prediction(as.numeric(predict),as.numeric(response))
  return(give_me_the_sen_spc(perf,print = FALSE))
}

ROCdata_save<-function(origin=NULL,perf,mark="none"){
  sroclong<-data.frame(TPF = perf[["TP"]], FPF = perf[["FP"]],
                       c = perf[["cut.values"]],
                       groups = rep(mark,
                                    length(perf[["FP"]])
                                    )
                       )
  sroclong<-rbind(origin,sroclong)
  sroclong
}

mutiple_time_ROC <- function(Test,y_surv,genesel) {
  requireNamespace(survivalROC) #version 1.0.3
  requireNamespace(prognosticROC) #version 0.7
  sroclong_all<-NULL
  for (n_time in c(3,5,7)) {
    ROC_res= survivalROC::survivalROC(Stime=y_surv[,1],
                         status=y_surv[,2],
                         marker =as.numeric(Test),
                         lambda = NULL,window = "asymmetric",
                         predict.time = n_time,method = "NNE",span =0.25*length(y_surv[,1])^(-0.20))#
    # AUC at 7 years: 0·670 (0·578–0·762)

    sroclong<-data.frame(TPF = ROC_res[["TP"]], FPF = ROC_res[["FP"]],
                         c = ROC_res[["cut.values"]],
                         groups = rep(paste("AUC at",n_time,"years:",round(ROC_res$AUC,2),collapse = " "),
                                     length(ROC_res[["FP"]])))
    sroclong_all<-rbind(sroclong_all,sroclong)
  }
  gg<-ggplot(sroclong_all, aes(x = FPF, y = TPF, label = c, color = groups))+
    coord_cartesian(ylim=c(0, 1.05),xlim=c(0, 1.05),expand=FALSE)+
    geom_roc(labels = FALSE, stat = "identity") + style_roc(theme = theme_light())+
    theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(),axis.line = element_line(colour = "black"))+
    theme(legend.position=c(0.8,0.35), legend.justification=c(1,1),legend.title = element_blank())+
    geom_abline(slope=1, colour="black")
  ggsave(filename = paste("Time ROC of ",genesel,".jpeg",sep=""),plot = gg,device ="jpeg" ,
         path = getwd(),dpi = 300,units = "in",width = 5, height = 4.5,
         limitsize=FALSE)
}




#'@title draw survival curve
#'@author Xin Yu
#'@description Draw a survival curve based on survminer package. This is a wrapper function of ggsurvplot.
#'@param exp expression level for variable.
#'@param living_days The survival time (days) for each individual.
#'@param living_events The survival event for each individual, 0 indicates alive and 1 indicates death.
#'Other choices are TRUE/FALSE (TRUE = death) or 1/2 (2=death). For interval censored data, the status indicator is 0=right censored, 1=event at time, 2=left censored, 3=interval censored.
#'@param write_name The name for pdf file which contains the result figure.
#'@param title_name The title for the result figure.
#'@param threshold Threshold used to indicate the high risk or low risk.
#'@param file If True, function will automatic generate a result pdf, otherwise it will
#'return a ggplot object. Default is FALSE.
#'
#'@return This function will generate a pdf file with 300dpi which compare survival curves using the Kaplan-Meier (KM) test.
#'@export
#'@examples 
#' data(demo_survival_data)
#' library(survival)
#' demo_set<-create_demo()
#' draw_survival_curve(demo_set[1,],
#'     living_days = demo_survival_data[,1],
#'     living_events =demo_survival_data[,2],
#'     write_name = "demo_data" )
#'
draw_survival_curve<-function(exp,living_days,living_events,write_name,title_name="",threshold=NA,file = FALSE){
  if (is.na(threshold)) threshold=median(exp)
  group_sur<-factor(1*(exp >threshold) ,
                    levels=0:1,
                    labels=c("low risk", "high risk"))
  pdata_gene <- list(time=living_days,
                     event=living_events,
                     gene=exp,
                     group_sur=group_sur)
  fit <- survival::survfit(survival::Surv(time, event) ~ group_sur, data = pdata_gene)
  data.survdiff <- survival::survdiff(survival::Surv(time, event) ~ group_sur, data = pdata_gene)
  p.val = 1 - stats::pchisq(data.survdiff$chisq, length(data.survdiff$n) - 1)
  if (p.val<0.00001){
    p_value="p value < 0.00001"
  }else{
    p_value=paste("p =",round(p.val,5),collapse = " ")
  }
  HR = (data.survdiff$obs[2]/data.survdiff$exp[2])/(data.survdiff$obs[1]/data.survdiff$exp[1])
  up95 = exp(log(HR) + qnorm(0.975)*sqrt(1/data.survdiff$exp[2]+1/data.survdiff$exp[1]))
  low95 = exp(log(HR) - qnorm(0.975)*sqrt(1/data.survdiff$exp[2]+1/data.survdiff$exp[1]))
  #require(ggfortify)
  gg<-survminer::ggsurvplot(fit,data = pdata_gene,  risk.table = TRUE,risk.table.y.text = TRUE,break.time.by=1,
                            fillLineSize=3,survSize=5,surv.linesize=5,pval = FALSE,
                            fontsize=2.5,
                            legend.labs = c("low risk","high risk"),xlab="Time (years)",
                        censor.alpha=0.5,palette=c("#00BFC4","#F8766D"),
                        title=paste(title_name ,"\n",p_value," HR ",round(HR,2),
                                         "(",round(up95,2) ,"-", round(low95,2),  ")",collapse = ""),
  ggtheme=ggplot2::theme(title =ggplot2::element_text(size=8, face='bold'),
          panel.grid.major =ggplot2::element_blank(), panel.grid.minor = ggplot2::element_blank(),
    panel.background = ggplot2::element_blank(),axis.line = ggplot2::element_line(colour = "black"),
    legend.position=c(1,1), legend.justification=c(1,1),
    legend.background = ggplot2::element_rect(fill = 'white', colour = 'white'))
  )
  if (file){
    ggplot2::ggsave(filename = paste("survival of ",write_name,".jpeg",sep=""),plot = print(gg),device ="jpeg" ,
                    path = getwd(),dpi = 300,units = "in",width = 5, height = 5,
                    limitsize=FALSE)
  }else{
    return(gg)
  }

}


is.integer0 <- function(x)
{
  is.integer(x) && length(x) == 0L
}

bs_ci <- function(data, indices, predict.time) { 
  d <- data[indices,] # allows boot to select sample
  d <- d[!is.na(d$survival),]
  surv.res = survivalROC::survivalROC.C(Stime = d$survival, 
                         status = d$survival_status, 
                         marker = d$marker, 
                         predict.time,
                         span = 0.25*NROW(d)^(-0.20))
  return(surv.res$AUC)
}

calculate_auc_ci <- function(survival, marker, predict_time,ci){
  if (ci){
    ci_res = boot::boot(data=data.frame(survival = survival[,1],
                                    survival_status = survival[,2],
                                    marker = marker),
                    statistic=bs_ci, R=1000, predict.time = predict_time)
    res = boot::boot.ci(ci_res,type="perc")
    return(list( AUC = res$t0, 
                 CI95 = paste(format(res$percent[,4], digits = 4),"-", format(res$percent[,5], digits = 4))
    ))
  }else{
    res = survivalROC::survivalROC.C(Stime = survival[,1], 
                               status = survival[,2], 
                               marker = marker, 
                               predict.time = predict_time,
                               span = 0.25*length(survival)^(-0.20))
    return(res)
  }

}

metricMCB.mean<-function(MCBset,MCB_matrix,Surv,data_set,show_bar=T){
  FunctionResults<-list()
  MCB_model_res<-NULL
  if (show_bar) {
    bar<-utils::txtProgressBar(min = 1,max = nrow(MCBset),char = "#",style = 3)
  }
  for (mcb in seq(nrow(MCBset))) {
    utils::setTxtProgressBar(bar, mcb)
    write_MCB<-rep(NA,5)
    #save the mcb number
    write_MCB[1]<-as.numeric(MCBset[mcb,'MCB_no'])
    write_MCB[2]<-MCBset[mcb,'CpGs']
    CpGs<-strsplit(MCBset[mcb,'CpGs']," ")[[1]]
    MCB_matrix[mcb,]<-colMeans(data_set[CpGs,])
    AUC_value<-survivalROC::survivalROC(Stime = Surv[,1],
                                        status = Surv[,2],
                                        marker = MCB_matrix[mcb,],
                                        predict.time = 5,
                                        method = "NNE",
                                        span =0.25*length(Surv)^(-0.20))$AUC
    write_MCB[3]<-AUC_value
    reg<-survival::survreg(Surv ~ MCB_matrix[mcb,])
    cindex<-survival::concordance(reg)
    write_MCB[4]<-cindex$concordance
    write_MCB[5]<-cindex$cvar
    MCB_model_res<-rbind(MCB_model_res,write_MCB)
  }
  colnames(MCB_model_res)<-c("MCB_no","CpGs","auc","C-index","C-index_SE")
  rownames(MCB_matrix)<-MCB_model_res[,'MCB_no']
  FunctionResults$MCB_matrix<-MCB_matrix
  FunctionResults$auc_results<-MCB_model_res
  return(FunctionResults)
}


ensemble_prediction.m<- function(ensemble_model,prediction_data) {
  prediction_data <- prediction_data[ensemble_model$cox$cox_model$CpGs,]
  if (nrow(prediction_data) != length(rownames(prediction_data))) {
    stop("ERROR: The predition data and the model have wrong dimensions!")
  }
  rank_svm <- stats::predict(ensemble_model$svm$svm_model, data.frame(t(prediction_data)))$predicted[1,]
  svm <- predict(ensemble_model$svm$hr_model,data.frame(rank_svm = rank_svm),type='lp')
  cox <-stats::predict(ensemble_model$cox$cox_model, data.frame(t(prediction_data)))
  enet <- stats::predict(ensemble_model$enet$enet_model, t(prediction_data), 
                         s = ensemble_model$enet$`corrected_lambda(min)`)
  mboost<-stats::predict(ensemble_model$mboost$mboost_model, t(prediction_data))[,1]
  data<-rbind(cox,
              svm,
              t(enet),
              mboost
  )
  rownames(data)<-c('cox','svm','enet','mboost')
  data<-t(data)
  data_f<-data.frame(cox = cox,svm=rank_svm,enet=as.numeric(enet),mboost = mboost)
  if (class(ensemble_model$stacking)[1] == "cv.glmnet")
    ensemble = as.numeric(stats::predict(ensemble_model$stacking, as.matrix(data_f)))
  else if (class(ensemble_model$stacking)[1] == "cph")
    ensemble = stats::predict(ensemble_model$stacking, data_f,type='lp')
  else stop("ERROR: The ensemble predition model is invaild !")
  return(t(cbind(data,ensemble)))
}

#inherent from survivalROC package to calculate ROC with small modification
survivalROC_C <- function (Stime, status, marker, entry = NULL, predict.time, 
                         cut.values = NULL, method = "NNE", lambda = NULL, span = NULL, 
                         window = "symmetric") 
{
  times = Stime
  x <- marker
  if (is.null(entry)) 
    entry <- rep(0, length(times))
  bad <- is.na(times) | is.na(status) | is.na(x) | is.na(entry)
  entry <- entry[!bad]
  times <- times[!bad]
  status <- status[!bad]
  x <- x[!bad]
  if (sum(bad) > 0) 
    cat(paste("\n", sum(bad), "records with missing values dropped. \n"))
  if (is.null(cut.values)) 
    cut.values <- unique(x)
  cut.values <- cut.values[order(cut.values)]
  ncuts <- length(cut.values)
  ooo <- order(times)
  times <- times[ooo]
  status <- status[ooo]
  x <- x[ooo]
  s0 <- 1
  unique.t0 <- unique(times)
  unique.t0 <- unique.t0[order(unique.t0)]
  n.times <- sum(unique.t0 <= predict.time)
  for (j in 1:n.times) {
    n <- sum(entry <= unique.t0[j] & times >= unique.t0[j])
    d <- sum((entry <= unique.t0[j]) & (times == unique.t0[j]) & 
               (status == 1))
    if (n > 0) 
      s0 <- s0 * (1 - d/n)
  }
  s.pooled <- s0
  roc.matrix <- matrix(NA, ncuts, 2)
  roc.matrix[ncuts, 1] <- 0
  roc.matrix[ncuts, 2] <- 1
  if (method == "KM") {
    for (c in 1:(ncuts - 1)) {
      s0 <- 1
      subset <- as.logical(x > cut.values[c])
      e0 <- entry[subset]
      t0 <- times[subset]
      c0 <- status[subset]
      if (!is.null(t0)) {
        unique.t0 <- unique(t0)
        unique.t0 <- unique.t0[order(unique.t0)]
        n.times <- sum(unique.t0 <= predict.time)
        if (n.times > 0) {
          for (j in 1:n.times) {
            n <- sum(e0 <= unique.t0[j] & t0 >= unique.t0[j])
            d <- sum((e0 <= unique.t0[j]) & (t0 == unique.t0[j]) & 
                       (c0 == 1))
            if (n > 0) 
              s0 <- s0 * (1 - d/n)
          }
        }
      }
      p0 <- mean(subset)
      roc.matrix[c, 1] <- (1 - s0) * p0/(1 - s.pooled)
      roc.matrix[c, 2] <- 1 - s0 * p0/s.pooled
    }
  }
  if (method == "NNE") {
    if (is.null(lambda) & is.null(span)) {
      cat("method = NNE requires either lambda or span! \n")
      stop(0)
    }
    x.unique <- unique(x)
    x.unique <- x.unique[order(x.unique)]
    S.t.x <- rep(0, length(x.unique))
    t.evaluate <- unique(times[status == 1])
    t.evaluate <- t.evaluate[order(t.evaluate)]
    t.evaluate <- t.evaluate[t.evaluate <= predict.time]
    for (j in 1:length(x.unique)) {
      if (!is.null(span)) {
        if (window == "symmetric") {
          ddd <- (x - x.unique[j])
          n <- length(x)
          ddd <- ddd[order(ddd)]
          index0 <- sum(ddd < 0) + 1
          index1 <- index0 + trunc(n * span + 0.5)
          if (index1 > n) 
            index1 <- n
          lambda <- ddd[index1]
          wt <- as.integer(((x - x.unique[j]) <= lambda) & 
                             ((x - x.unique[j]) >= 0))
          index0 <- sum(ddd <= 0)
          index2 <- index0 - trunc(n * span/2)
          if (index2 < 1) 
            index2 <- 1
          lambda <- abs(ddd[index2])
          set.index <- ((x - x.unique[j]) >= -lambda) & 
            ((x - x.unique[j]) <= 0)
          wt[set.index] <- 1
        }
        if (window == "asymmetric") {
          ddd <- (x - x.unique[j])
          n <- length(x)
          ddd <- ddd[order(ddd)]
          index0 <- sum(ddd < 0) + 1
          index <- index0 + trunc(n * span)
          if (index > n) 
            index <- n
          lambda <- ddd[index]
          wt <- as.integer(((x - x.unique[j]) <= lambda) & 
                             ((x - x.unique[j]) >= 0))
        }
      }
      else {
        wt <- exp(-(x - x.unique[j])^2/lambda^2)
      }
      s0 <- 1
      for (k in 1:length(t.evaluate)) {
        n <- sum(wt * (entry <= t.evaluate[k]) & (times >= 
                                                    t.evaluate[k]))
        d <- sum(wt * (entry <= t.evaluate[k]) & (times == 
                                                    t.evaluate[k]) * (status == 1))
        if (n > 0) 
          s0 <- s0 * (1 - d/n)
      }
      S.t.x[j] <- s0
    }
    S.all.x <- S.t.x[match(x, x.unique)]
    n <- length(times)
    S.marginal <- sum(S.all.x)/n
    for (c in 1:(ncuts - 1)) {
      p1 <- sum(x > cut.values[c])/n
      Sx <- sum(S.all.x[x > cut.values[c]])/n
      roc.matrix[c, 1] <- (p1 - Sx)/(1 - S.marginal)
      roc.matrix[c, 2] <- 1 - Sx/S.marginal
    }
  }
  sensitivity = roc.matrix[, 1]
  specificity = roc.matrix[, 2]
  x <- 1 - c(0, specificity)
  y <- c(1, sensitivity)
  n <- length(x)
  dx <- x[-n] - x[-1]
  mid.y <- (y[-n] + y[-1])/2
  area <- sum(dx * mid.y)
  list(cut.values = c(-Inf, cut.values), TP = y, FP = x, predict.time = predict.time, 
       Survival = s.pooled, AUC = area)
}

s_t_x<-function (Stime, status, marker, entry = NULL, predict.time) 
{
  span = 0.25*length(marker)^(-0.20)
  lambda = NULL
  times = Stime
  x <- marker
  if (is.null(entry)) 
    entry <- rep(0, length(times))
  bad <- is.na(times) | is.na(status) | is.na(x) | is.na(entry)
  entry <- entry[!bad]
  times <- times[!bad]
  status <- status[!bad]
  x <- x[!bad]
  if (sum(bad) > 0) 
    cat(paste("\n", sum(bad), "records with missing values dropped. \n"))
  cut.values <- unique(x)
  cut.values <- cut.values[order(cut.values)]
  ncuts <- length(cut.values)
  ooo <- order(times)
  times <- times[ooo]
  status <- status[ooo]
  x <- x[ooo]
  s0 <- 1
  unique.t0 <- unique(times)
  unique.t0 <- unique.t0[order(unique.t0)]
  n.times <- sum(unique.t0 <= predict.time)
  for (j in 1:n.times) {
    n <- sum(entry <= unique.t0[j] & times >= unique.t0[j])
    d <- sum((entry <= unique.t0[j]) & (times == unique.t0[j]) & 
               (status == 1))
    if (n > 0) 
      s0 <- s0 * (1 - d/n)
  }
  s.pooled <- s0
  roc.matrix <- matrix(NA, ncuts, 2)
  roc.matrix[ncuts, 1] <- 0
  roc.matrix[ncuts, 2] <- 1
  x.unique <- unique(x)
  x.unique <- x.unique[order(x.unique)]
  S.t.x <- rep(0, length(x.unique))
  t.evaluate <- unique(times[status == 1])
  t.evaluate <- t.evaluate[order(t.evaluate)]
  t.evaluate <- t.evaluate[t.evaluate <= predict.time]
  for (j in 1:length(x.unique)) {
    if (!is.null(span)) {
      ddd <- (x - x.unique[j])
      n <- length(x)
      ddd <- ddd[order(ddd)]
      index0 <- sum(ddd < 0) + 1
      index1 <- index0 + trunc(n * span + 0.5)
      if (index1 > n) 
        index1 <- n
      lambda <- ddd[index1]
      wt <- as.integer(((x - x.unique[j]) <= lambda) & 
                         ((x - x.unique[j]) >= 0))
      index0 <- sum(ddd <= 0)
      index2 <- index0 - trunc(n * span/2)
      if (index2 < 1) 
        index2 <- 1
      lambda <- abs(ddd[index2])
      set.index <- ((x - x.unique[j]) >= -lambda) & 
        ((x - x.unique[j]) <= 0)
      wt[set.index] <- 1
      
    }
    else {
      cat("requires a valid span! \n")
      stop(0)
    }
    s0 <- 1
    for (k in 1:length(t.evaluate)) {
      n <- sum(wt * (entry <= t.evaluate[k]) & (times >= 
                                                  t.evaluate[k]))
      d <- sum(wt * (entry <= t.evaluate[k]) & (times == 
                                                  t.evaluate[k]) * (status == 1))
      if (n > 0) 
        s0 <- s0 * (1 - d/n)
    }
    S.t.x[j] <- s0
  }
  S.all.x <- S.t.x[match(x, x.unique)]
  ks_res<-ks.test(x=S.all.x,y="exp",alternative = "t")
  list(S.all.x = S.all.x,ks_res = ks_res)
}

#error metrics -- Confusion Matrix
err_metric=function(CM,verbose = F)
{
  TP =CM[1,1]
  TN =CM[2,2]
  FP =CM[2,1]
  FN =CM[1,2]
  res = list()
  res$precision =TP / (TP + FP)
  res$recall_score =TP / (TP + FN)
  res$f1_score=2*((res$precision*res$recall_score)/(res$precision+res$recall_score))
  res$accuracy_model  = (TP+TN) / (TP+FN+FP+TN)
  res$False_positive_rate = (FP)/(FP+TN)
  res$False_negative_rate = (FN)/(FN+TP)
  if (verbose == T) {
    return(res)
  }else{
    print(paste("Precision value of the model: ",round(res$precision,2)))
    print(paste("Accuracy of the model: ",round(res$accuracy_model,2)))
    print(paste("Recall value of the model: ",round(res$recall_score,2)))
    print(paste("False Positive rate of the model: ",round(res$False_positive_rate,2)))
    print(paste("False Negative rate of the model: ",round(res$False_negative_rate,2)))
    print(paste("F1 score of the model: ",round(res$f1_score,2)))
  }
}

prepare_data_frame<-function(df_data){
  final_data<-NULL
  for (i in 1:ncol(df_data)) {
    one_data<-data.frame(samples= 1:nrow(df_data),
                         var_name=rep(colnames(df_data)[i],nrow(df_data)),
                         var_num=df_data[,i])
    final_data<-rbind(final_data,one_data)
  }
  final_data
}

get_the_loci <- function(CpGs,pos = TRUE){
  CpGs_t<-strsplit(CpGs, " ")[[1]]
  Illumina_Infinium_Human_Methylation_450K <- get450kAnno()
  met_cg_allgene <- Illumina_Infinium_Human_Methylation_450K[CpGs_t,]
  if (pos) met_cg_allgene[,'pos'] else met_cg_allgene
}

findAttractors<-function (methylation_matrix, annotation, class_vector, 
                          min.pwaysize = 5,base_method = c('Fstat','Tstat')[1], sec_method = c('ttest','kstest')[1] ,...) 
{
  
  dat.fr <- methylation_matrix
  all.probes <- rownames(dat.fr)
  dat.fr <- as.matrix(dat.fr)
  class.vector <- as.factor(class_vector)
  all_included <- strsplit(paste(annotation[,'CpGs'],collapse = ' ')," ")[[1]]
  probes.hits <- intersect(all.probes, all_included)
  dat.detect.w <- dat.fr[rownames(dat.fr) %in% probes.hits,]
  dat.detect.w <- as.matrix(dat.detect.w)
  incidence.matrix <- buildIncidenceMatrix(rownames(dat.detect.w), annotation)
  keep.mcb <- apply(incidence.matrix, 1, sum) >= min.pwaysize
  incidence.matrix <- incidence.matrix[keep.mcb,]
  new.order <- order(class.vector, colnames(dat.detect.w))
  dat.detect.w <- dat.detect.w[, new.order]
  class.vector <- class.vector[new.order]
  
  evalPway <- function(index, global,sec_method=c('ttest','kstest')[1]) {
    pway.vals <- global[index == 1]
    if (sec_method == 'ttest'){
      return(t.test(pway.vals, global)$p.value)
    }else if (sec_method == 'kstest'){
      return(ks.test(pway.vals, global,exact = F)$p.value)
    }
  }
  
  if (base_method=='Fstat'){
    fstat <- apply(dat.detect.w, 1, function(y, x) {
      anova(lm(y ~ x))[[4]][1]
    }, x = class.vector)
    fstat <- log(fstat, 2)
    t.pvals <- apply(incidence.matrix, 1, evalPway, global = fstat,sec_method=sec_method)
  }else if (base_method=='Tstat'){
    tstat <- apply(dat.detect.w, 1, function(y, x) {
      t.test(y ~ x)$statistic
    }, x = class.vector)
    tstat <- log(abs(tstat), 2)
    t.pvals <- apply(incidence.matrix, 1, evalPway, global = tstat,sec_method=sec_method)
  }else{
    stop('The method indicator base_method is invalid.')
  }
  t.pvals_adjust <- p.adjust(t.pvals, "BH")
  size <- apply(incidence.matrix, 1, sum)
  tab <- data.frame(MCBID = rownames(incidence.matrix),
                    annotation[sapply(X = rownames(incidence.matrix),FUN = function(x){strsplit(x, "_")[[1]][2]}),],
                    Pvalues = t.pvals,
                    AdjustedPvalues = t.pvals_adjust, 
                    NumberDetectedCpGs = size)
  #rownames(tab)<-sapply(X = rownames(tab),FUN = function(x){strsplit(x, "_")[[1]][2]})
  tab <- tab[order(t.pvals), ]
  return(tab)
}

buildIncidenceMatrix<-function (gene.ids, annotation) 
{
  CpGs_names<-annotation[,'CpGs']
  pway.genes <- strsplit(CpGs_names," ")
  convert.to.row <- function(x.genes, row.genes) {
    res <- integer(length(row.genes))
    res[row.genes %in% x.genes] <- 1
    return(res)
  }
  xmat <- t(sapply(lapply(pway.genes, convert.to.row, gene.ids), 
                   cbind))
  rownames(xmat) <- paste('MCB',annotation[,'MCB_no'],sep  = "_") 
  colnames(xmat) <- gene.ids
  return(xmat)
}

re_AnnotatedMMB<-function(MMB_list){
  Illumina_Infinium_Human_Methylation_450K <- get450kAnno()
  met_cg_allgene <- Illumina_Infinium_Human_Methylation_450K
  res<-rep(NA,7)
  names(res)<-c('chromosomes','start_site','end_site','CpGs','location','length','CpGs_num')
  res_all <-NULL
  for (MMB in MMB_list) {
    MMBannotated<-met_cg_allgene[MMB,]
    res['chromosomes'] <- unique(MMBannotated$chr)
    res['start_site']  <- min(MMBannotated$pos)
    res['end_site']  <- max(MMBannotated$pos)
    res['CpGs'] <-paste(MMBannotated$Name[order(MMBannotated$pos,decreasing = F)],collapse = " ")
    res['location'] <-paste(res['chromosomes'],":",res['start_site'],'-',res['chromosomes'],':',res['end_site'],sep = " ")
    res['length'] <- as.character(max(MMBannotated$pos)- min(MMBannotated$pos) )
    res['CpGs_num'] <- nrow(MMBannotated)
    res_all<-rbind(res_all,res)
  }
  rownames(res_all)<-1:nrow(res_all)
  res_all
}

bmeasures_otu<-
  function (x, y) 
  {
    dx <- dim(x)
    dy <- dim(y)
    if (dx[1] != dy[1]) {
      cat("Inputs have different number of rows \n")
    }
    else if (dx[2] != dy[2]) {
      cat("Inputs have different number of columns \n")
    }
    else {
      cat("")
    }
    input <- as.matrix(rbind(x, y))
    a <- 0
    b <- 0
    c <- 0
    d <- 0
    for (i in 1:dx[2]) {
      k <- input[1, i]
      l <- input[2, i]
      if ((k == 1) & (l == 1)) {
        a <- a + 1
      }
      else if ((k == 0) & (l == 1)) {
        b <- b + 1
      }
      else if ((k == 1) & (l == 0)) {
        c <- c + 1
      }
      else {
        d <- d + 1
      }
    }
    otu_out <- c(a = a, b = b, c = c, d = d)
    remove(a, b, c, d)
    return(otu_out)
  }

bmeasures <-
  function (x, y, method = 'eq_01') 
  {
    x<-matrix(x,1)
    y<-matrix(y,1)
    otuTable <- bmeasures_otu(x, y)
    a <- otuTable[1]
    b <- otuTable[2]
    c <- otuTable[3]
    d <- otuTable[4]
    n <- a + b + c + d
    switch(method, eq_01 = {
      coef <- a/(a + b + c)
    }, eq_02 = {
      coef <- a/(2 * a + b + c)
    }, eq_03 = {
      coef <- (2 * a)/(2 * a + b + c)
    }, eq_04 = {
      coef <- (3 * a)/(3 * a + b + c)
    }, eq_05 = {
      coef <- (2 * a)/((a + b) + (a + c))
    }, eq_06 = {
      coef <- a/(a + (2 * b) + (2 * c))
    }, eq_07 = {
      coef <- (a + d)/n
    }, eq_08 = {
      coef <- (2 * (a + d))/((2 * a) + b + c + (2 * d))
    }, eq_09 = {
      coef <- (a + d)/(a + (2 * (b + c)) + d)
    }, eq_10 = {
      coef <- (a + (0.5 * d))/n
    }, eq_11 = {
      coef <- (a + d)/(a + (0.5 * (b + c)) + d)
    }, eq_12 = {
      coef <- a
    }, eq_13 = {
      coef <- a + d
    }, eq_14 = {
      coef <- a/n
    }, eq_15 = {
      coef <- b + c
    }, eq_16 = {
      coef <- sqrt(b + c)
    }, eq_17 = {
      coef <- sqrt((b + c)^2)
    }, eq_18 = {
      coef <- sqrt((b + c)^2)
    }, eq_19 = {
      coef <- b + c
    }, eq_20 = {
      coef <- (b + c)/n
    }, eq_21 = {
      coef <- b + c
    }, eq_22 = {
      coef <- b + c
    }, eq_23 = {
      coef <- (b + c)/(4 * n)
    }, eq_24 = {
      coef <- ((b + c)^2)/(n^2)
    }, eq_25 = {
      coef <- ((n * (b + c)) - ((b - c)^2))/(n^2)
    }, eq_26 = {
      coef <- (4 * b * c)/(n^2)
    }, eq_27 = {
      coef <- (b + c)/(2 * a + b + c)
    }, eq_28 = {
      coef <- (b + c)/(2 * a + b + c)
    }, eq_29 = {
      coef <- 2 * (sqrt(1 - (a/sqrt((a + b) * (a + c)))))
    }, eq_30 = {
      coef <- (sqrt(2 * (1 - (a/(sqrt((a + b) * (a + c)))))))
    }, eq_31 = {
      coef <- a/(sqrt((a + b) * (a + c)))
    }, eq_32 = {
      coef <- log(a) - log(n) - log((a + b)/n) - log((a + c)/n)
    }, eq_33 = {
      coef <- a/(sqrt((a + b) * (a + c)))
    }, eq_34 = {
      coef <- (n * a)/((a + b) * (a + c))
    }, eq_35 = {
      coef <- (n * ((a - 0.5)^2))/((a + b) * (a + c))
    }, eq_36 = {
      coef <- (a^2)/((a + b) * (a + c))
    }, eq_37 = {
      coef <- a/(0.5 * (a * b + a * c) + (b * c))
    }, eq_38 = {
      coef <- a/(((a + b) * (a + c))^0.5)
    }, eq_39 = {
      coef <- ((a^2) - (b * c))/((a + b) * (a + c))
    }, eq_40 = {
      coef <- ((n * a) - (a + b) * (a + c))/((n * a) + (a + 
                                                          b) * (a + c))
    }, eq_41 = {
      coef <- ((a/2) * (2 * a + b + c))/((a + b) * (a + c))
    }, eq_42 = {
      coef <- (a/2) * ((1/(a + b)) + (1/(a + c)))
    }, eq_43 = {
      coef <- (a/(a + b)) + (a/(a + c))
    }, eq_44 = {
      coef <- ((a * d) - (b * c))/(sqrt(n * (a + b) * (a + 
                                                         c)))
    }, eq_45 = {
      coef <- a/(min((a + b), (a + c)))
    }, eq_46 = {
      coef <- a/(max((a + b), (a + c)))
    }, eq_47 = {
      coef <- (a/sqrt((a + b) * (a + c))) - (max((a + b), (a + 
                                                             c))/2)
    }, eq_48 = {
      coef <- ((n * a) - ((a + b) * (a + c)))/((n * min(a + 
                                                          b, a + c)) - (a + b) * (a + c))
    }, eq_49 = {
      coef <- 0.25 * ((a/(a + b)) + (a/(a + c)) + (d/(b + d)) + 
                        (d/(c + d)))
    }, eq_50 = {
      coef <- (a + d)/(sqrt((a + b) * (a + c) * (b + d) * (c + 
                                                             d)))
    }, eq_51 = {
      x2 <- (n * (((a * d) - (b * c))^2))/((a + b) * (a + c) * 
                                             (c + d) * (b + d))
      coef <- x2
    }, eq_52 = {
      x2 <- (n * (((a * d) - (b * c))^2))/((a + b) * (a + c) * 
                                             (c + d) * (b + d))
      coef <- sqrt(x2/(n + x2))
    }, eq_53 = {
      p <- ((a * d) - (b * c))/(sqrt((a + b) * (a + c) * (b + 
                                                            d) * (c + d)))
      coef <- sqrt(p/(n + p))
    }, eq_54 = {
      p <- ((a * d) - (b * c))/(sqrt((a + b) * (a + c) * (b + 
                                                            d) * (c + d)))
      coef <- p
    }, eq_55 = {
      coef <- cos((pi * sqrt(b * c))/(sqrt(a * d) + sqrt(b * 
                                                           c)))
    }, eq_56 = {
      coef <- (a + d)/(b + c)
    }, eq_57 = {
      coef <- (a * d)/(((a + b) * (a + c) * (b + d) * (c + 
                                                         d))^0.5)
    }, eq_58 = {
      coef <- (sqrt(2) * (a * d - b * c))/(sqrt(((a * d - b * 
                                                    c)^2) - (a + b) * (a + c) * (b + d) * (c + d)))
    }, eq_59 = {
      coef <- log10((n * (((abs(a * d - b * c)) - n/2)^2))/((a + 
                                                               b) * (a + c) * (b + d) * (c + d)))
    }, eq_60 = {
      coef <- (a * d)/(sqrt((a + b) * (a + c) * (b + d) * (c + 
                                                             d)))
    }, eq_61 = {
      coef <- (a * d - b * c)/(a * d + b * c)
    }, eq_62 = {
      coef <- ((2 * b * c)/(a * d + b * c))
    }, eq_63 = {
      coef <- (sqrt(a * d) - sqrt(b * c))/(sqrt(a * d) + sqrt(b * 
                                                                c))
    }, eq_64 = {
      coef <- a/(b + c)
    }, eq_65 = {
      coef <- a/((a + b) + (a + c) - a)
    }, eq_66 = {
      coef <- (a * d - b * c)/(n^2)
    }, eq_67 = {
      coef <- ((a + d) - (b + c))/n
    }, eq_68 = {
      coef <- (4 * (a * d - b * c))/(((a + d)^2) + ((b + c)^2))
    }, eq_69 = {
      sig <- max(a, b) + max(c, d) + max(a, c) + max(b, d)
      sigt <- max(a + c, b + d) + max(a + b, c + d)
      coef <- (sig - sigt)/(2 * n - sigt)
    }, eq_70 = {
      sig <- max(a, b) + max(c, d) + max(a, c) + max(b, d)
      sigt <- max(a + c, b + d) + max(a + b, c + d)
      coef <- (sig - sigt)/(2 * n)
    }, eq_71 = {
      coef <- (sqrt(a * d) + a)/(sqrt(a * d) + a + b + c)
    }, eq_72 = {
      coef <- (sqrt(a * d) + a - (b + c))/(sqrt(a * d) + a + 
                                             b + c)
    }, eq_73 = {
      coef <- (a * b + b * c)/((a * b) + (2 * b * c) + (c * 
                                                          d))
    }, eq_74 = {
      coef <- ((n^2) * (n * a - (a + b) * (a + c)))/((a + b) * 
                                                       (a + c) * (b + d) * (c + d))
    }, eq_75 = {
      coef <- (a * (c + d))/(c * (a + b))
    }, eq_76 = {
      coef <- abs((a * (c + d))/(c * (a + b)))
    }, eq_77 = {
      coef <- log(1 + a)/log(1 + n)
    }, eq_78 = {
      coef <- log(1 + a)/log(1 + a + b + c)
    }, eq_79 = {
      coef <- (log(1 + a * d) - log(1 + b * c))/log(1 + (n^2)/4)
    }, {
      cat("No desired equation. Please check it again.")
    })
    if (is.na(coef)) coef = 1
    result <- c(a, b, c, d, coef)
    result <- t(as.matrix(result))
    colnames(result) <- c("a", "b", "c", "d", 
                          "coef")
    return(coef)
  }

vector_to_matrix<-function(x){
  list_data<-x[grep('X',names(x))]
  matrix_data <-unlist((lapply(list_data,function(x){strsplit(x," ")[[1]]})))
  mat<-as.numeric_matrix(matrix(matrix_data,nrow = length(matrix_data)/14,ncol = 14))
  rownames(mat)<-strsplit(as.character(x['CpGs']),' ')[[1]]
  mat
}



# end load
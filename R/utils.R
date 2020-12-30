# env parameter
# platform       x86_64-redhat-linux-gnu
# arch           x86_64
# os             linux-gnu
# system         x86_64, linux-gnu
# status
# major          3
# minor          5.0
# year           2018
# month          04
# day            23
# svn rev        74626
# language       R
# version.string R version 3.5.0 (2018-04-23)
# nickname       Joy in Playing

# load private functions

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
#inherent from minifi package to reduce the dependencies
# getAnnotation <- function(object, what = "everything", lociNames = NULL,
#                           orderByLocation = FALSE, dropNonMapping = FALSE) {
#   .annoGet <- function(what, envir) {
#     pointer <- get(what, envir = envir)
#     if (is(pointer, "DataFrame")) return(pointer)
#     get(pointer$what, envir = as.environment(pointer$envir))
#   }
#   getAnnotationObject <- function(object) {
#     if (is.character(object)) {
#       if (!require(object, character.only = TRUE)) {
#         stop(sprintf("cannot load annotation package %s", object))
#       }
#       object <- get(object)
#     }
#     if (!is(object, "IlluminaMethylationAnnotation")) {
#       stop("Could not locate annotation object for 'object' of class",
#            class(object))
#     }
#     object
#   }
#   .availableAnnotation <- function(object) {
#     object <- getAnnotationObject(object)
#     allAnnoNames <- ls(object@data)
#     annoClasses <- sub("\\..*", "", allAnnoNames)
#     annoClassesChoices <- sub(".*\\.", "", allAnnoNames)
#     annoClassesChoices[grep("\\.", allAnnoNames, invert = TRUE)] <- ""
#     annoClassesChoices <- split(annoClassesChoices, annoClasses)
#     annoClasses <- unique(annoClasses)
#     defaults <- object@defaults
#     out <- list(
#       names = allAnnoNames,
#       classes = annoClasses,
#       classChoices = annoClassesChoices,
#       defaults = defaults)
#     out
#   }
#   # Processing of arguments and check
#   annoObject <- getAnnotationObject(object)
#   available <- .availableAnnotation(annoObject)
#   if ("everything" %in% what)  what <- available$defaults
#   if (!(all(what %in% available$names))) {
#     stop("the value of argument 'what' is not part of the annotation ",
#          "package or 'everything'")
#   }
#   n_choices <- vapply(
#     available$classes,
#     function(cl) length(grep(cl, what)),
#     integer(1))
#   if (any(n_choices > 1)) stop("only one choice per allowable class")
#   if (!any(grepl("^Locations", what)) &&
#       (orderByLocation || dropNonMapping)) {
#     stop("To use 'orderbyLocation' or 'dropNonMapping' specify Locations ",
#          "as part of 'what'")
#   }
#   # TODO: Ensure same order always
#   # Old code for inspiration
#   # bestOrder <- c("Locations", "Manifest", "IlluminaSNPs", "Annotation")
#   # what <- bestOrder[bestOrder %in% what]
#   
#   out <- do.call(cbind, lapply(what, function(wh) {
#     .annoGet(wh, envir = annoObject@data)
#   }))
#   
#   if (!is.null(lociNames)) {
#     lociNames <- lociNames[lociNames %in% rownames(out)]
#   }
#   if (is(object, "MethylSet") || is(object, "RatioSet") ||
#       is(object, "GenomicMethylSet") || is(object, "GenomicRatioSet")) {
#     rNames <- rownames(object)
#     if (is.null(lociNames)) {
#       lociNames <- rNames[rNames %in% rownames(out)]
#     } else {
#       lociNames <- rNames[rNames %in% lociNames]
#     }
#   }
#   if (!is.null(lociNames)) out <- out[match(lociNames, rownames(out)),]
#   out
# }

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
                         lambda = NULL,
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
                                        lambda = NULL,
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
                         lambda = NULL,
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
    cindex<-survival::survConcordance(Surv ~ MCB_matrix[mcb,])
    write_MCB[4]<-cindex$concordance
    write_MCB[5]<-cindex$std.err
    MCB_model_res<-rbind(MCB_model_res,write_MCB)
  }
  colnames(MCB_model_res)<-c("MCB_no","CpGs","auc","C-index","C-index_SE")
  rownames(MCB_matrix)<-MCB_model_res[,'MCB_no']
  FunctionResults$MCB_matrix<-MCB_matrix
  FunctionResults$auc_results<-MCB_model_res
  return(FunctionResults)
}


ensemble_prediction.m<- function(ensemble_model,predition_data) {
  predition_data<-predition_data[rownames(predition_data) %in% ensemble_model$cox$cox_model$CpGs,]
  if (nrow(predition_data)!=length(rownames(predition_data))) {
    stop("ERROR: The predition data and the model have wrong dimensions!")
  }
  svm<-stats::predict(ensemble_model$svm$svm_model, data.frame(t(predition_data)))$predicted
  cox<-stats::predict(ensemble_model$cox$cox_model, data.frame(t(predition_data)))
  enet<-stats::predict(ensemble_model$enet$enet_model,t(predition_data),s=ensemble_model$enet$`corrected_lambda(min)`)
  coxboost<-stats::predict(ensemble_model$coxboost$coxboost_model, t(predition_data))[,1]
  data<-rbind(cox,
              svm,
              t(enet),
              coxboost
  )
  rownames(data)<-c('cox','svm','enet','coxboost')
  data<-t(data)
  data_f<-as.data.frame(data)
  ensemble = stats::predict(ensemble_model$stacking, data_f)
  return(t(cbind(data,ensemble)))
}
# end load
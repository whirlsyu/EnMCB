#' @title  Identification of methylation correlated blocks
#'
#' @description This function is used to partition the genome into blocks of tightly co-methylated CpG sites,
#' Methylation correlated blocks. This function calculates Pearson correlation coefficients r^2 between
#' the beta values of any two CpGs r^2 < CorrelationThreshold was used to identify boundaries between any two
#' adjacent markers indicating uncorrelated methylation. Markers not separated by a boundary were combined into MCB. Pearson correlation coefficients between
#' two adjacent CpGs were calculated.
#'
#' @details Currently, only illumina 450k platform is supported, the methylation profile need to convert into matrix format.
#'
#' @param MethylationProfile Methylation matrix is used in the analysis.
#' @param CorrelationThreshold coef correlation threshold is used for define boundaries.
#' @param method method used for calculation of correlation, should be one of "pearson","spearman","kendall". Defualt is "pearson".
#' @param PositionGap CpG Gap between any two CpGs positioned CpG sites less than 1000 bp (default) will be calculated.
#' @param platform This parameter indicates the platform used to produce the methlyation profile.
#'
#' @author Xin Yu
#' @return
#' Object of class \code{list} with elements:
#'  \tabular{ll}{
#'    \code{MCBsites} \tab Character set contains all CpG sites in MCBs. \cr
#'    \code{MCBinformation} \tab Matrix contains the information of results. \cr
#'  }
#' @examples
#' #Not run! Remove # to run the example.
#' # require(EnMCB)
#' # data(demo_set)
#' #import exprs function
#' # library(SummarizedExperiment)
#' #import the demo TCGA data with 10000+ CpGs site and 455 samples
#' # res<-IdentifyMCB(assays(demo_set)[[1]])
#' # demo_MCBinformation<-res$MCBinformation
#'
#'
#' @export
#'
#' @references
#' Xin Yu et al. 2019 Predicting disease progression in lung adenocarcinoma patients based on methylation correlated blocks using ensemble machine learning classifiers (under review)
#'
#'
IdentifyMCB<-function(
  MethylationProfile,
  method=c("pearson","spearman","kendall")[1],
  CorrelationThreshold = 0.8,
  PositionGap = 1000,
  platform = "Illumina Methylation 450K"

){
  if (!(method %in% c("pearson","spearman","kendall"))) {
    stop(paste("Correlation method should be one of pearson, spearman and kendall."))
  }
  cat("Start calculating the correlation, this may take a while...\n")
  FunctionResults<-list()
  utils::data(Illumina_Infinium_Human_Methylation_450K,package = "EnMCB")
  met_cg_allgene<-Illumina_Infinium_Human_Methylation_450K[rownames(MethylationProfile),]
  chromosomes<-unique(met_cg_allgene[,'Chromosome'])
  chromosomes<-chromosomes[order(chromosomes)]
  res=NULL
  correlation_res<-NULL
  cat("(or you can try to use IdentifyMCB_parallel function instead)\n")
  if (length(chromosomes)>1) {
    bar<-utils::txtProgressBar(min = 1,max = length(chromosomes),char = "#",style = 3)
  }
  for(chr_no in 1:length(chromosomes)){
    if (length(chromosomes)>1) {
      utils::setTxtProgressBar(bar, chr_no)
    }
    chr_id<-chromosomes[chr_no]
    met_x<-MethylationProfile
    met_matrix<-met_x[met_cg_allgene[,'Chromosome'] %in% chr_id,]
    ann_matrix<-met_cg_allgene[met_cg_allgene[,'Chromosome'] %in% chr_id,]
    met_matrix<-met_matrix[order(as.numeric(ann_matrix[,'Start']),decreasing = F),]
    ann_matrix<-ann_matrix[order(as.numeric(ann_matrix[,'Start']),decreasing = F),]
    res<-NULL
    total<-nrow(met_matrix)
    for (i in 1:total) {
      # To investigate whether this indeed is evident in our data, we calculated Pearson
      # correlation coefficients r2 between beta values of any two CpGs positioned within
      # one kilobase (or indicated by PositionGap) of one another
      if (i+1<=total){
        if(as.numeric(ann_matrix[i+1,'Start'])-as.numeric(ann_matrix[i,'Start'])<PositionGap &
           ann_matrix[i+1,'Chromosome']==ann_matrix[i,'Chromosome']){
          res<-rbind(res,c(unlist(stats::cor.test(met_matrix[i,],met_matrix[i+1,],method = method))[1:5]))
        }else{
          res<-rbind(res,matrix(c(0,0,1,0,0),1,5))
        }
      }else{
        res<-rbind(res,matrix(c(0,0,1,0,0),1,5))
      }
    }
    rownames(res)<-rownames(met_matrix)
    correlation_res<-rbind(correlation_res,res)
  }
  cat("\n")
  cat("Now gathering the results, please wait ...\n")
  met_cg_allgene<-met_cg_allgene[rownames(correlation_res),]
  MCB_flag<-rep("boundary",nrow(correlation_res))
  MCB_flag[as.numeric(correlation_res[,'estimate.cor'])>CorrelationThreshold]<-"MCB"
  # This CpGs and next one are included.
  MCBsites<-union(grep("MCB",MCB_flag),grep("MCB",MCB_flag)+1)
  MCBsites<-MCBsites[order(MCBsites)]
  FunctionResults$MCBsites<-rownames(MethylationProfile)[MCBsites]

  MCB<-rep(NA,times=7)
  names(MCB)<-c("MCB_no","start","end","CpGs","location","chromosomes","length")
  MCB_block=F
  MCB_no=1
  total_res<-NULL
  for (i in 1:length(MCB_flag)) {
    flag<-MCB_flag[i]
    if (MCB_block==F&flag=="MCB"){
      MCB["start"] <- i
      MCB["MCB_no"] <- MCB_no
      MCB_no=MCB_no+1
      MCB_block = T
    }
    if (MCB_block==T&flag=="boundary"){
      MCB["end"] <- i
      MCB["CpGs"]<- paste(rownames(correlation_res)[as.numeric(MCB["start"]):(as.numeric(MCB["end"]))],collapse = " ")
      MCB["location"]<-paste(met_cg_allgene[as.numeric(MCB["start"]),'Chromosome'],":",
                                   met_cg_allgene[as.numeric(MCB["start"]),'Start'],"-",
                                   met_cg_allgene[as.numeric(MCB["end"]),'Chromosome'],":",
                                   met_cg_allgene[as.numeric(MCB["end"]),'End'],collapse = "")
      MCB["chromosomes"]<-met_cg_allgene[as.numeric(MCB["start"]),'Chromosome']
      MCB["length"]<-as.numeric(met_cg_allgene[as.numeric(MCB["end"]),'End'])-
        as.numeric(met_cg_allgene[as.numeric(MCB["start"]),'Start'])
      total_res<-rbind(total_res,MCB)
      MCB_block = F
    }
  }
  total_res<-cbind(total_res,CpGs_num=as.numeric(total_res[,'end'])-as.numeric(total_res[,'start'])+1)
  rownames(total_res)<-total_res[,'MCB_no']
  #Some of 'cgxxxxx' code point to same CpG, those cg-code are removed.
  total_res<-total_res[as.numeric(total_res[,'length'])>1,]
  cat("Statistics (",nrow(total_res)," MCBs in total):\n")
  for (chr_set in unique(total_res[,'chromosomes'])) {
    cat(chr_set,": ")
    cat("total MCBs:",length(total_res[total_res[,'chromosomes'] %in% chr_set,'chromosomes'])," ")
    cat("Mean Length:",mean(as.numeric(total_res[total_res[,'chromosomes'] %in% chr_set,'length']))," ")
    cat("(Range: ",range(as.numeric(total_res[total_res[,'chromosomes'] %in% chr_set,'length'])),")\n")
  }
  FunctionResults$MCBinformation<-total_res
  return(FunctionResults)
}

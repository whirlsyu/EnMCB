#' @title  Identification of methylation correlated blocks
#'
#' @description This function is used to partition the genome into blocks of tightly co-methylated CpG sites, \cr
#' Methylation correlated blocks. This function calculates Pearson correlation coefficients between \cr
#' the beta values of any two CpGs < CorrelationThreshold was used to identify boundaries between any two \cr
#' adjacent markers indicating uncorrelated methylation. Markers not separated by a boundary were combined into MCB. Pearson correlation coefficients between \cr
#' two adjacent CpGs were calculated.
#'
#' @details Currently, only illumina 450k platform is supported, the methylation profile need to convert into matrix format.
#'
#' @param MethylationProfile Methylation matrix is used in the analysis.
#' @param CorrelationThreshold coef correlation threshold is used for define boundaries.
#' @param method method used for calculation of correlation, \cr
#' should be one of "pearson","spearman","kendall". Defualt is "pearson".
#' @param PositionGap CpG Gap between any two CpGs positioned CpG sites less than 1000 bp (default) will be calculated.
#' @param platform This parameter indicates the platform used to produce the methlyation profile. 
#' You can use your own annotation file.
#' @param verbose True as default, which will print the block information for each chromosome. 
#'
#' @author Xin Yu
#' @return
#' Object of class \code{list} with elements:
#'  \tabular{ll}{
#'    \code{MCBsites} \tab Character set contains all CpG sites in MCBs. \cr
#'    \code{MCBinformation} \tab Matrix contains the information of results. \cr
#'  }
#' @examples
#' data('demo_data',package = "EnMCB")
#' 
#' #import the demo TCGA data with 10000+ CpGs site and 455 samples
#' #remove # to run
#' res<-IdentifyMCB(demo_data$realdata)
#' demo_MCBinformation<-res$MCBinformation
#'
#'
#' @export
#'
#' @references
#' Xin Yu, De-Xin Kong, EnMCB: an R/bioconductor package for predicting disease progression based on methylation correlated blocks using ensemble models, Bioinformatics, 2021, btab415
#'
IdentifyMCB<-function(
  MethylationProfile,
  method=c("pearson","spearman","kendall")[1],
  CorrelationThreshold = 0.8,
  PositionGap = 1000,
  platform = "Illumina Methylation 450K",
  verbose = T
){
  if (!(method %in% c("pearson","spearman","kendall"))) {
    stop(paste("Correlation method should be one of pearson, spearman and kendall."))
  }
  cat("Start calculating the correlation, this may take a while...\n")
  FunctionResults<-list()
  if (platform == "Illumina Methylation 450K"){
    Illumina_Infinium_Human_Methylation_450K<-get450kAnno()
    Illumina_Infinium_Human_Methylation_450K<-Illumina_Infinium_Human_Methylation_450K[!is.na(Illumina_Infinium_Human_Methylation_450K[,'pos']),]
    
    intersect_cpg<-intersect(rownames(Illumina_Infinium_Human_Methylation_450K),rownames(MethylationProfile))
    
    met_cg_allgene<-Illumina_Infinium_Human_Methylation_450K[intersect_cpg,]
  }else{
    met_cg_allgene = platform
  }

  MethylationProfile<-MethylationProfile[intersect_cpg,]
  
  chromosomes<-unique(met_cg_allgene[,'chr'])
  chromosomes<-chromosomes[order(chromosomes)]
  res=NULL
  correlation_res<-NULL
  cat("(or you can try to use IdentifyMCB_parallel function instead)\n")
  if (length(chromosomes)>1) {
    bar<-utils::txtProgressBar(min = 1,max = length(chromosomes),char = "#",style = 3)
  }
  for(chr_no in seq_along(chromosomes)){
    if (length(chromosomes)>1) {
      utils::setTxtProgressBar(bar, chr_no)
    }
    chr_id<-chromosomes[chr_no]
    met_x<-MethylationProfile
    if (sum(met_cg_allgene[,'chr'] %in% chr_id)<=2) {
      next
    }
    met_matrix<-met_x[met_cg_allgene[,'chr'] %in% chr_id,]
    ann_matrix<-met_cg_allgene[met_cg_allgene[,'chr'] %in% chr_id,]
    met_matrix<-met_matrix[order(as.numeric(ann_matrix[,'pos']),decreasing = FALSE),]
    ann_matrix<-ann_matrix[order(as.numeric(ann_matrix[,'pos']),decreasing = FALSE),]
    res<-NULL
    total<-nrow(met_matrix)
    for (i in seq_len(total)) {
      # To investigate whether this indeed is evident in our data, we calculated Pearson
      # correlation coefficients between beta values of any two CpGs positioned within
      # one kilobase (or indicated by PositionGap) of one another
      if (i+1<=total){
        if(as.numeric(ann_matrix[i+1,'pos'])-as.numeric(ann_matrix[i,'pos'])<PositionGap &
           ann_matrix[i+1,'chr']==ann_matrix[i,'chr']){
          res<-rbind(res,c(unlist(stats::cor.test(met_matrix[i,],met_matrix[i+1,],method = method))[1:5]))
        }else{
          if (method == "pearson") res<-rbind(res,matrix(c(0,0,1,0,0),1,5))
          else if (method == "spearman") res<-rbind(res,matrix(c(0,1,0,0,0),1,5))
          else if (method == "kendall") res<-rbind(res,matrix(c(0,1,0,0,0),1,5))
        }
      }else{
        if (method == "pearson") res<-rbind(res,matrix(c(0,0,1,0,0),1,5))
        else if (method == "spearman") res<-rbind(res,matrix(c(0,1,0,0,0),1,5))
        else if (method == "kendall") res<-rbind(res,matrix(c(0,1,0,0,0),1,5))
      }
    }
    rownames(res)<-rownames(met_matrix)
    correlation_res<-rbind(correlation_res,res)
  }
  if (verbose) {
    cat("\n")
    cat("Now gathering the results, please wait ...\n")
  }
  met_cg_allgene<-met_cg_allgene[rownames(correlation_res),]
  MCB_flag<-rep("boundary",nrow(correlation_res))
  if (method == "pearson") MCB_flag[as.numeric(correlation_res[,'estimate.cor'])>CorrelationThreshold]<-"MCB"
  else if (method == "spearman") MCB_flag[as.numeric(correlation_res[,'estimate.rho'])>CorrelationThreshold]<-"MCB"
  else if (method == "kendall") MCB_flag[as.numeric(correlation_res[,'estimate.tau'])>CorrelationThreshold]<-"MCB"
  # This CpGs and next one are included.
  MCBsites<-union(grep("MCB",MCB_flag),grep("MCB",MCB_flag)+1)
  MCBsites<-MCBsites[order(MCBsites)]
  FunctionResults$MCBsites<-rownames(MethylationProfile)[MCBsites]

  MCB<-rep(NA,times=10)
  names(MCB)<-c("MCB_no","start","end","CpGs","location","chromosomes",
                "length","MCB_Gene","Feature_Type","CGI_Coordinate")
  MCB_block=FALSE
  MCB_no=1
  total_res<-NULL
  for (i in seq_along(MCB_flag)) {
    flag<-MCB_flag[i]
    if (MCB_block==FALSE&flag=="MCB"){
      MCB["start"] <- i
      MCB["MCB_no"] <- MCB_no
      MCB_no=MCB_no+1
      MCB_block = TRUE
    }
    if (MCB_block==TRUE&flag=="boundary"){
      MCB["end"] <- i
      CpG_names<-rownames(correlation_res)[as.numeric(MCB["start"]):(as.numeric(MCB["end"]))]
      MCB["CpGs"]<- paste(CpG_names,collapse = " ")
      MCB["location"]<-paste(met_cg_allgene[as.numeric(MCB["start"]),'chr'],":",
                                   met_cg_allgene[as.numeric(MCB["start"]),'pos'],"-",
                                   met_cg_allgene[as.numeric(MCB["end"]),'chr'],":",
                                   met_cg_allgene[as.numeric(MCB["end"]),'pos'],collapse = "")
      MCB["chromosomes"]<-met_cg_allgene[as.numeric(MCB["start"]),'chr']
      MCB["MCB_Gene"]<-paste(unique(strsplit(paste(met_cg_allgene[CpG_names,'UCSC_RefGene_Name'],collapse = ";"),";")[[1]]),
                               collapse = " ")
      MCB["Feature_Type"]<-paste(unique(strsplit(paste(met_cg_allgene[CpG_names,'Relation_to_Island'],collapse = ";"),";")[[1]]),
                                   collapse = " ")
      MCB["CGI_Coordinate"]<-paste(unique(met_cg_allgene[CpG_names,'Islands_Name']),
                                     collapse = ";")
      MCB["length"]<-as.numeric(met_cg_allgene[as.numeric(MCB["end"]),'pos'])-
        as.numeric(met_cg_allgene[as.numeric(MCB["start"]),'pos'])
      total_res<-rbind(total_res,MCB)
      MCB_block = FALSE
    }
  }
  total_res<-cbind(total_res,CpGs_num=as.numeric(total_res[,'end'])-as.numeric(total_res[,'start'])+1)
  rownames(total_res)<-total_res[,'MCB_no']
  #Some of 'cgxxxxx' code point to same CpG, those cg-code are removed.
  total_res<-total_res[as.numeric(total_res[,'length'])>1,]
  if (verbose) {
    cat("Statistics (",nrow(total_res)," MCBs in total):\n")
    for (chr_set in unique(total_res[,'chromosomes'])) {
      cat(chr_set,": ")
      cat("total MCBs:",length(total_res[total_res[,'chromosomes'] %in% chr_set,'chromosomes'])," ")
      cat("Mean Length:",mean(as.numeric(total_res[total_res[,'chromosomes'] %in% chr_set,'length']))," ")
      cat("(Range: ",range(as.numeric(total_res[total_res[,'chromosomes'] %in% chr_set,'length'])),")\n")
    }
  }
  FunctionResults$MCBinformation<-total_res
  return(FunctionResults)
}

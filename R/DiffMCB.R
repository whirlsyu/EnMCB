#' @title  Differential expressed methylation correlated blocks
#'
#' @description This function is used to find the Methylation correlated blocks that differentially expressed between groups based on the attractor framework. 
#' This function calculates attractors of all the MCBs among the groups and find the attractor MCBs. \cr
#'
#' @details Currently, only illumina 450k platform is supported. \cr
#' If you want to use other platform, please provide the annotation file with CpG's chromosome and loci. \cr
#' The methylation profile need to convert into matrix format.
#'
#' @param methylation_matrix methylation profile matrix.
#' @param class_vector class vectors that indicated the groups.
#' @param mcb_matrix dataframe or matrix results returned by IdentifyMCB function.
#' @param min.cpgsize threshold for minimum CpGs must included in the individual MCBs.
#' @param pVals_num p value threshold for the test.
#' @param base_method base method used for calculation of differentially methylated regions,
#' should be one of 'Fstat','Tstat','eBayes'. Defualt is Fstat.
#' @param sec_method secondly method in attractor framework, should be one of 'kstest','ttest'. Defualt is ttest.
#' @param ... other parameters pass to the function.
#'
#' @author Xin Yu
#' @return
#' Object of class \code{list} with elements:
#'  \tabular{ll}{
#'    \code{global} \tab Character set contains statistical value for all CpG sites in MCBs. \cr
#'    \code{tab} \tab Matrix contains the information of results. \cr
#'  }
#' @examples
#' data('demo_data', package = "EnMCB")
#' data('demo_survival_data', package = "EnMCB")
#' data('demo_MCBinformation', package = "EnMCB")
#' #Using survival censoring as group label just for demo, 
#' #this may replace with disease and control group in real use.
#' diffMCB_results <- DiffMCB(demo_data$realdata,demo_survival_data[,2], 
#'                            demo_MCBinformation,
#'                            pVals_num = 1)
#'
#' @export
#'
#' @references
#' Xin Yu, De-Xin Kong, EnMCB: an R/bioconductor package for predicting disease progression based on methylation correlated blocks using ensemble models, Bioinformatics, 2021, btab415
#'
#'
DiffMCB<-function(
  methylation_matrix,
  class_vector, 
  mcb_matrix = NULL, 
  min.cpgsize = 5,
  pVals_num = 0.05,
  base_method = c('Fstat','Tstat','eBayes')[1], 
  sec_method = c('ttest','kstest')[1],
  ...
){
  if (is.matrix(methylation_matrix)) {
    dat.fr <- methylation_matrix
  } else {
    warning("The input parameter methylation_matrix must be a matrix. \n
            Now converting it into matrix.\n")
  }
  all.probes <- rownames(dat.fr)
  dat.fr <- as.matrix(dat.fr)
  class.vector <- as.factor(class_vector)
  annotation <- mcb_matrix
  all_included <- strsplit(paste(annotation[,'CpGs'],collapse = ' ')," ")[[1]]
  probes.hits <- intersect(all.probes, all_included)
  dat.detect.w <- dat.fr[rownames(dat.fr) %in% probes.hits,]
  dat.detect.w <- as.matrix(dat.detect.w)
  incidence.matrix <- buildIncidenceMatrix(rownames(dat.detect.w), annotation)
  keep.mcb <- apply(incidence.matrix, 1, sum) >= min.cpgsize
  incidence.matrix <- incidence.matrix[keep.mcb,]
  new.order <- order(class.vector, colnames(dat.detect.w))
  dat.detect.w <- dat.detect.w[, new.order]
  class.vector <- class.vector[new.order]
  
  evalMCB <- function(index, global,sec_method=c('ttest','kstest')[1]) {
    pway.vals <- global[index == 1]
    if (sec_method == 'ttest'){
      res<-t.test(pway.vals, global)
    }else if (sec_method == 'kstest'){
      res<-ks.test(pway.vals, global[index == 0])
    }
    return(c(res$statistic, res$p.value))
  }
  res<-list()
  if (base_method=='Fstat'){
    fstat <- apply(dat.detect.w, 1, function(y, x) {
      anova(lm(y ~ x))[[4]][1]
    }, x = class.vector)
    res$global <- fstat
    fstat <- log(fstat, 2)
    t.pvals <- apply(incidence.matrix, 1, evalMCB, global = fstat,sec_method=sec_method)
  }else if (base_method=='Tstat'){
    tstat <- apply(dat.detect.w, 1, function(y, x) {
      t.test(y ~ x)$statistic
    }, x = class.vector)
    res$global <- tstat
    tstat <- log(abs(tstat), 2)
    t.pvals <- apply(incidence.matrix, 1, evalMCB, global = tstat,sec_method=sec_method)
  }else if(base_method=='eBayes'){
    design<-model.matrix(~class.vector)
    fitlim<-limma::lmFit(dat.detect.w,design)
    fiteb<-limma::eBayes(fitlim)
    DEG<-limma::topTable(fiteb,coef=length(unique(class.vector)),n=nrow(dat.detect.w),p.value=1)
    DEG<-DEG[rownames(dat.detect.w),]
    tstat<-DEG$t
    res$global <- tstat
    tstat <- log(abs(tstat), 2)
    t.pvals <- apply(incidence.matrix, 1, evalMCB, global = tstat, sec_method=sec_method)
  }else{
    stop('The method indicator base_method is invalid.')
  }
  t.pvals_adjust <- p.adjust(as.numeric(t.pvals[2,]), "BH")
  ##cat(dim(t.pvals))
  
  size <- apply(incidence.matrix, 1, sum)
  tab <- data.frame(MCBID = rownames(incidence.matrix),
                    annotation[sapply(X = rownames(incidence.matrix),FUN = function(x){strsplit(x, "_")[[1]][2]}),],
                    statistic = as.numeric(t.pvals[1,]),
                    Pvalues = as.numeric(t.pvals[2,]),
                    AdjustedPvalues = t.pvals_adjust, 
                    NumberDetectedCpGs = size)
  #rownames(tab)<-sapply(X = rownames(tab),FUN = function(x){strsplit(x, "_")[[1]][2]})
  tab <- tab[order(t.pvals[2,]), ]
  tab <- subset(tab, AdjustedPvalues<=pVals_num)
  res$tab<-tab
  return(res)
}

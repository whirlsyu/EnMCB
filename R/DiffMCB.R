#' @title  Differential expressed methylation correlated blocks
#'
#' @description This function is used to find the Methylation correlated blocks that differentially expressed between groups. 
#' This function calculates attractors of all the MCBs among the groups and find the attractor MCBs. \cr
#'
#' @details Currently, only illumina 450k platform is supported, the methylation profile need to convert into matrix format.
#'
#' @param MCBs Methylation correlated blocks list.
#' @param method method used for calculation of differential expression, \cr
#' should be one of "attractors","t-test". Defualt is "attractors".
#' @param p_value p value threshold for the test.
#' @param min_CpGs threshold for minimum CpGs must included in the individual MCBs.
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
#' data('demo_data',package = "EnMCB")
#' 
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
  annotation, 
  min.cpgsize = 5,
  base_method = c('Fstat','Tstat')[1], 
  sec_method = c('ttest','kstest')[1],
  ...
){
  dat.fr <- methylation_matrix
  all.probes <- rownames(dat.fr)
  dat.fr <- as.matrix(dat.fr)
  class.vector <- as.factor(class_vector)
  
  Illumina_Infinium_Human_Methylation_450K<-get450kAnno()
  Illumina_Infinium_Human_Methylation_450K<-Illumina_Infinium_Human_Methylation_450K[!is.na(Illumina_Infinium_Human_Methylation_450K[,'pos']),]
  intersect_cpg<-intersect(rownames(Illumina_Infinium_Human_Methylation_450K),rownames(methylation_matrix))
  annotation<-Illumina_Infinium_Human_Methylation_450K[intersect_cpg,]
  
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
  }else if(base_method=='LIMMA'){
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
  res$tab<-tab
  return(res)
}

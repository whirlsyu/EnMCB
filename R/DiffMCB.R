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
  MCBs,
  method=c("attractors")[1],
  p_value = 0.05,
  min_CpGs = 5,
  platform = "Illumina Methylation 450K"
){
    union_cpgs<-list()
    i = 1
    for (MCB in MCBs) {
      nrow_MCB<-nrow(MCB)
      if (nrow_MCB==0){
        stop("compare MCB: the MCB info matrix doesn't include any data, please check the data.")
      }
      for (j in seq(nrow_MCB)) {
        CpG_list<-strsplit(paste(MCB[j, "CpGs"], collapse = " "), " ")[[1]]
        if (length(CpG_list)>=min_CpGs){
          union_cpgs[[i]]<-CpG_list
          i = i + 1
        }
      }
    }
    i = rep(1:length(union_cpgs), lengths(union_cpgs)) 
    j = factor(unlist(union_cpgs))
    tab = Matrix::sparseMatrix(i = i, j = as.integer(j), x = TRUE, dimnames = list(NULL, levels(j)))
    connects = Matrix::tcrossprod(tab, boolArith = TRUE)
    # 'graph_from_adjacency_matrix' seems to not work with the "connects" object directly. 
    # An alternative to coercing "connects" here would be to build it as 'tcrossprod(tab) > 0'
    group = igraph::clusters(igraph::graph_from_adjacency_matrix(as(connects, "lsCMatrix")))$membership
    #group
    #[1] 1 2 2 2 3 3
    CpGs_MMB<-tapply(union_cpgs, group, function(x) sort(unique(unlist(x))))
    ref_MMB<-re_AnnotatedMMB(CpGs_MMB)
    patterns<-matrix(NA,nrow = nrow(ref_MMB),ncol = length(MCBs))
    bscore<-rep(NA,nrow(ref_MMB))
    bar<-utils::txtProgressBar(min = 1,max = nrow(ref_MMB),char = "#",style = 3)
    for (i in seq(nrow(ref_MMB))) {
      utils::setTxtProgressBar(bar, i)
      CpG_list<-strsplit(paste(ref_MMB[i,'CpGs'], collapse = " "), " ")[[1]]
      pattern_codes<-NULL
      for (j in seq(length(MCBs))) {
        pattern_code<-NULL
        for (CpG in CpG_list) {
          
          if (is.integer0(grep(CpG,MCBs[[j]][,'CpGs']))){
            pattern_code<-c(pattern_code,0)
          }else{
            pattern_code<-c(pattern_code,1)
          }
        }
        patterns[i,j]<-paste(pattern_code,collapse = ' ')
        pattern_codes[[j]]<-pattern_code
      }
      bscore[i] <- mean(outer(pattern_codes,pattern_codes,Vectorize(bmeasures))
                        [upper.tri( matrix(NA,
                                           nrow = length(pattern_codes),
                                           ncol = length(pattern_codes))  ) ])
    }
    MCB_no<-seq(nrow(ref_MMB))
    data.frame(ref_MMB,patterns,bscore)
}

#'@title Preprocess the Beta value matrix
#'
#'@description This process is optional for the pipeline.
#'This function pre-process the Beta matrix and transform the Beta value into M value.
#'
#'@param met methylation matrix for CpGs. Rows are the CpG names, columns are samples.
#'@param Mvalue Boolean value, TRUE for the M transformation.
#'@param constant_offset the constant offset used in the M transformation formula.
#'@param remove_na Boolean value, if TRUE ,CpGs with NA values will be removed.
#'@param remove_percentage If precentage of NA value exceed the threshold(percentage), the whole CpG probe will be removed. Otherwise, the NA values are replaced with rowmeans.
#'
#'@usage pre_process_methylation(met,Mvalue,constant_offset,remove_na,remove_percentage)
#'
#'@export
#'
#'@return Object of class \code{matrix}.
#'
pre_process_methylation<-function(met,Mvalue=TRUE,constant_offset=0,remove_na=TRUE,remove_percentage=30){
  threshold=remove_percentage/100
  flag_na<-rep(TRUE,nrow(met))
  names(flag_na)<-rownames(met)
  n<-ncol(met)
  if (remove_na) {
    for (eachrow in rownames(met)) {
      if (sum(is.na(met[eachrow,])|met[eachrow,]==1|met[eachrow,]==0)>n*threshold){
        # whole row will be removed.
        flag_na[eachrow]<-FALSE
      }else{
        if (sum(is.na(met[eachrow,])|met[eachrow,]==1|met[eachrow,]==0)>0){
          met[eachrow,is.na(met[eachrow,])|met[eachrow,]==1|met[eachrow,]==0]<-mean(met[eachrow,(!is.na(met[eachrow,]))&met[eachrow,]!=1&met[eachrow,]!=0])
        }
      }
    }
    met<-met[flag_na,]
  }
  if (Mvalue){
    met<-log2((met+constant_offset)/(1-met+constant_offset))
  }
  met
}

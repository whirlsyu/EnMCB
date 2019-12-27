make_eset<-function(matrix,pData){
  requireNamespace("affycoretools")

  colnames<-rownames(pData)

  commdata<-intersect(colnames,colnames(matrix))

  pData<-pData[rownames(pData) %in% commdata,]

  matrix<-matrix[,colnames(matrix) %in% commdata]

  if (nrow(pData)!=ncol(matrix)){stop("matrix colnames do not match the pData rownames")}

  #pData<-pData[colnames(matrix),]
  matrix<-matrix[,rownames(pData)]

  metadata<-data.frame(labelDescription=colnames(pData),stringsAsFactors = FALSE)

  adf<-new("AnnotatedDataFrame",data=pData,varMetadata=metadata)

  eSet<-new("ExpressionSet",exprs=as.matrix(matrix),phenoData=adf)

  return(eSet)
}

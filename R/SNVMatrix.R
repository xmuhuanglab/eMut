# @import dplyr
#' SNV matrix construction for each sample
#'
#' contruct mutation matrix for one samples from Monopogen outputs
#' @param inputDir direction of somatic mutations of Monopogen output
#' @export
SNVMatrix<-function(inputDir)
{
  files<-list.files(path=inputDir,pattern="*.SNV_mat.RDS")
  if(length(files)>0){
    files<-paste0(inputDir,"/",files)
    mats<-lapply(files,function(x){
      data<-readRDS(x)
      data<-data[,19:ncol(data)] %>% as.matrix()
    })
    SNVs<-sapply(mats,function(x){ row.names(x)}) %>% unlist() %>% as.vector() %>% unique()
    cells<-sapply(mats,function(x){ colnames(x)}) %>% unlist() %>% as.vector() %>% unique()
    mat<-matrix("0/0",nrow=length(SNVs),ncol=length(cells),dimnames=list(SNVs,cells))
    for(i in 1:length(mats)){
      temp<-mats[[i]]
      mat[row.names(temp),colnames(temp)]<-temp
    }
    return(mat)
  }
}

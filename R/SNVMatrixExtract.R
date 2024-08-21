# @import dplyr
#' SNV matrix construction for each sample
#'
#' contruct mutation matrix for one samples from Monopogen outputs
#' @param inputDirs direction of somatic mutations of Monopogen output(named strings or named strings vector)
#' @export
SNVMatrixExtract<-function(inputDirs)
{
  if(is.null(names(inputDirs))){
    stop("please input the direction vector with names")
    }

  mat.list<-lapply(1:length(inputDirs),function(idx){
    if(!dir.exists(inputDirs[idx])){
      stop(paste("Direction not exist:",inputDirs[idx]))
    }

    ## extraction of mutation-cell matrix for each sample
    files<-list.files(path=inputDirs[idx],pattern="*.SNV_mat.RDS")
    if(length(files)>0){
      files<-paste0(inputDirs[idx],"/",files)
      mats<-lapply(files,function(x){
        data<-readRDS(x)
        data<-data[,19:ncol(data)] %>% as.matrix()
      })
    
      ## combination of mutation-cell matrix for each sample
      SNVs<-sapply(mats,function(x){ row.names(x)}) %>% unlist() %>% as.vector() %>% unique()
      cells<-sapply(mats,function(x){ colnames(x)}) %>% unlist() %>% as.vector() %>% unique()
      mat<-matrix("0/0",nrow=length(SNVs),ncol=length(cells),dimnames=list(SNVs,cells))
      for(i in 1:length(mats)){
        temp<-mats[[i]]
        mat[row.names(temp),colnames(temp)]<-temp
      }
      # sample#cell_barcode as identifier to avoid the same cell barcodes in different samples
      colnames(mat)<-paste(names(inputDirs)[idx],colnames(mat),sep="#")  
      return(mat)
    }
  })
  names(mat.list)<-names(inputDirs)

  ##  combined mat.list of each sample to a larger matrix
  SNVs<-sapply(mat.list,function(x){ row.names(x)}) %>% unlist() %>% unique()
  cells<-sapply(mat.list,function(x){ colnames(x)}) %>% unlist() %>% unique()
  mat<-matrix("0/0",nrow=length(SNVs),ncol=length(cells),dimnames=list(SNVs,cells))
  for(i in names(mat.list)){
    tmp<-mat.list[[i]]
    if(!is.null(tmp)){
      mat[row.names(tmp),colnames(tmp)]<-tmp
      }
    }
  return(mat)
}

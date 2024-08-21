# @import dplyr
# @import data.table
#' SNV extract 
#'
#' extract SNV from Monopogen outputs
#' @param inputDirs direction of somatic mutations of Monopogen output(named strings or named strings vector)
#' @export

SNVExtract<- function(inputDirs)
{
  require(plotly)
  if(is.null(names(inputDirs))){
    stop("please input the direction vector with names")
  }
  SNVDF.list<-lapply(1:length(inputDirs),function(idx){
    if(!dir.exists(inputDirs[idx])){
      stop(paste("Direction not exist:",inputDirs[idx]))
    }
    ###  extraction mutation
    files<-list.files(path=inputDirs[idx],pattern="*putativeSNVs.csv")
    if(length(files)>0){
    files<-paste0(inputDirs[idx],"/",files)
    SNV.df<-lapply(files,function(x){
      data<-fread(x,sep=",")
    }) %>% rbindlist() %>% distinct(.,chr,pos,Ref_allele,Alt_allele,.keep_all = TRUE)
    SNV.df$mut<-paste(SNV.df$chr,SNV.df$pos,SNV.df$Ref_allele,SNV.df$Alt_allele,sep=";")

    ### add sample information to the result
    if(!is.null(SNV.df)){
       SNV.df$sample<-rep(names(inputDirs)[idx],nrow(SNV.df))
       SNV.df  }
    }
  })

  ###  combination of mutations
  SNV.df<-rbindlist(SNVDF.list)
  SNV.df<-SNV.df[,c(1:4,13,5:12,14)]
  names(SNV.df)[1:5]<-c("CHROM","POS","REF","ALT","ID")
  return(SNV.df)
}

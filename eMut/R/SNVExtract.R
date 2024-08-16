# @import dplyr
# @import data.table
#' SNV extract
#'
#' extract SNV from Monopogen outputs
#' @param inputDir direction of somatic mutations of Monopogen output
#' @export
SNVExtract<- function(inputDir)
{
  files<-list.files(path=inputDir,pattern="*putativeSNVs.csv")
  if(length(files)>0){
    files<-paste0(inputDir,"/",files)
    SNV.df<-lapply(files,function(x){
      data<-fread(x,sep=",")
    }) %>% rbindlist() %>% distinct(.,chr,pos,Ref_allele,Alt_allele,.keep_all = TRUE)

    SNV.df$mut<-paste(SNV.df$chr,SNV.df$pos,SNV.df$Ref_allele,SNV.df$Alt_allele,sep=";")
    return(SNV.df)
  }
}

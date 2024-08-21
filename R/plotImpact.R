# @import ggplot2
# @import ComplexHeatmap
# @import dplyr
#' To map the effects of mutations on the expression of their target genes, including the spectrum of mutations and the expression of target genes significantly affected by them (up- or down-regulated).
#'
#' SNV Imputation
#' @param SNV.df Combined mutations for samples (or cells)
#' @param tpm gene expression matrix
#' @param result the comparsion of target gene expression between mutated samples(or cells) and wilde type [obtained from mutTargetExp()]
#' @param p.cutOff the cutoff of p-value for selected mutations
#' @export
plotImpact<-function(SNV.df,
                     tpm,
                     result,
                     p.cutOff=0.05){
  require(ggplot2)
  require(ComplexHeatmap)
  require(dplyr)
  
  print("select mutations:")
  ###    select mutations 
  SNV<-unique(SNV.df$ID)
  SNV.sam<-matrix(0,length(unique(SNV)),length(unique(SNV.df$sample)),dimnames=list(unique(SNV),unique(SNV.df$sample)))
  for(i in 1:nrow(SNV.df)){
    SNV.sam[SNV.df$ID[i],SNV.df$sample[i]]<-1
    }
  SNV.flt<-rbind(SNV.sam[row.names(SNV.sam) %in% result$mutation[result$adjsP<p.cutOff & result$foldchange>1],],
               SNV.sam[row.names(SNV.sam) %in% result$mutation[result$adjsP<p.cutOff & result$foldchange<1],])

  mutation.flt<-c(result$mutation[result$adjsP<p.cutOff & result$foldchange>1],
                  result$mutation[result$adjsP<p.cutOff & result$foldchange<1])

  tpm.flt<-tpm[result$gene[match(row.names(SNV.flt),result$mutation)],]
  tpm.flt<-apply(tpm.flt,1,function(x){scale(x)}) %>% t()
  colnames(tpm.flt)<-colnames(tpm)
  
  print("plot mutation profile:")
  ###    The heatmap of mutation profile
  p1<-Heatmap(SNV.flt, show_column_names =TRUE, show_row_names = FALSE,
        cluster_rows = FALSE, cluster_columns = FALSE, show_column_dend = F,
        raster_quality = 2,use_raster = TRUE,
        col = c("white","red"),
        border = TRUE,
        row_split=rep(c("up-regulated","down-regulated"),c(length(which(result$adjsP<p.cutOff & result$foldchange>1)),length(which(result$adjsP<p.cutOff & result$foldchange<1)))),
        heatmap_legend_param = list(title="SNV"))
  
  print("plot gene expression profile:")
  ###    The heatmap of matched target gene expression profile of above mutation profile 
  RNA.sample<-sapply(strsplit(colnames(tpm),split = "_"),function(x){x[[1]]}) %>% unlist()
  col.order<-lapply(unique(SNV.df$sample),function(x){
    which(RNA.sample==x)
  }) %>% unlist()
  p2<-Heatmap(tpm.flt[,col.order], show_column_names =TRUE, show_row_names = FALSE,
        cluster_rows = FALSE, cluster_columns = FALSE, show_column_dend = F,
        raster_quality = 2,use_raster = TRUE,
        border = TRUE,
        row_split=rep(c("up-regulated","down-regulated"),c(length(which(result$adjsP<p.cutOff & result$foldchange>1)),length(which(result$adjsP<p.cutOff & result$foldchange<1)))),
        heatmap_legend_param = list(title="scaled TPM"))
        
  return(p1+p2)
}

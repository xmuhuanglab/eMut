# @import ggplot2
# @import ArchR
# @import SCAVENGE
#' The umap of SNV Imputation
#'
#' SNV Imputation
#' @param SNV.df The combined SNV information from SNVExtract() (data.frame)
#' @param mut.show The index of mutations wanted to show (vector)
#' @export
plotSNVSamples<-function(SNV.df,
                         mut.show=NULL){

  require(ggplot2)
  require(ComplexHeatmap)

  ##  contruct SNV-sample matrix
  SNV<-unique(SNV.df$ID)
  SNV.sam<-matrix(0,length(unique(SNV)),length(unique(SNV.df$sample)),dimnames=list(unique(SNV),unique(SNV.df$sample)))
  for(i in 1:nrow(SNV.df)){
    SNV.sam[SNV.df$ID[i],SNV.df$sample[i]]<-1
  }

  ###  heatmap of SNV-sample matrix
  p<-Heatmap(SNV.sam, show_column_names =TRUE, show_row_names = FALSE,
        cluster_rows = FALSE, cluster_columns = FALSE, show_column_dend = F,
        raster_quality = 2,use_raster = TRUE,
        col = c("white","red"),
        border = TRUE,
        row_title = paste0(nrow(SNV.sam)," SNVs"),
        heatmap_legend_param = list(title="SNV"))
  
  if(!is.null(mut.show)){
    p<-p+rowAnnotation(link = anno_mark(at =mut.show, 
                                  labels = row.names(SNV.sam)[mut.show], labels_gp = gpar(fontsize = 9)))
  }
  return(p)
}

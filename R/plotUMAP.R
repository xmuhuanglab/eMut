# @import ggplot2
# @import ArchR
# @import SCAVENGE
#' The umap of SNV Imputation
#'
#' SNV Imputation
#' @param UMAP.loc The location of UMAP
#' @param TRS.df The result of SNV imputation of SNVImputation()
#' @export
plotUMAP<-function(UMAP.loc,
                   TRS.df){
  require(ggplot2)

  ###    add imputation of mutated cells to UMAP.loc
  UMAP.loc$imputed<-TRS.df$true_cell_top_idx[match(row.names(UMAP.loc),row.names(TRS.df))]
  UMAP.loc$seed<-TRS.df$seed[match(row.names(UMAP.loc),row.names(TRS.df))]

  ###    UMAP plot of raw mutated cells
  p1 <- ggplot()+geom_point(UMAP.loc[UMAP.loc$seed==F,],mapping=aes(x=UMAP1,y=UMAP2),size=0.5,color="gray80")+
    geom_point(UMAP.loc[UMAP.loc$seed==T,],mapping=aes(x=UMAP1,y=UMAP2),size=2,color="red3")+
    theme_bw()+
    theme(panel.grid =element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          panel.border = element_blank())+
    labs(x=NULL,y=NULL,title="Raw mutated cells")

  ###    UMAP plot of imputed mutated cells
  p2 <- ggplot()+geom_point(UMAP.loc[UMAP.loc$imputed==F,],mapping=aes(x=UMAP1,y=UMAP2),size=0.5,color="gray80")+
    geom_point(UMAP.loc[UMAP.loc$imputed==T,],mapping=aes(x=UMAP1,y=UMAP2),size=0.5,color="red3")+
    theme_bw()+
    theme(panel.grid =element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          panel.border = element_blank())+
    labs(x=NULL,y=NULL,title="Imputed mutated cells")
  p<-ggAlignPlots(p1,p2,type = "h")
  return(p)
}

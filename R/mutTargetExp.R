#' The impact of muttaion in target gene expression (single-cell level or sample-level)
#'
#' Cell type enrichment of mutant cells(raw or imputed)
#' @param SNV.df Combined mutations for samples (or cells)
#' @param SNV.anno mutation annotated file
#' @param tpm gene expression matrix
#' @export
mutTargetExp<-function(SNV.df,
                       SNV.anno,
                       tpm){
  ###
  library(pbapply)
  ###
  result<-pblapply(1:nrow(SNV.anno),function(x){
    mut<-SNV.anno[x,1]
    gene<-SNV.anno$Gene.Name[x]
    mut.sample<-SNV.df$sample[SNV.df$ID==mut] %>% unique()
    RNA.sample<-sapply(strsplit(colnames(tpm),split = "_"),function(x){x[[1]]}) %>% unlist()
    if(gene %in% row.names(tpm) & length(mut.sample)>0){
      mut.exp<-tpm[gene,RNA.sample %in% mut.sample] %>% t() %>% as.vector()
      nonMut.exp<-tpm[gene,!RNA.sample %in% mut.sample] %>% t() %>% as.vector()
      p<-wilcox.test(mut.exp,nonMut.exp)$p.value
      foldchange<-mean(mut.exp)/mean(nonMut.exp)   
      data.frame(gene=gene,mutation=mut,mutSamples=paste(mut.sample,collapse = ";"),pValue=p,foldchange=foldchange)
      }
    }) %>% rbindlist()
  result$adjsP<-p.adjust(result$pValue,method="fdr")
  return(result)
}


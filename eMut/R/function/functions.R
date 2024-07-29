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
    colnames(mat)<-paste(sam,colnames(mat),sep="#")
    return(mat)
  }    
}


#' SNV Imputation
#'
#' SNV Imputation
#' @param countMatrix direction of somatic mutations of Monopogen output
#' @param knnGraph direction of somatic mutations of Monopogen output
#' @param SNVMatrix direction of somatic mutations of Monopogen output
#' @param mutations direction of somatic mutations of Monopogen output
#' @param numk direction of somatic mutations of Monopogen output
#' @param queryCell_cutoff direction of somatic mutations of Monopogen output
#' @param ncors direction of somatic mutations of Monopogen output
#' @export
SNVImputation<-function(countMatrix,
                        knnGraph=NULL,
                        SNVMatrix,
                        mutations=row.names(SNVMatrix),
                        numk=30,
                        queryCell_cutoff=5,
                        ncors=10)
{
  ###  
  library(SCAVENGE)
  library(pbapply)
  
  ###  1) Constructing KNN Neighborhood Graphs
  if(is.null(knnGraph)){
    print("Constructing KNN Neighborhood Graphs:")
    tfidf.mat <- tfidf(bmat=countMatrix, 
                       mat_binary=TRUE,
                       TF=TRUE, 
                       log_TF=TRUE)
    
    lsi.mat <- do_lsi(mat = tfidf.mat, dims = 30)
    mutualknn <- getmutualknn(lsimat = lsi.mat, num_k = 30)
    cells<-sapply(strsplit(rownames(mutualknn),split="_"),function(x){x[1]}) %>% unlist()
  }else{
    print("KNN Neighborhood Graphs:")
    mutualknn <- knnGraph
    cells<-sapply(strsplit(rownames(mutualknn),split="_"),function(x){x[1]}) %>% unlist()
  }
  
  ### 2) imputation
  print(paste0("Imputation for",length(mutations),"mutations"))
  TRS.list<-pblapply(mutations,function(x){
    ##  mutant cells as seed cells 
    seed.cells<-colnames(SNVMatrix)[which(SNVMatrix[x,]=="0/1" | SNVMatrix[x,]=="1/1")]
    query.cells<-rownames(mutualknn)[match(seed.cells,cells)] %>% na.omit()
    
    if(length(query.cells)>=queryCell_cutoff){
      np_score <- randomWalk_sparse(intM = mutualknn, 
                                    queryCells = query.cells, 
                                    gamma = 0.05)
      ####  scale
      omit_idx <- np_score==0
      sum(omit_idx)
      mutualknn <- mutualknn[!omit_idx, !omit_idx]
      np_score <- np_score[!omit_idx]
      TRS <- capOutlierQuantile(x = np_score, 
                                q_ceiling = 0.95) |> max_min_scale()
      
      TRS_mat <- data.frame(np_score, TRS_score=TRS)
      TRS_mat$seed<-rep(FALSE,nrow(TRS_mat))
      TRS_mat.cells<-sapply(strsplit(rownames(TRS_mat),split="_"),function(x){x[1]}) %>% unlist()
      
      TRS_mat[match(seed.cells,TRS_mat.cells),"seed"]<-TRUE
      
      #######   permutated seed cells
      set.seed(123)
      mono_permu <- get_sigcell_simple(knn_sparse_mat=mutualknn, 
                                       seed_idx=TRS_mat$seed, 
                                       topseed_npscore=TRS_mat$np_score, 
                                       permutation_times=1000, # Increase to >=1000 in practice
                                       true_cell_significance=0.05, 
                                       rda_output=FALSE, 
                                       mycores=ncors,# Increase in practice
                                       rw_gamma=0.05)
      TRS_mat1 <- data.frame(TRS_mat, mono_permu)
      return(TRS_mat1)
    }else {
      return(NULL)}
  })
  names(TRS.list)<-mutations
  return(TRS.list)
}


#' SNV Imputation umap
#'
#' SNV Imputation
#' @param UMAP.loc The location of UMAP
#' @param TRS.df The result of SNV imputation of SNVImputation
#' @export
UMAPImputation<-function(UMAPLoc,
                         TRS.df){
  
  ###    add imputation of mutant cells to UMAP.loc
  UMAPLoc$imputed<-TRS.df$true_cell_top_idx[match(row.names(UMAPLoc),row.names(TRS.df))]
  UMAPLoc$seed<-TRS.df$seed[match(row.names(UMAPLoc),row.names(TRS.df))]
  
  ###    UMAP plot of raw mutant cells
  p1 <- ggplot()+geom_point(UMAPLoc[UMAPLoc$seed==F,],mapping=aes(x=UMAP1,y=UMAP2),size=0.5,color="gray80")+
    geom_point(UMAPLoc[UMAPLoc$seed==T,],mapping=aes(x=UMAP1,y=UMAP2),size=2,color="red3")+
    theme_bw()+
    theme(panel.grid =element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          panel.border = element_blank())+
    labs(x=NULL,y=NULL,title="raw mutant cells")
  
  ###    UMAP plot of imputed mutant cells
  p2 <- ggplot()+geom_point(UMAPLoc[UMAPLoc$imputed==F,],mapping=aes(x=UMAP1,y=UMAP2),size=0.5,color="gray80")+
    geom_point(UMAPLoc[UMAPLoc$imputed==T,],mapping=aes(x=UMAP1,y=UMAP2),size=0.5,color="red3")+
    theme_bw()+
    theme(panel.grid =element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          panel.border = element_blank())+
    labs(x=NULL,y=NULL,title="imputed mutant cells")
  p<-ggAlignPlots(p1,p2,type = "h")
  return(p)
}


#' Cell type enrichment
#'
#' Cell type enrichment of mutant cells(raw or imputed)
#' @param cells barcode of all cells
#' @param cellTypes The cell types of each cell
#' @param mutCells barcode of mutant cells
#' @param mutCells_ciuoff This threshold is used to determine how many mutant cells to exceed before cell type 
#'                        enrichment analysis is performed. 
#' @param cellType_cutoff This threshold to select cell types with a certain number of cells to consider for inclusion in the analysis 
#' @export
cellTypeEnrich<-function(cells,
                         cellTypes,
                         mutCells,
                         mutCells_cutoff=5,
                         cellType_cutoff=50){
  ###
  library(data.table)
  ###
  Types.mum<-table(cellTypes)
  Types<-names(Types.mum)[which(Types.mum>=cellType_cutoff)]
  if(length(mutCells)>mutCells_cutoff){
    nonmut.cells<-setdiff(cells,mutCells)
    p.value<-rep(NA,length(Types))
    mut.num<-rep(NA,length(Types))
    for(i in Types){
      types.cells<-cells[cellTypes==i]
      inter.cells<-intersect(mutCells,types.cells)
      p<-phyper(length(inter.cells)-1,length(types.cells),
                length(cells)-length(types.cells),length(mutCells),lower.tail=F)
      p.value[match(i,Types)]<-p
      mut.num[match(i,Types)]<-length(inter.cells)
      
    }
    result<-data.frame(p=p.value,types=Types,mut.num=mut.num)
    return(result)
  }
}




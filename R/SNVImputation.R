# @import Seurat
# @import Signac
# @import SCAVENGE
# @import dplyr
# @import pbapply
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
  require(SCAVENGE)
  require(pbapply)

  ###  1) Constructing KNN Neighborhood Graphs
  if(is.null(knnGraph)){
    print("Constructing KNN Neighborhood Graphs:")
    tfidf.mat <- tfidf(bmat=countMatrix,
                       mat_binary=TRUE,
                       TF=TRUE,
                       log_TF=TRUE)

    lsi.mat <- do_lsi(mat = tfidf.mat, dims = numk)
    mutualknn <- getmutualknn(lsimat = lsi.mat, num_k = numk)
    cells<-sapply(strsplit(rownames(mutualknn),split="_"),function(x){x[1]}) %>% unlist()
  }else{
    print("KNN Neighborhood Graphs:")
    mutualknn <- knnGraph
    cells<-sapply(strsplit(rownames(mutualknn),split="_"),function(x){x[1]}) %>% unlist()
  }
  
  ### 2) imputation
  if(is.null(row.names(SNVMatrix))){
    stop("please add mutation identifiers to the rownames of SNVMatrix")
  }
  print(paste("Imputation for",length(mutations),"mutations"))
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

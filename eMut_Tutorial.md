
# Tutorial for eMut

## Introduction
The eMut, an integrated pipeline for detecting, imputing, and characterizing non-coding mutations in cis-regulatory elements(CREs) with functional consequences at the single-cell level.

![image](https://github.com/xmuhuanglab/eMut/blob/main/Figures/eMut_workflow.png)


## workflow
Briefly, eMut consists of two main modules: mutation detection and functional interpretation. <br />
Step1. Mutation detection: eMut detects mutations in each cell by implementing methods such as Monopogen or GATK by single-cell chromatin accessibility data. <br />
Step2. Mutation imputation (optional): Given the sparse of scATAC-seq data, we further imputed candidate mutated cells by network propagation using mutated cells (seed cells) in cell-cell similarity graph. <br />
Step3. Functional interpretation: <br />
1) recognize cell type-specific or lineage-specific mutations; 
2) identify hypermutated CREs with significant excess of mutations to characterize potentially important enhancers; 
3) predict the effects of mutations on transcription factor motifs (loss or gain);
4) compare target gene expression changes between mutated cells (or samples) and wild-type. 

### Step1. Mutation detection
Input: <br>
Bam-format files of scATAC-seq data<br>
As an example, the [Monopogen](https://github.com/KChen-lab/Monopogen) implementation of somatic mutation prediction for a single sample (without matched normal samples), as well as the simultaneous detection of somatic and germline mutations based on [GATK Mutect2](https://github.com/broadinstitute/gatk), are shown here.
```
/MutationDetection/1.run_GATK.py
/MutationDetection/2.mutation_annotation.py
/MutationDetection/Mutation_detection.sh
```

### Step2. Mutation imputation (optional)
Considering dropout events due to sparsity of single-cell technical, we refer to [SCAVENGE](https://github.com/sankaranlab/SCAVENGE) with its default parameters to infer potential mutated cells. Briefly, a M-kNN graph was constructed based on scATAC-seq data to represent cell-cell similarity. For a given mutation, the mutated cells (as seed cells) were projected onto the M-kNN graph. Through network propagation of these seed nodes, relevant cells were identified as potential mutated cells.<br>
Input: <br>
(1) mutation profile: mutation-by-cell matrix (The matrix can be organized from Monopogen results introduced as step1); <br> 
(2) scATAC-seq data: peak-by-cell matrix or knnGraph (obatined from ArchR or signac object); <br>

```r
library(graphics)
library(ggforce)
library(scales)
library(ggpubr)
library(ggplot2)
library(pbapply)
library(SCAVENGE)
library(eMut)
library(Seurat)
library(Signac)

load("./mutualknn30.Rdata")      # m-knn graph   (see detials in SCAVENGE)
load("./SNVMat.Rdata")           # mutation-by-cell matrix from the Step1. Mutation detection 

TRS.list<-SNVImputation(countMatrix=NULL,            # peak-by-cell matrix(obatined from ArchR or signac object)
                        knnGraph=mutualknn30,        # knnGraph 
                        SNVMatrix=mat,               # mutation-by-cell matrix
                        mutations=row.names(mat),    # mutations
                        numk=30,                     # Number of nearest neighbors 
                        queryCell_cutoff=5,          # Minimum number of mutant cells, greater than this threshold for subsequent analysis
                        ncors=5)                     # Number of threads
```

### Step3. functional interpertaion
#### (1) cell type enrichment
Input: <br>
(1) mutation profile: mutation-by-cell matrix (raw or imputed); <br>
(2) scATAC-seq data: ArchR or signac object; <br>
Here is an example demonstration with an ArchR object .
```r
library(ArchR)
library(eMut)

load("./SNVMat.Rdata")
proj<- loadArchRProject(path = "./ArchR", force = FALSE, showLogo = FALSE)
cells<-row.names(getCellColData(proj))
cellTypes<-proj$NamedClust
```
#####  cell type enrichment for mutation profile (raw)
```r
mut.type<-pblapply(row.names(mat),function(x){
    mutCells<-colnames(mat)[which(mat[x,]=="0/1" | mat[x,]=="1/1")]  # mutated cells
    result<-cellTypeEnrich(cells,                  # The names of all cells
                         cellTypes,                # The cell type
                         mutCells,                 # Mutated cells
                         mutCells_cutoff=50,       # Muttaion with a number of mutated cells greater than this threshold were subsequently analyzed
                         cellType_cutoff=50)       # Cell type with a number of mutated cells greater than this threshold were subsequently analyzed
    if(!is.null(result)){
      result$mut<-rep(x,nrow(result))
      result  
    }
}) %>% rbindlist()
mut.type$p<-p.adjust(mut.type$p)
```

##### cell type enrichment for mutation profile (imputed)
```r
mut.type<-pblapply(names(TRS.list),function(x){
    mutCells<-row.names(TRS.list[[x]])[TRS.list[[x]]$true_cell_top_idx==TRUE] #imputed mutated cells
    result<-cellTypeEnrich(cells,
                         cellTypes,
                         mutCells,
                         mutCells_cutoff=50,
                         cellType_cutoff=50)
    if(!is.null(result)){
      result$mut<-rep(x,nrow(result))
      result  
    }
}) %>% rbindlist()
mut.type$p<-p.adjust(mut.type$p)
```

####  (2) hyperMutated CREs
Among all accessible regions, we adapted the [ActiveDriverWGS](https://github.com/reimandlab/ActiveDriverWGSR) method with modification to identify hypermuated CREs based on scATAC-seq data. Specifically, we changed adjacent flanking genomic regions to flanking accessible regions (±500 kbps) for training the model of expected mutations, we identified hypermuated CREs (observed excess expected mutations) in each sample. 
Input: <br>
(1) mutation profile: annotated mutation file(VCf/maf format); <br>
(2) peak file: The genomic location of chromatin accessible region; <br>

```r
library(ActiveDriverWGS)
library(tidyr)
library(BSgenome.Hsapiens.UCSC.hg38)
library(eMut)

##  load  mutation profile
mut.df<-data.table::fread(
    file = "./AML10.maf",
    sep = "\t",stringsAsFactors = FALSE,verbose = FALSE,data.table = TRUE,
    showProgress = TRUE,header = TRUE,fill = TRUE,
    skip =1,quote = "")
mut.df<-mut.df[mut.df$FILTER=="PASS",]
mut.df<-mut.df[,c("Chromosome","Start_Position","End_Position",
                             "Reference_Allele","Tumor_Seq_Allele2","Tumor_Sample_Barcode")]
colnames(mut.df)<-c("chr","pos1","pos2","ref","alt","patient")
mut.df<-mut.df[mut.df$chr!="chrM",]

##  load peak
peaks.df<-read.table("./peaks.bed",header=F,quote="",sep="\t")
colnames(peaks.df)<-c("chr","start","end","id")
openRegions<-GenomicRanges::GRanges(seqnames=peaks.df$chr, 
                                    IRanges::IRanges(peaks.df$start,peaks.df$end))
##  identify hypermutated CREs
hyperMut<-ActiveDriverCRE(mutations = mut.df,              # mutations
                         elements = peaks.df,              # Open region as elements 
                         ref_genome = "hg38",              # Reference genome
                         mc.cores=4,                       # Number of threads
                         window_size=1000000,              # The window size of flanking region as background 
                         detect_depleted_mutations=FALSE,  # The detection of hypomutated region
                         openRegions = openRegions,        # all open regions (union of all peaks)
                         recovery.dir="./tmp/")  # Temporary file paths to quickly recover results

```
####  (3) prediction of TF binding motif change
To explore the impact of mutations in their located enhancer, [motifbreakR](https://github.com/Simon-Coetzee/motifBreakR) was applied to predict TF motif disruptions (loss or gain) for a large number of single-nucleotide variants using several different sources of TF motifs (e.g. JASPAR and ENCODE). In the predicted results, "strong" effect motif change will be considered as the potential impact of mutations.
Input: <br>
(1) mutation profile: mutation file(VCf/maf format); <br>
Function parameters are consistent with motifbreakR
```r
library(motifbreakR)
library(MotifDb)
library(BSgenome.Hsapiens.UCSC.hg38)
library(eMut)

##   change mutation format
load("./SNV.Rdata")
mut<-data.frame(Chromosome=SNV.df$CHROM,Start_Position=SNV.df$POS-1,End_Position=SNV.df$POS,
                names=gsub(";",":",SNV.df$ID),score=rep(0,nrow(SNV.df)),strand=rep("+",nrow(SNV.df)))
mut<-mut %>% distinct(names,.keep_all = TRUE)
write.table(mut,file="./Monopogen/summary/SNVsForMotifBreakR.bed",sep="\t",
            col.names = F,row.names = F,quote=F)

##   motif change prediction
data(motifbreakR_motif)
ENCODE<-subset (motifbreakR_motif, dataSource=="ENCODE-motif" & organism=="Hsapiens")
snps.mb.frombed <- snps.from.file(file = "./Monopogen/summary/SNVsForMotifBreakR.bed", 
                                  search.genome = BSgenome.Hsapiens.UCSC.hg38,
                                  format = "bed")

SNV.motifs <- motifbreakR(snpList = snps.mb.frombed, filterp = TRUE,
                       pwmList = ENCODE,
                       threshold = 1e-4,
                       method = "ic",
                       bkg = c(A=0.25, C=0.25, G=0.25, T=0.25),
                       BPPARAM = BiocParallel::MulticoreParam(5))
```

####  (4) Comparsion of target gene expression between mutated samples(cells) and wild-type
Input: <br>
(1) mutation profile: Combined mutations for samples (or cells) <br>
(2) mutation annotated file: The data needs to contain information about the mutation and its corresponding target gene (nearest neighbor gene or mutation located enhancer linked gene); <br>
(3) gene expression matrix : gene-by-sample matrix or gene-by cell matrix; <br>
Since we don't have matching single-cell data for mutated cells, here is an example of transcriptome data for the mutated samples.
```r
load("./SNV.Rdata")  # mutation-sample information(from step1.mutation detection)
load("./SNV_anno.Rdata") ## mutation-gene information （homer annotation：annotatePeaks.pl SNV.bed hg38 > SNV_anno.txt）
load("./tpm.Rdata")      ## gene-sample information(3 sample replicates per patient)

result<-pblapply(1:nrow(SNV.anno),function(x){
  mut<-SNV.anno[x,1]
  gene<-SNV.anno$Gene.Name[x]
  mut.sample<-SNV.df$sample[SNV.df$ID==mut] %>% unique()
  RNA.sample<-sapply(strsplit(colnames(tpm),split = "_"),function(x){x[[1]]}) %>% unlist()
  mut.exp<-tpm[gene,RNA.sample %in% mut.sample] %>% t() %>% as.vector()
  nonMut.exp<-tpm[gene,!RNA.sample %in% mut.sample] %>% t() %>% as.vector()
  p<-wilcox.test(mut.exp,nonMut.exp)$p.value
  foldchange<-mean(mut.exp)/mean(nonMut.exp)   
  data.frame(gene=gene,mutation=mut,mutSamples=paste(mut.sample,collapse = ";"),pValue=p,foldchange=foldchange)
}) %>% rbindlist()
result$adjsP<-p.adjust(result$pValue,method="fdr")
```





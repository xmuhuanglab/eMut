# eMut
## Introduction
The eMut, an integrated pipeline for detecting, imputing, and characterizing non-coding mutations in CREs with functional consequences at the single-cell level. 

## Workflow

Briefly, eMut consists of two main modules: mutation detection and functional interpretation. <br />
Step1. Mutation detection: eMut detects mutations in each cell by implementing methods such as Monopogen or GATK by single-cell chromatin accessibility data. <br />
Step2. Mutation imputation (optional): Given the sparse of scATAC-seq data, we further imputed candidate mutated cells by network propagation using mutated cells (seed cells) in cell-cell similarity graph. <br />
Step3. Functional interpretation: <br />
1) recognize cell type-specific or lineage-specific mutations; 
2) identify hypermutated CREs with significant excess of mutations to characterize potentially important enhancers; 
3) predict the effects of mutations on transcription factor motifs (loss or gain);
4) compare target gene expression changes between mutated cells (or samples) and wild-type. 
<hr>

![image](https://github.com/xmuhuanglab/eMut/blob/main/Figures/eMut_workflow.png)

## Installation
1.Create a new environment for eMut
```
conda env create -f env.yml
```
2.To install the eMut R package, you can run the following commands in environment of R:
```
devtools::install_github("xmuhuanglab/eMut")
```
3.Dependency packages installation in environment of R:
```
if(!requireNamespace("devtools", quietly = TRUE)){
    install.packages("devtools") 
} 
devtools::install_github("GreenleafLab/ArchR", ref="master", repos = BiocManager::repositories())
devtools::install_github('https://github.com/reimandlab/ActiveDriverWGSR')
```

## Tutorial
The tutorial can be found at [tutorial](https://github.com/xmuhuanglab/eMut/blob/main/eMut_Tutorial.ipynb)

## Contact
For any inquiries or assistance, please feel free to open an issue.

## How to cite eMut
Zhu, M., Zhu, J., Zhu, Z. et al. Discovery of candidate functional non-coding mutations in acute myeloid leukemia using single-cell chromatin accessibility sequencing. Commun Biol 8, 808 (2025). https://doi.org/10.1038/s42003-025-08257-8


#' Cell type enrichment
#'
#' Cell type enrichment of mutant cells(raw or imputed)
#' @param cells barcode of all cells
#' @param cellTypes The cell types of each cell
#' @param mutCells barcode of mutant cells
#' @param mutCells_cutoff This threshold is used to determine how many mutant cells to exceed before cell type
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

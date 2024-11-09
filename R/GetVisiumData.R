GetVisiumData=function(){
    library(Seurat)
    load("inst/extdata/VisiumData.rda")
    return(spat)
}

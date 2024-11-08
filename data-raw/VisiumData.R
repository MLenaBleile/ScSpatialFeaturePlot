library(Seurat)

VisiumData <- readRDS("VisiumData.rds")
devtools::use_data(VisiumData, overwrite = TRUE)

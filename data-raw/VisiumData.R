library(Seurat)

VisiumData <- readRDS("data-raw/VisiumData.rds")
devtools::use_data(VisiumData, overwrite = TRUE)

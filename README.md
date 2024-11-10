ScSpatial feature plot is a function for plotting transcript-level gene expression data processed by Seurat. The function extends the capabilities of the Seurat function SpatialFeaturePlot by allowing the user to plot multiple genes on the same colour scale, thereby facilitating comparison of the magnitude of gene expression across the slide. ScSpatialFeaturePlot also allows the user to rotate the plot, and produces plots faster than SpatialFeaturePlot if you need to subset the data.

ScSpatialFeaturePlot can be installed using the following command:

```
devtools::install_github("MLenaBleile/ScSpatialFeaturePlot")
```

A demo of ScSpatialFeaturePlot is available at the following link:

https://marylenableile.com/ScSpatialFeaturePlot-Vignette.html

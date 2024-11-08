library(ggplot2)
library(Seurat)
ScSpatialFeaturePlot = function(spat, features,slot="counts",facet_labeller=NULL,
                                logs=F,pt.size.factor=1.6, cells=Cells(spat),LegLabel="Count", 
                                uq=.999, lq=0, flip=1, ncol=NULL, nrow=NULL, combine=T){
  if(!combine){
    plotlist=list()
    for(one.feature in features){
      plotlist[[one.feature]]= ScSpatialFeaturePlot(spat, one.feature,slot,facet_labeller,
                                                    logs,pt.size.factor, cells,LegLabel, 
                                                    uq, lq, flip, ncol=1, nrow=1, combine=T)+theme(legend.position="top")
    }
    return(plotlist)
  }
  fnames=features
  MinMaxq = function(X,uq=.99, lq=0.05){MinMax(X, min=quantile(X,lq, na.rm=T), max=quantile(X,uq, na.rm=T))}
  rawcoords=GetTissueCoordinates(spat)[cells,]
  if(flip==1){
    rawcoords$x = -rawcoords$imagecol
    rawcoords$y = rawcoords$imagerow
  }else if(flip==2){
    rawcoords$x = rawcoords$imagerow
    rawcoords$y = rawcoords$imagecol
  }else{
    rawcoords$x = rawcoords$imagecol
    rawcoords$y = -rawcoords$imagerow
  }
  all.plotdf=data.frame()
  spat$intfeat=NA
  for(fname in fnames){
    if(fname %in% colnames(spat@meta.data)){
      spat$intfeat[cells] = MinMaxq(spat[[fname]][cells,1], uq=uq, lq=lq)
                             }else{
        feat = FetchData(spat, slot=slot, vars=fname)[cells,1]
        mmfeat=MinMaxq(feat, uq=uq, lq=lq)
        names(mmfeat)=cells
        spat$intfeat=mmfeat
      }
    plotdf= data.frame(x=rawcoords$x, 
                       y= rawcoords$y,
                       feature=spat$intfeat[rownames(rawcoords)],
                       fname = rep(fname, nrow(rawcoords)))
                       #ftype = rep(names(fname), nrow(rawcoords)),
                       #fcols = color.gradient(spat$intfeat[rownames(rawcoords)],colors=red_grad, colsteps=10000),
                       #cells=factor(Cells(spat)))
    all.plotdf = rbind(all.plotdf, plotdf)
  }

  if(logs){
    all.plotdf$feature=log(all.plotdf$feature+1)

  }
  p1=ggplot(all.plotdf)+geom_point(mapping=aes(x=x, y=y, color=feature), size=pt.size.factor)+
    #scale_color_manual(values=plotdf$scfred)+
    scale_color_viridis_c(option="turbo")+
    theme_bw()+
    
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank())+labs(x="", y="",color=LegLabel )+facet_wrap(~fname, ncol=ncol, nrow=nrow)
  return(p1)
}
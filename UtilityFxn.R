#############################################################
#  Functions for single-cell analysis                       #
#  Lin Zhang linz@mednet.ucla.edu, linzhang26@g.ucla.edu    #
#  linzhangtuesday@gmail.com                                #
#############################################################

# loading packages, source scripts and setting up parallelization
Prepare <- function(verbose=TRUE){
  
  ## clean up environment
  rm(list=ls())
  gc()
  
  ## load packages
  library(dplyr)
  library(Seurat)
  library(patchwork)
  library(philentropy)
  library(useful)
  library(future)
  library(RANN)
  library(tibble)
  library(SeuratData)
  library(harmony)
  library(SeuratWrappers)
  library(cidr)
  library(SingleCellExperiment)
  library(cowplot)
  library(svglite)
  library(networkD3)
  library(tidyverse)
  library(htmlwidgets)
  library(TFBSTools)
  library(clipr)
  library(scater)
  library(ggplot2)
  library(ggstatsplot)
  library(gridExtra)
  library(limma)
  library(reshape2)
  library(xgboost)
  library(mclust)
  library(grDevices)
  library(ggplot2)
  library(RColorBrewer)
  library(ggpubr)
  library(utils)
  library(writexl)
  library(R.utils)
  library(tictoc)
  library(sqldf)
  library(ape)
  library(colorout)
  library(tidyverse)
  library(groupdata2)
  library(stringr)
  library(data.table)
  library(ggrepel)
  library(ggpubr)
  library(scales)
  library(philentropy)
  library(gdata)
  library(clustree)
  library(magrittr)
  library(dplyr)
  
  ## source function script
  source('UtilityFxn.R', encoding = 'UTF-8')
  
  ## Enable Parallelization
  plan()
  plan("multiprocess", workers = 16)
  options(future.globals.maxSize = 8000*2048^2,future.seed = TRUE,future.rng.onMisue = "ignore")
  plan() %>% print()
  set.seed(1)

}

# normalization, find variable genes, run pca and pca number selection 
NormalHVFPCAJack <- function(data=P0Marmoset,
                            file="P0Marmoset",
                            nfeatures=2000,
                            Jack=TRUE,
                            Jdim=50,
                            save=FALSE,
                            regress=FALSE,
                            regressv=c("percent.mt","percent.ribo","nCount_RNA"),
                            verbose=TRUE)
{ if(verbose)(print("Normalize and FindVariableFeatures,default is scaling all genes"))
  if(verbose)(print("Jdim currently not support over 100,but can be modified"))
  
  data <- NormalizeData(object = data, normalization.method = "LogNormalize", scale.factor = 1e4)
  data <- FindVariableFeatures(object = data,selection.method = 'vst', nfeatures = nfeatures)
  
  if(verbose)(print("Scale and Dimensional reduction using Scale data and PCA"))

  all.genes <- rownames(data)
  
  if(regress){
  data <- ScaleData(object=data, vars.to.regress = regressv,features = all.genes)
  }else
  {data <- ScaleData(object=data,features = all.genes)}

  data <- RunPCA(object = data,npcs=100) 
  
  if(verbose)(print("ElbowPlot"))
  
  plotname = paste0("PLOT/",file,"elbow.png")
  png(plotname,h=300,w=300)
  p1 <- ElbowPlot(object = data,ndims = Jdim)
  print(p1)
  dev.off()
  
  if(verbose)(print("DimHeatmap"))
  plotname = paste0("PLOT/",file,"dimheatmap.png")
  png(plotname,h=2000,w=300)
  p2 <- DimHeatmap(data, dims = 1:Jdim, cells = 500, balanced = TRUE) 
  print(p2)
  dev.off()
  
  if(Jack){
  if(verbose)(print("Jackstraw"))
  data=JackStraw(data,num.replicate=100,dims=Jdim)
  data <- ScoreJackStraw(data,dims=1:Jdim)
  
  plotname = paste0("PLOT/",file,"jackstraw.png")
  png(plotname,h=1000,w=1500)
  p3 <- JackStrawPlot(data,dims=1:Jdim)
  print(p3)
  dev.off()
  }
  if(save){
    filename = paste0("RDS/",file,"nHVG",nfeatures,"normalizedHVF_PCA.rds")
    saveRDS(data,file = filename)
  }
  
  return(data)
}

# batch correction and clustering
ClusteringSeurat <- function(data=P0Marmoset,
                               dim = 12,
                               reduction = "harmony", 
                               file = "P0Marmoset",
                               resolution=3,
                               min.dist=0.5,
                               n.neighbors=30,
                               batch = "orig.ident",
                               featureplot=FALSE,
                               tsne=FALSE,
                               assay="RNA",
                               cell="BC",
                               meta1="orig.ident", # factor 1 in metadata columns for contingency table
                               meta2="seurat_clusters",  # factor 2 in metadata columns for contingency table
                               hd=1800, # dim for dotplot
                               wd=1700, # dim for dotplot
                               hdi=2400, # dim for dimplot
                               wdi=1700, # dim for dimplot
                               features = c("RBPMS","SLC17A6","RBPMS2","POU4F2","POU4F3","SIX6","SPP1","GUCY1A3","EOMES","KCNA1","CA8","TBR1","TPBG","THY1","NEFL","NEFM","SNCG","CHTF18","VSX2","TMEM215","NVSX2","OTX2","GRM6","PRKCA","TRPM1","GRIK1","ISL1","VSX1","NXPH2","NXPH1","CABP5","APOE","GLUL","CRABP1","CLU","SLC1A3","DKK3","CRYM","RLBP1","CRYAB","GNGT2","HCN1","GUCA1C","ARR3","RCVRN","GNAT2","OPN1SW","LOC100412476","PDE6H","GUCA1A","SAG","RHO","PDC","NRL","GNAT1","GNGT1","NR2E3","GNB1","TFAP2B","TFAP2A","TFAP2C","GAD1","GAD2","GLYT1","SLC6A9","C1QL1","C1QL2","PAX6","SLC32A1","MEIS2","TCF4","CHAT","FEZF1","RND3","SLC18A3","SLC5A7","SOX2","SST","CRH","NPY","EBF3","NEUROD6","NEUROD3","CALB2","RET","SEPT4","PTN","CHN1","PCDH11X","TMOD1","ONECUT2","LHX1","SLC12A7","ONECUT1","CALB1","MGP","MYL9","COL4A1","C1QA","C1QB","C1QC","HEXB","CTSS","P2RY12","TMEM119","B2M","CLDN5","IGFBP7","RGS5","GSN","FN1","ZFX","SRY","SOX2","S100B","GFAP","ALDH1L1","ALDOC","GLT1","AQP4","TTR"), #major retina class marker panel
                               max.iterH = 10,
                               max.iterC = 20,
                               scaled=FALSE, # whether scale in dotplot
                               CEX=0.5,
                               Pointsize=NULL,
                               plot=TRUE,
                               save=FALSE,
                               verbose=TRUE 
){if(verbose)(print("adjust pc # accordingly"))
  features = features %>% unique()
  
  if (reduction  == "mnn"){
    if(verbose)(print("Scale and Dimensional reduction using FastMNN"))
    data <- RunFastMNN(object.list = SplitObject(data, split.by = batch),reduction.name = reduction)
  }
  
  
  if (reduction  == "harmony"){
    if(verbose)(print("Scale and Dimensional reduction using Harmony"))
    if(verbose)(print("require run reduction=pca first,otherwise do ScaleData and RunPCA first"))
    data <- RunHarmony(data,group.by.vars = batch, dims.use=1:dim,max.iter.harmony = max.iterH, max.iter.cluster = max.iterC) 
  }

  if(verbose)(print("Find Neighbors and clustering"))
  data <- FindNeighbors(object = data, dims = 1:dim,nn.eps = 0.5,reduction = reduction)
  data <- FindClusters(object = data, resolution = resolution, n.start = 10)
  
  if(tsne){
  if(verbose)(print("Fit-SNE takes long and high memory"))
  data <- RunTSNE(data, dims = 1:dim,nthreads = 8,reduction = reduction,max_iter = 2000,tsne.method = "FIt-SNE",late_exag_coeff=1.5,learning_rate=8000) 
  if(plot){
    plotname = paste0("PLOT/",file,cell,reduction,"_res_",resolution,batch,"tsne.png")
    png(plotname,h=1000,w=1300)
    p1 <- DimplotclusterSeurat(object = data,file = file,reduction = "tsne",Resolution = resolution)
    print(p1)
    dev.off()
  }
  }
  
  data <- SetIdent(data,value=paste0("RNA_snn_res.",resolution))
  
  if(verbose)(print("UMAP"))
  reductioname = paste0(reduction,cell,"umap")
  reductionkey = paste0(reductioname,"_")
  data <- RunUMAP(object = data, reduction = reduction,reduction.name= reductioname,reduction.key = reductionkey,dims = 1:dim, min.dist = min.dist,n.neighbors =n.neighbors)

  if(plot)
  {plotname=paste0("PLOT/",file,cell,"_res_",resolution,"dim_",dim,batch,"_vln.png")
   png(plotname,h=3000,w=4000)
   p1 <- VlnPlot(object = data, features = c("nFeature_RNA","nCount_RNA", "percent.mt","percent.ribo","RBPMS"),pt.size = 0.001, ncol = 1)+NoLegend() 
   print(p1)
   dev.off()
   }
  
  if(plot){
    plotname = paste0("PLOT/",file,reductioname,"_res_",resolution,batch,"dim_",dim,".png")
    png(plotname,h=1500,w=1300)
    p2 <- DimplotclusterSeurat(object = data,file = file,reduction = reductioname,h = hdi,w = wdi,Resolution = resolution,point.size = Pointsize)
    print(p2)
    dev.off()
  }
  
  if(plot){
    plotname = paste0("PLOT/",file,cell,batch,"_res_",resolution,"dotplot.png")
    png(plotname,h=hd,w=wd)
    p1 <- DotPlot(data, features = features,cluster.idents = TRUE, scale = scaled) +coord_flip()+ RotatedAxis()
    print(p1)
    dev.off()
  }

  if(save){
  filename = paste0("RDS/",file,assay,reduction,"clustering.rds")
  saveRDS(data,file = filename)}
  
  return(data)
}

# eliminate certain low quality clusters
ElilminateSeurat <- function(object=P0Marmoset,
                             file="P0Marmoset",
                             filter=c(0),
                             filtercell="RGC",
                             num=1,
                             verbose=TRUE
){if(verbose){print(paste0("Eliminate cluster", filter-1))}
  cluster <- Idents(object) %>% levels() %>% as.vector()
  
  ## save filtered object
  filetername = paste0("RDS/test2P0/Filtered/",file,"_",num,"_",filtercell,".rds")
  filtered <- subset(object,ident=cluster[filter])
  saveRDS(filtered,file=filetername)
  
  ## filter from object
  object <- subset(object,ident=cluster[-filter])
  
  return(object)
}

# resolution optimization using clustering trees
ClusterTree <- function(object=P0,
                        filename="P0",
                        RES=seq(0,6,0.1),
                        hs=30, # height of plot
                        ws=15, # width of plot
                        verbose=TRUE
                        )
{ # delete previous resolution info
  clusters <- object@meta.data %>% dplyr::select(grep("^RNA_snn_res.",names(object@meta.data))) %>% names()
  object@meta.data <- object@meta.data[,-which(colnames(object@meta.data) %in% clusters)]
  
  object= FindClusters(object,resolution = RES,n.start = 10)
  
  image = clustree(object,prefix = "RNA_snn_res.")
  ggsave(file=paste0("PLOT/",filename,"clustree.svg"), plot=image, width=ws, height=hs)
  ggsave(file=paste0("PLOT/",filename,"clustree.png"), plot=image, width=ws, height=hs)
}

# plot dimplots colored by meta data columns for batch correction effect check
DimplotclusterSeurat <- function(object=P0Marmoset,
                                 reduction ="tsne",
                                 file = "P0Marmoset",
                                 Resolution=0.3,
                                 h=2500,
                                 w=2500,
                                 hs=15,
                                 ws=15,
                                 point.size =NULL,
                                 label.size=5,
                                 full=TRUE)
{
  p1 <- DimPlot(object,reduction = reduction,label = TRUE, pt.size =point.size,repel = TRUE,label.size = label.size) +ggtitle(paste0("Cell Type"))+theme(plot.title = element_text(size = 40, face = "bold",hjust = 0))
  
  # the objects input should have corresponding meta data columns
  if(full){
    p3 <- DimPlot(object,reduction = reduction,label.size = label.size, pt.size =point.size,group.by = "orig.ident") +ggtitle(paste0("Sample ID"))+theme(plot.title = element_text(size = 40, face = "bold",hjust = 0))
    p5 <- DimPlot(object,reduction = reduction,group.by = "Chemistry", pt.size =point.size,label.size = label.size,cols = c('v1'="orange",'v2' = 'blue', 'v3' = 'grey',"v3.1"="pink")) +ggtitle(paste0("Chemistry"))+theme(plot.title = element_text(size = 40, face = "bold",hjust = 0)) 
    p6 <- DimPlot(object,reduction = reduction,group.by = "TCD", pt.size =point.size,label.size = label.size,cols = c('nocd' = 'grey','cd73' ='orange', 'cd90' = 'blue',"cd90cd73"='red')) +ggtitle(paste0("Antibody ID"))+theme(plot.title = element_text(size = 40, face = "bold",hjust = 0)) 
    p4 <- DimPlot(object,reduction = reduction,pt.size =point.size,label.size = label.size,group.by = "Region") +ggtitle(paste0("Region ID"))+theme(plot.title = element_text(size = 40, face = "bold",hjust = 0)) 
    p2 <- DimPlot(object,reduction = reduction,label = TRUE,pt.size =point.size,repel = TRUE,label.size = label.size,group.by = "seurat_clusters") +ggtitle(paste0("Cluster ID"))+theme(plot.title = element_text(size = 40, face = "bold",hjust = 0)) + NoLegend()
    p7 <- DimPlot(object,reduction = reduction,pt.size =point.size,label.size = label.size,group.by = "stage",cols = c('P0' = 'green',"Adult"="orange")) +ggtitle(paste0("Stage ID"))+theme(plot.title = element_text(size = 40, face = "bold",hjust = 0)) 
    
    p8 <- ggarrange(p2,p3,p4,p5,p6,p7, ncol=2,nrow=3)+ ggtitle(paste0(file)) + theme(plot.title = element_text(size = 60, face = "bold",hjust = 0.5))
 
    picname = paste0("PLOT/",file,reduction,Resolution,"Dimplot.png")
    png(picname,h=h,w=w) 
    print(p8)
    dev.off()
    
    ggsave(file=paste0("PLOT/",file,reduction,Resolution,"Dimplot.svg"), plot=p8, width=ws, height=hs)
  }else{
    picname = paste0("PLOT/",file,reduction,Resolution,"Dimplot.png")
    png(picname,h=600,w=800)
    print(p1)
    dev.off()
    
    ggsave(file=paste0("PLOT/",file,reduction,Resolution,"Dimplot.svg"), plot=p1,width=ws, height=hs)
  }
}

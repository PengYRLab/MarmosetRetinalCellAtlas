#############################################################
#  Demo for iterative single-cell analysis                  #
#  using functions from UtilityFxn.R                        #
#  Lin Zhang  linzhang26@g.ucla.edu                         #
#  linz@mednet.ucla.edu, linzhangtuesday@gmail.com          #
#############################################################

# set up the analysis environment, please create a /PLOT /DOC /RDS folder in your working directory for plots, metadata documents and objects respectively
source('UtilityFxn.R', encoding = 'UTF-8')
Prepare()

# use RGC as an example to showcase condensed clustering pipeline with batch correction implemented
# this is an iterative procedure, please refer to dimplot, dotplot and vlnplot at each round to assess the cell quality 
RGC_F1 <- ReadRDS(file="RGC_F1.rds")

# Round 01 of clustering and filtering
RGC_F1 <- NormalHVFPCAJack(RGC_F1,file="RGC_F1",Jack=TRUE,Jdim=80,regress=TRUE)
## choose dim based on Jackstraw plot in PLOT/
RGC_F1 <- ClusteringSeurat(RGC_F1,file="RGC_F1_test",batch="orig.ident",dim=40,cell="RGC")
## check full marker dotplot, dimplot and vlnplot for low-quality cell clusters
##  eliminate low-quality cell clusters
RGC_F2 <- ElilminateSeurat(RGC_F1,filter=c(9,16),filtercell="lowqua",num=2)


# Round 02 of clustering and filtering
# After filtering, rerun normalization, PCA et al
RGC_F2 <- NormalHVFPCAJack(RGC_F2,file="RGC_F2",Jack=FALSE,regress=TRUE)
RGC_F2 <- ClusteringSeurat(RGC_F2,file="RG C_F2_test",batch="orig.ident",dim=40,cell="RGC")
## check full marker dotplot, dimplot and vlnplot for low-quality cell clusters
# eliminate low quality cell clusters
RGC_F3 <- ElilminateSeurat(RGC_F2,filter=c(29),filtercell="lowqua",num=3)


# Round 03 of clustering and filtering
RGC_F3 <- NormalHVFPCAJack(RGC_F3,file="RGC_F3",Jack=FALSE,regress=TRUE)
## use previous PC from Jackstraw, here use 50 for demo
RGC_F3 <- ClusteringSeurat(RGC_F3,file="RGC_F3_test",batch="orig.ident",dim=50 ,cell="RGC")
## eliminate low quality cell clusters
RGC_F4 <- ElilminateSeurat(RGC_F3,filter=c(16),filtercell="lowqua",num=4)


# Round 04 of clustering and filtering
RGC_F4 <- NormalHVFPCAJack(RGC_F4,file="RGC_F4",Jack=FALSE,regress=TRUE)
RGC_F4 <- ClusteringSeurat(RGC_F4,file="RGC_F4_test",batch="orig.ident",dim=50 ,cell="RGC")
# within cluster low-quality cell filtering by checking vlnplot of nFeature_RNA, percentMT, percent_Ribo
low.det.v2 <- WhichCells(RGC_F4,expression = nFeature_RNA > 5000 & (seurat_clusters %in% c(0,1,32)))
# check how many cells are being filtered
length(low.det.v2)
RGC_F5 <- subset(RGC_F4,cells = setdiff(WhichCells(RGC_F4),c(low.det.v2)))

# Round 05 of clustering and filtering
# Rerun starting from normalization 
RGC_F5 <- NormalHVFPCAJack(RGC_F5,file="RGC_F5",Jack=FALSE,regress=TRUE)
# now we are happy with the cell quality, let's check potential optimal resolution for the data set
ClusterTree(RGC_F5,file="RGC_F5_test")

# Final Round of clustering
# re-clustering with optimized resolution
RGC_F5 <- ClusteringSeurat(RGC_F4,file="RGC_F5_test",batch="orig.ident",dim=50,resolution=1.7,cell="RGC")

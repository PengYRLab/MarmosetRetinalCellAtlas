
#' @title  signed Gene Set Enrichnment Analyisis(sGSEA)
#' @description This function generates weighted Stouffer integrated gene signatures using both up-regulated and down-regulated genes
#' @param  seu.query The seurat object of the query data
#' @param  seu.ref The seurat object of the reference data
#' @param  sig.method "two.sided", "up" or "down"
#' @param  top.nn the number of top (up/down) changed features used
#' @author  Junqiang Wang
#' @export
sGSEA<-function(seu.query,
                seu.ref,
                top.nn=100,
                sig.method="two.sided",
                rank.transform=FALSE,
                q.lb=0.01,
                wt="log2FC",
                de.method="MAST",
                orderDEG.method="p_val",
                redo.sct=FALSE,
                vars.to.regress=NULL,
                # vars.to.regress = c("percent.mt", "percent.ribo", "nCount_RNA"),
                de.assay="SCT",
                de.slot="counts",
                signature.z.log2CPM=FALSE
                ){
  
  require(dplyr)
  require(tidyr)

seu.query@meta.data$sig.condition<-"1" %>% as.factor()
seu.ref@meta.data$sig.condition<-"-1" %>% as.factor()

seu<-merge(x=seu.query, y=seu.ref)


# get the signature

if(isTRUE(redo.sct)){
  
DefaultAssay(seu)<-"RNA"

seu<-SCTransform(seu, vars.to.regress = vars.to.regress, verbose = FALSE, method="glmGamPoi")

}


# use log2CPM

if(isTRUE(signature.z.log2CPM)){
  
dat<-seu@assays$RNA@counts %>% as.matrix()

dat<-Log2CPM(dat)

} else{
  
  dat<-seu@assays$SCT@data %>% as.matrix()
}

# 

dat.1<-dat[,which(seu@meta.data$sig.condition=="1")]

dat.2<-dat[,which(seu@meta.data$sig.condition=="-1")] 


# The 
if(isFALSE(rank.transform)){

mean.ref<-apply(dat.2, 1, mean)

sd.ref<-apply(dat.2, 1, sd)

z<-apply(dat.1, 2, FUN=function(x, mean.ref, sd.ref){
  x<-(x-mean.ref)
  x<-x/sd.ref
}, mean.ref=mean.ref, sd.ref=sd.ref)

z<-z %>% na.omit() %>% as.matrix()

}

# rank transformation
if(isTRUE(rank.transform)){
  
  q.mat<-matrix(data=NA, nrow=nrow(dat.1), ncol=ncol(dat.1))
  rownames(q.mat)<-rownames(dat.1)
  colnames(q.mat)<-colnames(dat.1)
  

  for(i in 1:nrow(dat.1)){
    
    percentile<-ecdf(dat.2[i,])
    
    q.mat[i,]<-percentile(dat.1[i,])
    
  }
  
  
  q.mat[q.mat<q.lb]<-q.lb
  
  q.mat[q.mat>1-q.lb]<-1-q.lb
  
  z<-qnorm(q.mat)

}

# The

# do DE analysis


Idents(seu)<-"sig.condition"

DefaultAssay(seu)<-"RNA"

markers.all <- FindAllMarkers(seu, 
                              assay=de.assay,
                              test.use=de.method,
                              only.pos = TRUE, 
                              min.pct = 0.25, 
                              logfc.threshold = 0.25,
                              slot=de.slot,
                              latent.vars = c("percent.mt", "percent.ribo", "nCount_RNA")
                            #  recorrect_umi = FALSE
                              )


if(orderDEG.method=="p_val"){

if(sig.method=="two.sided"){
  
  markers.all<-markers.all %>%
    arrange(p_val)
  
}else if(sig.method=="up"){
  
  markers.all<-markers.all %>%
    dplyr::filter(cluster==1) %>%
    arrange(p_val)
  
}else if(sig.method=="down"){
  
  markers.all<-markers.all %>%
    dplyr::filter(cluster==-1) %>%
    arrange(p_val)
  
}
  
}


if(orderDEG.method=="avg_log2FC"){
  
  if(sig.method=="two.sided"){
    
    markers.all<-markers.all %>%
      arrange(desc(avg_log2FC))
    
  }else if(sig.method=="up"){
    
    markers.all<-markers.all %>%
      dplyr::filter(cluster==1) %>%
      arrange(desc(avg_log2FC))
    
  }else if(sig.method=="down"){
    
    markers.all<-markers.all %>%
      dplyr::filter(cluster==-1) %>%
      arrange(desc(avg_log2FC))
    
  }
  
}





# filter out the markers with NA z

markers.all<-markers.all %>%
  dplyr::filter(gene %in% rownames(z))

tmp<-grep("^MT-|^RPL|^RPS", markers.all$gene)



if(length(tmp)>0){

markers.all<-markers.all[-tmp, ]

}


markers.all<-markers.all %>% as_tibble()

print(markers.all, n=100)


# remove mt, rpl, rps gene


top.nn<-min(top.nn, nrow(markers.all))

stat<-markers.all[1:top.nn,]

z<-z[match(stat$gene, rownames(z)),]

#nes

if(is.null(wt)==FALSE){
  
wt<-stat$avg_log2FC * as.numeric(unfactor(stat$cluster)) %>% as.matrix() %>% t()

} else if(is.null(wt)){
  
  wt<-as.numeric(unfactor(stat$cluster)) %>% as.matrix() %>% t()
  
}


print(wt)

dim(wt)

dim(z)

es <- wt %*% z

es.sd<-wt^2 %>% sum() %>% sqrt()

nes<-es/es.sd

rownames(nes)<-"nes"

p<-pnorm(nes, lower.tail = FALSE)

p.fdr<-p.adjust(p, method="fdr")

p.fdr.log10<-log10(p.fdr)*(-1)

zGSEA<-list(nes=nes, p.fdr=p.fdr, p.fdr.log10=p.fdr.log10)

return(zGSEA)

}





#' @title Log2CPM
#' @description This function convert count matrix to log2CPM matrix
#' @param count The count matrix
#' @return The log2CPM matrix
#' @author Junqiang Wang
#' @export
Log2CPM<-function(count){
  
  count <- t(t(count) / (colSums(count) / 1e6))
  count<-log2(count+1)
  
  count<-round(count, digits=2)
  
  return(count)
  
}




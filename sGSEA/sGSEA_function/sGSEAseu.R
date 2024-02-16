


#' @description This functon performs the sGSEA for all cell types
#' @author Junqiang Wang
#' @export
sGSEAseu<-function(seu, 
                   idents.seu=NULL,
                   idents.query=NULL,
                   idents.ref=NULL,
                   idents.class=NULL,
                   top.nn=100,
                   cells.n.cutoff=10,
                   wt.gsea=NULL,
                   rank.transform=FALSE,
                   redo.sct=FALSE,
                   signature.z.log2CPM=TRUE,
                   vars.to.regress=NULL,
                   de.method='MAST',
                   de.assay="RNA",
                   de.slot="counts"
                   ){


seu.0<-seu

seu.0<-SCTransform(seu.0, verbose = FALSE, method="glmGamPoi")

#
# class<-seu.0@meta.data$class %>% sort() %>% unique()
#

tmp<-seu.0@meta.data

class<-tmp[, match(idents.class, colnames(tmp))] %>% sort() %>% unique() %>% as.vector()


#

p.list<-vector(mode = "list", length=length(class))

names(p.list)<-class


for(i in 1:length(class)){
  
  message("processing cell type: ", class[i], "  id: ", i)
  
  Idents(seu.0)<-idents.class
  
  seu<-subset(seu.0, idents=class[i])
  
  # if not paired, pass
  # tmp<-seu@meta.data$Tissue %>% unique() %>% sort() %>% as.vector()
  # 
  # if(!identical(tmp, c("foveal", "periphery"))){
  #   
  #   message(class[i], " is not included...")
  #   
  #   p.list[[i]]<-NULL
  #   
  #   next
  #   
  # }else{
  #   
  #   message(class[i], " passed the 1-round filter ...")
  #   
  # }
  
  
  Idents(seu)<-idents.seu
  
  seu.query<-subset(seu, idents=idents.query)
  
  seu.ref<-subset(seu, idents=idents.ref)
  
  min.cells<-min(ncol(seu.ref), ncol(seu.query))
  
  if(min.cells<cells.n.cutoff){
    
    message(class[i], " is not included...")
    
    p.list[[i]]<-NULL
    
    next
    
  }else{
    
    message(class[i], " passed the 2-round filter ...")
    
  }
  
  
  # test
  stat<-sGSEA(seu.query=seu.query,
              seu.ref=seu.ref,
              rank.transform=rank.transform,
              vars.to.regress=vars.to.regress,
              sig.method="two.sided",
              wt=wt.gsea,
              de.method=de.method,
              orderDEG.method="p_val",
              redo.sct=redo.sct,
              de.assay=de.assay,
              de.slot=de.slot,
              signature.z.log2CPM= signature.z.log2CPM,
              top.nn=top.nn)
  
  
  stat
  
  df<-do.call(rbind, stat) %>% t() %>% data.frame()
  
  df$class<-rep(class[i], length=nrow(df))
  
  df$cell.id<-colnames(stat$nes)
  
  df<-df %>% data.frame()
  
  colnames(df)<-c("nes", "p.fdr", "p.fdr.log10", "class", "cell.id")
  
  p.list[[i]]<-df
  
}



df<-do.call(rbind, p.list)

df$class<-df$class %>% as.factor()

df<-df %>% as.data.frame()

return(df)

}





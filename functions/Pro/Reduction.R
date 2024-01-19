

#### PCA/UMAP OF INTEGRATED DATA

reduction_plot <- function(object=object_data,reduced=F){
  cat('\nLaunch reduction analysis....\n')
  ### create main folder for reduction plot
  if(reduced==F){
    plot_save_dir <- c(dir_create('PLOTs','All',''),dir_create('PLOTs',object$folder,'Group') )
  }else{
    plot_save_dir <- c(dir_create('reduced_PLOTs','All',''))
  }
  
  reduction_plot_batch <- function(object_,plot_save_dir){
    object_ <- ScaleData(object_)
    object_ <- RunPCA(object_, features = VariableFeatures(object = object_),npcs=UMAP_max_dims )
    object_ <- RunUMAP(object_,  dims = 1:  UMAP_max_dims)
    object_ <- FindNeighbors(object_, dims = 1:  UMAP_max_dims)
    object_ <- FindClusters(object_, resolution = UMAP_resolution)
    
    ########## PCA main plot
    cat('\nLaunch PCA plotting....\n')
    lapply(plot_save_dir,function(x) subdir_create(paste0(x,'/PCA')) )%>% invisible()
    pdf(paste0(plot_save_dir,'/PCA/',Project(object_),'_PCA.pdf'))
    print(DimPlot(object_, reduction = "pca"))
    dev.off()
    
    ######## UMAP main plot
    cat('\nLaunch UMAP plotting....\n')
    lapply(plot_save_dir,function(x) subdir_create(paste0(x,'/UMAP')) )%>% invisible()
    pdf(paste0(plot_save_dir,'/UMAP/',Project(object_),'_UMAP.pdf'))
    print(DimPlot(object_, reduction = "umap",label = 'TRUE'))
    dev.off()
    
    ######### within group sample grouped 
    if( length( unique(object_$group)) ==1 & length(unique(object_$orig.ident))>1  ){    
      #UMAP
      pdf(paste0(plot_save_dir,'/UMAP/',Project(object_),'_UMAP_sample_grouped.pdf'))
      print(DimPlot(object_, reduction = "umap",label = 'FALSE', group.by = 'orig.ident'),width=7)
      dev.off()
      #PCA
      pdf(paste0(plot_save_dir,'/PCA/',Project(object_),'_PCA_sample_grouped.pdf'))
      print(DimPlot(object_, reduction = "pca", group.by = 'orig.ident'),width=7*1.618 )
      dev.off()
    }
    ######### within group sample split 
    if( length( unique(object_$group)) ==1 & length(unique(object_$orig.ident))>1  ){
      #UMAP
      pdf(paste0(plot_save_dir,'/UMAP/',Project(object_),'_UMAP_sample_grouped.pdf'))
      print(DimPlot(object_, reduction = "umap",label = 'FALSE', split.by = 'orig.ident'),width=7*length(unique(object_$orig.ident)))
      dev.off()
      #PCA
      pdf(paste0(plot_save_dir,'/PCA/',Project(object_),'_PCA_sample_split.pdf'))
      print(DimPlot(object_, reduction = "pca", split.by = 'orig.ident'),width=7*length(unique(object_$orig.ident)) )
      dev.off()
    }
    ######## inter-group grouped 
    if( length( unique(object_$group)) >1  ){   
      #UMAP
      pdf(paste0(plot_save_dir,'/UMAP/',Project(object_),'_UMAP_group_grouped.pdf'),width=7)
      print(DimPlot(object_, reduction = "umap",label = 'FALSE', group.by = 'group'))
      dev.off()
      #PCA
      pdf(paste0(plot_save_dir,'/PCA/',Project(object_),'_PCA_group_grouped.pdf'))
      print(DimPlot(object_, reduction = "pca", group.by = 'group'),width=7*1.618 )
      dev.off()
    }
    ######## inter-group split 
    if( length( unique(object_$group)) >1  ){       
      pdf(paste0(plot_save_dir,'/UMAP/',Project(object_),'_UMAP_group_split.pdf'),width=7*length(unique(object_$group)))
      print(DimPlot(object_, reduction = "umap",label = 'FALSE', split.by = 'group'))
      dev.off()
      #PCA
      pdf(paste0(plot_save_dir,'/PCA/',Project(object_),'_PCA_group_split.pdf'))
      print(DimPlot(object_, reduction = "pca", split.by = 'group'),width=7*length(unique(object_$group)) )
      dev.off()
    }
    return(object_)
  }
  if(reduced==F){
    return(sapply(names(object$object),function(x) reduction_plot_batch(object$object[[x]],plot_save_dir[which(basename(plot_save_dir)==x)]) ,USE.NAMES=T))
  }else{
    return(sapply('All',function(x) reduction_plot_batch(object$reduced[[x]],plot_save_dir[which(basename(plot_save_dir)==x)]) ,USE.NAMES=T))
  }
}


#####  heatmap and other plottings

multi_plot <- function(object=object_data,reduced=F,extra_features=c()){
  cat('----------------------------------Start multi plotting------------------------------------')
  n=length(plot_feature)
  if(reduced==F){
    save_dir <- c(dir_create('PLOTs','All',''),dir_create('PLOTs',object$folder,'Group') )
  }else{
    save_dir <- c(dir_create('reduced_PLOTs','All','') )
  }
  
  heatmap_plot_batch <- function(object_=object_data$object$KP3,plot_save_dir){
    cat('\nLaunch heatmap plotting for dimensions ....\n')
    lapply(plot_save_dir,function(x) subdir_create(paste0(x,'/heatmap'))) %>% invisible()
    
    ############ variable heatmap
    pdf(paste0(plot_save_dir,'/heatmap/',Project(object_),'_cluster_variable_marker_top_100_heatmap.pdf'),height=12,width=0.6*length(unique(object_$seurat_clusters)  ))
    p<-DoHeatmap(object_,features = VariableFeatures(object_)[1:100],size = 4, angle = 90) +scale_fill_gradientn(colours = rev(colorRampPalette(brewer.pal(10, "RdYlBu"))(256) ))
    print(p)
    dev.off()
    ############# dim heatmap
    htmap <- function(dim_=1){
      if(dim_ < 50){
        pdf(paste0(plot_save_dir,'/heatmap/',Project(object_),'_dim_from_',dim_,'_to_',min(50,dim_ + dimension_plotting_intervals-1),'_heatmap.pdf'))
        p<-DimHeatmap(object = object_, 
                      dims = dim_ : min(50,dim_ + dimension_plotting_intervals-1), 
                      cells = n_cells_for_plot, 
                      balanced = TRUE)
        ##+scale_fill_gradientn(colours = rev(colorRampPalette(brewer.pal(10, "RdYlBu"))(256) ))
        print(p)
        dev.off()
        dim_ = dim_ + dimension_plotting_intervals
        return(htmap(dim_))}
    }
    htmap(dim_=1)
    
    cat('\nLaunch heatmap plotting for variable features....\n')
    
    
    ################ Elbowplot
    cat('\n---Launch Elbowplot plotting....')
    lapply(plot_save_dir,function(x) subdir_create(paste0(x,'/Elbowplot')) )%>% invisible()
    pdf(paste0(plot_save_dir,'/Elbowplot/',Project(object_),'_Elbowplot_PCA.pdf'))
    print(ElbowPlot(object_, ndims = ElbowPlot_ndims))
    dev.off()
    
    ######## switch assay to original data
    DefaultAssay(object_) <- 'RNA'
    
    ##################### individual vln plot
    
    cat('\n---Launch Vln plotting....')
    lapply(plot_save_dir,function(x) subdir_create(paste0(x,'/Vln_plot')) )%>% invisible()
    individual_vln_plot <- function(single_feature){
      cat(paste0('\nPlotting Featture : ',single_feature ))
      
      pdf(paste0(plot_save_dir,'/Vln_plot/',Project(object_),'_',single_feature,'_Vln_plot.pdf'),width = max(7,1.618*length(plot_feature)))
      print(VlnPlot(object_ , features = single_feature))
      dev.off()
      
      if( length(unique(object_$group)) >1  ){
        lapply(plot_save_dir,function(x) subdir_create(paste0(x,'/Vln_plot_split')) )%>% invisible() 
        pdf(paste0(plot_save_dir,'/Vln_plot_split/',Project(object_),'_',single_feature,'_Vln_plot_split.pdf'),width = max(7,2.2*length(plot_feature)))
        print(VlnPlot(object_ , features = single_feature,split.by = 'group',split.plot = F))
        dev.off()
      }
      
    }
    lapply(c("nFeature_RNA",plot_feature),individual_vln_plot)%>% invisible()
    
    ############### Feature plot individually for each signature
    cat('\n---Launch Feature plotting....')
    feature_plot <- function(data=object_,feature = plot_feature ){
      cat(paste0('\nplotting Feature : ',feature))
      lapply(plot_save_dir,function(x) subdir_create(paste0(x,'/feature_plot')) )%>% invisible()
      
      pdf(paste0(plot_save_dir,'/feature_plot/',Project(object_),'_',feature,'_feature_plot.pdf'),width = 7)
      print(FeaturePlot(object_, features = feature))
      dev.off()
      
      if( length(unique(object_$group)) >1  ){
        lapply(plot_save_dir,function(x) subdir_create(paste0(x,'/feature_plot_split')) )%>% invisible()
        pdf(paste0(plot_save_dir,'/feature_plot_split/',Project(object_),'_',feature,'_feature_plot_split.pdf'),width = 7*length(unique(object_$group)))
        print(FeaturePlot(object_, features = feature,split.by = 'group'))
        dev.off()
        
      }
    }
    lapply(c("nFeature_RNA",plot_feature), function(x) feature_plot(feature=x) )%>% invisible()
    
    ######### total feature plot
    pdf(paste0(plot_save_dir,'/feature_plot/',Project(object_),'_total_feature_plot.pdf'),width = 3*as.integer(n^0.5),height = 3*as.integer(n^0.5))
    print(FeaturePlot(object_, features = plot_feature))
    dev.off()
    
    ## dot plot
    cat('\n---Launch dot plotting....')
    lapply(plot_save_dir,function(x) subdir_create(paste0(x,'/dot_plot')) )%>% invisible()
    pdf(paste0(plot_save_dir,'/dot_plot/',Project(object_),'_dot_plot.pdf'))
    print(DotPlot(object_, features = plot_feature,scale = F) + RotatedAxis())
    dev.off()
    
    ######### dot plot for all
    if( length(unique(object_$group)) >1  ){ 
      cat('\n---Launch grouped_dot plotting....\n')
      lapply(plot_save_dir,function(x) subdir_create(paste0(x,'/dot_plot')) )%>% invisible()
      pdf(paste0(plot_save_dir,'/dot_plot/',Project(object_),'_dot_plot_group_grouped.pdf'))
      print(DotPlot(object_, features = plot_feature, group.by = "group",scale = F) + RotatedAxis())
      dev.off()
      ############# split dot plot
      pdf(paste0(plot_save_dir,'/dot_plot/',Project(object_),'_dot_plot_group_split.pdf'),width = 1*length(plot_feature),height = 1*length(unique(object_$seurat_clusters)))
      print(DotPlot(object_, features = plot_feature, split.by  = "group",scale = F) + RotatedAxis())
      dev.off()
    }
    return(object_)
  }
  if(reduced==F){
    sapply(names(object$object),function(x) heatmap_plot_batch(object$object[[x]],save_dir[which(basename(save_dir)==x)]))
  }else{
    plot_feature <- c(plot_feature,extra_features)
    sapply('All',function(x) heatmap_plot_batch(object$reduced[[x]],save_dir[which(basename(save_dir)==x)]))
  }
  cat('\n------------------multi Plotting finished----------------------\n')
}
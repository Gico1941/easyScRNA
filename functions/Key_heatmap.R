######################### heatmap plot

key_heatmap <- function(object_ = reduced_data, 
                        key_file='test.xlsx',
                        group_level=c('WT','KP3'),
                        slot='scale.data',
                        row_cluster=F,
                        col_cluster=F,
                        aggregate='Cell'){  ###### aggregate = cell or group
  
  
  DefaultAssay(object_) <- 'RNA'
  
  #convert cluster number to numeric
  cluster_annotaion <- readxl::read_excel(key_file,sheet = 1)
  
  if(cluster_annotaion$Cluster[1]=='All'){
    cluster_annotaion$Cluster = list(unique(object_$seurat_clusters))
  }else{
    cluster_annotaion$Cluster <- lapply(cluster_annotaion$Cluster,function(x) if(length(grep(',',x))>0 ){unlist(strsplit(x,','))%>% as.numeric()}else{x} ) 
  }
  
  #import gene info for heatmap
  genes_for_plot <- readxl::read_excel(key_file,sheet = 2)
  genes_for_plot <- genes_for_plot[which(genes_for_plot$Genes%in%rownames(object_)),]
  genes_for_plot <- genes_for_plot[order(genes_for_plot$Order),]
  
  object_ <- object_ [genes_for_plot$Genes,]
  #### within cluster reordering by group/condition/treatment   
  if(length(group_level) == 0){
    group_level <- unique( object_$group)
  }
  
  extrac_matrix <- function(cluster_info=cluster_annotaion[1,]){
    matrix <- object_[,object_$seurat_clusters%in%(cluster_info$Cluster%>%unlist())]
    matrix$cell_name <- lapply(matrix$seurat_clusters,function(x) cluster_info$`Cell annotation`[which(grepl(x,cluster_info$Cluster))] ) %>%unlist()
    Project(matrix) <- cluster_info$`Cell annotation`
    return(matrix)
  }
  matrix_list <- lapply(1:nrow(cluster_annotaion),function(x) extrac_matrix(cluster_annotaion[x,]))
  
  

  
  
  #####################
  


  group_ordered=c()
  cell_ordered = c()
  plot_matrix_ordered=c()
  cell_name = c()
  groups_and_annotation = c()
  for (mt in matrix_list){
  
      #mt <- ScaleData(mt,features = unlist(genes_for_plot$Genes) )
      plot_matrix <- as.matrix(GetAssayData(object = mt, slot = slot))
      groups <- mt$group
      groups_and_annotation <-  c(groups_and_annotation,paste0(mt$cell_name,' ',mt$group)[order(factor(groups,levels=group_level))])
      group_ordered <-  c(group_ordered,groups[order(factor(groups,levels=group_level))])
      cell_name <- c(cell_name,rep(Project(mt),ncol(plot_matrix)))
      
     if(length(plot_matrix_ordered)==0){
      plot_matrix_ordered <- plot_matrix[,order(factor(groups,levels=group_level))]
      }else{
        plot_matrix_ordered <-cbind(plot_matrix_ordered,plot_matrix[,order(factor(groups,levels=group_level))])
  }}

  
    ha <- HeatmapAnnotation(
      `Cell Type` = cell_name,
      #empty = anno_empty(border = FALSE, height=unit(0.5, "mm")),
      Group =  group_ordered ,
     col = list(Group = c(group = "blue", "KP3" = "red")),
      empty = anno_empty(border = FALSE, height=unit(1.5, "mm"))
    )
    
    # Combine the heatmap and the annotation

    ########## row -cluster
    if(row_cluster==T){
      row_split = NULL
      ifrowcluster='_row_Clustered'
    }else{
      row_split = row_split = genes_for_plot$Order[order(factor(genes_for_plot$Genes,levels=rownames(plot_matrix_ordered)))]
      ifrowcluster='_row_Manually_ordered'
    }
    
    ########## col -cluster
    if(col_cluster==T){
      ifcolcluster='_col_Clustered'
      column_split=NULL
    }else{
      if(aggregate=='Cell'){
        column_split = cell_name
      }else{
        column_split = group_ordered
      }
        
      ifcolcluster='_col_Manually_ordered'
    }
#colnames(plot_matrix_ordered)
    #unique(colnames(plot_matrix_ordered))
    
   p <- Heatmap(plot_matrix_ordered, name = " Scaled Normalized Expression", 
                 row_names_gp = gpar(fontsize = 5),
                 column_names_gp =gpar(fontsize = 6) ,
                 #cluster_rows = color_branches(row_dend, k = 5),
                 #cluster_columns = color_branches(col_dend, k = 2),
                 col = circlize::colorRamp2(c(-2,0,2),c("#4575B4","#FFFFBF","#D73027")),         #brewer.pal(10, "RdYlBu")
                 row_split = row_split,
                 column_split = column_split,
                 top_annotation = ha,show_column_names = F,cluster_rows =row_cluster,
                 cluster_columns = col_cluster,heatmap_legend_param = list(direction = "horizontal"))
  
  save_dir <- dir_create('reduced_PLOTs','All','')
  subdir_create(paste0(save_dir,'/heatmap') )
  pdf(paste0(save_dir,'/heatmap/',aggregate,slot,'_heatmap_',Project(object_),ifrowcluster,ifcolcluster,key_file,'.pdf'),width = max(7,0.002*ncol(plot_matrix_ordered)) )
  draw(p,auto_adjust = FALSE)
  ####### extract expression matrix and meta data
  dev.off()
  
  
  return('Heatmap complete')
  #show_col(hue_pal(h = c(270, 360))(9))
  
}

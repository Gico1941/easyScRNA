######################### heatmap plot

key_heatmap <- function(object_ = reduced_data, 
                        key_file='DEG_heatmap_key.xlsx',
                        group_level=c('WT','KP3'),
                        cell_level=c(),
                        slot='scale.data',
                        row_cluster=F,
                        col_cluster=F){
  
  
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
  genes_for_plot <- genes_for_plot[order(genes_for_plot$Order),]
  
  #### filter plot objects and add cell type annotation
  object_ <- object_[unlist(genes_for_plot$Genes),object_$seurat_clusters %in% unlist(cluster_annotaion$Cluster) ]
  object_$cell_name <- lapply(object_$seurat_clusters,function(x) cluster_annotaion$`Cell annotation`[which(grepl(x,cluster_annotaion$Cluster))[1]] ) %>%unlist()
  
  
  
  #### within cluster reordering by group/condition/treatment   
  if(length(group_level) == 0){
    group_level <- unique( object_$group)
  }
  if(length(cell_level) == 0 ){
    cell_level <- unique( object_$cell_name)
  }
  
  
  ####### extract expression matrix and meta data
  object_ <- ScaleData(object_,features = unlist(genes_for_plot$Genes) )
  plot_matrix <- as.matrix(GetAssayData(object = object_, slot = slot))
  
  groups <-  object_$group
  cell_names <- object_$cell_name
  
  group_ordered <-  groups[order(factor(cell_names,levels = cell_level),factor(groups,levels=group_level))]
  cell_name_ordered <- cell_names[order(factor(cell_names,levels = cell_level),factor(groups,levels=group_level))]
  plot_matrix_ordered <- plot_matrix[genes_for_plot$Genes,order(factor(cell_names,levels = cell_level),factor(groups,levels=group_level))]
  
  ha <- HeatmapAnnotation(
    `Cell annotation`= cell_name_ordered,
    #empty = anno_empty(border = FALSE, height=unit(0.5, "mm")),
    Group =   group_ordered,
    empty = anno_empty(border = FALSE, height=unit(1.5, "mm"))
  )
  
  # Combine the heatmap and the annotation
  
  
  save_dir <- dir_create('reduced_PLOTs','All','')
  subdir_create(paste0(save_dir,'/heatmap') )
  ########## row -cluster
  if(row_cluster==T){
    row_split = NULL
    ifrowcluster='_row_Clustered'
  }else{
    row_split = genes_for_plot$Order
    ifrowcluster='_row_Manually_ordered'
  }
  
  ########## col -cluster
  if(col_cluster==T){
    ifcolcluster='_col_Clustered'
    column_split=NULL
  }else{
    if(length(unique(cell_name_ordered))>1){
      column_split = cell_name_ordered
    }else{
      column_split = group_ordered
    }
    ifcolcluster='_col_Manually_ordered'
  }
  
  pdf(paste0(save_dir,'/heatmap/',slot,'_heatmap_',Project(object_),ifrowcluster,ifcolcluster,key_file,'.pdf') )
  
  draw(Heatmap(plot_matrix_ordered, name = "Scaled Normalized Expression", 
               row_names_gp = gpar(fontsize = 5),
               column_names_gp =gpar(fontsize = 6) ,
               #cluster_rows = color_branches(row_dend, k = 5),
               #cluster_columns = color_branches(col_dend, k = 2),
               col = circlize::colorRamp2(c(-2,0,2),c("#4575B4","#FFFFBF","#D73027")),         #brewer.pal(10, "RdYlBu")
               row_split = row_split,
               column_split = column_split,
               top_annotation = ha,show_column_names = F,cluster_rows =row_cluster,
               cluster_columns = col_cluster,heatmap_legend_param = list(direction = "horizontal")
  ),heatmap_legend_side = "bottom")
  
  dev.off()
  
  
  return('Heatmap complete')
  #show_col(hue_pal(h = c(270, 360))(9))
  
}
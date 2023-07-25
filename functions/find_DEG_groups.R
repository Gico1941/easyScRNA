

####################################### find_DEG_bewteen_groups (in specific clusters)

find_DEG_bewteen_groups <- function(object_reduced,subset_name='All_cluster',subset_cluster=c(),control_group=1,variable_group=2,logfc.threshold =DEG.logfc.threshold,min.pct=DEG.min.pct,save_folder='DEG'){
  object_reduced@active.assay = 'RNA'  
  if(length(subset_cluster)!=0 ){
    if(subset_cluster[1] > 0){
      print(paste0('Looking for DEG in clusters : ',paste0(subset_cluster,collapse = ',') ))
      object_reduced <- subset(reduced_data,idents = subset_cluster,invert = FALSE)
    }
    if( subset_cluster[1] < 0){
      clusters= unique(levels(reduced_data$seurat_clusters))
      print(paste0('Looking for DEG in clusters : ',paste0(clusters[clusters %in%abs(subset_cluster)==F],collapse = ',')))
      object_reduced <- subset(reduced_data,idents = unique(reduced_data$seurat_clusters)[unique(reduced_data$seurat_clusters)%in%abs(subset_cluster)==F],invert = FALSE)
    }}else{
      print('Looking for DEG across all clusters ')
    }
  
  dir <- dir_create(save_folder,subset_name,'')
  if(is.numeric(control_group[1])==T){
    control_group=unique(object_reduced$group)[control_group]
  }
  
  if(is.numeric(variable_group[1])==T){
    variable_group=unique(object_reduced$group)[variable_group]
  }
  
  print(paste0('set : ', control_group,' as control') )
  print(paste0('set : ', variable_group,' as variable') )
  
  DEG <- FindMarkers(object = object_reduced ,ident.1=variable_group,ident.2=control_group,group.by='group',logfc.threshold =logfc.threshold ,min.pct=min.pct)
  
  write.csv(DEG,paste0(dir,'/',subset_name,'_group_',paste0(variable_group,collapse='_'),'_vs_group_',paste0(control_group,collapse='_'),'.csv'))
  #return(DEG)
}

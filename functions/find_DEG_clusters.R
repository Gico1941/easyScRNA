

################################################## find_DEG_bewteen_clusters (in specific groups)

find_DEG_bewteen_clusters <- function(object_reduced=reduced_data,subset_name='All_group',subset_group=c(),control_cluster=1,variable_cluster=2,logfc.threshold =DEG.logfc.threshold,min.pct=DEG.min.pct,save_folder='DEG'){
  object_reduced@active.assay = 'RNA'  
  if(length(subset_group)!=0 ){
    print(paste0('Looking for DEG in groups : ',paste0(subset_group,collapse = ',') ))
    object_reduced <- object_reduced[,object_reduced@meta.data[["group"]] %in% subset_group]
  }else{
    print(paste0('Looking for DEG across all groups '))
  }
  
  all_cluster <- object_reduced$seurat_clusters%>%unique()
  all_cluster <- all_cluster[order(all_cluster)]
  
  variable_cluster_name <- variable_cluster
  control_cluster_name <- control_cluster
  
  if(control_cluster%>%length() == 0 & variable_cluster%>%length() !=0){
    control_cluster= all_cluster[all_cluster%in%variable_cluster==F]
    control_cluster_name <-'others'
  }
  if(control_cluster%>%length() !=0 & variable_cluster%>%length() == 0){
    variable_cluster=all_cluster[all_cluster%in%control_cluster==F]
    variable_cluster_name <- 'others'
  }
  
  
  dir <- dir_create(save_folder,subset_name,'')
  
  print(paste0('set cluster : ', paste0(control_cluster,collapse=',') ,' as control') )
  print(paste0('set cluster : ', paste0(variable_cluster,collapse=','),' as variable') )
  
  DEG <- FindMarkers(object = object_reduced ,ident.1=variable_cluster,ident.2=control_cluster,group.by='seurat_clusters',logfc.threshold =logfc.threshold ,min.pct=min.pct)
  
  write.csv(DEG,paste0(dir,'/',subset_name,'_cluster_',paste0(variable_cluster_name,collapse='_')  ,'_vs_cluster_',paste0(control_cluster_name ,collapse='_') ,'.csv'))
  #return(DEG)
}
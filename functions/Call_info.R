####################################### call info
call_info <- function(data){
  cat(paste0('\n\n-----------------------------------useful info for find DEG--------------------------------- :\n' ,
             'groups in dataset :',paste0(unique(data$group),collapse = ','),'\n','clusters in dataset :',paste0(levels(data$seurat_clusters ),collapse = ' ')))
}

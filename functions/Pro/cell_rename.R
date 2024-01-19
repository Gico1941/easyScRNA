
######## rename cell extension with sample info for inter-group integration

cell_extesion_rename <- function(object=object_data){
  pro_rename <- function(index){
    ##### rename cells for further integration
    
    project <- object$object[[index]]
    project <- RenameCells(project,new.names = unlist(lapply( rownames(project@meta.data), function(x) paste0(strsplit(x,'-')[[1]][1],'-',index) ))) # or index here can be replace with : object$folder[,level][index]
    
    
    ##### add group info to meta data
    
    cat(paste0('\nAdding group info for :', object$folder$Group[index] ,'...'))
    
    project <- AddMetaData(project,rep(object$folder$Group[index],length(project@meta.data$orig.ident)),col.name = 'group')
    
    return(project)
    
  }
  ####
  cat('\nRenaming cells....')
  ### replace extension with sample name
  object$object <- sapply(object$folder$Sample ,function(x) pro_rename(which(object$folder$Sample ==x)),USE.NAMES =T )
  cat('\nrenaming complete....')
  return(object)
}
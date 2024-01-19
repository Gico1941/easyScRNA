
####### file_extraction

file_extraction <- function(parent='Raw'){
  cat('reading files...')
  dirs <- list.dirs(parent,recursive=T)
  
  dirs_matrix=data.frame(dir=dirs,level=lapply(dirs,function(x) length(unlist(gregexec('/',x)))) %>% unlist())
  
  folder_dirs <- dirs_matrix[which(dirs_matrix$level==max(dirs_matrix$level)),]
  folder_dirs[,'Group'] <- unlist(lapply(folder_dirs$dir,function(x) strsplit(x,'/')[[1]][2] ))
  folder_dirs[,'Sample'] <- unlist(lapply(folder_dirs$dir,function(x) strsplit(x,'/')[[1]][3] ))
  
  return(folder_dirs)
  
}

########## import data and create seurat object

Object_import <- function(fd,min.features=200,min.cells=3){
  cat('\nCreating Seurat object...')
  
  object_list<- lapply(c(1:nrow(fd)),function(x) CreateSeuratObject(counts = Read10X(fd$dir[x]), 
                                                                    project = fd$Sample[x], 
                                                                    min.cells = min.cells, 
                                                                    min.features = min.features))
  
  
  return(list(object=object_list,folder=fd))
}

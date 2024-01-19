############## dir_create  level = Group or Sample
dir_create <- function(key,folder,level){
  if(level!=''){
    save_dir <- paste0(key,'/',unique(folder[,level]))
    lapply(save_dir,function(x) if(dir.exists(x)!=T){
      dir.create(x,recursive = T)
      
    })
  }else{
    save_dir <- paste0(key,'/',folder)
    if(dir.exists(save_dir)!=T){
      dir.create(save_dir,recursive = T)
    } 
  }
  
  return(save_dir)}
################## subdir for plot folder organization

subdir_create <- function(subdir){
  if(dir.exists(subdir)==F){
    dir.create(subdir,recursive = T)
  }
  return('Creat subdir for plot saving...')}

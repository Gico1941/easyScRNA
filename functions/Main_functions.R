############## main funcions
standard_processing <- function(read=F,save=T){
  if(read!=T){
    #### read folders
    fd <- file_extraction(raw_file_folder)
    
    ##### import data
    object_data <-  Object_import(fd,min.features=min.features,min.cells=min.cells)
    
    ##### rename cell extension with group name   ## will barcode of one group do not overlap across all samples?
    object_data <- cell_extesion_rename(object_data)
    
    ##### replace object with QC filtering
    object_data$object <- SimpleQC(object_data$object,object_data$folder)
    
    ##### normalization and variable determination
    object_data$object <- normalization_and_variable_finding(object_data$object)
    
    ##### group integration  ##### also required for integration-free dataset
    object_data$object <- group_integration(object_data)

    #### PCA and plot
    object_data$object <- reduction_plot(object_data)
    
    ##### batch_plot heatmap
    multi_plot(object_data)%>%invisible()
    

    if(save==T){
      if(dir.exists('RDS')!=T){
        dir.create(('RDS'),recursive = T)
      }
      print('Saving current work to .rds file ...')
      saveRDS(object_data,'RDS/Raw.rds')
    }
  }else{
    object_data<-readRDS('RDS/Raw.rds')
  }
  return(object_data)
}


processing_reduced_data<-function(object_data,read=F,save=T,idents_ = idents_you_want_to_keep){
  
  if(read==F){
    ######## exclude undesired clusers or skip this line
    object_data$reduced$All  <- subset(x=object_data$object$All , idents = idents_) 
    
    ####### re-find variable
    print('Set default assay as : RNA')
    DefaultAssay(object_data$reduced$All) <- 'RNA'
    object_data$reduced <- normalization_and_variable_finding(object_data$reduced)
    
    ###### re-clustering 
    object_data$reduced <- reduction_plot(object_data,reduced=T)
    ###### plot
    multi_plot(object_data,reduced=T,extra_features = extra_features)
    
    reduced_data <- object_data$reduced$All
    

    if(save==T){
      if(dir.exists('RDS')!=T){
        dir.create(('RDS'),recursive = T)
      }
      print('Saving current work to .rds file...')
      saveRDS(reduced_data,'RDS/reduced_all.rds')
    }
  }else{
    reduced_data <- readRDS('RDS/reduced_all.rds')
  }
  
  return(reduced_data)}

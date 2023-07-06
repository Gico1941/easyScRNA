
########### data_normalization

normalization_and_variable_finding <- function(object = object_data$object){
  
  nm_vf_batch <- function(object){
    cat(paste0('Processing batch :',Project(object),'.....\n'))
    object_normalized <- NormalizeData(object)
    cat('find_variables...\n')
    object_normalized <- FindVariableFeatures(object_normalized, selection.method = selection_method)
    #length(x = VariableFeatures(object = object))
    return(object_normalized)
  }
  cat('\nLaunch normalization........\n')
  return(lapply(object,nm_vf_batch))
  
}

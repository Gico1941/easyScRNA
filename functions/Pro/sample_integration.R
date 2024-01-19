

########  sample integration by condition/group/treatment

group_integration <- function(object_data){
  cat('\nIntegrating data...')
  
  ####### rename list name of all object
  
  grouped_object <- lapply(unique(object_data$folder$Group),function(x) object_data$object[grep(x,object_data$folder$Group)] )
  
  #### determine group that will have multiple samples and will be integrated
  #multi_sample_group_folder <- object_data$folder[which(unlist(lapply(object_data$folder$Group,function(x) grep(x,object_data$folder$Group) %>% length() > 1 ))) ,]
  #integration_plot_save_dir_list <-  dir_create('Integrated',multi_sample_group_folder,'Group')
  
  
  integration_batch <- function(list_of_object_for_integration,rds_save_dir=F,group_name='All'){
    cat(paste0('\nIntegrating following data: ',paste0(unlist(lapply(list_of_object_for_integration,Project)),collapse = ' *and* '),'\n'))
    features <- SelectIntegrationFeatures(object.list = list_of_object_for_integration)
    integrates.all.anchors <- FindIntegrationAnchors(object.list = list_of_object_for_integration, dims = 1:100,anchor.features = features)
    integrates.all <- IntegrateData(anchorset = integrates.all.anchors, dims = 1:100)
    
    #if(file.exists(paste0(rds_save_dir,'/',group_name,'_integrated.rds'))==F){
    #  saveRDS(integrates.all, file = paste0(rds_save_dir,'/',group_name,'_integrated.rds'))
    #}
    
    Project(integrates.all) <- group_name
    return(integrates.all)
  }
  
  integration_condition <- function(x){
    if(grep(x,object_data$folder$Group) %>% length() > 1)
    {
      integration_batch(grouped_object[[which(unique(object_data$folder$Group)==x)]],
                        integration_plot_save_dir_list[which(unlist(lapply(integration_plot_save_dir_list,basename))==x)],
                        x)
    }
    else{
      Project(grouped_object[[which(unique(object_data$folder$Group)==x)]][[1]] ) <- x
      
      
      return(grouped_object[[which(unique(object_data$folder$Group)==x)]][[1]] )
    } 
  }
  
  ##### all groups merge into one 
  #all.save.dir <- dir_create('Integrated','All','') 
  all.save.dir <- F
  object.all <- integration_batch(object_data$object ,all.save.dir,'All')
  
  #### grouped_data :
  #grouped objects: are grouped based on their group identity/condition/treatment 
  #individual:origin sample named objects
  #all:all conditions are integrated together
  
  
  return(append(list(All=object.all),sapply(unique(object_data$folder$Group),integration_condition,USE.NAMES=T) )) }  

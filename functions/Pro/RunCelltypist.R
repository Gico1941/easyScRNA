################### list all avaliable models
Model_list <- function(model_path){
  available_models <- lapply(list.files(model_path,pattern = '.pkl'),function(x) gsub('.pkl','',x))
  cat(paste0('\n-----------Installed Celltypist models : --------\n\n ', paste0(available_models,collapse = ' \n ' )))
}

############################## download models for cell annotating
Model_download<-function(model_path = 'CelltypistModel'){
  print('Download models :')
  subdir_create(model_path) %>% invisible()
  celltypist$models$models_path <- model_path
  celltypist$models$download_if_required()
}

#################### model_import
Model_import <- function(models){
  models <- paste0(models,'.pkl')
  model_imported<- sapply(models,function(x) celltypist$models$Model$load(model=x) )
  return(model_imported)
}

Runcelltypist <- function(dt = seura_obj,model=c(),majority_voting = F,test_model=F,assay=NULL){
  if(length(model) ==0){
    return()
  }
  md_list <- Model_import(model)
  
  ########## tests 
  if(test_model==T){
    pbmc.data <- Read10X('test_dataset/filtered_gene_bc_matrices/hg19/')
    pbmc <- CreateSeuratObject(counts = pbmc.data, project = "pbmc3k", min.cells = 3, min.features = 200)
    pbmc <- NormalizeData(pbmc)
    pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
    pbmc <- ScaleData(pbmc)
    pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
    pbmc <- FindNeighbors(pbmc, dims = 1:10)
    pbmc <- FindClusters(pbmc, resolution = 0.5)
    pbmc <- RunUMAP(pbmc, dims = 1:10)
    #DimPlot(pbmc, reduction = "umap")
    dt<-pbmc
  }
  
  if(length(assay) == 0){
    assay <- DefaultAssay(dt)
  }
  
  ############
  dt.data = scanpy$AnnData(X = numpy$array(as.matrix(t(as.matrix(GetAssayData(dt,slot = 'data'))))),
                         obs = pandas$DataFrame(dt@meta.data),
                         var = pandas$DataFrame(data.frame(gene = rownames(GetAssayData(dt,slot = 'data')),
                                                           row.names = rownames(GetAssayData(dt,slot = 'data'))))
  )
  # Add cluster information to adata
  dt.data$obs$seurat_clusters <- dt$seurat_clusters
  
  predictions = celltypist$annotate(dt.data, model = md_list[[1]], majority_voting = majority_voting,over_clustering='seurat_clusters')  # set majority_voting for simplified_annotation
  if(majority_voting==T){
    dt  = AddMetaData(dt, predictions$predicted_labels['majority_voting'], col.name ="typist_prediction") 
  }else{
    dt  = AddMetaData(dt, predictions$predicted_labels['predicted_labels'], col.name ="typist_prediction") 
  }
  
  
  
  #SetIdent(dt,value = "typist_prediction")
  #DimPlot(dt,group.by='typist_prediction')
  return(dt)
}

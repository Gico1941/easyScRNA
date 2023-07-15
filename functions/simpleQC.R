
######### simple QC and data_filtering

SimpleQC <- function(object_list=object_data$object,folder=object_data$folder){
  
  QC_plot_save_dir_list <- dir_create('QC',folder,'Sample')
  
  QC_batch <- function(object,plot_save_dir){
    cat(paste0('\nPerforming QC on ',Project(object)))
    pdf(paste0(plot_save_dir,'/',Project(object),'.pdf'))
    mito.features <- grep(pattern = "^mt-", x = rownames(x = object), value = TRUE)
    percent.mito <- Matrix::colSums(x = GetAssayData(object = object, slot = 'counts')[mito.features, ]) / Matrix::colSums(x = GetAssayData(object = object, slot = 'counts'))
    object[['percent.mito']] <- percent.mito
    object.QC <- subset(x=object, subset = nFeature_RNA > nFeature_RNA_threshold_min & nFeature_RNA < nFeature_RNA_threshold_max & percent.mito < percent.mito_threshold)
    print(VlnPlot(object = object.QC, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3, pt.size = 0.1 ))
    dev.off()
    return(object.QC)
  }
  cat('\nLaunch QC........')
  return(lapply(1:length(object_list), function(x) QC_batch(object_list[[x]],QC_plot_save_dir_list[[x]])))
  cat('QC finished')
}
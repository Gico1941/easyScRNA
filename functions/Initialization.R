remotes <- c('Seurat')

CRAN <- c('remotes','tidyverse','devtools','cowplot','ggrepel','ggtext','cluster','scales','circlize','RColorBrewer','dendextend','BiocManager','Seurat','installr')

devtools <- c("thomasp85/patchwork","jokergoo/ComplexHeatmap")

biocManager <- c('enrichplot')


CRAN_install <- function(x){if(require(x,character.only = T)==F){
  cat(paste0('Installing package : ',x,'\n'))
  install.packages(x,character.only = T)
  require(x,character.only = T)
}}

invisible(lapply(CRAN,CRAN_install))


library(devtools)
devtools_install <-  function(x){if(require(basename(x),character.only = T)==F){
  cat(paste0('Installing package : ',x,'\n'))
  devtools::install_github(x,character.only = T)
  require(basename(x),character.only = T)
}}

invisible(lapply(devtools,devtools_install))


biocManager_install <- function(x){if(require(x,character.only = T)==F){
  cat(paste0('Installing package : ',x,'\n'))
  BiocManager::install(x,character.only = T)
  require(x,character.only = T)
}}

invisible(lapply(biocManager,biocManager_install))


##### load all R functions

source_all <- function(function_path='functions'){
  utilities <- list.files(function_path,pattern = '.R',full.names = T,recursive = F)
  utilities <- utilities[-grep('Initialization.R',utilities)]
  invisible(lapply(utilities,source) )
}

source_all() 


######### global variables may be ueful :
###### import parameters
min.cells = 3
min.features = 200

#########  QC and filtering parameters
nFeature_RNA_threshold_max = 20000
percent.mito_threshold = 0.1
nFeature_RNA_threshold_min = 200

####### normalization parameters
selection_method = "vst"
n_of_variant_feature = 2000

#### heatmap dimension
dimension_plotting_intervals = 10   # n of plots in same pdf
n_cells_for_plot = 200               # n of cols in heatmap
ElbowPlot_ndims = 50

#### UMAP resolution 
UMAP_resolution <- 0.6
UMAP_max_dims <- 50

### DEG parameter
DEG.logfc.threshold = 0.1
DEG.min.pct=0.1

### rnk convert parameter
rnk.padj_hold = 0.1
rnk.log2fc_hold =0
#### volcano plot parameters

vol_plot.padj_hold= 0.01
vol_plot.log2fc_hold = 0.3


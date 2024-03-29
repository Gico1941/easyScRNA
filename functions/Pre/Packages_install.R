#remotes <- c('Seurat')
remotes <- c('')

CRAN <- c("downloader",'remotes','tidyverse','devtools','cowplot','ggrepel','ggtext','cluster','scales','circlize','RColorBrewer','dendextend','BiocManager','Seurat','installr','reticulate')

devtools <- c("thomasp85/patchwork","jokergoo/ComplexHeatmap")

biocManager <- c('enrichplot')
############### install R packages

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


remote_install <- function(x){if(require(x,character.only = T)==F){
  cat(paste0('Installing package : ',x,'\n'))
  remotes::install_github(x,character.only = T)
  require(x,character.only = T)
}}

if (remotes != c('')){
  invisible(lapply(remotes,remote_install))
}


Seurat_install <- function(x){if(require(x,character.only = T)==F){
  cat(paste0('Installing package : ',x,'\n'))
  # Install the remotes package 
  require('remotes')
  # Replace Seurat with '4.3.0' 
  remotes::install_version(package = 'Seurat', version = package_version('4.3.0'))
  
  require(x,character.only = T)
}}

Seurat_install('Seurat')


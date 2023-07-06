
GSEA_download <- function(url='https://data.broadinstitute.org/gsea-msigdb/gsea/software/desktop/4.3/GSEA_4.3.2.zip'){
  if(dir.exists(GSEA_path)!=T ){
    print(paste0('Downloading ',basename(url)) )
    downloader::download(url,basename(url))
    unzip(basename(url))
    file.remove(basename(url))
  }    
  return('GSEA downloaded')

}



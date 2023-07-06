library(installr)

version=11

java_install <- function(){
  if(dir.exists(paste0(GSEA_path,'/jdk-',version) )==F){
    installr::install.java(version = version,
                           path=GSEA_path)
    print('Insalling Java')
  }
  return('Java already installed...')
}


 


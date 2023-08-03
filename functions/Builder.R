
easyBuild <- function(force=F){
  
  conda_install_ <- function(){
    if(system('conda -v') ==F | force == T )
      installr::install.conda()
  }
  #############download microsoft visual studio for celltypist installation
  celltypist_initialization <- function(){
    if(length(grep('vs_BuildTools.exe',list.files())) ==0 | force == T){
      print('Download visual builder')
      download.file('https://aka.ms/vs/17/release/vs_BuildTools.exe','vs_BuildTools.exe', mode = "wb")
      wd <- gsub('/','\\\\',getwd())
      system('cmd.exe',input =  paste0(wd,'\\vs_buildtools.exe --norestart --passive --downloadThenInstall --includeRecommended --add Microsoft.VisualStudio.Workload.NativeDesktop --add Microsoft.VisualStudio.Workload.VCTools --add Microsoft.VisualStudio.Workload.MSBuildTools'))
      }else{
      print('Visual builder installed...')
      }}

  ############ create conda environment and install celltypist
  Celltypist_install <- function(packs){
    if(force==T){
      conda_remove('celltypist')
    }
    if('celltypist'%in%conda_list()$name==F | force==T){
      print('Install celltypist...')
      conda_create('celltypist')
    }
      installed_packs <- py_list_packages()$package
      for( i in packs){
        if(i %in% installed_packs ==F){
          conda_install('celltypist',packages = i,channel=c('anaconda','conda-forge','r','bioconda'))
        }}
    
      print('CellTypist installed...')
    }
  
  ########### install java and GSEA command line
  
  GSEA_download <- function(url='https://data.broadinstitute.org/gsea-msigdb/gsea/software/desktop/4.3/GSEA_4.3.2.zip',
                            GSEA_path='GSEA_4.3.2'){
    if(dir.exists(GSEA_path)!=T | force == T  ){
      print(paste0('Downloading GSEA: ',basename(url)) )
      downloader::download(url,basename(url))
      unzip(basename(url))
      file.remove(basename(url))
    }    
    print('GSEA installed...')
  }
 
  java_install <- function(GSEA_path='GSEA_4.3.2',version=11){
    if(dir.exists(paste0(GSEA_path,'/jdk-',version) )==F | force==T){
      print('Insalling Java')
      installr::install.java(version = version,
                             path=GSEA_path)
    }
    print('Java installed...')
  }
  conda_install_()
  celltypist_initialization()
  Celltypist_install(packs = c('celltypist','scanpy','pandas','numpy','leidenalg'))
  #system('conda install -c bioconda -c conda-forge celltypist')
  GSEA_download()
  java_install()
  }



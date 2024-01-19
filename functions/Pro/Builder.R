
easyBuild <- function(force=F){
  
  conda_install_ <- function(){
    if(file.exists(reticulate::miniconda_path()) == F | force == T ){
      reticulate::install_miniconda()
    }
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
    if(force==T & 'celltypist' %in% conda_list()$name == T){
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
  
  
  
  install.java_new <- function (version = 11, page_with_download_url = "https://jdk.java.net/java-se-ri/", 
                            path = "C:/java") 
  {
    
    URL <- 'https://download.java.net/openjdk/jdk11/ri/openjdk-11+28_windows-x64_bin.zip'
    filename <- file.path(tempdir(), file.name.from.url(URL))
    download.file(URL, destfile = filename, mode = "wb")
    if (grepl("zip$", URL)) {
      unzip(zipfile = filename, exdir = path)
    }
    if (grepl("gz$", URL)) {
      untar(tarfile = filename, exdir = path)
    }
    path_list <- list.dirs(path)
    home_path <- grep("jdk-[0-9]+$", path_list, value = T)
    home_path <- grep(version, home_path, value = T)
    profiled <- paste0("Sys.setenv(JAVA_HOME='", home_path, "')")
    if (!file.exists("~/.Rprofile")) {
      file.create("~/.Rprofile")
    }
    pre <- readLines("~/.Rprofile", warn = F)
    if (any(grepl("JAVA_HOME", pre))) {
      pre <- pre[-grep("JAVA_HOME", pre)]
    }
    profiled <- c(pre, profiled)
    write(profiled, "~/.Rprofile")
    Sys.setenv(JAVA_HOME = home_path)
  }
  
  java_install <- function(GSEA_path='GSEA_4.3.2',version=11){
    
    if(dir.exists(paste0(GSEA_path,'/jdk-',version) )==F | force==T){
      print('Insalling Java')
      install.java_new(version = version,
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












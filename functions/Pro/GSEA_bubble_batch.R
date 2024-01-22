

#############################
GSEA_bubble <- function(GSEA_folder='GSEA',GSEA_fdr_hold=0.5,fdr_top=20){
  files <- list.files(GSEA_folder,pattern='gsea_report_for',recursive = T,full.names = T)
  files <- data.frame(file_name=grep('.tsv',files,value = T))
  files$ident <- lapply(files$file_name,function(x) str_sub(x,-17,-1)) %>% unlist() %>% invisible()
  save_dirs <- files$file_name %>% dirname() %>% dirname()
  
  batch <- function(iden=files$ident[3],fl = files$file_name[3]){
    print(paste0('Launch GSEA plot for ' ,dirname(files$file_name[files$ident==iden]) %>% basename() %>% unique()) )
    data <- lapply(files$file_name[files$ident==iden],function(x) read.table(x,sep="\t",quot="",header = T)%>%invisible())
    
    plot <- function(dt=data[[1]]){
      if(nrow(dt)==0){
        return('Ship empty data')
      }
      if(dt$NES[1]<0){
        group='Suppressed'
      }else{
        group='Activated'
      }
      
      dt <-  dt[dt$`FDR.q.val`<GSEA_fdr_hold,]
      dt <- dt[order(dt$`FDR.q.val`),]
      if(nrow(dt) > fdr_top){
        dt <- dt[1:20,]
      }
      
      dt$NAME <- gsub('_',' ',dt$NAME)
      dt$NAME <- str_wrap(dt$NAME, width = 40,  indent = 2,whitespace_only = T)
      dt$NES <- abs(dt$NES)
      dt$NAME <- factor(dt$NAME,levels = dt$NAME[order(dt$NES)])

          ggplot(dt,aes(x=NES,y=NAME,size=SIZE,color=`FDR.q.val`))+
            geom_point()+
            theme_bw()+
            scale_color_continuous(high='red',low='blue')+
            ylab(NULL)+
            xlab('Absolute NES')+
            ggtitle(paste0(basename(fl),'_GSEA_enrichment_',group))+
            theme(axis.text.y=element_text(size=10))
 
      ggsave(paste0(strsplit(dirname(fl),'.Gsea')[[1]][1],'_',group,'.pdf'),height = max(7,0.5*nrow(dt)))
    }
    lapply(1:length(data),function(x) plot(data[[x]]) )
    
  }
  lapply(1:nrow(files),function(x) try(batch(files$ident[x],files$file_name[x])))
  return('GSEA plot finished ....')
}



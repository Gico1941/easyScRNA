
###################### volcano plot for all DEG csv

volcano <- function(DEG_path='DEG',p_hold = vol_plot.padj_hold,log2_fc_hold = vol_plot.log2fc_hold,tops=10,hightlight_top_by='p_val_adj',hightlight_by_keys=T){
  
  if(hightlight_top_by %in% c('p_val_adj','avg_log2FC')==F){
    return('Plot failed, avaliable hightlight_top_by options : p_val_adj,avg_log2FC')
  }
  
  folders <- list.dirs(DEG_path,recursive = F)
 
  
  
  
  batch <- function(path=folders[1]){
    csv_file <- list.files(path,recursive=T,pattern = '.csv',full.names = T)
    hightlight_keys <- list.files(path,recursive=T,pattern = 'highlight_key.txt',full.names = T)
    
    
    
    data <- read.csv(csv_file,row.names = 1)
    data <- data[data$p_val_adj<p_hold,]
    
    data$significance <- 'Stable'
    data$significance[data$avg_log2FC>log2_fc_hold] <- 'Up'
    data$significance[data$avg_log2FC< -log2_fc_hold] <- 'Down'
    data$gene <- rownames(data)
    
    Colors <- c(Up='red',Stable='grey',Down='blue')
    
    
    if(length(hightlight_keys) != 0 & hightlight_by_keys==T){
      hightlight_keys <- read_tsv(hightlight_keys,col_names = F) %>% unlist()
      data_highlight <- data[hightlight_keys,]
      
      hightlight_by = 'keys'
    }else{
      data_highlight <- rbind(top_n(data[data$significance=='Up',],
                                    min(nrow(data[data$significance=='Up',]),tops),
                                    data[data$significance=='Up',][,hightlight_top_by]),
                              
                              top_n(data[data$significance=='Down',],
                                    max(-nrow(data[data$significance=='Down',]),-tops),
                                    data[data$significance=='Down',][,hightlight_top_by])) 
      hightlight_by = hightlight_top_by
    }
  
    
    plot<-ggplot(data,aes(y=-log10(p_val_adj),x=avg_log2FC,color=significance))+
      theme_classic()+
      geom_point(data=data_highlight,aes(x=avg_log2FC,y=-log10(p_val_adj)),colour='#f3e529',size=1.75,fill='white',alpha=0.5)+
      geom_point(size=1,alpha=0.5)+
      
      scale_color_manual(values=Colors)+
      geom_vline(xintercept =log2_fc_hold,linetype='dotted',alpha=0.5 )+
      geom_vline(xintercept =-log2_fc_hold,linetype='dotted',alpha=0.5)+
      geom_hline(yintercept =-log10(p_hold ),linetype='dotted',alpha=0.5 )+
      geom_text_repel(data=data_highlight,aes(x=avg_log2FC,y=-log10(p_val_adj),color=significance,label = gene),position='jitter',show.legend = FALSE,size=1)
    
    ggsave(paste0(strsplit(csv_file,'.csv')[[1]][1],'_',hightlight_by,'.pdf') ,plot,width = 7,height=7)
    
    return('Volcano plot Complete')}
  
  lapply(folders,batch)
}

###################### volcano plot for all DEG csv

volcano <- function(DEG_path='DEG',p_hold = vol_plot.padj_hold,log2_fc_hold = vol_plot.log2fc_hold,tops=10,highlight_top_by='p_val_adj',highlight_by_keys=T){
  
  if(highlight_top_by %in% c('p_val_adj','avg_log2FC')==F){
    return('Plot failed, avaliable highlight_top_by options : p_val_adj,avg_log2FC')
  }
  
  folders <- list.dirs(DEG_path,recursive = F)
 
  
  
  
  batch <- function(path=folders[2]){
    csv_file <- list.files(path,recursive=T,pattern = '.csv',full.names = T)
    highlight_keys <- list.files(path,recursive=T,pattern = 'highlight_key.txt',full.names = T)
    
    cat(paste0('Plotting: ',basename(csv_file),'\n'))
    
    data <- read.csv(csv_file,row.names = 1)
    data <- data[data$p_val_adj<p_hold,]
    
    data$significance <- 'Stable'
    data$significance[data$avg_log2FC>log2_fc_hold] <- 'Up'
    data$significance[data$avg_log2FC< -log2_fc_hold] <- 'Down'
    data$gene <- rownames(data)
    
    Colors <- c(Up='red',Stable='grey',Down='blue')
    
    
    if(length(highlight_keys) != 0 & highlight_by_keys==T){
      highlight_keys <- read.table(highlight_keys,col.names = F) %>% unlist()
      
      data_highlight <- data[highlight_keys,] %>% na.omit()
      
      if(nrow(data_highlight)==0){
        return('All of highlight_key genes were not found in DEG')
      }
      highlight_by = 'keys'
      
    }else{
      if(length(highlight_keys) == 0 & highlight_by_keys ==T){
        cat('"highlight_key.txt" file not found, use tops instead\n')
      }
      data_highlight <- rbind(top_n(data[data$significance=='Up',],tops,
                                    data[data$significance=='Up',][,highlight_top_by]),
                              
                              top_n(data[data$significance=='Down',],-tops,
                                    data[data$significance=='Down',][,highlight_top_by])) 
      highlight_by = highlight_top_by
    }

    
    plot<-ggplot(data,aes(y=-log10(p_val_adj),x=avg_log2FC,color=significance))+
      theme_classic()+
      
      geom_point(size=1,alpha=0.5)+
      scale_color_manual(values=Colors)+
      geom_vline(xintercept =log2_fc_hold,linetype='dotted',alpha=0.5 )+
      geom_vline(xintercept =-log2_fc_hold,linetype='dotted',alpha=0.5)+
      geom_hline(yintercept =-log10(p_hold ),linetype='dotted',alpha=0.5 )+
      geom_point(data=data_highlight,aes(x=avg_log2FC,y=-log10(p_val_adj)),colour='#f3e529',size=0.4,fill='white',alpha=0.7)+
      geom_text_repel(data=data_highlight,aes(x=avg_log2FC,y=-log10(p_val_adj),color='black',label = gene),position='jitter',show.legend = FALSE,size=1)
      
    
    ggsave(paste0(strsplit(csv_file,'.csv')[[1]][1],'_',highlight_by,'.pdf') ,plot,width = 7,height=7)
    
    
    
    return('Volcano plot Complete\n')}
  
  for (i in folders){cat(paste0(batch(i),'\n'))}
}
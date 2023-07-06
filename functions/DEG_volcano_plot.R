
###################### volcano plot for all DEG csv

volcano <- function(DEG_path='DEG',p_hold = vol_plot.padj_hold,log2_fc_hold = vol_plot.log2fc_hold,tops=10){
  
  files <- c(list.files(DEG_path,recursive=T,pattern = '.csv',full.names = T))
  batch <- function(csv_file=files[1]){
    data <- read.csv(csv_file,row.names = 1)
    data <- data[data$p_val_adj<p_hold,]
    
    data$significance <- 'Stable'
    data$significance[data$avg_log2FC>log2_fc_hold] <- 'Up'
    data$significance[data$avg_log2FC< -log2_fc_hold] <- 'Down'
    data$gene <- rownames(data)
    
    Colors <- c(Up='red',Stable='grey',Down='blue')
    data_top <- rbind( top_n(data[data$significance=='Up',],tops,data[data$significance=='Up',]$avg_log2FC)
                       ,top_n(data[data$significance=='Down',],-tops,data[data$significance=='Down',]$avg_log2FC)) 
    
    
    plot<-ggplot(data,aes(y=-log10(p_val_adj),x=avg_log2FC,color=significance))+
      theme_classic()+
      geom_point(data=data_top,aes(x=avg_log2FC,y=-log10(p_val_adj)),colour='#f3e529',size=1.75,fill='white',alpha=0.5)+
      geom_point(size=1,alpha=0.5)+
      
      scale_color_manual(values=Colors)+
      geom_vline(xintercept =log2_fc_hold,linetype='dotted',alpha=0.5 )+
      geom_vline(xintercept =-log2_fc_hold,linetype='dotted',alpha=0.5)+
      geom_hline(yintercept =-log10(p_hold ),linetype='dotted',alpha=0.5 )+
      geom_text_repel(data=data_top,aes(x=avg_log2FC,y=-log10(p_val_adj),color=significance,label = gene),position='jitter',show.legend = FALSE)
    
    ggsave(paste0(strsplit(csv_file,'.csv')[[1]][1],'.pdf') ,plot,width = 7,height=7)
    
    return('Volcano plot Complete')}
  
  lapply(files,batch)
}
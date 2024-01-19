

###############################################  convert DEG csv to .rnk for further GSEA analysis

DEG2RNK <- function(DEG_path='DEG',p_hold=0.1,log2fc_hold=0,remove_mt_header=F){
  files <- c(list.files(DEG_path,recursive=T,pattern = '.csv',full.names = T))
  batch <- function(csv_file=files[1]){
    data <- read.csv(csv_file,row.names = 1)
    if(remove_mt_header==T){
      rownames(data) <- lapply(rownames(data),function(x) gsub('mt-','',x)) %>% unlist()  #################### 
      
    }
    
    #data <- data[data$p_val_adj<p_hold,]
    up <- data[data$avg_log2FC>log2fc_hold,]
    down <- data[data$avg_log2FC< (-log2fc_hold),]
    all <- rbind(up,down)
    if(nrow(up)!=0){
      write_tsv(cbind(data.frame(rownames(up)),up[,'avg_log2FC']),paste0(strsplit(csv_file,'.csv')[[1]][1],'_up.rnk') ,col_names = F)
    }
    if(nrow(down)!=0){
      write_tsv(cbind(data.frame(rownames(down)),down[,'avg_log2FC']),paste0(strsplit(csv_file,'.csv')[[1]][1],'_down.rnk'),col_names = F)
    }
    write_tsv(cbind(data.frame(rownames(all)),all[,'avg_log2FC']),paste0(strsplit(csv_file,'.csv')[[1]][1],'_all.rnk'),col_names = F)
    return('DEG to RNK conversion Complete')}
  lapply(files,batch)
}
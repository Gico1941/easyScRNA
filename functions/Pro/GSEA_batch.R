


############################################   GSEA batch pipeline # reauires GSEA command line package

GSEA_batch <- function(
    DEG_path='DEG',
    species = 'mouse',
    gene_sets =c(`hallmark gene sets`='mh.all.v2023.1.Mm.symbols.gmt',
                 `positional gene sets`='m1.all.v2023.1.Mm.symbols.gmt',
                 `curated gene sets`='m2.all.v2023.1.Mm.symbols.gmt',
                 `regulatory target gene sets`='m3.all.v2023.1.Mm.symbols.gmt',
                 `ontology gene sets`='m5.all.v2023.1.Mm.symbols.gmt',
                 `cell type signature gene sets`='m8.all.v2023.1.Mm.symbols.gmt'),
    symbol_chip='Mouse_Gene_Symbol_Remapping_MSigDB.v2023.1.Mm.chip',
    out_dir='GSEA',
    GSEA_plots_number=30,
    collapse='Remap_Only'
){
  
  rnks <- list.files(DEG_path,recursive=T,pattern = '_all.rnk',full.names = T,include.dirs=T)
  save_dirs <- lapply(rnks, function(x) dir_create(out_dir,strsplit(x,'/')[[1]][2],''))
  rnks <- paste0(getwd(),'\\\\',rnks)
  
  batch <- function(rnk=rnks[1],gene_set=gene_sets[4],save_dir=save_dirs[1],name){
    
    rnk <- gsub('/','\\\\',rnk,fixed = TRUE)
    
    
    save_dir <-  paste0(getwd(),'\\\\',save_dir)
    save_dir <- gsub('/','\\\\',save_dir,fixed = TRUE)
    
    fileConn <- file("functions/GSEA_batch.bat")
    command <- readLines(fileConn)
    close(fileConn)
    
    #Replace Values
    command <- gsub("_rnk_file_", rnk, command)
    command <- gsub("_chip_", symbol_chip, command)
    command <- gsub("_gene_set_", gene_set, command)
    command <- gsub("_save_dir_", save_dir, command)
    command <- gsub("_the_label_", paste0(basename(save_dir),'--',gsub(' ','_',name)), command,fixed = T)
    command <- gsub("_cd_",GSEA_path, command)
    command <- gsub("_numberplot_",GSEA_plots_number, command)
    command <- gsub("_collapse_",collapse, command)
    command <- gsub("_species_",species, command)
    print(paste0('Performing enrichment analysis of ---  ',basename(rnk),'  --- with ---  ' ,name,' : ',gene_set,' .....' ))
    system("cmd.exe",input=command,show.output.on.console = F)
    Sys.sleep(1)
  }
  for(i in 1:length(rnks)){
    for(name in names(gene_sets) ){
      batch(rnks[i],gene_sets[name],save_dirs[i],name)}}
}

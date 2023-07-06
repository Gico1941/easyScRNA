
source('functions/Initialization.R')

##-------------------------------------------------------------------------- configurations
### GSEA_install folder

GSEA_path = 'GSEA_4.3.2'

###### import parameters
raw_file_folder='Raw'
min.cells = 3
min.features = 200

#########  QC and filtering parameters
nFeature_RNA_threshold_max = 20000
percent.mito_threshold = 0.1
nFeature_RNA_threshold_min = 200

####### normalization parameters
selection_method = "vst"
n_of_variant_feature = 2000

#### heatmap dimension
dimension_plotting_intervals = 10   # n of plots in same pdf
n_cells_for_plot = 200               # n of cols in heatmap
ElbowPlot_ndims = 50

#### UMAP resolution 
UMAP_resolution <- 0.6
UMAP_max_dims <- 50

## feature plot 
plot_feature = c('Slc17a6','Snap25','Syn1','Phox2b','Cd74','Cldn5','Apoe','Prdm12','Trpv1')

##### additional feature to plot after reduction

extra_features = c('Gpr65','Kcng1')

### DEG parameter
DEG.logfc.threshold = 0.1
DEG.min.pct=0.1

### rnk convert parameter
rnk.padj_hold = 0.1
rnk.log2fc_hold =0
#### volcano plot parameters

vol_plot.padj_hold= 0.01
vol_plot.log2fc_hold = 0.3



#----------------------------------------------------- make sure all folder are at the same dir level

object_data <- standard_processing(read=F,save=T)


#object_data <- standard_processing(read=T)  ######## read processed file for next step - processing_reduced_data

#-------------------------------------------------------------------stop here and select cluster you want to keep for further analysis 

idents_you_want_to_keep = c(4,7,8,10,11,15,16,22,24,26,27)

reduced_data <-processing_reduced_data(object_data,read=F,save=T,idents_ = idents_you_want_to_keep)  ### process from here

#reduced_data <-processing_reduced_data(read=T)      # or read previously processed rds file for next step - find_DEG

call_info(reduced_data)

################## if you want to keep all cluster, skip processing_reduced_data 
#reduced_data <- object_data$object$All


#--------------------------------------------------------------- DIY DEG discovery part 


###### find DEGs examples  / supports multiple cluster and group input        no white space in subset_name

######### DEG within all groups in all clusters
find_DEG_bewteen_groups(reduced_data,subset_name='All_cluster',control_group = 'WT',variable_group = 'KP3',logfc.threshold=0.1,min.pct=0.1)

######### DEG between two groups (WT and KP3) in kcng1+ (7 seventh) clusters
find_DEG_bewteen_groups(reduced_data,subset_cluster = 6,subset_name = 'kcng1+',control_group = 'WT',variable_group = 'KP3',save_folder='DEG')  

find_DEG_bewteen_groups(reduced_data,subset_cluster = 8,subset_name = 'Gpr65+',control_group = 'WT',variable_group = 'KP3',save_folder='DEG')  

######## DEG between two groups (WT and KP3) in all clusters but  (7,8 clusters) !!
#find_DEG_bewteen_groups(reduced_data,subset_cluster = c(-7,-8),subset_name = 'cluster 7&8 excluded',control_group = 'WT',variable_group = 'KP3',save_folder='DEG')  

########## cluster6,7 vs cluster5 in all groups
#find_DEG_bewteen_clusters(reduced_data,subset_name='All_group',control_cluster=5,variable_cluster=c(6,7))

########### cluster 7 vs other clusters in WT   #### set one cluster(variable / control)as c()  to use all other clusters
find_DEG_bewteen_clusters(reduced_data,subset_name='WT_only_kcng1+',control_cluster=c(),subset_group='WT',variable_cluster=6)

########### cluster 7 vs other clusters in KP3   #### set one cluster(variable / control) as c() to use all other clusters
find_DEG_bewteen_clusters(reduced_data,subset_name='KP3_only_kcng1+',control_cluster=c(),subset_group='KP3',variable_cluster=6)

########### cluster 7 vs other clusters in WT   #### set one cluster(variable / control) as c()to c() to use all other clusters
find_DEG_bewteen_clusters(reduced_data,subset_name='WT_only_Gpr65+',control_cluster=c(),subset_group='WT',variable_cluster=8)

########### cluster 7 vs other clusters in KP3   #### set one cluster(variable / control) as c() to use all other clusters
find_DEG_bewteen_clusters(reduced_data,subset_name='KP3_only_Gpr65+',control_cluster=c(),subset_group='KP3',variable_cluster=8)

########### cluster 7 vs other clusters in WT   #### set one cluster(variable / control) as c() to use all other clusters
find_DEG_bewteen_clusters(reduced_data,subset_name='WT_only_Gpr65+_kcng1+',control_cluster=6,subset_group='WT',variable_cluster=8)

########### cluster 7 vs other clusters in KP3   #### set one cluster to c() to use all other clusters
find_DEG_bewteen_clusters(reduced_data,subset_name='KP3_only_Gpr65+_kcng1+',control_cluster=6,subset_group='KP3',variable_cluster=8)

#----------------------------------------------------------------- other plots  

############# DEG csv result to GSEA rnk
DEG2RNK(DEG_path='DEG',p_hold=rnk.padj_hold,log2fc_hold=rnk.log2fc_hold,remove_mt_header=F)

############# volcano plot
volcano(DEG_path='DEG',p_hold = vol_plot.padj_hold,log2_fc_hold = 0.4,tops=10)

############ heat map plot   level is the order of plotting,for detailed configuration please refer to key_file
key_heatmap(reduced_data,key_file='DEG_heatmap_GOI_KW_All.xlsx',group_level=c('WT','KP3'),cell_level=c(),slot='scale.data',row_cluster=T)
key_heatmap(reduced_data,key_file='DEG_heatmap_GOI_KW_All.xlsx',group_level=c('WT','KP3'),cell_level=c(),slot='scale.data',row_cluster=F)

key_heatmap(reduced_data,key_file='DEG_heatmap_GOI_KW.xlsx',group_level=c('WT','KP3'),cell_level=c(),slot='scale.data',row_cluster=T)
key_heatmap(reduced_data,key_file='DEG_heatmap_GOI_KW.xlsx',group_level=c('WT','KP3'),cell_level=c(),slot='scale.data',row_cluster=F)

############ GSEA auto
GSEA_batch(
  DEG_path='DEG',
  gene_sets =c(`hallmark gene sets`='mh.all.v2023.1.Mm.symbols.gmt',
               `positional gene sets`='m1.all.v2023.1.Mm.symbols.gmt',
               `curated gene sets`='m2.all.v2023.1.Mm.symbols.gmt',
               `regulatory target gene sets`='m3.all.v2023.1.Mm.symbols.gmt',
               `ontology gene sets`='m5.all.v2023.1.Mm.symbols.gmt',
               `cell type signature gene sets`='m8.all.v2023.1.Mm.symbols.gmt'),
  symbol_chip='Mouse_Gene_Symbol_Remapping_MSigDB.v2023.1.Mm.chip',
  out_dir='GSEA'
)

############## bubble plot of GSEA
GSEA_bubble(GSEA_folder='GSEA',GSEA_fdr_hold=0.5,fdr_top=20)

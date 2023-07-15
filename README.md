# easyScRNA: An easy-to-use single cell analysis pipeline
An automatic, batch-supported scRNA analysis pipeline integrates following functions: Recursive file import and group integration; UMAP dimension reduction and visualization; Multi-plots includes feature plot, violin plot and variable feature heatmap; cluster marker (differentially expressed genes) discovery and visualization; GSEA (Gene Set Enrichment Analysis) and visualization.

## Workflow overview :
<img src="https://github.com/Gico1941/easyScRNA/assets/127346166/05dba2f4-a2f6-4703-879d-5ca79883cece" width="800" />


## Prerequisites :
R (version 4.3.0)

R studio


### What will be automatically built in pipeline ?

Microsoft C++ Build Tools (https://visualstudio.microsoft.com/visual-cpp-build-tools/)

CellTypist (https://github.com/Teichlab/celltypist)

GSEA v4.3.2 for the command line (https://www.gsea-msigdb.org/gsea/downloads.jsp)

Java, Conda (with python) and other R packages ....



## Manual :

###  Installation
1.	Download code and extract files (The “easyScRNA-main” folder will be the work folder).
2.	Create your own R script for running the pipeline. Working with the “example.R” is recommended (with example commands and parameter settings).
3.	Download all necessary R packages plus other softwares and load them with command:
    ```
    source('functions/Main_functions.R')
    ```
### Data orgranization
#### important : for GSEA command line integration, white space is not allowed in any path in work folder
1. To process raw 10X matrix files, a specific file sctructure is required (Since a recursive file reading process is introduced, the program supports multiple samples and groups): 
<img src="https://github.com/Gico1941/easyScRNA/assets/127346166/dfffceeb-0dbc-4c30-8f2e-a75b126d3a2d" width="500" />

                        
2. A heatmap_key.xlsx file is required before plot grouped heatmap ("Cell annotation" is manually specified, "Cluster" is the seurat cluster number which will be presented in UMAP/violin plot, use a ',' to seperate multiple cluster numbers ):
   
sheet 1 :
| Cell annotation  | Cluster |
| ------------- | -------- |
| Cell A | 7       |
| Cell B  | 8,9,10  | 
| .... | .... |

sheet 2 :
| Genes  | Order |
| ------------- | -------- |
| Gene A | 1       |
| Gene B  | 1  | 
| Gene C  | 3  | 
| Gene D  | 2 | 
| .... | .... |

order represents the plotting order of each gene and corresponding cluster


### Data processing
1. Specify the raw file folders :
  ```
  raw_file_folder='Raw'
  ```
and features will be presented in featureplot, violinplot and dotplot :
  ```
  e.g for gene feature plots :  plot_feature = c('Syn1','Cd74','Apoe')
  ```
specify other required numeric parameters (as shown in example.R):
  ```
  nFeature_RNA_threshold_max = 20000
  percent.mito_threshold = 0.1
  nFeature_RNA_threshold_min = 200
  .....
  ```
2. Run the command to start standard process automatically (if save=T is specified, group-inttegratetd and dimentional reduced data will be saved in RDS/Raw.rds
   if read=T is specified, it will skip the standard processing and read the RDS/Raw.rds for downstream analysis):
  ```
  object_data <- standard_processing(save=T)
  ```
  or
  ```
  object_data <- standard_processing(read=T)
  ```
3. object_data will contain two types of seurat object : 1. "All" : Object integrated from all groups . 2. Object of each individual group.
4. Plots will be saved in PLOTs/, QC/ and RDS/
   
### cluster filtering
1. continue with the "standard processing", specify additional features that you may want to plot :
   ```
   extra_features = c('Cldn5')
   ```
or skip it by specify :
   ```
   extra_features = c()
   ```
2. To exclude unwanted cell populations, set :
  ```
  idents_you_want_to_keep = c(4,7,9,10,15,16)
  ```
cluster numbers are corresponding cell cluster numbers shown in UMAP plot 

3.Run the command to re-clustre remaining cells and generate plots
  ```
  reduced_data <-processing_reduced_data(object_data,,save=T,idents_ = idents_you_want_to_keep) 
  ```
4. results will be saved under folder "reduced_PLOTs"
   
5. Auto annotation with CellTypist
```
Model_download(model_path = 'CelltypistModel')           ######### download models 
Model_list(model_path = 'CelltypistModel')             ####### list avalible models

predicted_data <- Runcelltypist(dt = reduced_data,model='Immune_All_Low',majority_voting = T)    ### choose one model from avaliable model list  e.g. 'Immune_All_low'
                   
d1 <- DimPlot(predicted_datadas)
d2 <- DimPlot(predicted_datadas,group.by = "typist_prediction")
d1+d2
```




### DEG discovery
Pipeline contains two modes of DEG discovery :
1. Find DEG between different groups but in same cell populatiton (find DEG in cluster 8 between untreated group and treated group):
   ```
   find_DEG_bewteen_groups(reduced_data,subset_cluster = 8,subset_name = 'Cd74+',control_group = 'UNTREATED',variable_group = 'TREATED',save_folder='DEG')
   ```
subset_cluster is the corresponding cluster number of "reduced_data" (reduced_PLOTs); subset_name is customized label for that population; Control_group and variable_group should be exact match to the origin group names (same as group folder names).

for multiple input :
   ```
   find_DEG_bewteen_groups(reduced_data,subset_cluster = c(4,5,8),subset_name = 'Cd74+',control_group = c('UNTREATED','WATER'),variable_group = C('CHEMICAL_A','CHEMICAL_B'),save_folder='DEG')
   ```
2. Find DEG between different population but in same group (find DEG between cluster 6 and 8 in untreated group):
   ```
   find_DEG_bewteen_clusters(reduced_data,subset_name='Subset_A',control_cluster=6,subset_group='UNTREATED',variable_cluster=8)
   ```
subset_group name should be exact match to group name (group folder name)
   
   ```
   find_DEG_bewteen_clusters(reduced_data,subset_name='Subset_A',control_cluster=c(),subset_group='UNTREATED',variable_cluster=8)
   ```
specify one cluster as variable while set control cluster as "c()" to use all other clusters as control.

   ```
   find_DEG_bewteen_clusters(reduced_data,subset_name='Subset_A',control_cluster=8,subset_group='UNTREATED',variable_cluster=c())
   ```
find DEG between cluster 8 and all other clusters in water and untreated merged group
   ```
   find_DEG_bewteen_clusters(reduced_data,subset_name='Subset_A',control_cluster=8,subset_group=c('WATER','UNTREATED'),variable_cluster=c())
   ```

3. For downstream enrichment analysis, convert DEG .csv file to .rnk file with :
   ```
   DEG2RNK(DEG_path='DEG',p_hold=0.1,log2fc_hold=0)
   ```

Generate volcano plot for all DEG .csv files :
```
volcano(DEG_path='DEG',p_hold = 0.01,log2_fc_hold = 0.4,tops=10)
```

tops is the number of top DEGs that will be highlighted

### heatmap visualization
To visualize gene expression difference, generate heatmap with heatmap_key.xlsx and reduced_data. Plot the row by manual order :
```
key_heatmap(reduced_data,key_file='DEG_heatmap_key.xlsx',group_level=c('WT','KP3'),cell_level=c(6,7,8),slot='scale.data',row_cluster=F)

```
Or Plot the row by cluster order :
```
key_heatmap(reduced_data,key_file='DEG_heatmap_key.xlsx',group_level=c('WT','KP3'),cell_level=c(6,7,8),slot='scale.data',row_cluster=F)

```

group_level / cell_level represents the plotting orders and should be exactly match to group folder name / cluster numbers provided in heatmap_key.xlsx

### Gene set enrichment and visualization
To perform batch Gene set enrichment analysis with GSEA software and MsigDB, run (if host is mouse):
```
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
```
or replace the symbol_chip and gene_sets for human data analysis
```
    gene_sets =c(`hallmark gene sets`='h.all.v2023.1.Hs.symbols.gmt',
                 `positional gene sets`='c1.all.v2023.1.Hs.symbols.gmt',
                 `curated gene sets`='c2.all.v2023.1.Hs.symbols.gmt',
                 `regulatory target gene sets`='c3.all.v2023.1.Hs.symbols.gmt',
                 `ontology gene sets`='c6.all.v2023.1.Hs.symbols.gmt',
                 `cell type signature gene sets`='c8.all.v2023.1.Hs.symbols.gmt')
    symbol_chip='Human_Gene_Symbol_with_Remapping_MSigDB.v2023.1.Hs.chip'
```

visualization with bubble plot :
```
GSEA_bubble(GSEA_folder='GSEA',GSEA_fdr_hold=0.5,fdr_top=20)
```
specify fdr_top to plot the top enriched pathways ranked by fdr



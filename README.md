# easyScRNA: An integrated, easy-to-use pipeline for standard single cell sequencing analysis
An automatic, batch-supported scRNA analysis pipeline integrates following functions: Recursive file import and group integration; UMAP dimension reduction and visualization; Multi-plots includes feature plot, violin plot and variable feature heatmap; cluster marker (differentially expressed genes) discovery and visualization; GSEA (Gene Set Enrichment Analysis) and visualization.

## Workflow overview :
<img src="https://github.com/Gico1941/easyScRNA/assets/127346166/05dba2f4-a2f6-4703-879d-5ca79883cece" width="800" />


## Prerequisites :
R (version 4.3.0)

R studio


### What will be automatically built in pipeline ?
Seurat 4.3.0 (https://satijalab.org/seurat)

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
#### Example results (use GSE166635 as example)
#### QC of GSM5076750
<img src="https://github.com/Gico1941/easyScRNA/assets/127346166/ba54faaa-d8d4-4815-b207-2a36fba467ea" width="300" />

#### UMAP of GSM5076750
<img src="https://github.com/Gico1941/easyScRNA/assets/127346166/7ef0801e-ff22-4446-a28a-c7e7d50d40bb" width="300" />

#### heatmap of top100 variable genes
<img src="https://github.com/Gico1941/easyScRNA/assets/127346166/6fe24aa5-6c27-490f-8948-173f326beb37" width="300" />

#### Split violin plot for certain feature
<img src="https://github.com/Gico1941/easyScRNA/assets/127346166/6ddd1c5e-af53-4d50-b10c-4316564f25bc" width="300" />

#### Dot plot, Feature plot, PCA and other plots

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




### DEG (differentially expressed genes) discovery
Pipeline contains two modes of DEG discovery :
1. Find DEG between different groups but for same cell populatiton (find DEG in cluster 6 between HCC1 group and HCC2 group):
   ```
   find_DEG_between_groups(reduced_data,subset_cluster = 6,subset_name = 'cluster6',control_group = 'HCC1',variable_group = 'HCC2',save_folder='DEG') 
   ```
   
![1705878888644](https://github.com/Gico1941/easyScRNA/assets/127346166/17347345-9562-4c5d-bad3-9b5de9086062)


subset_cluster is the corresponding cluster number of "reduced_data" (reduced_PLOTs); subset_name is customized label for that population; Control_group and variable_group should be exact match to the origin group names (same as group folder names).

For multiple subset input :
   ```
   find_DEG_bewteen_groups(reduced_data,subset_cluster = c(4,5,8),subset_name = 'subsetA',control_group = c('UNTREATED','WATER'),variable_group = C('CHEMICAL_A','CHEMICAL_B'),save_folder='DEG')
   ```
2. Find DEG between different population but in same group (find DEG between cluster 6 and 8 in untreated group):
   ```
   find_DEG_bewteen_clusters(reduced_data,subset_name='Subset_A',control_cluster=6,subset_group='UNTREATED',variable_cluster=8)
   ```

subset_group name should be exact match to group name (group folder name), specify one cluster as variable while set control cluster as "c()" to use all other clusters as control.
   
   ```
   find_DEG_between_clusters(reduced_data,subset_name='HCC1_only_cluster6_vs_others',control_cluster=c(),subset_group='HCC1',variable_cluster=6)
   ```
<img src="https://github.com/Gico1941/easyScRNA/assets/127346166/635fa2bc-db01-4430-8d0a-e44f7ebc22ed" width="800" />

#### Or
   ```
   find_DEG_bewteen_clusters(reduced_data,subset_name='Subset_A',control_cluster=8,subset_group='UNTREATED',variable_cluster=c())
   ```
find DEG between cluster 8 and all other clusters in water and untreated merged group
   ```
   find_DEG_bewteen_clusters(reduced_data,subset_name='Subset_A',control_cluster=8,subset_group=c('WATER','UNTREATED'),variable_cluster=c())
   ```

3. For downstream GSEA (Gene set enrichment analysis), convert DEG .csv file to .rnk file with the command :
   ```
   DEG2RNK(DEG_path='DEG',p_hold=0.1,log2fc_hold=0)
   ```

That will add three .rnk files in addition to the DEG .csv file


![1705879279416](https://github.com/Gico1941/easyScRNA/assets/127346166/ba35ec72-2e6f-4051-befc-1d47abee2450)


Generate volcano plot for all DEG .csv files (tops: number of top dots you want to annotate and highlight):
```
volcano(DEG_path='DEG',p_hold = vol_plot.padj_hold,log2_fc_hold = 0.1,tops=10,highlight_top='avg_log2FC',highlight_by_keys=F,height=7,width=4)
```
<img src="https://github.com/Gico1941/easyScRNA/assets/127346166/2a562cfc-b6e2-4b47-9658-4408af58a4b8" width="600" />


To highlight specific genes in volcanoplot: create a "highlight_key.txt" under the subfolders of DEG, e.g:

highlight_key.txt :
```
Cd8a
Pdcd1
Cd3e
... 
```
```
volcano(DEG_path='DEG',p_hold = vol_plot.padj_hold,log2_fc_hold = 0.1,tops=10,highlight_by_keys=T,width=4,height=4)

```
For example (highlight the CD55,IL7R,FABP1 in previous volcano plot)

<img src="https://github.com/Gico1941/easyScRNA/assets/127346166/0ba547b2-bf8f-4f5c-bb95-9efa5b20b5f9" width="600" />


### heatmap visualization
To visualize gene expression difference, generate heatmap with heatmap_key.xlsx and reduced_data :
```
key_heatmap(object_ = reduced_data, 
                        key_file='test.xlsx',
                        group_level=c('WT','KP3'),
                        slot='scale.data',
                        row_cluster=F,
                        col_cluster=F,
                        aggregate='Cell') # aggregate by "cell" or "group"

```
### heatmap example 1  | Aggregate by Cell:
sheet 1 :
| Cell annotation  | Cluster |
| ------------- | -------- |
| cell set1 | 6       |
| cell set2  | 7,8,9  | 
| cell set3 | 5 |


sheet 2 :
| Genes  | Order |
| ------------- | -------- |
| CD55 | 1       |
| IL7R  | 2  | 
| FABP1  | 3  | 


```
key_heatmap(object_ = reduced_data, 
                        key_file='DEG_heatmap_key_example.xlsx',
                        group_level=c('HCC1','HCC2'),
                        slot='scale.data',
                        row_cluster=F,
                        col_cluster=F,
                        aggregate='Cell',
            group_color = c('HCC1'='red','HCC2'='blue'))

```

<img src="https://github.com/Gico1941/easyScRNA/assets/127346166/b29f0a62-5677-4bb7-b2ee-38571eab044d" width="600" />


### heatmap example 2  | Aggregate by Group:


sheet 1 :
| Cell annotation  | Cluster |
| ------------- | -------- |
| cell set1 | 6       |
| cell set2  | 7,8,9  | 
| cell set3 | 5 |


sheet 2 :
| Genes  | Order |
| ------------- | -------- |
| CD55 | 1       |
| IL7R  | 2  | 
| FABP1  | 3  | 


```
key_heatmap(object_ = reduced_data, 
            key_file='DEG_heatmap_key_example.xlsx',
            group_level=c('HCC1','HCC2'),
            slot='scale.data',
            row_cluster=F,
            col_cluster=F,
            aggregate='Group',
            group_color = c('HCC1'='red','HCC2'='blue'))
```

<img src="https://github.com/Gico1941/easyScRNA/assets/127346166/fbadff9b-4921-4c33-a63c-dccf7f117f2f" width="600" />



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
  out_dir='GSEA',
  GSEA_plots_number=30
)   #GSEA_plots_number=30 : max lines displayed in GSEA result (default : top 20)

```
or replace the symbol_chip and gene_sets for human data analysis
```
GSEA_batch(
  DEG_path='DEG',
  gene_sets =c(`hallmark gene sets`='h.all.v2023.1.Hs.symbols.gmt',
                 `positional gene sets`='c1.all.v2023.1.Hs.symbols.gmt',
                 `curated gene sets`='c2.all.v2023.1.Hs.symbols.gmt',
                 `regulatory target gene sets`='c3.all.v2023.1.Hs.symbols.gmt',
                 `ontology gene sets`='c6.all.v2023.1.Hs.symbols.gmt',
                 `cell type signature gene sets`='c8.all.v2023.1.Hs.symbols.gmt')
  symbol_chip='Human_Gene_Symbol_with_Remapping_MSigDB.v2023.1.Hs.chip',
  out_dir='GSEA',
  GSEA_plots_number=30
)   #GSEA_plots_number=30 : max lines displayed in GSEA result (default : top 20)
```

visualization with bubble plot :
```
GSEA_bubble(GSEA_folder='GSEA',GSEA_fdr_hold=0.5,fdr_top=20)
```
specify fdr_top to plot the top enriched pathways ranked by fdr



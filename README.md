##  Project: Study of Diabetes in MES13 cells
### Study ID: 
### Scientist: Wenji Li
### Data Analysis: Davit Sargsyan 
### Created: 09/11/2017 

---    

## Table of Contents
[Files](#files)
[Daily Logs](#logs)  
[Results](#results)   

## To Do (01/22/2018)
1. use edgeR to get p-vals and select (Done, 01/24/2018)

## Files<a name="files"></a>
FastQ file legend (both, Methyl-seq and RNA-seq):    
~/kidney.ro1/mes13/docs/legend.xlsx    

1. Original FastQ files are located here (4Tb internal hard drive):    
*/media/administrator/datastorage/FastQ_2017/Methyl_seq/July/pools/Kong_MouseMethyl_pool6.zip*    
*/media/administrator/datastorage/FastQ_2017/RNA_seq/BGI_RNA*   
   
2. All alignemnt files and documents are located here (1Tb internal SSD):    
a. Methyl-seq processed by Davit:    
*/media/administrator/datastorage/Processed_BAM_Files/David_MethylSeq_Processed/Davit/BAM Files*   
b. Methyl-seq processed by Renyi:  
*/media/administrator/datastorage/Processed_BAM_Files/Wenji_MES13_MethylSeq_Processed/Renyi_12292017/Results*, *combined_WJ_anno.xlsx*   
c. RNA-seq processed by Davit:      
*/media/administrator/datastorage/Processed_BAM_Files/David_BGI_RNAseq_Processed/Davit*   
NOTE: there are more, processed by Davit, Renyi and John, see *David_MethylSeq_Processed* and *David_BGI_RNAseq_Processed* folders.    
NOTE: some of the files were too large to push to GitHub; saved on the Lab230 machine in *data* folders.    
d. RNA-seq processed by Renyi    
*/media/administrator/datastorage/Processed_BAM_Files/Wenji_MES13_RNASeq_Processed/Renyi_12292017/Results_Hisat2*, *mes13_featurecounts_Dec2017_david.csv* (raw counts) and *mes13_fpkm_Dec2017_david.csv* (FPKM data converted from counts by DESeq2)

## Daily Logs<a name="logs"></a>
### 01/24/2018
* Added script for RNA-seq analysis with *edgeR*

### 01/16/2018
* Methyl-seq from RUCDR rerun is imported to R; first analysis.

### 11/06/2017
* Added SD estimation with difference vs. mean method(*mes13_rnaseq_lg_vs_hg_diff_vs_mean_pval_v1.R* and *mes13_rnaseq_hg_vs_mic_diff_vs_mean_pval_v1.R*).        
* Edited this (*README.md*) file updated FastQ and BAM file locations as the files were moved to the 4Tb internal storage (*datastorage*).   

### 09/27/2017
* Recreatedfigures 9 and 10 for better resolution
* Lists of genes from inflamation/oxidative stress pathways

### 09/23/2017
* Added db/db methyl-seq analysis

### 09/20/2017
* Added reverse mapping to TNF signaling pathway etc.

### 09/18/2017
* Added network plots
* Pathway lists explored with KEGG.db and reactome.db R databases

### 09/14/2017
* Improved differentially methylated gene analysis
* Added pathway analysis (*ReactomePA* package)

### 09/12/2017
* Added log2 difference analysis of % methilation
* Analysis script Version 2: collapse CpGs by regions and gene

### 09/11/2017
* Added analysis script for % methylation

## Results <a name="results"></a>
a. RO1 porposal results:   
~/kidney.ro1/mes13/docs/results_09142017.pptx
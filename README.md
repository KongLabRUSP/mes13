##  Project: Study of Diabetes in MES13 cells
### Study ID: 
### Scientist: Wenji Li, David, Renyi
### Data Analysis: Davit Sargsyan 
### Created: 09/11/2017 

---    

## Table of Contents
[file legend](#leg)
[Daily Logs](#logs)  
[Results](#results)   
[Files](#files)

## File Legend<a name="files"></a>
### Scripts
**mes13_rnaseq_DEGseq_MITC_v1.R**: current (05/12/2018) version of RNA-seq data analysis script for LG, HG and MITC. Used for David's paper and for Kidney RO1.        
**mes13_methylseq_DSS_MITC_v2.R**: current (05/12/2018) version of DNA methyl-seq data analysis script for LG, HG and MITC. Used for David's paper and for Kidney RO1. This version DOES NOT actually use package DSS even though it is in the title.   
    
**mes13_rnaseq_DEGseq_TIIA_v2.R**: current (05/12/2018) version of RNA-seq data analysis script for LG, HG and TIIA. Used for Wenji's paper.    
**mes13_methylseq_DSS_TIIA_v1.R**: current (05/12/2018) version of DNA methyl-seq data analysis script for LG, HG and MITC. Used for wenji's paper. This version DOES NOT actually use package DSS even though it is in the title.

### Raw Data
FastQ file legend (both, Methyl-seq and RNA-seq):    
~/mes13/docs/legend.xlsx    

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
### To Do (01/22/2018)
1. use edgeR to get p-vals and select (Done, 01/24/2018)

### 05/12/2018
* Methyl-seq analysis of LG, HG and MITC for RO1 (**mes13_methylseq_DSS_MITC_v2.R**)

### 04/13/2018
* Added DEGseq script for RNA-seq TIIA (Wenji)

### 03/15/2018
* Analysis plan (Davit)    
a. We are splitting data to LG, HG and MIC (David) and LG, HG and TII (Wenji), and rerunning DEGseq analysis of RNA-seq    
b. We are splitting Methyl-seq data same way and analyzing using John's DMRfinder to  cluster CpG followed by a new method to get p-values for data with no reduplicates.    

### 02/22/2018
* Finished *DEGseq* analysis of MES13 data, all treatments vs. HG

### 01/29/2018
* Added script for RNA-seq analysis with *DEGseq*

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
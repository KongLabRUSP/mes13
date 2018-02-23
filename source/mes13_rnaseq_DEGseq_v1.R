# |----------------------------------------------------------------------------------|
# | Project: Study of Diabetes in MES13 cells                                        |
# | Script: RNA-seq data analysis and visualization using edgeR                      |
# | Author: Davit Sargsyan                                                           |
# | Created: 01/29/2018                                                              |
# | Modified:                                                                        |
# |----------------------------------------------------------------------------------|
# sink(file = "tmp/log_mes13_rnaseq_DEGseq_v1.R")
# Source: 
# https://bioconductor.org/packages/release/bioc/html/DEGseq.html

# source("https://bioconductor.org/biocLite.R")
# biocLite("DEGseq")

# Header----
require(data.table)
require(ggplot2)
require(DEGseq)
require(knitr)

# MES13 data----
# Treatment legend----
trt.names <- c("LG",
               "HG",
               "MIC 1.5uM",
               "TIIA-5 uM",
               "FX 1uM",
               "Gen-10 uM",
               "Ber 6uM")

# Load data----
# Question to Renyi: how was the data processed and annotated?
# If by John's DMR Finder, what were the settings?
# dt1 <- fread("data/Renyi_RNAseq_12292017/mes13_fpkm_Dec2017_david.csv")
dt1 <- fread("data/Renyi_RNAseq_12292017/mes13_featurecounts_Dec2017_david.csv",
             skip = 1)
dt1

# Keep only gene IDs and counts
dt1 <- dt1[, c("Geneid",
               "WJ1.dedup.bam",
               "WJ2.dedup.bam",
               "WJ3.dedup.bam",
               "WJ4.dedup.bam",
               "WJ5.dedup.bam",
               "WJ6.dedup.bam",
               "WJ7.dedup.bam")]
dt1
names(dt1) <- c("gene", 
                paste("WJ",
                      1:7,
                      sep = ""))
dt1

# Remove genes with low counts----
summary(dt1[, -1])
tmp <- rowSums(dt1[, -1])
# Remove if total across 5 samples is no more than 20
dt1 <- droplevels(subset(dt1,
                         tmp > 20))
dt1
# 13,954 genes left, down from 24,421 genes

# DEGseq----
write.table(dt1,
            file = "tmp/dt1.txt",
            sep = "\t",
            row.name = FALSE)

ndx <- c(2, 4:8)

for(i in ndx){
  DEGexp2(geneExpFile1 = "tmp/dt1.txt", 
          geneCol1 = 1, 
          expCol1 = 3, 
          groupLabel1 = "HG",
          
          geneExpFile2 = "tmp/dt1.txt", 
          geneCol2 = 1, 
          expCol2 = i,
          groupLabel2 = trt.names[i - 1],
          
          foldChange = 2,
          qValue = 0.1,
          thresholdKind = 5, 
          rawCount = TRUE,
          normalMethod = "none",
          method = "MARS",
          outputDir = "tmp")
  
  # Read in the output file
  dt2 <- fread("tmp/output_score.txt")
  dt2
  
  # Write as CSV----
  write.csv(dt2,
            file = paste("tmp/mes13_rnaseq_DEGseq_HG-",
                         trt.names[i - 1],
                         ".csv",
                         sep = ""),
            row.names = FALSE)
}

# Significant genes----
tmp <- dt2[dt2$`Signature(q-value(Storey et al. 2003) < 0.1)`, ]
tmp
write.csv(tmp,
          file = "tmp/mes13_rnaseq_DEGseq2_signif_out.csv",
          row.names = FALSE)

# sink()
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

# DEGset Example-----
# Read files
geneExpFile <- system.file("extdata", "GeneExpExample5000.txt", package="DEGseq")
geneExpMatrix1 <- readGeneExp(file=geneExpFile, geneCol=1, valCol=c(7,9,12,15,18))
geneExpMatrix2 <- readGeneExp(file=geneExpFile, geneCol=1, valCol=c(8,10,11,13,16))
head(geneExpMatrix1)
head(geneExpMatrix2)

# Compare treatments with replicates
DEGexp(geneExpMatrix1=geneExpMatrix1, expCol1=2, groupLabel1="kidneyR1L1", 
       geneExpMatrix2=geneExpMatrix2, expCol2=2, groupLabel2="liverR1L2",
       replicateExpMatrix1=geneExpMatrix1, expColR1=3, replicateLabel1="kidneyR1L3", 
       replicateExpMatrix2=geneExpMatrix1, expColR2=4, replicateLabel2="kidneyR1L7",
       method="MATR", 
       outputDir = "tmp")

# Compare two samples: R1L1Kidney & R1L3Kidney
geneExpMatrix1 <- readGeneExp(file=geneExpFile, geneCol=1, valCol=c(7))
geneExpMatrix2 <- readGeneExp(file=geneExpFile, geneCol=1, valCol=c(9))
head(geneExpMatrix1)
head(geneExpMatrix2)

DEGexp(geneExpMatrix1=geneExpMatrix1, expCol1=2, groupLabel1="kidneyR1L1", 
       geneExpMatrix2=geneExpMatrix2, expCol2=2, groupLabel2="liverR1L2",
       method="CTR",
       outputDir = "tmp")

# Calculate differences and q-values etc.----
DEGexp2(geneExpFile1 = geneExpFile, 
        geneCol1 = 1, 
        expCol1 = 7, 
        groupLabel1="R1L1Kidney",
        
        geneExpFile2 = geneExpFile, 
        geneCol2 = 1, 
        expCol2 = 9,
        groupLabel2="R1L3Kidney",
        
        outputDir = "C:/git_local/mes13/tmp")

# Read in teh output file----
dt1 <- fread("tmp/output_score.txt")
dt1

# Clean all----
rm(list = ls())

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
dt1 <- fread("data/Renyi_RNAseq_12292017/mes13_featurecounts_Dec2017_david.csv")
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

# # Remove genes with low counts----
# summary(dt1[, -1])
# tmp <- rowSums(dt1[, -1])
# # Remove if total across 5 samples is no more than 100
# dt1 <- droplevels(subset(dt1,
#                          tmp > 100))
# dt1

# DEGseq----C:/git_local/mes13/
write.table(dt1,
            file = "tmp/dt1.txt",
            sep = "\t",
            row.name = FALSE)

DEGexp2(geneExpFile1 = "tmp/dt1.txt", 
        geneCol1 = 1, 
        expCol1 = 2, 
        groupLabel1="LG",
        
        geneExpFile2 = "tmp/dt1.txt", 
        geneCol2 = 1, 
        expCol2 = 3,
        groupLabel2="HG",
        
        foldChange = 2,
        qValue = 0.1,
        thresholdKind = 5, 
        rawCount = TRUE,
        normalMethod = "none",
        method = "MARS",
        outputDir = "tmp")

# Read in teh output file----
dt2 <- fread("tmp/output_score.txt")
dt2

# Write as CSV----
write.csv(dt2,
          file = "tmp/mes13_rnaseq_DEGseq2_out.csv",
          row.names = FALSE)

# Significant genes----
tmp <- dt2[dt2$`Signature(q-value(Storey et al. 2003) < 0.1)`, ]
tmp
write.csv(tmp,
          file = "tmp/mes13_rnaseq_DEGseq2_signif_out.csv",
          row.names = FALSE)

# sink()
# |----------------------------------------------------------------------------------|
# | Project: Study of Diabetes in MES13 cells                                        |
# | Script: RNA-seq data analysis and visualization using edgeR, TIIA only (Wenji)   |
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
# dt1 <- fread("data/Renyi_RNAseq_12292017/mes13_fpkm_Dec2017_david.csv")
dt1 <- fread("data/rna_seq/Renyi_12292017/mes13_featurecounts_Dec2017_david.csv",
             skip = 1)
dt1

# CHECK
dt1[Geneid == "Tnfrsf25",]
dt1[Geneid == "Fcer1g",]

# Keep only gene IDs and counts
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
colnames(dt1) <- c("gene",
                trt.names)

dt1

# Remove genes with low counts----
summary(dt1[, -1])
tmp <- rowSums(dt1[, -1])
# Remove if total across 3 samples is no more than 10
dt1 <- droplevels(subset(dt1,
                         tmp > 20))
dt1
# 13,954 genes left, down from 24,421 genes

# Leave the 3 treatments only----
dt1 <- subset(dt1,
              select = c(1:3, 5))
dt1

# DEGseq----
# a. (HG - LG)----
DEGexp(geneExpMatrix1 = dt1,
       geneCol1 = 1, 
       expCol1 = 3, 
       groupLabel1 = colnames(dt1)[3],
       
       geneExpMatrix2 = dt1,
       geneCol2 = 1, 
       expCol2 = 2,
       groupLabel2 = colnames(dt1)[2],
       
       foldChange = 2,
       qValue = 0.1,
       thresholdKind = 5, 
       rawCount = TRUE,
       normalMethod = "none",
       method = "MARS",
       outputDir = "tmp")

hg_lg <- fread("tmp/output_score.txt")
hg_lg
hg_lg[hg_lg$`Signature(q-value(Storey et al. 2003) < 0.1)`,]

# Write as CSV----
write.csv(hg_lg,
          file = "tmp/mes13_rnaseq_DEGseq_HG-LG.csv",
          row.names = FALSE)

# MA Plot----
hg_lg[, mu := (log2(value1) + log2(value2))/2]
hg_lg[, diff := log2(value1) - log2(value2)]

tiff(filename = "tmp/mes13_rnaseq_DEGseq_HG-LG_maplot.tiff",
     height = 6,
     width = 6,
     units = 'in',
     res = 300,
     compression = "lzw+p")
plot(hg_lg$diff ~ hg_lg$mu,
     pch = ".",
     xlab = "Mean",
     ylab = "Difference",
     main = "MES13 Gene Expression, HG-LG, FDR < 0.1")
points(hg_lg$diff[hg_lg$`q-value(Storey et al. 2003)` < 0.3 & hg_lg$diff > 0] ~ hg_lg$mu[hg_lg$`q-value(Storey et al. 2003)` < 0.3 & hg_lg$diff > 0] ,
       pch = "x",
       col = "green")
points(hg_lg$diff[hg_lg$`q-value(Storey et al. 2003)` < 0.3 & hg_lg$diff < 0] ~ hg_lg$mu[hg_lg$`q-value(Storey et al. 2003)` < 0.3 & hg_lg$diff < 0] ,
       pch = "x",
       col = "red")
abline(h = c(-0.3, 0.3),
       lty = 2)
graphics.off()

# b. (TIIA - HG)----
DEGexp(geneExpMatrix1 = dt1,
       geneCol1 = 1, 
       expCol1 = 4, 
       groupLabel1 = colnames(dt1)[4],
       
       geneExpMatrix2 = dt1,
       geneCol2 = 1, 
       expCol2 = 3,
       groupLabel2 = colnames(dt1)[3],
       
       foldChange = 2,
       qValue = 0.1,
       thresholdKind = 5, 
       rawCount = TRUE,
       normalMethod = "none",
       method = "MARS",
       outputDir = "tmp")

tiia_hg <- fread("tmp/output_score.txt")
tiia_hg
tiia_hg[tiia_hg$`Signature(q-value(Storey et al. 2003) < 0.1)`, ]

# Write as CSV----
write.csv(tiia_hg,
          file = "tmp/mes13_rnaseq_DEGseq_TIIA-HG.csv",
          row.names = FALSE)

# MA Plot----
tiia_hg[, mu := (log2(value1) + log2(value2))/2]
tiia_hg[, diff := log2(value1) - log2(value2)]

tiff(filename = "tmp/mes13_rnaseq_DEGseq_TIIA-HG_maplot.tiff",
     height = 6,
     width = 6,
     units = 'in',
     res = 300,
     compression = "lzw+p")
plot(tiia_hg$diff ~ tiia_hg$mu,
     pch = ".",
     xlab = "Mean",
     ylab = "Difference",
     main = "MES13 Gene Expression, TIIA-HG, FDR < 0.1")
points(tiia_hg$diff[tiia_hg$`q-value(Storey et al. 2003)` < 0.2 & tiia_hg$diff > 0] ~ tiia_hg$mu[tiia_hg$`q-value(Storey et al. 2003)` < 0.2 & tiia_hg$diff > 0] ,
       pch = "x",
       col = "green")
points(tiia_hg$diff[tiia_hg$`q-value(Storey et al. 2003)` < 0.2 & tiia_hg$diff < 0] ~ tiia_hg$mu[tiia_hg$`q-value(Storey et al. 2003)` < 0.2 & tiia_hg$diff < 0] ,
       pch = "x",
       col = "red")
abline(h = c(-0.3, 0.3),
       lty = 2)
graphics.off()

# Venn diagram----
g1 <- hg_lg[`q-value(Storey et al. 2003)` < 0.5 & 
              `log2(Fold_change) normalized` > 0.3,]$GeneNames
# 263 genes
g2 <- hg_lg[`q-value(Storey et al. 2003)` < 0.5 & 
              `log2(Fold_change) normalized` < -0.3,]$GeneNames
# 207 genes

g3 <- tiia_hg[`q-value(Storey et al. 2003)` < 0.5 & 
                `log2(Fold_change) normalized` > 0.3,]$GeneNames
# 1,120 genes
g4 <- tiia_hg[`q-value(Storey et al. 2003)` < 0.5 & 
                `log2(Fold_change) normalized` < -0.3,]$GeneNames
# 1,393 genes

up.dn <- g1[g1 %in% g4]
# 124 genes
dn.up <- g2[g2 %in% g3]
# 89 genes

# Heatmap----
ll <- unique(c(up.dn,
               dn.up))

t1 <- merge(hg_lg[hg_lg$GeneNames %in% ll, 
                  c("GeneNames",
                    "log2(Fold_change) normalized")],
            tiia_hg[tiia_hg$GeneNames %in% ll, 
                    c("GeneNames",
                      "log2(Fold_change) normalized")],
            by = "GeneNames")
colnames(t1) <- c("Gene",
                  "HG-LG",
                  "TIIA-HG")
t1 <- t1[order(t1$`HG-LG`,
               decreasing = TRUE), ]
t1
write.csv(t1,
          file = "tmp/mes13_tiia_genes_q-0.5_log2-0.3.csv",
          row.names = FALSE)

ll <- melt.data.table(data = t1,
                      id.vars = 1,
                      measure.vars = 2:3,
                      variable.name = "Comparison",
                      value.name = "Gene Expression Diff")
ll$Comparison <- factor(ll$Comparison,
                        levels = c("TIIA-HG", 
                                   "HG-LG"))
lvls <- ll[ll$Comparison == "HG-LG", ]
ll$Gene <- factor(ll$Gene,
                  levels = lvls$Gene[order(lvls$`Gene Expression Diff`)])
# Keep top 100 genes for the plot----
gene.keep <- unique(ll$Gene[order(abs(ll$`Gene Expression Diff`)) < 101])
ll <- droplevels(subset(ll,
                        Gene %in% gene.keep))
ll

p1 <- ggplot(data = ll) +
  coord_polar("y",
              start = 0,
              direction = -1) +
  geom_tile(aes(x =  as.numeric(Comparison),
                y = Gene, 
                fill = `Gene Expression Diff`),
            color = "white") +
  geom_text(data = ll[Comparison == "HG-LG", ],
            aes(x = rep(1.75,
                        nlevels(Gene)),
                y = Gene,
                label = unique(Gene),
                angle = 90 + seq(from = 0,
                            to = 360,
                            length.out = nlevels(Gene))[as.numeric(Gene)]),
            hjust = 0) +
  scale_fill_gradient2(low = "red", 
                       high = "green", 
                       mid = "grey", 
                       midpoint = 0, 
                       name = "Gene Expr Diff") +
  scale_x_continuous(limits = c(0, 
                                max(as.numeric(ll$Comparison)) + 0.5),
                     expand = c(0, 0)) + 
  scale_y_discrete("",
                   expand = c(0, 0)) +
  ggtitle("Changes in Gene Expression
          Fold-Change > 0.3 and q-Value < 0.5") + 
  theme(plot.title = element_text(hjust = 0.5),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank())
p1

tiff(filename = "tmp/mes13_rnaseq_DEGseq_heatmap.tiff",
     height = 10,
     width = 10,
     units = 'in',
     res = 300,
     compression = "lzw+p")
print(p1)
graphics.off()

# sessionInfo()
# sink()
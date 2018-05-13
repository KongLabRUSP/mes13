# |----------------------------------------------------------------------------------|
# | Project:     Study of Diabetes in MES13 cells, LG, HG and MITC                   |
# | Script:      Methyl-seq data analysis and visualization for RO1 grant            |
# | Coordinator: Wenji, David, Renyi Wu                                              |
# | Author:      Davit Sargsyan                                                      |
# | Created:     05/12/2018                                                          |
# | Modified:                                                                        |
# |----------------------------------------------------------------------------------|
# sink(file = "tmp/log_mes13_methylseq_DSS_MITC_v1.txt")
date()

# NOTE: several packages, e.g. Rcpp, MASS, etc., might be deleted manually and reinstalled
# Workflow: https://www.bioconductor.org/help/workflows/rnaseqGene/
# source("https://bioconductor.org/biocLite.R")
# biocLite("TxDb.Mmusculus.UCSC.mm10.knownGene")
# biocLite("ChIPseeker")
# biocLite("org.Mm.eg.db")

require(data.table)
require(ggplot2)
require(knitr)

require(ChIPseeker)
require(TxDb.Mmusculus.UCSC.mm10.knownGene)

# Load and view raw counts (no annoation)----
dt01 <- fread("data/Renyi_Methylseq_12292017/combined_WJ_anno.csv")
dt01

# NOTE: there are 14 rows representing mitochondrial DNA
unique(dt01$chr)
dt01[dt01$chr == "chrM",]

# Annotate----
peakAnno1 <- annotatePeak(peak = "data/Renyi_Methylseq_12292017/combined_WJ_anno.csv", 
                          tssRegion = c(-3000, 3000), 
                          TxDb = TxDb.Mmusculus.UCSC.mm10.knownGene,
                          annoDb = "org.Mm.eg.db")

head(peakAnno1@detailGenomicAnnotation)
t1 <- peakAnno1@annoStat
t1$Feature <- factor(t1$Feature,
                     levels = as.character(t1$Feature[order(t1$Frequency,
                                                            decreasing = TRUE)]))
t1
p1 <- ggplot(t1,
             aes(x = "",
                 y = Frequency,
                 fill = Feature)) +
  geom_bar(width = 1, 
           stat = "identity",
           color = "black") +
  coord_polar("y",
              start = 0,
              direction = -1) +
  scale_x_discrete("") +
  ggtitle("Annotation by Region (%)")
p1 

tiff(filename = "tmp/mes13_anno_by_reg.tiff",
     height = 5,
     width = 5,
     units = 'in',
     res = 300,
     compression = "lzw+p")
print(p1)
graphics.off()

# Make data table----
dt1 <- data.table(start = peakAnno1@anno@ranges@start,
                  as.data.frame(peakAnno1@anno@elementMetadata@listData))

# Remove unmapped regions
dt1 <- dt1[!is.na(dt1$SYMBOL == "NA"), ]
# Removed 12 rows

# Subset data: LG, HG and MITC----
dt1 <- data.table(gene = dt1$SYMBOL,
                  anno = dt1$annotation,
                  geneId = dt1$geneId,
                  chr = dt1$geneChr,
                  pos = dt1$start,
                  reg = NA,
                  dt1[, CpG:WJ03.X],
                  geneName = dt1$GENENAME)
dt1

# Remove genes with all NAs
dt1 <- droplevels(dt1[rowSums(dt1[, WJ01.N:WJ03.X],
                              na.rm = TRUE) > 0, ])

# Regions----
kable(data.table(table(substr(dt1$anno, 1, 9))))
  # |V1        |     N|
  # |:---------|-----:|
  # |3' UTR    |  4273|
  # |5' UTR    |   745|
  # |Distal In | 58992|
  # |Downstrea |  2572|
  # |Exon (uc0 | 11272|
  # |Intron (u | 52225|
  # |Promoter  | 80765|

# Separate Promoter, Body and Downstream----
dt1$reg <- as.character(dt1$anno)

# a. Promoter: up to 3kb upstream
dt1$reg[substr(dt1$anno, 
               1,
               8) == "Promoter"] <- "Promoter"

# b. Body: exons and introns
dt1$reg[substr(dt1$anno, 
               1, 
               4) %in% c("Exon",
                         "Intr")] <- "Body"

# c. Downstream: Distal Intergenic and  Downstream
dt1$reg[substr(dt1$anno, 
               1, 
               4) %in% c("Dist",
                         "Down")] <- "Downstream"

dt1$reg <- factor(dt1$reg,
                  levels = c("Promoter",
                             "5' UTR",
                             "Body",
                             "3' UTR",
                             "Downstream"))
kable(data.table(table(dt1$reg)))
  # |V1         |     N|
  # |:----------|-----:|
  # |Promoter   | 80765|
  # |5' UTR     |   745|
  # |Body       | 63497|
  # |3' UTR     |  4273|
  # |Downstream | 61564|

# CpG distribution and coverage----
p2 <- ggplot(dt1,
             aes(x = CpG)) +
  facet_wrap(~ reg,
             scale = "free_y") +
  geom_histogram(color = "black",
                 fill = "grey",
                 binwidth = 5) +
  scale_x_continuous(name = "Number of CpG",
                     breaks = c(3, 
                                seq(from = 5,
                                    to = 60,
                                    by = 5))) +
  coord_cartesian(xlim=c(3, 60)) +
  scale_y_continuous(name = "Counts") +
  ggtitle("Distribution of DMR by Number of CpG and Region")
p2

tiff(filename = "tmp/mes13_CpG_by_reg_hist.tiff",
     height = 8,
     width = 12,
     units = 'in',
     res = 300,
     compression = "lzw+p")
print(p2)
graphics.off()

# Percent methylation----
tmp <- as.matrix(dt1[, WJ01.N:WJ03.X])
head(tmp)

dtN <- tmp[, seq(1,
                 ncol(tmp) - 1, 
                 by = 2)]
head(dtN)

dtX <- tmp[, seq(2,
                 ncol(tmp), 
                 by = 2)]
head(dtX)

# Add 0.5 to all NAs and zeros in meth. hits
# NOTE: if there were no hits (N = NA or 0), the pct will be NA anyway
dtX <- apply(dtX,
             2,
             function(a) {
               a[is.na(a)] <- a[a == 0] <- 0.5
               return(a)
             })
head(dtX)

pct <- dtX/dtN
colnames(pct) <- substr(colnames(pct),
                        1,
                        nchar(colnames(pct)) - 2)

dt.pct <- data.table(dt1[, gene:CpG],
                  pct)
dt.pct

# Hits per CpG average (i.e. vertical coverage)----
t1 <- apply(dtN,
            2,
            function(a) {
              return(round(a/dt.pct$CpG,
                           1))
            })
t1

mu <- list()
for (i in 1:ncol(t1)) {
  x1 <- aggregate(x = t1[, i],
                  FUN = mean,
                  by = list(dt.pct$reg),
                  na.rm = TRUE)
  x2 <- aggregate(x = t1[, i],
                  FUN = mean,
                  by = list(dt1$reg),
                  na.rm = TRUE)
  x3 <- merge(x1, x2, by = "Group.1")
  mu[[i]] <- data.table(reg = x3[, 1],
                        mu = (x3[, 2] + x3[, 3])/2)
}
names(mu) <- unique(substr(colnames(t1),
                           1,
                           4))
mu

# Average methylation per region per treatment/time
mumth <- list()
for (i in 1:ncol(pct)) {
  x1 <- aggregate(x = c(pct[, i],
                        pct[, i]),
                  FUN = mean,
                  by = list(rep(dt1$reg, 2)),
                  na.rm = TRUE)
  
  x2 <- aggregate(x = c(pct[, i],
                        pct[, i]),
                  FUN = sd,
                  by = list(rep(dt1$reg, 2)),
                  na.rm = TRUE)
  mumth[[i]] <- data.table(rep(colnames(pct)[[i]],
                               5),
                           merge(x1, 
                                 x2,
                                 by = "Group.1"))
  
  colnames(mumth[[i]]) <- c("trt",
                            "reg",
                            "mu",
                            "std")
}
mumth <- do.call("rbind",
                 mumth)
mumth

mumth$trt <- factor(mumth$trt,
                    levels = c("WJ02",
                               "WJ01",
                               "WJ03"),
                    labels = c("HG",
                               "LG",
                               "MITC"))

mumth$`Methylation (%)` <- mumth$mu

p1 <- ggplot(mumth,
             aes(x = reg,
                 y = 100*`Methylation (%)`,
                 group = trt,
                 fill = trt)) +
  geom_bar(position = position_dodge(),
           stat="identity",
           color = "black") +
  scale_x_discrete("Region") +
  scale_y_continuous("Methylation (%)",
                     limits = c(0, 100)) +
  scale_fill_discrete("Treatment") +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 45,
                                   hjust = 1))
p1
tiff(filename = "tmp/mes13_avg_methyl_by_reg.tiff",
     height = 6,
     width = 7,
     units = 'in',
     res = 300,
     compression = "lzw+p")
print(p1 + ggtitle("Percent of Methylated CpG by Region"))
graphics.off()

# Gene level analysis----
trt.names <- c("LG",
               "HG",
               "MITC")
# Collapse by gene
out <- list()
for (i in 7:13) {
  out[[i - 6]] <- aggregate(dt1[, i, with = FALSE],
                            by = list(dt1$gene,
                                      dt1$reg),
                            FUN = sum,
                            na.rm = TRUE)
}

dt.gene <- data.table(Reduce(merge, out))
dt.gene 

# Calculate percent methylation in each sample----
dt.gene <- data.table(gene = dt.gene$Group.1,
                      reg = dt.gene$Group.2,
                      CpG = dt.gene$CpG,
                      pct = dt.gene[, c(5, 7, 9)]/
                                    dt.gene[, c(4, 6, 8)])
colnames(dt.gene)[4:6] <- trt.names
dt.gene

# Promoters only----
dt.pr <- droplevels(subset(dt.gene,
                           reg == "Promoter"))
dt.pr

# Ratios----
dt.pr$`diff(lg,hg)` <- dt.pr$LG - dt.pr$HG
dt.pr$`mean(lg,hg)` <- (dt.pr$LG + dt.pr$HG)/2

dt.pr$`diff(mitc,hg)` <- dt.pr$MITC - dt.pr$HG
dt.pr$`mean(mitc,hg)` <- (dt.pr$MITC + dt.pr$HG)/2

tiff(filename = "tmp/diff_vs_mean_lg_hg.tiff",
     height = 5,
     width = 5,
     units = 'in',
     res = 300,
     compression = "lzw+p")
plot(dt.pr$`diff(lg,hg)` ~ dt.pr$`mean(lg,hg)`,
     xlab = "Mean(LG,HG)",
     ylab = "Diff(LG,HG)",
     pch = ".")
abline(h = c(-0.1, 0.1),
       lty = 2,
       col = "red",
       lw = 2)
graphics.off()

tiff(filename = "tmp/diff_vs_mean_mitc_hg.tiff",
     height = 5,
     width = 5,
     units = 'in',
     res = 300,
     compression = "lzw+p")
plot(dt.pr$`diff(mitc,hg)` ~ dt.pr$`mean(mitc,hg)`,
     xlab = "Mean(MIC,HG)",
     ylab = "Diff(MIC,HG)",
     pch = ".")
abline(h = c(-0.1, 0.1),
       lty = 2,
       col = "red",
       lw = 2)
graphics.off()

# Rank genes by the LG vs. HG differences,
dt.sorted <- subset(dt.pr,
                    !is.nan(`diff(lg,hg)`))
dt.sorted <- dt.sorted[order(dt.sorted$`diff(lg,hg)`,
                             decreasing = TRUE), ]
write.csv(dt.sorted,
          file = "tmp/mes13-mict_prom_meth_by_gene.csv",
          row.names = FALSE)

# Top 20 positive differences----
tmp <- dt.sorted[1:20, ]
tmp

dtl <- melt.data.table(data = tmp,
                       id.vars = "gene",
                       measure.vars = 4:6)
dtl

p2 <- ggplot(data = dtl) +
  geom_tile(aes(x =  variable,
                y = gene,
                fill = value),
            color = "black") +
  scale_fill_gradient2(high = "red",
                       limit = c(0, 1),
                       name = "Methylation(%)") +
  scale_x_discrete("Treatment",
                   expand = c(0, 0)) +
  scale_y_discrete("Gene",
                   expand = c(0, 0)) +
  theme(axis.text.x = element_text(angle = 30,
                                   hjust = 1),
        plot.title = element_text(hjust = 0.5))
p2

tiff(filename = "tmp/top20_pos_lg.vs.hg_prom.tiff",
     height = 6,
     width = 6,
     units = 'in',
     res = 300,
     compression = "lzw+p")
print(p2 + ggtitle("MES13: Top 20 Genes With Largest\nPositive Differences in LG vs. HG in Promoter"))
graphics.off()

# Top 20 positive differences----
tmp <- dt.sorted[(nrow(dt.sorted) - 19):nrow(dt.sorted), ]
tmp

dtl <- melt.data.table(data = tmp,
                       id.vars = "gene",
                       measure.vars = 4:6)
dtl

p3 <- ggplot(data = dtl) +
  geom_tile(aes(x =  variable,
                y = gene,
                fill = value),
            color = "black") +
  scale_fill_gradient2(high = "red",
                       limit = c(0, 1),
                       name = "Methylation(%)") +
  scale_x_discrete("Treatment",
                   expand = c(0, 0)) +
  scale_y_discrete("Gene",
                   expand = c(0, 0)) +
  theme(axis.text.x = element_text(angle = 30,
                                   hjust = 1),
        plot.title = element_text(hjust = 0.5))
p3

tiff(filename = "tmp/top20_neg_lg.vs.hg_prom.tiff",
     height = 6,
     width = 6,
     units = 'in',
     res = 300,
     compression = "lzw+p")
print(p2 + ggtitle("MES13: Top 20 Genes With Largest\nNegative Differences in LG vs. HG in Promoter"))
graphics.off()

# Figure 5 of RO1----
tiff(filename = "tmp/figure5_ro1.tiff",
     height = 5,
     width = 12,
     units = 'in',
     res = 600,
     compression = "lzw+p")
gridExtra::grid.arrange(p2, p3, p1,
                        nrow = 1)
graphics.off()

# DNA vs. RNA----
pctMeth <- data.table(gene = dt1$gene,
                      anno = dt1$anno,
                      pct)
pctMeth$`HG-LG DNA` <- 100*(pctMeth$WJ02 - pctMeth$WJ01)
pctMeth$`MITC-HG DNA` <- 100*(pctMeth$WJ03 - pctMeth$WJ02)
pctMeth

# Load RNA DiffExp----
# NOTE: produced by mes13_rnaseq_DEGseq_MITC_v2.R script on 04/30/2018!
expRNA <- fread("data/mes13_mitc_genes_q-0.5_log2-0.3.csv")
colnames(expRNA)[1] <- "gene"

rna_dna <- merge(pctMeth,
                 expRNA,
                 by = "gene")
rna_dna

# Separate regions---
rna_dna$reg <- as.character(rna_dna$anno)
rna_dna$reg[substr(rna_dna$anno,
                   1,
                   8) == "Promoter"] <- "Promoter"
rna_dna$reg[substr(rna_dna$anno,
                   1,
                   4) %in% c("Exon",
                             "Intr")] <- "Body"
rna_dna$reg[substr(rna_dna$anno,
                   1,
                   4) %in% c("Dist",
                             "Down")] <- "Downstream"
rna_dna$reg <- factor(rna_dna$reg,
                      levels = c("Promoter",
                                 "5' UTR",
                                 "Body",
                                 "3' UTR",
                                 "Downstream"))
kable(data.table(table(dt1$reg)))
  # |V1         |     N|
  # |:----------|-----:|
  # |Promoter   | 80765|
  # |5' UTR     |   745|
  # |Body       | 63497|
  # |3' UTR     |  4273|
  # |Downstream | 61564|

g1 <- rna_dna[rna_dna$`HG-LG DNA` >= 10 &
                rna_dna$`HG-LG` <= -0.3 &
                reg == "Promoter"]
g1
length(unique(g1$gene))

g2 <- rna_dna[rna_dna$`HG-LG DNA` <= -10 &
                rna_dna$`HG-LG` >= 0.3 &
                reg == "Promoter"]
g2
length(unique(g2$gene))

g3 <- rna_dna[rna_dna$`MITC-HG DNA` >= 10 &
                rna_dna$`MITC-HG` <= -0.3 &
                reg == "Promoter"]
g3
length(unique(g3$gene))

g4 <- rna_dna[rna_dna$`MITC-HG DNA` <= -10 &
                rna_dna$`MITC-HG` >= 0.3 &
                reg == "Promoter"]
g4
length(unique(g4$gene))

write.csv(g1,
          file = "tmp/rna.up_dna.dn_hg-lg.csv",
          row.names = FALSE)
write.csv(g2,
          file = "tmp/rna.dn_dna.up_hg-lg.csv",
          row.names = FALSE)
write.csv(g3,
          file = "tmp/rna.up_dna.dn_mitc-hg.csv",
          row.names = FALSE)
write.csv(g4,
          file = "tmp/rna.dn_dna.up_mitc-lg.csv",
          row.names = FALSE)

p1 <- ggplot(data = rna_dna,
             aes(x = `HG-LG DNA`,
                 y = `HG-LG`,
                 fill = reg)) +
  geom_point(alpha = 0.7,
             size = 2,
             shape = 21) +
  geom_text(data = unique(rna_dna[gene %in% unique(g1$gene) &
                                    reg == "Promoter",
                                  c("gene",
                                    "reg",
                                    "HG-LG")]),
            aes(x = 40,
                y = `HG-LG`,
                label = gene),
            color = "blue",
            size = 2) +
  geom_text(data = unique(rna_dna[gene %in% unique(g2$gene) &
                                    reg == "Promoter",
                                  c("gene",
                                    "reg",
                                    "HG-LG")]),
            aes(x = -40,
                y = `HG-LG`,
                label = gene),
            color = "blue",
            size = 2) +
  geom_hline(yintercept = c(-0.3, 0.3),
             linetype = "dashed") +
  geom_vline(xintercept = c(-10, 10),
             linetype = "dashed") +
  scale_x_continuous("DNA Methylation Difference(%)",
                     breaks = seq(-30, 30, 10)) +
  scale_y_continuous("RNA Expression Difference (log2)") +
  ggtitle("HG - LG") +
  scale_fill_manual("Region",
                    values = c("Promoter" = "green",
                               "5' UTR" = "white",
                               "Body" = "blue",
                               "3' UTR" = "grey",
                               "Downstream" = "red")) +
  theme(plot.title = element_text(hjust = 0.5))
p1
tiff(filename = "tmp/starburst_hg-lg.tiff",
     height = 10,
     width = 10,
     units = 'in',
     res = 300,
     compression = "lzw+p")
print(p1)
graphics.off()

p2 <- ggplot(data = rna_dna,
             aes(x = `MITC-HG DNA`,
                 y = `MITC-HG`,
                 fill = reg)) +
  geom_point(alpha = 0.7,
             size = 2,
             shape = 21) +
  geom_text(data = unique(rna_dna[gene %in% unique(g3$gene) &
                                    reg == "Promoter",
                                  c("gene",
                                    "reg",
                                    "MITC-HG")]),
            aes(x = 40,
                y = `MITC-HG`,
                label = gene),
            color = "blue",
            size = 2) +
  geom_text(data = unique(rna_dna[gene %in% unique(g4$gene) &
                                    reg == "Promoter",
                                  c("gene",
                                    "reg",
                                    "MITC-HG")]),
            aes(x = -40,
                y = `MITC-HG`,
                label = gene),
            color = "blue",
            size = 2) +
  geom_hline(yintercept = c(-0.3, 0.3),
             linetype = "dashed") +
  geom_vline(xintercept = c(-10, 10),
             linetype = "dashed") +
  scale_x_continuous("DNA Methylation Difference(%)",
                     breaks = seq(-30, 30, 10)) +
  scale_y_continuous("RNA Expression Difference (log2)") +
  ggtitle("MITC-HG") +
  scale_fill_manual("Region",
                    values = c("Promoter" = "green",
                               "5' UTR" = "white",
                               "Body" = "blue",
                               "3' UTR" = "grey",
                               "Downstream" = "red")) +
  theme(plot.title = element_text(hjust = 0.5))
p2
tiff(filename = "tmp/starburst_mitc-hg.tiff",
     height = 10,
     width = 10,
     units = 'in',
     res = 300,
     compression = "lzw+p")
print(p2)
graphics.off()

# sessionInfo()
# sink()
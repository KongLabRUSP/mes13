# |----------------------------------------------------------------------------------|
# | Project: Skin UVB SKH1 mouse model treated with UA/SFN                           |
# | Script: Methyl-seq data analysis and visualization using DSS                     |
# | Coordinator: Ran Yin, Renyi Wu                                                   |
# | Author: Davit Sargsyan                                                           |
# | Created: 03/17/2018                                                              |
# | Modified:04/05/2018, DS: changed hitmaps to donut plots; added more comparisons  |
# |          09/04/2018, DS: TIIA methyl-seq sample is WJ5, not WJ4 as previously    |
# |----------------------------------------------------------------------------------|
# sink(file = "tmp/log_skin_uvb_dna_v2.txt")
date()

# NOTE: several packages, e.g. Rcpp, MASS, etc., might be deleted manually and reinstalled
# Workflow: https://www.bioconductor.org/help/workflows/rnaseqGene/
# source("https://bioconductor.org/biocLite.R")
# biocLite("TxDb.Mmusculus.UCSC.mm10.knownGene")
# biocLite("ChIPseeker")
# biocLite("DO.db")
# biocLite("GenomicRanges")
# biocLite("org.Mm.eg.db")
# biocLite("DSS")
# biocLite("bsseq")

require(data.table)
require(ggplot2)
require(knitr)

require(ChIPseeker)
require(TxDb.Mmusculus.UCSC.mm10.knownGene)
require(DSS)

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

# Subset data: LG, HG and TIIA----
dt1 <- data.table(gene = dt1$SYMBOL,
                  anno = dt1$annotation,
                  geneId = dt1$geneId,
                  chr = dt1$geneChr,
                  pos = dt1$start,
                  distanceToTSS = dt1$distanceToTSS,
                  reg = NA,
                  dt1[, CpG:WJ02.X],
                  dt1[, WJ05.N:WJ05.X],
                  geneName = dt1$GENENAME)
dt1

# # Dispersion Shrinkage for Sequencing data (DSS)----
# # This is based on Wald test for beta-binomial distribution.
# # Source: https://www.bioconductor.org/packages/release/bioc/vignettes/DSS/inst/doc/DSS.pdf
# # The DM detection procedure implemented in DSS is based on a rigorous Wald test for betabinomial
# # distributions. The test statistics depend on the biological variations (characterized
# # by dispersion parameter) as well as the sequencing depth. An important part of the algorithm
# # is the estimation of dispersion parameter, which is achieved through a shrinkage estimator
# # based on a Bayesian hierarchical model [1]. An advantage of DSS is that the test can be
# # performed even when there is no biological replicates. That’s because by smoothing, the
# # neighboring CpG sites can be viewed as “pseudo-replicates", and the dispersion can still be
# # estimated with reasonable precision.

# Regions----
kable(data.table(table(substr(dt1$anno, 1, 9))))
  # |V1        |     N|
  # |:---------|-----:|
  # |3' UTR    |  4396|
  # |5' UTR    |   758|
  # |Distal In | 61324|
  # |Downstrea |  2683|
  # |Exon (uc0 | 11546|
  # |Intron (u | 54267|
  # |Promoter  | 82137|

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
  # |Promoter   | 82137|
  # |5' UTR     |   758|
  # |Body       | 65813|
  # |3' UTR     |  4396|
  # |Downstream | 64007|


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
tmp <- as.matrix(dt1[, WJ01.N:WJ05.X])
head(tmp)

# Remove rows with all NAs
ndx.keep <- rowSums(is.na(tmp)) < 6
sum(ndx.keep)
# 211,940 out of 217,111

dt1 <- dt1[ndx.keep, ]
tmp <- tmp[ndx.keep, ]

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
dtX

pct <- dtX/dtN
colnames(pct) <- substr(colnames(pct),
                        1,
                        nchar(colnames(pct)) - 2)
head(pct)

# Remove rows with all zeros----
dim(pct[rowSums(pct) == 0, ])
dim(pct[is.na(rowSums(pct)), ])
dim(pct)

ndx.keep <- rowSums(pct) != 0 & !is.na(rowSums(pct))
pct <- pct[ndx.keep, ]
dt1 <- dt1[ndx.keep, ]
dtN <- dtN[ndx.keep, ]
dtX <- dtX[ndx.keep, ]
dim(dtX)
# 187,879  remaine

# Hits per CpG average (i.e. vertical coverage)----
t1 <- apply(dtN,
            2,
            function(a) {
              return(round(a/dt1$CpG,
                           1))
            })

mu <- list()
for (i in 1:ncol(t1)) {
  x1 <- aggregate(x = t1[, i],
                  FUN = mean,
                  by = list(dt1$reg))
  x2 <- aggregate(x = t1[, i],
                  FUN = mean,
                  by = list(dt1$reg))
  x3 <- merge(x1, x2, by = "Group.1")
  mu[[i]] <- data.table(reg = x3[, 1],
                        mu = (x3[, 2] + x3[, 3])/2)
}
names(mu) <- unique(substr(colnames(t1),
                           1,
                           4))
mu
# more than 10

# Average methylation per region per treatment/time
mumth <- list()
for (i in 1:ncol(pct)) {
  x1 <- aggregate(x = c(pct[, i],
                        pct[, i]),
                  FUN = mean,
                  by = list(rep(dt1$reg, 2)))
  
  x2 <- aggregate(x = c(pct[, i],
                        pct[, i]),
                  FUN = sd,
                  by = list(rep(dt1$reg, 2)))
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
                               "WJ05"),
                    labels = c("HG",
                               "LG",
                               "TIIA"))

mumth$`Methylation (%)` <- 100*mumth$mu

p1 <- ggplot(mumth,
             aes(x = reg,
                 y = `Methylation (%)`,
                 group = trt,
                 fill = trt)) +
  geom_bar(position = position_dodge(),
           stat="identity",
           color = "black") +
  scale_x_discrete("Region") +
  scale_y_continuous(limits = c(0, 100)) +
  ggtitle("Percent of Methylated CpG by Region") +
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
print(p1)
graphics.off()

p2 <- ggplot(mumth,
             aes(x = reg,
                 y = mu,
                 group = trt,
                 fill = trt)) +
  geom_errorbar(aes(ymax = mu + std,
                    ymin = mu),
                width = 0.5,
                position = position_dodge(0.9)) +
  geom_point(position = position_dodge(0.9),
             size = 1) +
  geom_bar(position = position_dodge(0.9),
           stat="identity",
           color = "black") +
  scale_x_discrete("Region") +
  scale_y_continuous() +
  ggtitle("Proportion of Methylated CpG by Region") +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 45,
                                   hjust = 1))
p2
tiff(filename = "tmp/mes13_avg_sd_methyl_by_reg.tiff",
     height = 6,
     width = 7,
     units = 'in',
     res = 300,
     compression = "lzw+p")
print(p2)
graphics.off()

# DNA vs. RNA----
pctMeth <- data.table(gene = dt1$gene,
                      anno = dt1$anno,
                      pct)
pctMeth$`HG-LG DNA` <- 100*(pctMeth$WJ02 - pctMeth$WJ01)
pctMeth$`TIIA-HG DNA` <- 100*(pctMeth$WJ05 - pctMeth$WJ02)
pctMeth

# Load RNA DiffExp----
# NOTE: produced by mes13_rnaseq_DEGseq_TIIA_v2.R script on 04/30/2018!
# expRNA <- fread("data/mes13_tiia_genes_q-0.5_log2-0.3.csv")
expRNA <- fread("data/mes13_tiia_rnaseq_degseq_genes_q-0.5_log2-0.3.csv")
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
  # |Promoter   | 75408|
  # |5' UTR     |   695|
  # |Body       | 55203|
  # |3' UTR     |  3808|
  # |Downstream | 52765|

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

g3 <- rna_dna[rna_dna$`TIIA-HG DNA` >= 10 & 
                rna_dna$`TIIA-HG` <= -0.3 &
                reg == "Promoter"]
g3
length(unique(g3$gene))

g4 <- rna_dna[rna_dna$`TIIA-HG DNA` <= -10 & 
                rna_dna$`TIIA-HG` >= 0.3 &
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
          file = "tmp/rna.up_dna.dn_tiia-hg.csv",
          row.names = FALSE)
write.csv(g4,
          file = "tmp/rna.dn_dna.up_tiia-lg.csv",
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
             aes(x = `TIIA-HG DNA`,
                 y = `TIIA-HG`,
                 fill = reg)) +
  geom_point(alpha = 0.7,
             size = 2,
             shape = 21) +
  geom_text(data = unique(rna_dna[gene %in% unique(g3$gene) &
                                    reg == "Promoter", 
                                  c("gene",
                                    "reg",
                                    "TIIA-HG")]),
            aes(x = 40,
                y = `TIIA-HG`,
                label = gene),
            color = "blue",
            size = 2) +
  geom_text(data = unique(rna_dna[gene %in% unique(g4$gene) &
                                    reg == "Promoter", 
                                  c("gene",
                                    "reg",
                                    "TIIA-HG")]),
            aes(x = -40,
                y = `TIIA-HG`,
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
  ggtitle("TIIA-HG") +
  scale_fill_manual("Region",
                    values = c("Promoter" = "green",
                               "5' UTR" = "white",
                               "Body" = "blue",
                               "3' UTR" = "grey",
                               "Downstream" = "red")) +
  theme(plot.title = element_text(hjust = 0.5))
p2
tiff(filename = "tmp/starburst_tiia-hg.tiff",
     height = 10,
     width = 10,
     units = 'in',
     res = 300,
     compression = "lzw+p")
print(p2)
graphics.off()

# Detailes on genes with significant expression differences----
# Genes with RNA down/DNA up----
gene.keep <- unique(g3$gene)
gene.keep

gene.keep <- c("Fgl2",
               "Gulo/(GLO)",
               "Kcnip2/KChIP2",
               "Nmu",
               "Bmp8b")

gene.keep <- c("Fgl2",
               "Gulo",
               "GLO",
               "Kcnip2",
               "KChIP2",
               "Nmu",
               "Bmp8b")


gene.keep

# DNA----
dna <- data.table(gene = dt1$gene,
                  CpG = dt1$CpG,
                  annotation = as.character(dt1$anno),
                  distanceToTSS = dt1$distanceToTSS,
                  pct)

dna <- dna[gene %in% gene.keep, ]
dna
class(dna$distanceToTSS)

# Differences
dna$`HG-LG` <- 100*(dna$WJ02 - dna$WJ01)
dna$`TIIA-HG` <- 100*(dna$WJ05 - dna$WJ02)

dna$reg <- "5 to 10"
dna$reg[dna$CpG > 10] <- "11 to 20"
dna$reg[dna$CpG > 20] <- ">20"
dna$reg <- factor(dna$reg,
                  levels = c("5 to 10",
                             "11 to 20",
                             ">20"))

dna

setkey(dna,
        gene,
        distanceToTSS)
dna[, distRank := rev(rank(distanceToTSS)),
    by = gene]

# Long data
dt3 <- melt.data.table(data = dna,
                       id.vars = c("gene",
                                   "reg",
                                   "annotation",
                                   "distanceToTSS",
                                   "distRank"),
                       measure.vars = c("HG-LG",
                                        "TIIA-HG"),
                       variable.name = "Treatment",
                       value.name = "DNA")
dt3$Treatment <- as.character(dt3$Treatment)

dt3$annotation[substr(dt3$annotation, 1, 4) == "Exon"] <- "Exon"
dt3$annotation[substr(dt3$annotation, 1, 6) == "Intron"] <- "Intron"
dt3$annotation[substr(dt3$annotation, 1, 8) == "Promoter"] <- "Promoter"
dt3$annotation[substr(dt3$annotation, 1, 4) == "Down"] <- "Downstream"
dt3$annotation <- factor(dt3$annotation)

# RNA
# RNA data Long format----
dt2 <- melt.data.table(data = expRNA,
                       id.vars = "gene",
                       measure.vars = c("HG-LG",
                                        "TIIA-HG"),
                       variable.name = "Treatment",
                       value.name = "RNA")
dt2$Treatment <- as.character(dt2$Treatment)
dt2

# Merge DNA with RNA----
dt3 <- merge(dt3,
             dt2,
             by = c("gene",
                    "Treatment"))
dt3

# Isolate genes----
for (i in 1:length(unique(dna$gene))) {
  gX <- unique(dt3$gene)[i]
  dna.gX <- dt3[dt3$gene %in% gX, ]
  dna.gX$y0 <- 0
  
  dna.gX$Treatment <- paste(dna.gX$Treatment,
                            " (RNA = ",
                            round(dna.gX$RNA, 3),
                            ")",
                            sep = "")
  
  p1 <- ggplot(dna.gX,
               aes(x = distRank,
                   y = DNA)) +
    facet_wrap(.~ Treatment,
               scales = "free_y",
               ncol = 1) +
    geom_rect(aes(xmin = -Inf,
                  xmax = Inf,
                  ymin = -Inf,
                  ymax = -10),
              fill = "red",
              alpha = 0.1) +
    geom_rect(aes(xmin = -Inf,
                  xmax = Inf,
                  ymin = 10,
                  ymax = Inf),
              fill = "green",
              alpha = 0.1) +
    geom_hline(yintercept = 0) +
    geom_hline(yintercept = c(10,
                              -10),
               linetype = "dashed") +
    geom_segment(aes(x = distRank,
                     y = y0,
                     xend = distRank,
                     yend = DNA)) + 
    geom_point(aes(x = distRank,
                   y = DNA,
                   fill = annotation,
                   size = reg),
               shape = 21) +
    # geom_rect(aes(xmin = -Inf,
    #               xmax = Inf,
    #               ymin = -10,
    #               ymax = 10),
    #           fill = "white",
    #           alpha = 0.1) +
    ggtitle(paste("Gene:",
                  gX)) +
    scale_x_continuous("Distance from TSS",
                       breaks = dna.gX$distRank,
                       labels = dna.gX$distanceToTSS) +
    scale_y_continuous("% Methylation") +
    scale_fill_manual("Region",
                      values = c("Distal Intergenic" = "purple",
                                 "Exon" = "blue",
                                 "Intron" = "white",
                                 "Promoter" = "brown",
                                 "3' UTR" = "black",
                                 "5' UTR" = "yellow",
                                 "Downstream" = "orange")) +
    scale_size_manual("Number of CpG-s",
                      values = c("5 to 10" = 5,
                                 "11 to 20" = 6,
                                 ">20" = 7)) +
    guides(fill = guide_legend(override.aes = list(size = 7))) +
    theme(plot.title = element_text(hjust = 0.5),
          #legend.position = "top",
          axis.text.x = element_text(angle = 45,
                                     hjust = 1))
  p1
  tiff(filename = paste("tmp/",
                        gX,
                        ".tiff",
                        sep = ""),
       height = 8,
       width = 8,
       units = 'in',
       res = 300,
       compression = "lzw+p")
  print(p1)
  graphics.off()
}

# sessionInfo()
# sink()
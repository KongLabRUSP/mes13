# |----------------------------------------------------------------------------------|
# | Project: Study of Diabetes in MES13 cells                                        |
# | Script: Methyl-seq data analysis and visualization                               |
# | Author: Davit Sargsyan                                                           |
# | Created: 09/12/2017                                                              |
# | Created: 01/15/2018                                                              |
# |----------------------------------------------------------------------------------|
# Header----
# Save consol output to a log file
# sink(file = "tmp/log_mes13_methylseq_data_analysis_v4.txt")

# Source: file:///C:/R/R-3.3.2/library/ChIPseeker/doc/ChIPseeker.html
# Source: http://bioconductor.org/packages/release/data/annotation/html/org.Mm.eg.db.html

# # NOTE: FOR LINUX, RUN THIS IN THE TERMINAL AS SUDO!
# # Source: https://support.bioconductor.org/p/70093/
# sudo R
# source("http://bioconductor.org/biocLite.R")
# biocLite("TxDb.Mmusculus.UCSC.mm9.knownGene",
#          suppressUpdates = TRUE)
# biocLite("org.Mm.eg.db")
# biocLite("ReactomePA")

require(data.table)
require(ggplot2)
# require(ChIPseeker)
# require(TxDb.Mmusculus.UCSC.mm9.knownGene)
require(knitr)
# require(ReactomePA)

# Treatment legend----
trt.names <- c("LG",
               "HG",
               "MIC 1.5uM",
               "FX 1uM",
               "Ber 6uM")

# Load data----
# Question to Renyi: how was the data processed and annotated?
# If by John's DMR Finder, what were the settings?
dt1 <- fread("data/Renyi_Methylseq_12292017/combined_WJ_anno.csv")
dt1[dt1$gene == unique(dt1$gene)[2], ]

# Remove unmapped regions
dt1 <- dt1[!is.na(dt1$gene == "NA"), ]
# None missing
length(unique(dt1$gene))
# 21,353 genes
summary(dt1)
dt1

# Subset data----
dt1 <- data.table(gene = dt1$SYMBOL,
                  anno = dt1$annotation,
                  geneId = dt1$geneId,
                  reg = NA,
                  CpG = dt1$CpG,
                  dt1[, 2:11], 
                  geneName = dt1$GENENAME)
kable(data.table(table(substr(dt1$feature,
                              1, 
                              9))))
  # |V1        |      N|
  # |:---------|------:|
  # |3' UTR    |  3,751|
  # |5' UTR    |    465|
  # |Distal In | 66,373|
  # |Downstrea |  2,827|
  # |Exon (NM_ | 11,662|
  # |Exon (NR_ |    539|
  # |Intron (N | 54,519|
  # |Promoter  | 77,041|

# Separate Promoter, Body and Downstream; remove everything else
dt1$reg <- dt1$feature
# a. Promoter: up to 3kb upstream
dt1$reg[substr(dt1$feature, 
               1,
               8) == "Promoter"] <- "Promoter"

# b. Body: exons and introns
dt1$reg[substr(dt1$feature, 
               1, 
               4) %in% c("Exon",
                         "Intr")] <- "Body"

# c. Downstream: Distal Intergenic and  Downstream
dt1$reg[substr(dt1$feature, 
               1, 
               4) %in% c("Dist",
                         "Down")] <- "Downstream"
dt1$reg <- factor(dt1$reg,
                  levels = c("5' UTR",
                             "Promoter",
                             "Body",
                             "Downstream",
                             "3' UTR"))

kable(data.table(table(dt1$reg)))
  # |V1         |      N|
  # |:----------|------:|
  # |5' UTR     |    465|
  # |Promoter   | 77,041|
  # |Body       | 66,720|
  # |Downstream | 69,200|
  # |3' UTR     |  3,751|

# Part I: total methylation----
# Aggregate data by region
out <- list()
for (i in 4:14) {
  out[[i - 3]] <- aggregate(dt1[, i, with = FALSE],
                            by = list(dt1$reg),
                            FUN = sum,
                            na.rm = TRUE)
}
dt.reg <- data.table(Reduce(merge, out))

# Change row order
dt.reg <- dt.reg[c(2, 5, 3, 4, 1), ]
dt.reg

# Hits per CpG average (i.e. vertical coverage)
t1 <- dt.reg[, c(1, 2, 3, 5, 7, 9, 11)]
t1[, 3:7] <- do.call("cbind",
                     lapply(t1[, 3:7],
                            function(a) {
                              return(round(a/t1[, 2],
                                           1))
                            }))
colnames(t1) <- c("Gene Region",
                  "Total CpG Count",
                  trt.names)
kable(t1)
  # |Gene Region | Total CpG Count|   LG|   HG| MIC 1.5uM| FX 1uM| Ber 6uM|
  # |:-----------|---------------:|----:|----:|---------:|------:|-------:|
  # |5' UTR      |            6423|  9.0|  9.6|      11.1|   11.0|    11.6|
  # |Promoter    |         1229957|  8.8|  9.4|      10.8|   10.8|    11.5|
  # |Body        |          472884| 10.5| 11.2|      13.0|   13.1|    13.5|
  # |Downstream  |          489164| 10.8| 11.5|      13.3|   13.5|    13.9|
  # |3' UTR      |           34952|  9.4| 10.1|      11.8|   11.7|    12.2|

# Compare to the previous data ("bad" samples):
  # |Gene Region | Total CpG Count|  LG|  HG| MIC 1.5uM| FX 1uM| Ber 6uM|
  # |:-----------|---------------:|---:|---:|---------:|------:|-------:|
  # |Promoter    |           22856| 0.6| 1.3|       1.2|    0.9|     1.2|
  # |Body        |            5153| 1.0| 1.9|       1.6|    1.5|     1.6|
  # |Downstream  |            6458| 2.2| 2.9|       3.0|    4.1|     3.2|
write.csv(t1,
          file = "tmp/t1.csv",
          row.names = FALSE)

# Calculate percent methylation in each sample----
dt.reg <- data.table(Region = dt.reg$Group.1,
                     CpG = dt.reg$CpG,
                     pct = round(100*dt.reg[, c(4, 6, 8, 10, 12)]/
                                   dt.reg[, c(3, 5, 7, 9, 11)],
                                 1))
dt.reg

# Melt data
dt.reg.l <- melt.data.table(data = dt.reg,
                            id.vars = 1:2,
                            measure.vars = 3:7,
                            variable.name = "Treatment",
                            value.name = "Methylation(%)")

dt.reg.l$Treatment <- factor(dt.reg.l$Treatment)
levels(dt.reg.l$Treatment) <- trt.names
dt.reg.l

# Plot
p1 <- ggplot(dt.reg.l,
             aes(x = Region,
                 y = `Methylation(%)`,
                 group = Treatment,
                 fill = Treatment)) +
  geom_bar(position = position_dodge(),
           stat="identity",
           color = "black") +
  scale_x_discrete("Region") +
  scale_y_continuous(limits = c(0, 100)) +
  ggtitle("Total Methylation (%)") +
  theme(plot.title = element_text(hjust = 0.5))
p1
tiff(filename = "tmp/global_methyl_by_region.tiff",
     height = 5,
     width = 5,
     units = 'in',
     res = 300,
     compression = "lzw+p")
print(p1)
graphics.off()

# Part II: gene level analysis----
# Collapse by gene
out <- list()
for (i in 4:14) {
  out[[i - 3]] <- aggregate(dt1[, i, with = FALSE],
                            by = list(dt1$gene,
                                      dt1$reg),
                            FUN = sum,
                            na.rm = TRUE)
}

dt.gene <- data.table(Reduce(merge, out))

# Calculate percent methylation in each sample----
dt.gene <- data.table(gene = dt.gene$Group.1,
                      reg = dt.gene$Group.2,
                      CpG = dt.gene$CpG,
                      pct = round(100*dt.gene[, c(5, 7, 9, 11, 13)]/
                                    dt.gene[, c(4, 6, 8, 10, 12)],
                                  1))
colnames(dt.gene)[4:8] <- trt.names
dt.gene

# Promoters only----
dt.pr <- droplevels(subset(dt.gene,
                           reg == "Promoter"))
dt.pr

hist(dt.pr$LG, 100)
hist(log2(dt.pr$LG), 100)
hist(dt.pr$HG, 100)

tiff(filename = "tmp/prom_methyl_by_trt.tiff",
     height = 10,
     width = 10,
     units = 'in',
     res = 300,
     compression = "lzw+p")
plot(dt.pr[, 4:8])
graphics.off()

# Ratios----
dt.pr$`diff(lg,hg)` <- dt.pr$LG - dt.pr$HG
dt.pr$`mean(lg,hg)` <- (dt.pr$LG + dt.pr$HG)/2

dt.pr$`diff(mic,hg)` <- dt.pr$`MIC 1.5uM` - dt.pr$HG
dt.pr$`mean(mic,hg)` <- (dt.pr$`MIC 1.5uM` + dt.pr$HG)/2

dt.pr$`diff(fx,hg)` <- dt.pr$`FX 1uM` - dt.pr$HG
dt.pr$`mean(fx,hg)` <- (dt.pr$`FX 1uM` + dt.pr$HG)/2

dt.pr$`diff(ber,hg)` <- dt.pr$`Ber 6uM` - dt.pr$HG
dt.pr$`mean(ber,hg)` <- (dt.pr$`Ber 6uM` + dt.pr$HG)/2

tiff(filename = "tmp/diff_vs_mean_lg_hg.tiff",
     height = 5,
     width = 5,
     units = 'in',
     res = 300,
     compression = "lzw+p")
plot(dt.pr$`diff(lg,hg)` ~ dt.pr$`mean(lg,hg)`,
     xlab = "Mean(LG,HG)",
     ylab = "Diff(LG,HG)")
abline(h = c(-20, 20),
       lty = 2)
graphics.off()

tiff(filename = "tmp/diff_vs_mean_mic_hg.tiff",
     height = 5,
     width = 5,
     units = 'in',
     res = 300,
     compression = "lzw+p")
plot(dt.pr$`diff(mic,hg)` ~ dt.pr$`mean(mic,hg)`,
     xlab = "Mean(MIC,HG)",
     ylab = "Diff(MIC,HG)")
abline(h = c(-20, 20),
       lty = 2)
graphics.off()

tiff(filename = "tmp/diff_vs_mean_fx_hg.tiff",
     height = 5,
     width = 5,
     units = 'in',
     res = 300,
     compression = "lzw+p")
plot(dt.pr$`diff(fx,hg)` ~ dt.pr$`mean(fx,hg)`,
     xlab = "Mean(FX,HG)",
     ylab = "Diff(FX,HG)")
abline(h = c(-20, 20),
       lty = 2)
graphics.off()

tiff(filename = "tmp/diff_vs_mean_ber_hg.tiff",
     height = 5,
     width = 5,
     units = 'in',
     res = 300,
     compression = "lzw+p")
plot(dt.pr$`diff(ber,hg)` ~ dt.pr$`mean(ber,hg)`,
     xlab = "Mean(Ber,HG)",
     ylab = "Diff(Ber,HG)")
abline(h = c(-20, 20),
       lty = 2)
graphics.off()

# Scale and center
dt.pr$scale.lg.hg <- scale(dt.pr$`diff(lg,hg)`)

hist(dt.pr$`diff(lg,hg)`, 100)
hist(dt.pr$scale.lg.hg, 100)
hist(abs(dt.pr$scale.lg.hg), 100)

# Rank genes by the absolute differences,
# plot top20 genes
tmp <- subset(dt.pr,
              !is.nan(scale.lg.hg))
tmp <- tmp[order(abs(tmp$scale.lg.hg),
                 decreasing = TRUE), ]
tmp <- tmp[1:20, ]
tmp

dtl <- melt.data.table(data = tmp,
                       id.vars = "gene",
                       measure.vars = 4:8)
dtl

p2 <- ggplot(data = dtl) +
  geom_tile(aes(x =  variable,
                y = gene,
                fill = value),
            color = "black") +
  scale_fill_gradient2(high = "red",
                       limit = c(0, 100),
                       name = "Methylation(%)") +
  scale_x_discrete("Treatment",
                   expand = c(0, 0)) +
  scale_y_discrete("Gene",
                   expand = c(0, 0)) +
  ggtitle("MES13: Top 20 Genes\nWith Largest Differences in LG vs. HG Promoters") +
  theme(axis.text.x = element_text(angle = 30,
                                   hjust = 1),
        plot.title = element_text(hjust = 0.5))
p2

tiff(filename = "tmp/top20_lg.vs.hg_all.tiff",
     height = 6,
     width = 6,
     units = 'in',
     res = 300,
     compression = "lzw+p")
print(p2)
graphics.off()

# Differences
tmp[, c(9, 11, 13, 15)]
dtl <- melt.data.table(data = tmp,
                       id.vars = "gene",
                       measure.vars = c(9, 11, 13, 15))
dtl
summary(dtl)

p3 <- ggplot(data = dtl) +
  geom_tile(aes(x =  variable,
                y = gene,
                fill = value),
            color = "black") +
  scale_fill_gradient2(high = "green",
                       mid = "black",
                       low = "red",
                       #limit = c(0, 100),
                       name = "Methylation Diff(%)") +
  scale_x_discrete("Treatment",
                   expand = c(0, 0)) +
  scale_y_discrete("Gene",
                   expand = c(0, 0)) +
  ggtitle("MES13: Top 20 Genes\nWith Largest Differences in LG vs. HG Promoters") +
  theme(axis.text.x = element_text(angle = 30,
                                   hjust = 1),
        plot.title = element_text(hjust = 0.5))
p3

tiff(filename = "tmp/top20_lg.vs.hg_diff.tiff",
     height = 6,
     width = 6,
     units = 'in',
     res = 300,
     compression = "lzw+p")
print(p3)
graphics.off()

# Largest change between LG vs. HG and MIC vs. HG and in opposite direction----
dt.pr$lg.hg_mic.hg <- dt.pr$`diff(lg,hg)` - dt.pr$`diff(mic,hg)`
tmp <- subset(dt.pr,
              !is.nan(lg.hg_mic.hg))
tmp <- tmp[order(abs(tmp$lg.hg_mic.hg),
                 decreasing = TRUE), ]
tmp <- tmp[abs(tmp$`diff(lg,hg)`) > 5 & 
             abs(tmp$`diff(mic,hg)`) > 5 &
             tmp$`diff(lg,hg)`*tmp$`diff(mic,hg)` < 0, ]
tmp <- tmp[1:40, ]

# Order by LG vs. HG
tmp$gene <- factor(tmp$gene,
                   levels = tmp$gene[order(tmp$`diff(lg,hg)`)])

# Differences
tmp[, c(9, 11)]
dtl <- melt.data.table(data = tmp,
                       id.vars = "gene",
                       measure.vars = c(9, 11))
dtl
summary(dtl)

p4 <- ggplot(data = dtl) +
  geom_tile(aes(x =  variable,
                y = gene,
                fill = value),
            color = "black") +
  scale_fill_gradient2(high = "green",
                       mid = "black",
                       low = "red",
                       #limit = c(0, 100),
                       name = "Methylation Diff(%)") +
  scale_x_discrete("Treatment",
                   expand = c(0, 0)) +
  scale_y_discrete("Gene",
                   expand = c(0, 0)) +
  ggtitle("MES13: Top 40 Genes With Largest Differences in\n(HG vs. LG) and (MIC vs. HG) Promoters") +
  theme(axis.text.x = element_text(angle = 30,
                                   hjust = 1),
        plot.title = element_text(hjust = 0.5))
p4

tiff(filename = "tmp/top40_lg.vs.hg_mic.vs.hg.tiff",
     height = 8,
     width = 5,
     units = 'in',
     res = 300,
     compression = "lzw+p")
print(p4)
graphics.off()

# Largest change between HG vs. LG and FX vs. HG----
dt.pr$hg.lg_fx.hg <- dt.pr$`diff(hg,lg)` - dt.pr$`diff(fx,hg)`
tmp <- subset(dt.pr,
              !is.nan(hg.lg_fx.hg))
tmp <- tmp[order(abs(tmp$hg.lg_fx.hg),
                 decreasing = TRUE), ]
tmp <- tmp[1:40, ]
# Order by HG vs. LG
tmp$gene <- factor(tmp$gene,
                   levels = tmp$gene[order(tmp$`diff(hg,lg)`)])

# Differences
tmp[, c(9, 13)]
dtl <- melt.data.table(data = tmp,
                       id.vars = "gene",
                       measure.vars = c(9, 13))
dtl
summary(dtl)

p5 <- ggplot(data = dtl) +
  geom_tile(aes(x =  variable,
                y = gene,
                fill = value),
            color = "black") +
  scale_fill_gradient2(high = "green",
                       mid = "black",
                       low = "red",
                       #limit = c(0, 100),
                       name = "Methylation Diff(%)") +
  scale_x_discrete("Treatment",
                   expand = c(0, 0)) +
  scale_y_discrete("Gene",
                   expand = c(0, 0)) +
  ggtitle("MES13: Top 40 Genes With Largest Differences in\n(HG vs. LG) and (FX vs. HG) Promoters") +
  theme(axis.text.x = element_text(angle = 30,
                                   hjust = 1),
        plot.title = element_text(hjust = 0.5))
p5

tiff(filename = "tmp/top40_hg.vs.lg_fx.vs.hg.tiff",
     height = 8,
     width = 5,
     units = 'in',
     res = 300,
     compression = "lzw+p")
print(p5)
graphics.off()

# Largest change between HG vs. LG and Ber vs. HG----
dt.pr$hg.lg_ber.hg <- dt.pr$`diff(hg,lg)` - dt.pr$`diff(ber,hg)`
tmp <- subset(dt.pr,
              !is.nan(hg.lg_ber.hg))
tmp <- tmp[order(abs(tmp$hg.lg_ber.hg),
                 decreasing = TRUE), ]
tmp <- tmp[1:40, ]
# Order by HG vs. LG
tmp$gene <- factor(tmp$gene,
                   levels = tmp$gene[order(tmp$`diff(hg,lg)`)])

# Differences
tmp[, c(9, 15)]
dtl <- melt.data.table(data = tmp,
                       id.vars = "gene",
                       measure.vars = c(9, 15))
dtl
summary(dtl)

p6 <- ggplot(data = dtl) +
  geom_tile(aes(x =  variable,
                y = gene,
                fill = value),
            color = "black") +
  scale_fill_gradient2(high = "green",
                       mid = "black",
                       low = "red",
                       #limit = c(0, 100),
                       name = "Methylation Diff(%)") +
  scale_x_discrete("Treatment",
                   expand = c(0, 0)) +
  scale_y_discrete("Gene",
                   expand = c(0, 0)) +
  ggtitle("MES13: Top 40 Genes With Largest Differences in\n(HG vs. LG) and (Ber vs. HG) Promoters") +
  theme(axis.text.x = element_text(angle = 30,
                                   hjust = 1),
        plot.title = element_text(hjust = 0.5))
p6

tiff(filename = "tmp/top40_hg.vs.lg_ber.vs.hg.tiff",
     height = 8,
     width = 5,
     units = 'in',
     res = 300,
     compression = "lzw+p")
print(p6)
graphics.off()

# Save the dataset----
write.csv(dt.pr,
          file = "tmp/mes13_promoters_gene_level.csv",
          row.names = FALSE)

# Genes that changed significantly (>15%)----
# Down-up regulated
g1 <- dt.pr$gene[dt.pr$`diff(hg,lg)` > 15 &
                   dt.pr$`diff(mic,hg)` < -15 &
                   !is.na(dt.pr$hg.lg_mic.hg)]
g1

g2 <- dt.pr$gene[abs(dt.pr$hg.lg_fx.hg) > 20 &
                   !is.na(dt.pr$hg.lg_fx.hg)]
g2

g3 <- dt.pr$gene[abs(dt.pr$hg.lg_ber.hg) > 20 &
                   !is.na(dt.pr$hg.lg_ber.hg)]
g3

# # Part III: Pathway analysis with Reactome----
# # Get Entrezgene IDs
# geneID <- as.character(unique(dt1$geneId[as.character(dt1$gene) %in% gene.sorted]))
# 
# # Save Entrezgene IDs and gene names
# out <- data.table(geneID = geneID,
#                   geneName = gene.sorted)
# write.csv(out,
#           file = "mes13/tmp/gene.entrezID.csv",
#           row.names = FALSE)
# 
# # Get pathways----
# m1 <- enrichPathway(gene = geneID,
#                     pvalueCutoff = 1, 
#                     readable = T, 
#                     organism = "mouse")
# t2 <- as.data.table(m1)
# t2
# write.csv(t2,
#           file = "mes13/tmp/pathways.csv",
#           row.names = FALSE)
# 
# # Pathway barplot----
# tiff(filename = "mes13/tmp/pathway.tiff",
#      height = 3,
#      width = 10,
#      units = 'in',
#      res = 300,
#      compression = "lzw+p")
# barplot(m1)
# graphics.off()
# 
# # Part IV: get specific pathways----
# path.names <- do.call("rbind", 
#                       as.list(reactome.db::reactomePATHID2NAME))
# path.names <- data.table(reactome_id = rownames(path.names),
#                          path = path.names[, 1])
# path.names
# 
# # Keep mouse pathways only
# path.mm <- subset(path.names,
#                   substr(path, 1, 12) == "Mus musculus")
# path.mm$path <- gsub(pattern = "Mus musculus: ",
#                      replacement = "",
#                      x = path.mm$path)
# path.mm
# 
# ##  David's selection of pathways, 09/20/2017----
# # TNF signaling
# # Activation of NF-kappaB in B cells
# # Activation of the AP-1 family of trTNF signalinganscription factors
# # Biological oxidations
# # COX reactions
# # Cell-extracellular matrix interactions
# # Extracellular matrix organization
# # IL-6-type cytokine receptor ligand interactions
# # Interleukin-6 family signaling
# # Regulation of TNFR1 signaling
# # TNFR2 non-canonical NF-kB pathway
# # NF-kB is activated and signals survival
# # Inflammasomes
# # Oxidative Stress Induced Senescence
# # TGF-beta receptor signaling activates SMADs
# # TGF-beta receptor signaling in EMT (epithelial to mesenchymal transition)
# 
# # Separate TNF signaling pathway----
# path.list <- c("TNF signaling",
#                "Activation of NF-kappaB in B cells",
#                "Activation of the AP-1 family of trTNF signalinganscription factors",
#                "Biological oxidations",
#                "COX reactions",
#                "Cell-extracellular matrix interactions",
#                "Extracellular matrix organization",
#                "IL-6-type cytokine receptor ligand interactions",
#                "Interleukin-6 family signaling",
#                "Regulation of TNFR1 signaling",
#                "TNFR2 non-canonical NF-kB pathway",
#                "NF-kB is activated and signals survival",
#                "Inflammasomes",
#                "Oxidative Stress Induced Senescence",
#                "TGF-beta receptor signaling activates SMADs",
#                "TGF-beta receptor signaling in EMT (epithelial to mesenchymal transition)")
# path.mm[path %in% path.list]
# 
# # Map path IDs to Entrez gene IDs----
# path2gene <- as.list(reactome.db::reactomePATHID2EXTID)
# path2gene <- path2gene[names(path2gene) %in% path.mm$reactome_id]
# 
# # Get TNF signaling gene Entrez IDs----
# get.genes <- path2gene[names(path2gene) %in% path.mm$reactome_id[path.mm$path %in% path.list]]
# get.genes <- unique(do.call("c", get.genes))
# 
# dt2 <- merge(dt1[, 1:2],
#              dt2, 
#              by = "gene")
# genes.from.path <- unique(dt2[dt2$geneId %in% get.genes, ])
# write.csv(genes.from.path,
#           file = "mes13/tmp/mes13_meth_genes.from.path.csv")
# 
# # CHECKPOINT: do the genes I found map to pathways correctly?
# m1 <- enrichPathway(gene = as.character(genes.from.path$geneId),
#                     pvalueCutoff = 0.01, 
#                     readable = TRUE, 
#                     organism = "mouse")
# t1 <- as.data.table(m1)
# t1
# barplot(m1, 
#         showCategory = 10)
# 
# tiff(filename = "mes13/tmp/mes13_meth_genes.from.path_pathways.tiff",
#      height = 6,
#      width = 10,
#      units = 'in',
#      res = 300,
#      compression = "lzw+p")
# barplot(m1, 
#         showCategory = 20)
# graphics.off()
# 
# # Heatmap of the selected genes----
# dt3 <- melt.data.table(data = genes.from.path,
#                        id.vars = c(1, 3),
#                        measure.vars = 5:8,
#                        variable.name = "Treatment",
#                        value.name = "Methylation(%)")
# 
# dt3$Treatment <- factor(dt3$Treatment)
# dt3
# 
# p3 <- ggplot(data = dt3) +
#   facet_wrap(~ reg,
#              #scales = "free_y",
#              nrow = 1) +
#   geom_tile(aes(x =  Treatment,
#                 y = gene,
#                 fill = `Methylation(%)`),
#             color = "black") +
#   scale_fill_gradient2(high = "red",
#                        limit = c(0, 100),
#                        name = "Methylation(%)") +
#   scale_x_discrete(expand = c(0, 0)) +
#   scale_y_discrete("Gene",
#                    expand = c(0, 0)) +
#   # ggtitle("MES13: Top 30 Largest Negative Diff in Methylation \n High - Low Glucose") +
#   theme(axis.text.x = element_text(angle = 30,
#                                    hjust = 1),
#         # legend.position = "top",
#         plot.title = element_text(hjust = 0.5))
# p3
# 
# tiff(filename = "mes13/tmp/mes13_meth_genes.from.pathways_heatmap.tiff",
#      height = 15,
#      width = 6,
#      units = 'in',
#      res = 300,
#      compression = "lzw+p")
# print(p3)
# graphics.off()
# write.csv(genes.from.path,
#           file = "mes13/tmp/mes13_meth_genes.from.pathways_heatmap.csv")

# sink()
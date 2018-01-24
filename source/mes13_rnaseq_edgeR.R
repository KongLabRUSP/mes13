# |----------------------------------------------------------------------------------|
# | Project: Study of Diabetes in MES13 cells                                        |
# | Script: RNA-seq data analysis and visualization using edgeR                      |
# | Author: Davit Sargsyan                                                           |
# | Created: 01/24/2018                                                              |
# | Modified:                                                                        |
# |----------------------------------------------------------------------------------|
# Source: 
# https://web.stanford.edu/class/bios221/labs/
# https://web.stanford.edu/class/bios221/labs/rnaseq/lab_4_rnaseq.html

# Header----
require(data.table)
require(ggplot2)
require(edgeR)
require(knitr)

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

# Remove genes with low counts----
tmp <- rowSums(dt1[, -1])
dt1 <- droplevels(subset(dt1,
                         tmp > 100))
dt1

cpm(dt1[, -1])

# edgeR----
dt2 <- DGEList(counts = as.matrix(dt1[, -1]),
               group = trt.names)
dt2

# Data normalization----
dt2 <- calcNormFactors(dt2)
dt2

# Explore the data----
# a. Multidimensional scaling plot of distances between 
#    gene expression profiles. 
#    This plots samples on a two-dimensional scatterplot so 
#    that distances on the plot approximate the typical 
#    log2 fold changes between the samples.
plotMDS(dt2)
#    Options: method = bcv: ??
plotMDS(dt2,
        method = "bcv")

# PCA----
m.pca <- prcomp(t(dt2$counts),
                center = TRUE,
                scale = TRUE)
summary(m.pca)
plot(m.pca)

# Keep only the most important cytokines (Javier)----
# Select PC-s to pliot (PC1 & PC2)
choices <- 1:2

nobs.factor <- sqrt(nrow(m.pca$x) - 1)
d <- m.pca$sdev
u <- m.pca$x
v <- m.pca$rotation

# Scores
df.u <- data.frame(u[, choices])
# Add grouping variable
df.u$grp <- rownames(df.u)
df.u

# Directions
df.v <- as.data.frame(v[, choices])

# Annotate
df.v$Gene <- dt1$gene

# Separate top variables (largest arrows)
df.v$lgth <- sqrt(df.v$PC1^2 + df.v$PC2^2)
df.v <- df.v[order(df.v$lgth,
                   decreasing = TRUE), ]

df.v$Gene <- factor(df.v$Gene,
                        levels = unique(df.v$Gene))
df.v

p1 <- ggplot(df.v[1:200, ]) +
  geom_bar(aes(x = Gene,
               y = lgth),
           stat = "identity") +
  ggtitle("Genes") +
  scale_x_discrete("Gene") +
  scale_y_continuous("Axis Length") +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 45,
                                   hjust = 1))
# p1
tiff(filename = "tmp/mes13_methylseq_pca_axis.tiff",
     height = 5,
     width = 20,
     units = 'in',
     res = 300,
     compression = "lzw+p")
print(p1)
graphics.off()

# Axis labels
u.axis.labs <- paste(colnames(df.v)[1:2], 
                     sprintf('(%0.1f%% explained var.)', 
                             100*m.pca$sdev[choices]^2/sum(m.pca$sdev^2)))

# Which cytokines to display
# var.keep.ndx <- 1:10
# # Genes with largest PC1 and PC1
# var.keep.ndx <- unique(c(which(order(abs(df.v$PC1)) %in% 1:10),
#                          which(order(abs(df.v$PC2)) %in% 1:10)))
# NEW(DS  01/24/2018): genes that are most influential in PC1 direction
dff <- df.v$PC1/df.v$PC2
var.keep.ndx <- which(order(dff) %in% c(1:10,
                                        (length(dff)-10):length(dff)))


xx <- 10000

p2 <- ggplot(data = df.v[var.keep.ndx,], 
             aes(x = PC1,
                 y = PC2)) +
  coord_equal() +
  geom_segment(aes(x = 0,
                   y = 0, 
                   xend = 0.95*xx*PC1,
                   yend = 0.95*xx*PC2),
               arrow = arrow(length = unit(1/2, 'picas')), 
               color = "black",
               size = 0.5) +
  geom_text(aes(x = xx*PC1,
                y = xx*PC2,
                label = df.v$Gene[var.keep.ndx]),
            size = 3,
            angle = 30,
            hjust = 0.5) +
  geom_point(data = df.u,
             aes(fill = grp),
             shape = 21,
             size = 3,
             alpha = 0.8) +
  scale_x_continuous(u.axis.labs[1]) +
  scale_y_continuous(u.axis.labs[2]) +
  scale_fill_discrete(name = "Treatment",
                      labels = trt.names) +
  ggtitle("Biplot of genes with 20 largest and\n20 smallest PC1/PC2 ratios\nMES13 RNA-seq") +
  coord_equal(ratio = 0.75) +
  theme(plot.title = element_text(hjust = 0.5,
                                  size = 20))
# p2

tiff(filename = "tmp/biplot_cytokines.tiff",
     height = 8,
     width = 8,
     units = 'in',
     res = 300,
     compression = "lzw+p")
print(p2)
graphics.off()

# Estimating dispersion----
dt2 <- estimateCommonDisp(dt2, verbose=T)
dt2
dt2 <- estimateTagwiseDisp(dt2)
dt2
plotBCV(dt2)
# Not possible without replicates! Set it manually
dt2$common.dispersion <- 0.05
plotBCV(dt2)

# Compute genewise exact tests for differences in the means 
# between two groups of negative-binomially distributed counts.
t1 <- edgeR::exactTest(dt2,
                       pair = c(1, 2))
# t1 <- edgeR::exactTest(dt2,
#                        pair = c(2, 3))
t1
topTags(t1, 
        n = 10)

sum1 <- decideTestsDGE(t1, 
                       adjust.method = "none",
                       p.value = 0.05)
summary(sum1)

plotSmear(t1)
abline(h = c(-0.5, 0.5), 
       col = "blue",
       lty = 2)

# GLM----
design.mat <- matrix(0, nrow = 7, ncol = 2)
design.mat[1:2, 1] <- design.mat[2:3, 2] <- 1
design.mat

m1 <- glmFit(dt2, design.mat)
summary(m1)

# Contrasts: HG - LG
lrt12 <- glmLRT(fit, 
                contrast = c(-1, 1))
lrt12
topTags(lrt12, 
        n = 10)
sum2 <- decideTestsDGE(lrt12, 
                       adjust.method = "none",
                       p.value = 0.05)
summary(sum2)


de2tags12 <- rownames(dt2)[as.logical(sum2)]
plotSmear(lrt12, 
          de.tags=de2tags12)
abline(h = c(-1, 1), 
       col = "blue",
       lty = 2)

# Log2 differences----
dt1$`diff(wj2,wj1)` <- log2(dt2$counts[, 2]/dt2$counts[, 1])

# Log2 means----
dt1$`mean(wj2,wj1)` <- log2(dt2$counts[, 2]*dt2$counts[, 1])/2

# Lot log2 differences vs. means----
plot(dt1$`diff(wj2,wj1)` ~ dt1$`mean(wj2,wj1)`)
plot(abs(dt1$`diff(wj2,wj1)`) ~ dt1$`mean(wj2,wj1)`)
# |----------------------------------------------------------------------------------|
# | Project: Study of Diabetes in MES13 cells                                        |
# | Script: RNA-seq SD estimate using difference vs. mean, LG vs. HG only            |
# | Author: Davit Sargsyan                                                           |
# | Created: 12/06/2017                                                              |
# |----------------------------------------------------------------------------------|
# Header----
# Save consol output to a log file
# sink(file = "tmp/log_mes13_rnaseq_diff_vs_mean_pval_v1.txt")
require(data.table)
require(ggplot2)

# Load data----
dt1 <- fread("data/rna_seq/gene_exp.csv")
dt1

# Log2 of (values +1)
dt1$lv1 <- log2(dt1$value_1 + 1)
dt1$lv2 <- log2(dt1$value_2 + 1)

# Keep LG and HG only
dt2 <- droplevels(subset(dt1,
                         sample_1 == "LG" &
                           sample_2 == "HG",
                         select = colnames(dt1) %in% c("gene",
                                                       "sample_1",
                                                       "sample_2",
                                                       "value_1",
                                                       "value_2",
                                                       "lv1",
                                                       "lv2")))
dt2
dt1 <- dt2

# p-Values for single-replica samples----
# Differences vs. means, original scale----
dt2$diff <- dt2$value_1 - dt2$value_2
dt2$mu <- (dt2$value_1 + dt2$value_2)/2
tiff(filename = "tmp/mes13_diff_vs_mu_original.tiff",
     height = 8,
     width = 8,
     units = 'in',
     res = 300,
     compression = "lzw+p")
plot(dt2$diff ~ dt2$mu,
     xlab = "Means of Expressions",
     ylab = "Differences of Expressions",
     main = "MES13 RNA-Seq: LG vs. HG")
graphics.off()

# Differences vs. means, log2 scale----
dt2$diff <- dt2$lv1 - dt2$lv2
dt2$mu <- (dt2$lv1 + dt2$lv2)/2

# Plot differences vs. means, original scale---
tiff(filename = "tmp/mes13_diff_vs_mu_log2.tiff",
     height = 8,
     width = 8,
     units = 'in',
     res = 300,
     compression = "lzw+p")
plot(dt2$diff ~ dt2$mu,
     xlab = "Mean of Log2(Expression)",
     ylab = "Difference of Log2(Expression)",
     main = "MES13 RNA-Seq: LG vs. HG")
graphics.off()

# Remove all zeros and plot again----
dt2 <- dt2[value_1 != 0 & value_2 != 0, ]
dt2$diff <- dt2$value_1 - dt2$value_2
dt2$mu <- (dt2$value_1 + dt2$value_2)/2
dt2$diff <- dt2$lv1 - dt2$lv2
dt2$mu <- (dt2$lv1 + dt2$lv2)/2

tiff(filename = "tmp/mes13_diff_vs_mu_log2_no_zeros.tiff",
     height = 8,
     width = 8,
     units = 'in',
     res = 300,
     compression = "lzw+p")
plot(dt2$diff ~ dt2$mu,
     xlab = "Mean of Log2(Expression)",
     ylab = "Difference of Log2(Expression)",
     main = "MES13 RNA-Seq: LG vs. HG, No Zeros")
graphics.off()

# Checkpoint
dt2[diff < -6, ]

# Regularize by the above within margin of epsilon
epln <- 0.1
dt2$sd <- NA

for (i in 1:nrow(dt2)) {
  tmp <- subset(dt2,
                (mu <= dt2$mu[i] + epln) &
                  (mu >= dt2$mu[i] - epln))
  dt2$sd[i] <- sd(tmp$diff)
}

tiff(filename = "tmp/mes13_diff_vs_mu_log2_no_zeros_sd_hist.tiff",
     height = 5,
     width = 5,
     units = 'in',
     res = 300,
     compression = "lzw+p")
hist(dt2$sd, 100)
graphics.off()

dt2

# Smooth SD with LOESS----
m1 <- loess(sd ~ mu, 
            data = dt2)
summary(m1)

dt2$sd.fit <- predict(m1,
                      newdata = data.frame(mu = dt2$mu))
dt2 <- dt2[order(dt2$mu), ]
dt2

# Plot estimated and predicted SD----
tiff(filename = "tmp/mes13_diff_vs_mu_log2_no_zeros_sd_pred_and_fitted.tiff",
     height = 5,
     width = 6,
     units = 'in',
     res = 300,
     compression = "lzw+p")
plot(dt2$sd ~ dt2$mu,
     xlab = "Means",
     ylab = "SD",
     main = "MES13 RNA-Seq: LG vs. HG, No Zeros")
lines(dt2$sd.fit ~ dt2$mu,
      col = "red",
      lw = 2)
graphics.off()

# Assuming t follows normal distribution (it does not!)
dt2$t <- dt2$diff/dt2$sd.fit

tiff(filename = "tmp/mes13_diff_vs_mu_log2_no_zeros_qqplot.tiff",
     height = 5,
     width = 5,
     units = 'in',
     res = 300,
     compression = "lzw+p")
qqnorm(dt2$t)
abline(0, 1)
graphics.off()

tiff(filename = "tmp/mes13_diff_vs_mu_log2_no_zeros_stat_hist.tiff",
     height = 5,
     width = 5,
     units = 'in',
     res = 300,
     compression = "lzw+p")
hist(dt2$t, 50,
     main = "MES13 RNA-Seq Test Statistic's Histogram",
     xlab = "Test Statistic")
graphics.off()

dt2$p <- 2*pnorm(-abs(dt2$t))

tiff(filename = "tmp/mes13_diff_vs_mu_log2_no_zeros_pval.tiff",
     height = 5,
     width = 5,
     units = 'in',
     res = 300,
     compression = "lzw+p")
hist(dt2$p, 100,
     main = "MES13 RNA-Seq p_Value Histogram",
     xlab = "p-Value")
graphics.off()

# Save as CSV----
write.csv(dt2, file = "tmp/mes13_pvals.csv")

# Predict on all 24k+ genes----
dt1$diff <- dt1$lv1 - dt1$lv2
dt1$mu <- (dt1$lv1 + dt1$lv2)/2
dt1
dt1$sd.fit <- predict(m1,
                      newdata = data.frame(mu = dt1$mu))
dt1
sum(!is.na(dt1$sd.fit))
dt1$t <- dt1$diff/dt1$sd.fit

tiff(filename = "tmp/mes13_diff_vs_mu_log2_no_zeros_stat_hist_all.tiff",
     height = 5,
     width = 5,
     units = 'in',
     res = 300,
     compression = "lzw+p")
hist(dt1$t, 100,
     main = "MES13 RNA-Seq Test Statistic's Histogram, All Genes",
     xlab = "Test Statistic")
graphics.off()

dt1$p <- 2*pnorm(-abs(dt1$t))

tiff(filename = "tmp/mes13_diff_vs_mu_log2_no_zeros_pval_all_genes.tiff",
     height = 5,
     width = 5,
     units = 'in',
     res = 300,
     compression = "lzw+p")
hist(dt1$p, 100,
     main = "MES13 RNA-Seq p_Value Histogram, All Genes",
     xlab = "p-Value")
graphics.off()

# Plot estimated and predicted SD, all 24k+ genes----
tiff(filename = "tmp/mes13_diff_vs_mu_log2_all_genes_sd_pred_and_fitted_all_genes.tiff",
     height = 5,
     width = 6,
     units = 'in',
     res = 300,
     compression = "lzw+p")
plot(dt1$sd.fit ~ dt1$mu,
     xlab = "Means",
     ylab = "SD",
     main = "MES13 RNA-Seq: LG vs. HG, No Zeros",
     col = "red")
graphics.off()

# Save as CSV----
write.csv(dt1, file = "tmp/mes13_pvals_all_genes.csv")

# sink()
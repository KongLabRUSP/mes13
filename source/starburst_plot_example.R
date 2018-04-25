require(data.table)
require(ggplot2)

dt1 <- data.table(gene = paste("g", 1:1000, sep = ""),
                  diffDNA = runif(n = 1000, 
                                  min = -100,
                                  max = 100),
                  qDNA = runif(n = 1000,
                               min = 0,
                               max = 1),
                  log2RNA = rnorm(n = 1000,
                                  mean = 0,
                                  sd = 4),
                  qDNA = runif(n = 1000,
                               min = 0,
                               max = 1),
                  grp = 1)
dt1
dt1$grp[dt1$diffDNA <  -10 & dt1$log2RNA > 2] <- 2
dt1$grp[dt1$diffDNA > 10 & dt1$log2RNA < -2] <- 3
dt1$grp <- factor(dt1$grp)

ggplot(data = dt1,
       aes(x = diffDNA,
           y = log2RNA,
           colour = grp)) +
  geom_point(alpha = 0.5) +
  geom_hline(yintercept = c(-2, 2),
             linetype = "dotted") +
  geom_vline(xintercept = c(-10, 10),
             linetype = "dotted") +
  geom_text(data = dt1[1:5, ],
            aes(x = diffDNA,
                y = log2RNA,
                label = gene),
            color = "blue",
            size = 10) +
  scale_color_manual(values = levels(dt1$grp))

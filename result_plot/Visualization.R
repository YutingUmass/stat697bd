### Load FDR Data
RR_result_ratio_1_1 <- as.data.frame(read.csv("~/Downloads/RR_result_ratio_1_1.csv"))
RR_result_ratio_1_3 <- as.data.frame(read.csv("~/Downloads/RR_result_ratio_1_3.csv"))
RR_result_ratio_1_9 <- as.data.frame(read.csv("~/Downloads/RR_result_ratio_1_9.csv"))

library(ggpubr)
variation = c("Low", "Moderate", "Strong")
group = c('Chi-Squared', 'Mann-Whitney', 'Wilcoxon', 'Kolmogorov–Smirnov')
p1 <- ggplot(color = group) + geom_line(RR_result_ratio_1_1 ,mapping = aes(x = variation, y = RR_result_ratio_1_1[,2], group = group[1],color = 'Chi-Squared')) +
  geom_point(RR_result_ratio_1_1 ,mapping = aes(x = variation, y = RR_result_ratio_1_1[,2], group = group[1],color = 'Chi-Squared')) + 
  geom_line(RR_result_ratio_1_1 ,mapping = aes(x = variation, y = RR_result_ratio_1_1[,3], group = group[2], color = 'Mann-Whitney')) +
  geom_point(RR_result_ratio_1_1 ,mapping = aes(x = variation, y = RR_result_ratio_1_1[,3], group = group[2], color = 'Mann-Whitney')) +
  geom_line(RR_result_ratio_1_1 ,mapping = aes(x = variation, y = RR_result_ratio_1_1[,4], group = group[3], color = 'Wilcoxon')) +
  geom_point(RR_result_ratio_1_1 ,mapping = aes(x = variation, y = RR_result_ratio_1_1[,4], group = group[3], color = 'Wilcoxon')) +
  geom_line(RR_result_ratio_1_1 ,mapping = aes(x = variation, y = RR_result_ratio_1_1[,5], group = group[4], color = 'Kolmogorov–Smirnov')) +
  geom_point(RR_result_ratio_1_1 ,mapping = aes(x = variation, y = RR_result_ratio_1_1[,5], group = group[4], color = 'Kolmogorov–Smirnov')) +
  ylim(0,1) + labs(x = "Batch Effect", y = "Rejection Rate") + theme(plot.title = element_text(hjust = 0.5)) + ggtitle("r = 1")

p2 <- ggplot(color = group) + geom_line(RR_result_ratio_1_1 ,mapping = aes(x = variation, y = RR_result_ratio_1_3[,2], group = group[1],color = 'Chi-Squared')) +
  geom_point(RR_result_ratio_1_1 ,mapping = aes(x = variation, y = RR_result_ratio_1_3[,2], group = group[1],color = 'Chi-Squared')) + 
  geom_line(RR_result_ratio_1_1 ,mapping = aes(x = variation, y = RR_result_ratio_1_3[,3], group = group[2], color = 'Mann-Whitney')) +
  geom_point(RR_result_ratio_1_1 ,mapping = aes(x = variation, y = RR_result_ratio_1_3[,3], group = group[2], color = 'Mann-Whitney')) +
  geom_line(RR_result_ratio_1_1 ,mapping = aes(x = variation, y = RR_result_ratio_1_3[,4], group = group[3], color = 'Wilcoxon')) +
  geom_point(RR_result_ratio_1_1 ,mapping = aes(x = variation, y = RR_result_ratio_1_3[,4], group = group[3], color = 'Wilcoxon')) +
  geom_line(RR_result_ratio_1_1 ,mapping = aes(x = variation, y = RR_result_ratio_1_3[,5], group = group[4], color = 'Kolmogorov–Smirnov')) +
  geom_point(RR_result_ratio_1_1 ,mapping = aes(x = variation, y = RR_result_ratio_1_3[,5], group = group[4], color = 'Kolmogorov–Smirnov')) +
  ylim(0,1) + labs(x = "Batch Effect") + theme(axis.title.y=element_blank(),legend.position = "none",plot.title = element_text(hjust = 0.5)) + ggtitle("r = 1/3")

p3 <- ggplot(color = group) + geom_line(RR_result_ratio_1_1 ,mapping = aes(x = variation, y = RR_result_ratio_1_9[,2], group = group[1],color = 'Chi-Squared')) +
  geom_point(RR_result_ratio_1_1 ,mapping = aes(x = variation, y = RR_result_ratio_1_9[,2], group = group[1],color = 'Chi-Squared')) + 
  geom_line(RR_result_ratio_1_1 ,mapping = aes(x = variation, y = RR_result_ratio_1_9[,3], group = group[2], color = 'Mann-Whitney')) +
  geom_point(RR_result_ratio_1_1 ,mapping = aes(x = variation, y = RR_result_ratio_1_9[,3], group = group[2], color = 'Mann-Whitney')) +
  geom_line(RR_result_ratio_1_1 ,mapping = aes(x = variation, y = RR_result_ratio_1_9[,4], group = group[3], color = 'Wilcoxon')) +
  geom_point(RR_result_ratio_1_1 ,mapping = aes(x = variation, y = RR_result_ratio_1_9[,4], group = group[3], color = 'Wilcoxon')) +
  geom_line(RR_result_ratio_1_1 ,mapping = aes(x = variation, y = RR_result_ratio_1_9[,5], group = group[4], color = 'Kolmogorov–Smirnov')) +
  geom_point(RR_result_ratio_1_1 ,mapping = aes(x = variation, y = RR_result_ratio_1_9[,5], group = group[4], color = 'Kolmogorov–Smirnov')) +
  ylim(0,1) + labs(x = "Batch Effect") + theme(axis.title.y=element_blank(),legend.position = "none",plot.title = element_text(hjust = 0.5)) + ggtitle("r = 1/9")



ggarrange(p1, p2, p3, ncol=3, nrow=1, common.legend = TRUE, legend="bottom")



### PCA plot 
library(splatter)
library(Seurat)
library(scuttle)
install.packages("EMA")
library(scater)
sim1 <- splatSimulate(nGenes = 1000, batchCells = c(250, 250),
                      batch.facLoc = 0.01, batch.facScale = 0.01,
                      verbose = FALSE)
sim1 <- logNormCounts(sim1)
sim1 <- runPCA(sim1)
p1 <- plotPCA(sim1, colour_by = "Batch") + ggtitle("Low Batch Effect")


sim2 <- splatSimulate(nGenes = 1000, batchCells = c(250, 250),
                      batch.facLoc = 0.05, batch.facScale = 0.05,
                      verbose = FALSE)
sim2 <- logNormCounts(sim2)
sim2 <- runPCA(sim2)
p2 <- plotPCA(sim2, colour_by = "Batch") + ggtitle("Moderate Batch Effect")


sim3 <- splatSimulate(nGenes = 1000, batchCells = c(250, 250),
                      batch.facLoc = 0.1, batch.facScale = 0.1,
                      verbose = FALSE)
sim3 <- logNormCounts(sim3)
sim3 <- runPCA(sim3)
p3 <- plotPCA(sim3, colour_by = "Batch") + ggtitle("Strong Batch Effect")


ggarrange(p1, p2, p3, ncol=3, nrow=1, common.legend = TRUE, legend="bottom")


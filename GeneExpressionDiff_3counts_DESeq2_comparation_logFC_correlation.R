# Salmon, featureCounts, HTSeq-count 3中方法计算的gene-count； 用DESeq2做基因表达差异分析，比较分析结果
# Target: 3种方法DESeq2 差异基因比较(padj<0.05)
# 加载必要的包
library(dplyr)
library(readr)
library(ggVennDiagram) # 韦恩图
library(ggplot2)
# library(gridExtra) # 用于组合图形（可选）
# BiocManager::install("")
# library(ggvenn) # 更简洁的2个样品的韦恩图，不好看
# library(Vennerable) #韦恩图
#library(GGally)  # scaterplot
library(ggpubr) # scaterplot
library(tidyr) 
library(patchwork)   # 用于合并图片

# 设置工作目录
setwd("/Users/lixia/Data/data/Fractionation_seq")
# 创建结果保存目录
dir.create("3counts_camperation", showWarnings = FALSE)

### ————————————————————————————— Step 1: 读取3种方法DESeq2 差异基因较(padj<0.05) 数据———————————————--——-—————-------------— ###
#（1）读取result数据
HTSeqcount_dif_file <- "HTSeqcount_labE/DESeq2_results/HTSeqcount_DESeq2_res_ged_counts_geneName_TPM_Siggene.csv"
featureCounts_dif_file <- "featureCounts_labE/DESeq2_results/featureCount_DESeq2_res_ged_counts_geneName_TPM_Siggene.csv"
Salmon_dif_file <- "Salmon/Ensembl/DESeq2_results/Salmon_DESeq2_res_ged_counts_geneName_TPM_Siggene.csv"

HTSeqcount_dif <- as.data.frame(read_csv(HTSeqcount_dif_file, col_names = TRUE, locale = locale(encoding = "UTF-8")))
featureCounts_dif <- as.data.frame(read_csv(featureCounts_dif_file, col_names = TRUE, locale = locale(encoding = "UTF-8")))
Salmon_dif <- as.data.frame(read_csv(Salmon_dif_file, col_names = TRUE, locale = locale(encoding = "UTF-8")))

print(paste("HTSeqcount检测到的差异表达基因数量为:", nrow(HTSeqcount_dif)))
print(paste("featureCounts检测到的差异表达基因数量为:", nrow(featureCounts_dif)))
print(paste("Salmon检测到的差异表达基因数量为:", nrow(Salmon_dif)))

### ————————————————————————————— Step 2: 检查关键基因(XIST, GAPDH)———————————————--------------------——-—————-------------— ###

#（2）检查XIST，GAPDH的log2FoldChange
result1 <- HTSeqcount_dif$log2FoldChange[HTSeqcount_dif$gene_name == c("XIST")]
result2 <- featureCounts_dif$log2FoldChange[featureCounts_dif$gene_name == "XIST"]
result3 <- Salmon_dif$log2FoldChange[Salmon_dif$gene_name == "XIST"]

result4 <- HTSeqcount_dif$log2FoldChange[HTSeqcount_dif$gene_name == c("GAPDH")]
result5 <- featureCounts_dif$log2FoldChange[featureCounts_dif$gene_name == "GAPDH"]
result6 <- Salmon_dif$log2FoldChange[Salmon_dif$gene_name == "GAPDH"]

cat("XIST在HTSeqcount,featureCounts,Salmon中log2FoldChange分别是:", result1, result2, result3) 
cat("GAPDH在HTSeqcount,featureCounts,Salmon中log2FoldChange分别是:", result4, result5, result6)

### ————————————————————————————— Step3: 创建gene_id log2FC 的3种方法并集表——————————————---------------------——-—————---------— ###

# 从每个数据集中提取gene_id和log2FoldChange
htseq <- HTSeqcount_dif %>% 
  select(gene_id, log2FC_HTSeq = log2FoldChange)

feature <- featureCounts_dif %>% 
  select(gene_id, log2FC_featureCounts = log2FoldChange)

salmon <- Salmon_dif %>% 
  select(gene_id, log2FC_Salmon = log2FoldChange)

# 创建总表（全外连接）
combined <- full_join(htseq, feature, by = "gene_id") %>%
  full_join(salmon, by = "gene_id") %>%
  replace(is.na(.), 0)  # 将NA替换为0

### ————————————————————————————— Step4: sig_gene log2FC 的scaterplot (3种方法取并集)———————————————---------------------——-—————---------— ###

# 计算相关系数（Pearson和Spearman）
cor_pearson_union3 <- cor(combined[, c("log2FC_HTSeq", "log2FC_featureCounts", "log2FC_Salmon")], method = "pearson")
cor_spearman_union3 <- cor(combined[, c("log2FC_HTSeq", "log2FC_featureCounts", "log2FC_Salmon")], method = "spearman")

# 创建correlation coefficient 保存文件
content <- c(
  "correlation_information",
  "methods_person_(up)_sperman_(down)",
  "",
  "cor_pearson_union3", 
  "cor_spearman_union3",
  "",
  "HTSeq_feature_union_pearson",
  "HTSeq_Salmon_union_pearson",
  "feature_Salmon_union_spearman",
  "",
  "HTSeq_feature_union_spearman",
  "HTSeq_Salmon_union_spearman",
  "feature_Salmon_union_spearman",
  "",
  "cor_pearson_inter3",
  "cor_spearman_inter3",
  "",
  "HTSeq_feature_inter_pearson",
  "HTSeq_Salmon_inter_pearson",
  "feature_Salmon_inter_spearman"
)

# 写入文本文件
writeLines(content,"3counts_camperation/correlation_coefficient.txt")

#保存相关系数，追加到记录表
write.table(cor_pearson_union3, "3counts_camperation/correlation_coefficient.txt", append = TRUE, 
            col.names = T, row.names = T)
write.table(cor_spearman_union3, "3counts_camperation/correlation_coefficient.txt", append = TRUE, 
            col.names = T, row.names = T)

# 绘制两两散点图函数
create_scatter <- function(data, x_var, y_var) {
  pearson <- cor(data[[x_var]], data[[y_var]], method = "pearson")
  spearman <- cor(data[[x_var]], data[[y_var]], method = "spearman")
  
  ggplot(data, aes_string(x = x_var, y = y_var)) +
    geom_point(alpha = 0.5, size = 1) +
    geom_smooth(method = "lm", se = FALSE, color = "red") +
    labs(title = paste(x_var, "vs", y_var)) +
    annotate("text", x = Inf, y = Inf, 
             label = sprintf("Pearson: %.2f\nSpearman: %.2f", pearson, spearman),
             hjust = 1.1, vjust = 1.1, size = 5)
}

# 生成并保存散点图
p1 <- create_scatter(combined, "log2FC_HTSeq", "log2FC_featureCounts")
p2 <- create_scatter(combined, "log2FC_HTSeq", "log2FC_Salmon")
p3 <- create_scatter(combined, "log2FC_featureCounts", "log2FC_Salmon")

# 组合图形
logFC_scaterplot <- ggarrange(p1, p2, p3, ncol = 2, nrow = 2)
#combined_1_2_3 <- (p1+ p2) / (p3) 

# 保存图形（PNG格式）
ggsave("3counts_camperation/sig_gene_logFC_scaterplot_unino3.png", logFC_scaterplot, width = 8, height = 6, dpi = 300)

### ————————————————————————————— Step5: sig_gene log2FC 的scaterplot (3种方法两两之间取并集)———————————————---------------------——-—————---------— ###

HTSeq_feature_union <- combined %>%
  filter(!(log2FC_HTSeq == 0 & log2FC_featureCounts == 0))

HTSeq_Salmon_union <- combined %>%
  filter(!(log2FC_HTSeq == 0 & log2FC_Salmon == 0))

feature_Salmon_union <- combined %>%
  filter(!(log2FC_featureCounts == 0 & log2FC_Salmon == 0))

HTSeq_feature_union_pearson <- cor(HTSeq_feature_union[, c("log2FC_HTSeq", "log2FC_featureCounts")], method = "pearson")
HTSeq_Salmon_union_pearson <- cor(HTSeq_Salmon_union[, c("log2FC_HTSeq", "log2FC_Salmon")], method = "pearson")
feature_Salmon_union_pearson <- cor(feature_Salmon_union[, c("log2FC_Salmon", "log2FC_featureCounts")], method = "pearson")

HTSeq_feature_union_spearman <- cor(HTSeq_feature_union[, c("log2FC_HTSeq", "log2FC_featureCounts")], method = "spearman")
HTSeq_Salmon_union_spearman <- cor(HTSeq_Salmon_union[, c("log2FC_HTSeq", "log2FC_Salmon")], method = "spearman")
feature_Salmon_union_spearman <- cor(feature_Salmon_union[, c("log2FC_Salmon", "log2FC_featureCounts")], method = "spearman")

#保存相关系数，追加到记录表
write.table(HTSeq_feature_union_pearson, "3counts_camperation/correlation_coefficient.txt", append = TRUE, 
            col.names = T, row.names = T)
write.table(HTSeq_Salmon_union_pearson, "3counts_camperation/correlation_coefficient.txt", append = TRUE, 
            col.names = T, row.names = T)
write.table(feature_Salmon_union_pearson, "3counts_camperation/correlation_coefficient.txt", append = TRUE, 
            col.names = T, row.names = T)
write.table(HTSeq_feature_union_spearman, "3counts_camperation/correlation_coefficient.txt", append = TRUE, 
            col.names = T, row.names = T)
write.table(HTSeq_Salmon_union_spearman, "3counts_camperation/correlation_coefficient.txt", append = TRUE, 
            col.names = T, row.names = T)
write.table(feature_Salmon_union_spearman, "3counts_camperation/correlation_coefficient.txt", append = TRUE, 
            col.names = T, row.names = T)

# 绘制两两散点图函数
create_scatter <- function(data, x_var, y_var) {
  pearson <- cor(data[[x_var]], data[[y_var]], method = "pearson")
  spearman <- cor(data[[x_var]], data[[y_var]], method = "spearman")
  
  ggplot(data, aes_string(x = x_var, y = y_var)) +
    geom_point(alpha = 0.5, size = 1) +
    geom_smooth(method = "lm", se = FALSE, color = "red") +
    labs(title = paste(x_var, "vs", y_var)) +
    annotate("text", x = Inf, y = Inf, 
             label = sprintf("Pearson: %.2f\nSpearman: %.2f", pearson, spearman),
             hjust = 1.1, vjust = 1.1, size = 5)
}

# 生成并保存散点图
HTSeq_feature_union_cor <- create_scatter(HTSeq_feature_union, "log2FC_HTSeq", "log2FC_featureCounts")
HTSeq_Salmon_union_cor <- create_scatter(HTSeq_Salmon_union, "log2FC_HTSeq", "log2FC_Salmon")
feature_Salmon_union_cor <- create_scatter(feature_Salmon_union, "log2FC_featureCounts", "log2FC_Salmon")

# 组合图形
logFC_scaterplot_unino2 <- ggarrange(HTSeq_feature_union_cor, HTSeq_Salmon_union_cor, feature_Salmon_union_cor, ncol = 2, nrow = 2)
#combined_1_2_3 <- (p1+ p2) / (p3) 

# 保存图形（PNG格式）
ggsave("3counts_camperation/sig_gene_logFC_scaterplot_unino2.png", logFC_scaterplot_unino2, width = 8, height = 6, dpi = 300)

### ————————————————————————————— Step6: sig_gene log2FC 的scaterplot (3种方法取交集)———————————————---------------------——-—————---------— ###

combined_nonzero <- combined %>%
  filter(log2FC_HTSeq != 0 & 
           log2FC_featureCounts != 0 & 
           log2FC_Salmon != 0)

# 计算相关系数（Pearson和Spearman）
cor_pearson_inter3 <- cor(combined_nonzero[, c("log2FC_HTSeq", "log2FC_featureCounts", "log2FC_Salmon")], method = "pearson")
cor_spearman_inter3 <- cor(combined_nonzero[, c("log2FC_HTSeq", "log2FC_featureCounts", "log2FC_Salmon")], method = "spearman")

#保存相关系数，追加到记录表
write.table(cor_pearson_inter3, "3counts_camperation/correlation_coefficient.txt", append = TRUE, 
            col.names = T, row.names = T)
write.table(cor_spearman_inter3, "3counts_camperation/correlation_coefficient.txt", append = TRUE, 
            col.names = T, row.names = T)

# 绘制两两散点图函数
create_scatter <- function(data, x_var, y_var) {
  pearson <- cor(data[[x_var]], data[[y_var]], method = "pearson")
  spearman <- cor(data[[x_var]], data[[y_var]], method = "spearman")
  
  ggplot(data, aes_string(x = x_var, y = y_var)) +
    geom_point(alpha = 0.5, size = 1) +
    geom_smooth(method = "lm", se = FALSE, color = "red") +
    labs(title = paste(x_var, "vs", y_var)) +
    annotate("text", x = Inf, y = Inf, 
             label = sprintf("Pearson: %.2f\nSpearman: %.2f", pearson, spearman),
             hjust = 1.1, vjust = 1.1, size = 5)
}

# 生成并保存散点图
p1_inter3 <- create_scatter(combined_nonzero, "log2FC_HTSeq", "log2FC_featureCounts")
p2_inter3 <- create_scatter(combined_nonzero, "log2FC_HTSeq", "log2FC_Salmon")
p3_inter3 <- create_scatter(combined_nonzero, "log2FC_featureCounts", "log2FC_Salmon")

# 组合图形
logFC_scaterplot_inter3 <- ggarrange(p1_inter3, p2_inter3, p3_inter3, ncol = 2, nrow = 2)
#combined_1_2_3 <- (p1+ p2) / (p3) 

# 保存图形（PNG格式）
ggsave("3counts_camperation/sig_gene_logFC_scaterplot_intersection3.png", logFC_scaterplot_inter3, width = 8, height = 6, dpi = 300)

### ————————————————————————————— Step7: sig_gene log2FC 的scaterplot (3种方法两两之间取交集)———————————————---------------------——-—————---------— ###

HTSeq_feature_inter <- combined %>%
  filter((log2FC_HTSeq != 0 & log2FC_featureCounts != 0))

HTSeq_Salmon_inter <- combined %>%
  filter((log2FC_HTSeq != 0 & log2FC_Salmon != 0))

feature_Salmon_inter <- combined %>%
  filter((log2FC_featureCounts != 0 & log2FC_Salmon != 0))

HTSeq_feature_inter_pearson <- cor(HTSeq_feature_inter[, c("log2FC_HTSeq", "log2FC_featureCounts")], method = "pearson")
HTSeq_Salmon_inter_pearson <- cor(HTSeq_Salmon_inter[, c("log2FC_HTSeq", "log2FC_Salmon")], method = "pearson")
feature_Salmon_inter_pearson <- cor(feature_Salmon_inter[, c("log2FC_Salmon", "log2FC_featureCounts")], method = "pearson")

HTSeq_feature_inter_spearman <- cor(HTSeq_feature_inter[, c("log2FC_HTSeq", "log2FC_featureCounts")], method = "spearman")
HTSeq_Salmon_inter_spearman <- cor(HTSeq_Salmon_inter[, c("log2FC_HTSeq", "log2FC_Salmon")], method = "spearman")
feature_Salmon_inter_spearman <- cor(feature_Salmon_inter[, c("log2FC_Salmon", "log2FC_featureCounts")], method = "spearman")

#保存相关系数，追加到记录表
write.table(HTSeq_feature_inter_pearson, "3counts_camperation/correlation_coefficient.txt", append = TRUE, 
            col.names = T, row.names = T)
write.table(HTSeq_Salmon_inter_pearson, "3counts_camperation/correlation_coefficient.txt", append = TRUE, 
            col.names = T, row.names = T)
write.table(feature_Salmon_inter_pearson, "3counts_camperation/correlation_coefficient.txt", append = TRUE, 
            col.names = T, row.names = T)
write.table(HTSeq_feature_inter_spearman, "3counts_camperation/correlation_coefficient.txt", append = TRUE, 
            col.names = T, row.names = T)
write.table(HTSeq_Salmon_inter_spearman, "3counts_camperation/correlation_coefficient.txt", append = TRUE, 
            col.names = T, row.names = T)
write.table(feature_Salmon_inter_spearman, "3counts_camperation/correlation_coefficient.txt", append = TRUE, 
            col.names = T, row.names = T)

# 绘制两两散点图函数
create_scatter <- function(data, x_var, y_var) {
  pearson <- cor(data[[x_var]], data[[y_var]], method = "pearson")
  spearman <- cor(data[[x_var]], data[[y_var]], method = "spearman")
  
  ggplot(data, aes_string(x = x_var, y = y_var)) +
    geom_point(alpha = 0.5, size = 1) +
    geom_smooth(method = "lm", se = FALSE, color = "red") +
    labs(title = paste(x_var, "vs", y_var)) +
    annotate("text", x = Inf, y = Inf, 
             label = sprintf("Pearson: %.2f\nSpearman: %.2f", pearson, spearman),
             hjust = 1.1, vjust = 1.1, size = 5)
}

# 生成并保存散点图
HTSeq_feature_inter_cor <- create_scatter(HTSeq_feature_inter, "log2FC_HTSeq", "log2FC_featureCounts")
HTSeq_Salmon_inter_cor  <- create_scatter(HTSeq_Salmon_inter, "log2FC_HTSeq", "log2FC_Salmon")
feature_Salmon_inter_cor  <- create_scatter(feature_Salmon_inter, "log2FC_featureCounts", "log2FC_Salmon")

# 组合图形
logFC_scaterplot_inter2 <- ggarrange(HTSeq_feature_inter_cor, HTSeq_Salmon_inter_cor, feature_Salmon_inter_cor, ncol = 2, nrow = 2)
#combined_1_2_3 <- (p1+ p2) / (p3) 

# 保存图形（PNG格式）
ggsave("3counts_camperation/sig_gene_logFC_scaterplot_intersection2.png", logFC_scaterplot_inter2, width = 8, height = 6, dpi = 300)



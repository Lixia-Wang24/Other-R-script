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
limma_res_file <-   "HTSeqcount_labE/limma_results/limma_DEG_results.csv"

HTSeqcount_dif <- as.data.frame(read_csv(HTSeqcount_dif_file, col_names = TRUE, locale = locale(encoding = "UTF-8")))
limma_res_dif <- as.data.frame(read_csv(limma_res_file, col_names = TRUE, locale = locale(encoding = "UTF-8")))

print(paste("DESeq2检测到的差异表达基因数量为:", nrow(HTSeqcount_dif)))
print(paste("limma检测到的差异表达基因数量为:", nrow(limma_res_dif)))


### ————————————————————————————— Step3: 创建gene_id log2FC 的3种方法并集表——————————————---------------------——-—————---------— ###

# 从每个数据集中提取gene_id和log2FoldChange
DESeq2 <- HTSeqcount_dif %>% 
  select(gene_id, log2FC_DESeq2 = log2FoldChange)

limma <- limma_res_dif %>% 
  select(gene_id, log2FC_limma = logFC)

# 创建总表（全外连接）
combined <- full_join(DESeq2, limma, by = "gene_id") %>%
  replace(is.na(.), 0)  # 将NA替换为0


### ————————————————————————————— Step7: sig_gene log2FC 的scaterplot (3种方法两两之间取交集)———————————————---------------------——-—————---------— ###

DESeq2_limma_inter <- combined %>%
  filter((log2FC_DESeq2 != 0 & log2FC_limma != 0))

DESeq2_limma_inter_pearson <- cor(DESeq2_limma_inter[, 2:3], method = "pearson")

DESeq2_limma_inter_spearman <- cor(DESeq2_limma_inter[, 2:3], method = "spearman")


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
DESeq2_limma_inter_cor <- create_scatter(DESeq2_limma_inter, "log2FC_DESeq2", "log2FC_limma")

# 保存图形（PNG格式）
ggsave("3counts_camperation/DESeq2_limma_inter_cor_logFC_scaterplot.png", DESeq2_limma_inter_cor, width = 8, height = 6, dpi = 300)

# 结果解释
# 图中右上角出现两条近似平行的黑点线，是由于 低表达基因在计算 log2FC 时引入伪计数（pseudocount）导致的离散化效应。
# DESeq2 默认使用 log2(n + 1) 的伪计数策略（n 为原始计数）。
# limma 默认使用 log2(n + 0.125) 的伪计数（通过 prior.count 参数调整）。
# 不同伪计数会导致 log2FC 值离散化，形成多条平行线。
# fit <- glmQLFit(y, design, prior.count = 1)  修改伪计数为log2(n + 1)。


### ————————————————————————————— Step5: sig_gene log2FC 的scaterplot (3种方法两两之间取并集)———————————————---------------------——-—————---------— ###

DESeq2_limma_union <- combined %>%
  filter(!(log2FC_DESeq2 == 0 & log2FC_limma == 0))

DESeq2_limma_union_pearson <- cor(DESeq2_limma_union[, 2:3], method = "pearson")

DESeq2_limma_union_spearman <- cor(DESeq2_limma_union[, 2:3], method = "spearman")

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
DESeq2_limma_union_cor <- create_scatter(DESeq2_limma_union, "log2FC_DESeq2", "log2FC_limma")

# 保存图形（PNG格式）
ggsave("3counts_camperation/DESeq2_limma_union_cor_logFC_scaterplot.png", DESeq2_limma_union_cor, width = 8, height = 6, dpi = 300)


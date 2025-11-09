# Salmon, featureCounts, HTSeq-count 3中方法计算的gene-count； 用DESeq2做基因表达差异分析，比较分析结果
# 加载必要的包
library(dplyr)
library(readr)
library(ggVennDiagram) # 韦恩图
library(ggplot2)
# library(gridExtra) # 用于组合图形（可选）
# BiocManager::install("")
# library(ggvenn) # 更简洁的2个样品的韦恩图，不好看
# library(Vennerable) #韦恩图
library(GGally)  # scaterplot
library(tidyr) # pivot_longer()
library(pheatmap)


# 设置工作目录
setwd("/Users/lixia/Data/data/Fractionation_seq")
# 创建结果保存目录
dir.create("3counts_camperation", showWarnings = FALSE)

### ————————————————————————————— Step 1: 3种方法count数据 ——————————————————————————-----———-----------———-—————-------------— ###
#（1）读取TPM 数据
HTSeqcount_count_file <- "HTSeqcount_labE/prepare_counts/HTSeqcounts_counts_filter1.csv"
featureCounts_count_file <- "featureCounts_labE/prepare_counts/featureCounts_counts_filter1.csv"

HTSeqcount_count<- as.data.frame(read_csv(HTSeqcount_count_file, col_names = TRUE, locale = locale(encoding = "UTF-8")))
featureCounts_count <- as.data.frame(read_csv(featureCounts_count_file, col_names = TRUE, locale = locale(encoding = "UTF-8")))

# (2) 检查基因数量
print(paste("HTSeqcount检测到的表达基因数量为:", nrow(HTSeqcount_count)))
print(paste("featureCounts检测到的表达基因数量为:", nrow(featureCounts_count)))

head(HTSeqcount_count)
head(featureCounts_count)
dim(featureCounts_count)
dim(HTSeqcount_count)

怎样取HTSeqcount_count和featureCounts_count中gene_id并集对应的每个样品的counts(最后4列)，并集共8列；计算并集8列counts间的person和spearman correlation，并画scatterplot


### ————————————————————————————— Step 3: 3种方法计算的TPM scaterpot —————————————————————————------———-—————-------------— ###

# 步骤1: 获取基因ID并集
union_genes <- union(HTSeqcount_count$gene_id, featureCounts_count$gene_id)

# 步骤2: 创建合并计数矩阵
combined_counts <- matrix(0, nrow = length(union_genes), ncol = 8)
rownames(combined_counts) <- union_genes
colnames(combined_counts) <- c(
  paste0("HTSeq.", colnames(HTSeqcount_count)[2:5]),
  paste0("featureCounts.", colnames(featureCounts_count)[3:6])
)

# 步骤3: 填充HTSeq数据
idx_htseq <- match(union_genes, HTSeqcount_count$gene_id)
htseq_vals <- as.matrix(HTSeqcount_count[idx_htseq, 2:5])
htseq_vals[is.na(htseq_vals)] <- 0
combined_counts[, 1:4] <- htseq_vals

# 步骤4: 填充featureCounts数据（跳过gene_length列）
idx_feature <- match(union_genes, featureCounts_count$gene_id)
feature_vals <- as.matrix(featureCounts_count[idx_feature, 3:6])
feature_vals[is.na(feature_vals)] <- 0
combined_counts[, 5:8] <- feature_vals

# 步骤5: 计算相关系数
pearson_cor <- cor(combined_counts, method = "pearson")
spearman_cor <- cor(combined_counts, method = "spearman")

# (4) 保存数据
write.csv(as.data.frame(pearson_cor), file = "3counts_camperation/HTSeqcounts_featireCounts_pearson_cor.csv", row.names = T)
write.csv(as.data.frame(spearman_cor), file = "3counts_camperation/HTSeqcounts_featireCounts_spearman_cor.csv", row.names = T)

### ------------------------------------------ log10 count scatterplot -------------------------------------- ###
# 步骤6: 绘制散点图（以rep1胞质样本为例）
# 安装并加载必要包
library(gridExtra)

# 获取基因ID并集
union_genes <- union(HTSeqcount_count$gene_id, featureCounts_count$gene_id)

# 创建合并计数矩阵
combined_counts <- matrix(0, nrow = length(union_genes), ncol = 8)
rownames(combined_counts) <- union_genes
colnames(combined_counts) <- c(
  paste0("HTSeq.", colnames(HTSeqcount_count)[2:5]),
  paste0("featureCounts.", colnames(featureCounts_count)[3:6])
)

# 填充HTSeq数据
idx_htseq <- match(union_genes, HTSeqcount_count$gene_id)
htseq_vals <- as.matrix(HTSeqcount_count[idx_htseq, 2:5])
htseq_vals[is.na(htseq_vals)] <- 0
combined_counts[, 1:4] <- htseq_vals

# 填充featureCounts数据（跳过gene_length列）
idx_feature <- match(union_genes, featureCounts_count$gene_id)
feature_vals <- as.matrix(featureCounts_count[idx_feature, 3:6])
feature_vals[is.na(feature_vals)] <- 0
combined_counts[, 5:8] <- feature_vals

# 计算相关系数
pearson_cor <- cor(combined_counts, method = "pearson")
spearman_cor <- cor(combined_counts, method = "spearman")

# 修正后的绘图函数 - 确保相关系数正确显示
create_scatter <- function(htseq_col, feature_col) {
  df <- data.frame(
    HTSeq = combined_counts[, htseq_col],
    featureCounts = combined_counts[, feature_col]
  )
  
  # 计算对数转换后的值（用于绘图）
  df$log_HTSeq <- log10(df$HTSeq + 1)
  df$log_featureCounts <- log10(df$featureCounts + 1)
  
  # 获取相关系数
  pearson_val <- pearson_cor[htseq_col, feature_col]
  spearman_val <- spearman_cor[htseq_col, feature_col]
  
  # 提取样本名称（去除前缀）
  sample_name <- gsub("HTSeq\\.|featureCounts\\.", "", htseq_col)
  
  ggplot(df, aes(x = log_HTSeq, y = log_featureCounts)) +
    geom_point(alpha = 0.3, size = 1, color = "steelblue") +
    geom_smooth(method = "lm", color = "red", se = FALSE) +
    labs(
      title = paste("Sample:", sample_name),
      x = "HTSeq-count (log10(count+1))",
      y = "featureCounts (log10(count+1))"
    ) +
    theme_bw() +
    # 添加相关系数标注 - 修正参数传递
    annotate("text", 
             x = min(df$log_HTSeq, na.rm = TRUE), 
             y = max(df$log_featureCounts, na.rm = TRUE),
             label = sprintf("Pearson r = %.3f\nSpearman ρ = %.3f", pearson_val, spearman_val),
             hjust = 0, vjust = 1, size = 4, 
             color = "darkred", fontface = "bold") +
    # 添加对角线参考线
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray50")
}

# 创建所有样本的散点图
plot_list <- list()
samples <- c("K562.LabE.polyA.cyt.rep1", "K562.LabE.polyA.cyt.rep2",
             "K562.LabE.polyA.nuc.rep1", "K562.LabE.polyA.nuc.rep2")

for (i in 1:4) {
  plot_list[[i]] <- create_scatter(
    htseq_col = paste0("HTSeq.", samples[i]),
    feature_col = paste0("featureCounts.", samples[i])
  )
}

# 组合并显示图形
log10count_sactterplot <- grid.arrange(
  grobs = plot_list,
  ncol = 2,
  top = "HTSeq-count vs featureCounts Comparison"
)

# 保存图形（PNG格式）
ggsave("3counts_camperation/HTSeq_featureCounts_log10count_sactterplot.png", log10count_sactterplot, width = 8, height = 8, dpi = 300)

### ------------------------------------------ raw count scatterplot ---------------------------------------------------------------- ###










# 保存图形（PNG格式）
ggsave("3counts_camperation/HTSeq_featureCounts_rawcount_sactterplot.png", log10count_sactterplot, width = 8, height = 8, dpi = 300)



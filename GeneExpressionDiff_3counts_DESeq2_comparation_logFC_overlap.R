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

### ————————————————————————————— Step 3: Venn图检测Sig_gene overlap———————————————---------------------——-—————---------— ###

# 提取每个数据表的 gene_id 列表
HTSeq_dif_genes <- HTSeqcount_dif$gene_id
featureCounts_dif_genes <- featureCounts_dif$gene_id
Salmon_dif_genes <- Salmon_dif$gene_id

# 创建命名列表
gene_dif_list3 <- list(
  HTSeqcount = HTSeq_dif_genes,
  featureCounts = featureCounts_dif_genes,
  Salmon = Salmon_dif_genes
)

# 绘制韦恩图
venn_plot_dif <- ggVennDiagram(
  gene_dif_list3,
  label_alpha = 0,  # 去除标签背景
  edge_size = 0.5   # 调整边框粗细
) +
  scale_fill_gradient(low = "white", high = "red") +  # 设置填充颜色
  labs(title = "Significant genes Overlap") +                   # 添加标题
  theme(plot.title = element_text(hjust = 0.5))       # 标题居中

# 显示图形
print(venn_plot_dif)

# 可选：保存图形（PNG格式）
ggsave("3counts_camperation/Sig_gene_overlap_venn.png", venn_plot_dif, width = 8, height = 6, dpi = 300)

### ————————————————————————————— Step 4: Venn图检测up_gene overlap———————————————---------------------——-—————---------— ###

# 提取每个数据表的 gene_id 列表
HTSeq_up_genes <- HTSeqcount_dif$gene_id[HTSeqcount_dif$log2FoldChange >1]
featureCounts_up_genes <- featureCounts_dif$gene_id[featureCounts_dif$log2FoldChange >1]
Salmon_up_genes <- Salmon_dif$gene_id[Salmon_dif$log2FoldChange >1]

# 创建命名列表
gene_up_list3 <- list(
  HTSeqcount = HTSeq_up_genes,
  featureCounts = featureCounts_up_genes,
  Salmon = Salmon_up_genes
)

# 绘制韦恩图
venn_plot_up <- ggVennDiagram(
  gene_up_list3,
  label_alpha = 0,  # 去除标签背景
  edge_size = 0.5   # 调整边框粗细
) +
  scale_fill_gradient(low = "white", high = "red") +  # 设置填充颜色
  labs(title = "Significant genes Overlap") +                   # 添加标题
  theme(plot.title = element_text(hjust = 0.5))       # 标题居中

# 显示图形
print(venn_plot_up)

# 可选：保存图形（PNG格式）
ggsave("3counts_camperation/up_gene_overlap_venn.png", venn_plot_up, width = 8, height = 6, dpi = 300)

### ————————————————————————————— Step 5: Venn图检测down_gene overlap———————————————---------------------——-—————---------— ###

# 提取每个数据表的 gene_id 列表
HTSeq_dn_genes <- HTSeqcount_dif$gene_id[HTSeqcount_dif$log2FoldChange < -1]
featureCounts_dn_genes <- featureCounts_dif$gene_id[featureCounts_dif$log2FoldChange < -1]
Salmon_dn_genes <- Salmon_dif$gene_id[Salmon_dif$log2FoldChange < -1]

# 创建命名列表
gene_dn_list3 <- list(
  HTSeqcount = HTSeq_dn_genes,
  featureCounts = featureCounts_dn_genes,
  Salmon = Salmon_dn_genes
)

# 绘制韦恩图
venn_plot_dn <- ggVennDiagram(
  gene_dn_list3,
  label_alpha = 0,  # 去除标签背景
  edge_size = 0.5   # 调整边框粗细
) +
  scale_fill_gradient(low = "white", high = "red") +  # 设置填充颜色
  labs(title = "Significant genes Overlap") +                   # 添加标题
  theme(plot.title = element_text(hjust = 0.5))       # 标题居中

# 显示图形
print(venn_plot_dn)

# 保存图形（PNG格式）
ggsave("3counts_camperation/dn_gene_overlap_venn.png", venn_plot_dn, width = 8, height = 6, dpi = 300)
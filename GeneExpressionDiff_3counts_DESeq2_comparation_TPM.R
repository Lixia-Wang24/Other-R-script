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

### ————————————————————————————— Step 1: 3种方法TPM数据 ——————————————————————————-----———-----------———-—————-------------— ###
#（1）读取TPM 数据
HTSeqcount_TPM_file <- "HTSeqcount_labE/prepare_counts/HTSeqcounts_gene_TPM_geneName_filter1.csv"
featureCounts_TPM_file <- "featureCounts_labE/prepare_counts/featureCounts_gene_TPM_geneName_filter1.csv"
Salmon_TPM_file <- "Salmon/Ensembl/prepare_counts/Salmon_gene_name_TPM_filter1.csv"

HTSeqcount_TPM_file <- as.data.frame(read_csv(HTSeqcount_TPM_file, col_names = TRUE, locale = locale(encoding = "UTF-8")))
featureCounts_TPM_file <- as.data.frame(read_csv(featureCounts_TPM_file, col_names = TRUE, locale = locale(encoding = "UTF-8")))
Salmon_TPM_file <- as.data.frame(read_csv(Salmon_TPM_file, col_names = TRUE, locale = locale(encoding = "UTF-8")))

# (2) 检查基因数量
print(paste("HTSeqcount检测到的表达基因数量为:", nrow(HTSeqcount_TPM_file)))
print(paste("featureCounts检测到的表达基因数量为:", nrow(featureCounts_TPM_file)))
print(paste("Salmon检测到的表达基因数量为:", nrow(Salmon_TPM_file)))

### ————————————————————————————— Step 2: 3种方法检测到的表达基因overlap —————————————————————————------———-—————-------------— ###

# 提取每个数据表的 gene_id 列表
HTSeq_genes <- HTSeqcount_TPM_file$gene_id
featureCounts_genes <- featureCounts_TPM_file$gene_id
Salmon_genes <- Salmon_TPM_file$gene_id

# 创建命名列表
gene_list3 <- list(
  HTSeqcount = HTSeq_genes,
  featureCounts = featureCounts_genes,
  Salmon = Salmon_genes
)

# 绘制韦恩图
venn_plot_TPM <- ggVennDiagram(
  gene_list3,
  label_alpha = 0,  # 去除标签背景
  edge_size = 0.5   # 调整边框粗细
) +
  scale_fill_gradient(low = "white", high = "red") +  # 设置填充颜色
  labs(title = "Gene TPM Overlap") +                   # 添加标题
  theme(plot.title = element_text(hjust = 0.5))       # 标题居中

# 显示图形
print(venn_plot_TPM)

# 可选：保存图形（PNG格式）
ggsave("3counts_camperation/TPM_gene_overlap_venn.png", venn_plot_TPM, width = 8, height = 6, dpi = 300)

### ————————————————————————————— Step 3: 3种方法计算的TPM scaterpot —————————————————————————------———-—————-------------— ###

HTSeqcount_TPM_file <- as.data.frame(read_csv(HTSeqcount_TPM_file, col_names = TRUE, locale = locale(encoding = "UTF-8")))
featureCounts_TPM_file <- as.data.frame(read_csv(featureCounts_TPM_file, col_names = TRUE, locale = locale(encoding = "UTF-8")))
Salmon_TPM_file <- as.data.frame(read_csv(Salmon_TPM_file, col_names = TRUE, locale = locale(encoding = "UTF-8")))


# (1) 提取共同基因ID
common_genes <- intersect(
  intersect(HTSeqcount_TPM_file$gene_id, featureCounts_TPM_file$gene_id),
  Salmon_TPM_file$gene_id
)

# 3. 为每个数据集创建统一格式的长表
prepare_long_data <- function(df, tool_name) {
  df %>%
    filter(gene_id %in% common_genes) %>%
    pivot_longer(
      cols = starts_with("TPM_") | starts_with("K562."),
      names_to = "sample",
      values_to = "TPM"
    ) %>%
    mutate(
      sample = gsub("TPM_|K562\\.", "", sample),  # 统一样本名格式
      tool = tool_name
    ) %>%
    select(gene_id, sample, TPM, tool)
}

# 4. 合并三个数据集
combined_tpm <- bind_rows(
  prepare_long_data(HTSeqcount_TPM_file, "HTSeqcount"),
  prepare_long_data(featureCounts_TPM_file, "featureCounts"),
  prepare_long_data(Salmon_TPM_file, "Salmon")
)

# 5. 转换为宽格式（每个工具+样本作为一列）
wide_tpm <- combined_tpm %>%
  unite(col = "tool_sample", tool, sample) %>%
  pivot_wider(
    names_from = tool_sample,
    values_from = TPM
  ) %>%
  select(-gene_id)  # 移除基因ID列

# 6. 计算相关系数矩阵
cor_matrix_spearman <- cor(wide_tpm, method = "spearman", use = "complete.obs")
cor_matrix_pearson <- cor(wide_tpm, method = "pearson", use = "complete.obs")

# 7. 绘制热图展示相关性

Spearman_heatmap <- pheatmap(cor_matrix_spearman, 
         main = "Spearman Correlation of TPM Values",
         display_numbers = TRUE)

Pearson_heatmap <- pheatmap(cor_matrix_pearson, 
         main = "Pearson Correlation of TPM Values",
         display_numbers = TRUE)

# 保存图形（PNG格式）
ggsave("3counts_camperation/gene_TPM_Spearman_heatmap.png", Spearman_heatmap, width = 8, height = 6, dpi = 300)
ggsave("3counts_camperation/gene_TPM_Pearson_heatmap.png", Pearson_heatmap, width = 8, height = 6, dpi = 300)


# 8. 绘制散点图矩阵（选择代表性样本避免过载）
# 重点比较：相同样本不同工具 vs 相同工具不同重复
selected_cols <- c(
  "HTSeqcount_LabE.polyA.cyt.rep1",
  "featureCounts_LabE.polyA.cyt.rep1",
  "Salmon_LabE.polyA.cyt.rep1",
  "HTSeqcount_LabE.polyA.nuc.rep1",
  "featureCounts_LabE.polyA.nuc.rep1",
  "Salmon_LabE.polyA.nuc.rep1"
)

# 开始画图
TPM_spearman <- ggpairs(
  data = log2(wide_tpm[, selected_cols] + 1),  # log2(TPM+1) 转换
  title = "TPM Comparison (log2 scale)",
  lower = list(continuous = wrap("smooth", alpha = 0.3)),
  diag = list(continuous = "densityDiag"),
  upper = list(continuous = wrap("cor", method = "spearman", size = 4))
) + 
  theme_minimal() +
theme(
  # 修改轴标签（坐标轴上的数字）
  axis.text = element_text(size = 10),  # 坐标轴刻度文字大小
  # 修改变量名标签（对角线上的名称）
  strip.text = element_text(
    size = 5,  # 增大标签字体
    face = "bold",  # 加粗
    margin = margin(5, 0, 5, 0)  # 增加上下边距
  ),
  # 修改标题
  plot.title = element_text(
    size = 14,  # 标题字体大小
    face = "bold",  # 加粗
    hjust = 0.5  # 居中
  )
)

# 保存图形（PNG格式）
ggsave("3counts_camperation/gene_TPM_scaterplot_spearman.png", TPM_spearman, width = 8, height = 8, dpi = 300)


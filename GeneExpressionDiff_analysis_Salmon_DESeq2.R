# 加载必要的包
library(dplyr)
library(readr)
library(tximport)
library(DESeq2)

# 设置工作目录
setwd("/Users/lixia/Data/data/Fractionation_seq/Salmon/Ensembl")
# 创建结果保存目录
dir.create("DESeq2_results", showWarnings = FALSE)

### ---------------------------- Step1: 读取并查看Salmoun transcript-level count数据 ------------------------------------------------------- ###

#(1) 备样本信息
samples <- c("K562.LabE.polyA.cyt.rep1", "K562.LabE.polyA.cyt.rep2", 
             "K562.LabE.polyA.nuc.rep1", "K562.LabE.polyA.nuc.rep2")
# 创建样本信息表
sample_info <- data.frame(
  sample = samples,
  condition = c("cytoplasm", "cytoplasm", "nucleus", "nucleus"),
  replicate = c("rep1", "rep2", "rep1", "rep2"),
  row.names = samples,
  stringsAsFactors = FALSE)

# (2) 准备salmon文件路径
salmon_files <- file.path(samples, "quant.sf")
names(salmon_files) <- samples  #需要

# 检查文件是否存在
if(!all(file.exists(salmon_files))) {
  stop("部分quant.sf文件不存在，请检查文件路径")
}

### ---------------------------- Step2: 读取Salmoun gene-level 数据 -------------------------------------------------------------------- ###

# (3) 导入tx2gene映射表
tx2gene_file <- "/Users/lixia/Data/database/ref_RNA/Ensembl/GeneID_GeneName_Transcript_mapping_table/gene_id_name_tx_id_map_Ensembl_rtracklayer.csv"
if(!file.exists(tx2gene_file)) {
  stop("tx2gene.tsv文件不存在，请检查路径")
}

# gene_id, gene_name, transcript_id
tx2gene_ano <- as.data.frame(read_csv(tx2gene_file, col_names = TRUE, locale = locale(encoding = "UTF-8")))
# transcript_id, gene_id,
tx2gene <- tx2gene_ano[,c("transcript_id", "gene_id")]

# (4) 导入salmon定量数据并汇总到基因水平
txi <- tximport(salmon_files, 
                type = "salmon", 
                tx2gene = tx2gene,
                countsFromAbundance = "lengthScaledTPM") # DESeq2推荐参数

# gene_count <- as.data.frame(txi$counts) 
# gene_count$gene_id <- row.names(gene_count)

### ---------------------------- Step3: DESeq2基因表达差异分析------------------------------------------------------------------------- ###
# (1) 设置参考水平（细胞质作为参考）
sample_info$condition <- factor(sample_info$condition, levels = c("cytoplasm", "nucleus"))

# (2) 创建DESeq2对象
dds <- DESeqDataSetFromTximport(txi, colData = sample_info, design = ~ condition)

# (3) 过滤低表达基因（至少在2个样本中计数>=10）
# keep <- rowSums(counts(dds) >= 10) >= 2
# dds <- dds[keep,]
# print(paste("过滤后保留基因数量:", nrow(dds)))

# (4) 运行DESeq2分析
dds <- DESeq(dds)

# (5) 查看结果名称
resultsNames(dds)

### ————————————————————————————— Step4: 提取差异表达结果、gene counts  ——————————————————————————————------------------------———-—————— ###

# (1) 提取表达差异 nucleus vs cytoplasm
res <- results(dds, contrast = c("condition", "nucleus", "cytoplasm"))
res$gene_id <- rownames(res) #必须标记gene_id，因为后续merge会改变res(dds)的rowname，重新按照numeric排行号
res_ged <- as.data.frame(res)

# (2) 提取gene counts: 分为3类，raw_counts；rlog，vst。rlog最精准，此处选rlog。
rld <- rlogTransformation(dds, blind = FALSE)  #FALSE考虑实验设计信息
rlog_matrix <- as.data.frame(assay(rld)) #assay(rld)是matrix,不能给rowname gene_id
rlog_matrix$gene_id <- rownames(rlog_matrix)
colnames(rlog_matrix) <- c("counts.K562.LabE.polyA.cyt.rep1", "counts.K562.LabE.polyA.cyt.rep2",
                           "counts.K562.LabE.polyA.nuc.rep1", "counts.K562.LabE.polyA.nuc.rep2",
                           "gene_id")

# (3) 合并差异表达结果和gene counts
res_ged_counts <- merge(res_ged, rlog_matrix, by="gene_id", sort=FALSE) #FALSE表示合并后的数据框不会按照by进行排序
rownames(res_ged_counts) <- res_ged_counts$gene_id ##画图需要

# (4) 保存差异分析结果
write.csv(res_ged_counts , file = "DESeq2_results/Salmon_DESeq2_res_ged_counts.csv", row.names = FALSE)

### ————————————————————————————— Step5: 注释gene_name, TPM  ——————————————————————————————---------------------------------———-—————— ###

# (1) 读取注释数据
# gene_id, gene_name, transcript_id
anno_file <- "/Users/lixia/Data/data/Fractionation_seq/Salmon/Ensembl/prepare_counts/Salmon_gene_name_TPM_filter1.csv"
tx2gene_ano <- read_csv(anno_file, col_names = TRUE, locale = locale(encoding = "UTF-8"))
colnames(tx2gene_ano) <- c("gene_id", "gene_name", "TPM_K562.LabE.polyA.cyt.rep1", "TPM_K562.LabE.polyA.cyt.rep2", 
                           "TPM_K562.LabE.polyA.nuc.rep1", "TPM_K562.LabE.polyA.nuc.rep2")

# (2) 通过gene_id 注释gene_name；若无对应gene_name,则用gene_id代替gene_name
res_ged_counts_geneName_TPM <- res_ged_counts %>%
  left_join(
    tx2gene_ano %>% distinct(gene_id, .keep_all = TRUE), 
    by = "gene_id"
  ) %>%
  mutate(gene_name = coalesce(gene_name, gene_id)) %>% # 若gene_name有值，保留原值；若 gene_name缺失，则使用gene_id
  relocate(gene_id, gene_name) 

# (3) 检查多少gene_id注释了gene_name
success_anno <- sum(res_ged_counts_geneName_TPM$gene_name != res_ged_counts_geneName_TPM$gene_id, na.rm = TRUE)
print(paste(nrow(res_ged_counts),"个gene注释了", sum(success_anno), "个gene_name"))

### ————————————————————————————— Step 6: 保存差异分析和注释结果 ——————————————————————————————---------------------------------———-—————— ###

# (1) 按照padj值排序,headmap画图需要
res_ged_counts_geneName_TPM_Allgene <- res_ged_counts_geneName_TPM[order(res_ged_counts_geneName_TPM$padj),]

# (2) 保存所有基因文件
write.csv(res_ged_counts_geneName_TPM_Allgene, 
          file = "DESeq2_results/Salmon_DESeq2_res_ged_counts_geneName_TPM_Allgene.csv",
          row.names = FALSE)

# (3) 保存显著差异基因文件
res_ged_counts_geneName_TPM_Siggene <- subset(res_ged_counts_geneName_TPM, padj < 0.05 & abs(log2FoldChange) > 1)
write.csv(as.data.frame(res_ged_counts_geneName_TPM_Siggene), 
          file = "DESeq2_results/Salmon_DESeq2_res_ged_counts_geneName_TPM_Siggene.csv",
          row.names = FALSE)

# (4) 查看差异基因
print(paste("显著差异表达基因数量:", nrow(res_ged_counts_geneName_TPM_Siggene)))
print(paste("上调基因数量:", sum(res_ged_counts_geneName_TPM_Siggene$log2FoldChange > 1)))
print(paste("下调基因数量:", sum(res_ged_counts_geneName_TPM_Siggene$log2FoldChange < -1)))

###------------------------------------ Step 7: 数据质量控制图  -------------------------------------------------------------------------- ###
library(ggplot2)
library(RColorBrewer)
library(pheatmap)
library(ggrepel)

# 12. 数据质量控制图
rld <- rlogTransformation(dds, blind = FALSE)

# 12.1 样本距离热图
sampleDists <- dist(t(assay(rld)))
sampleDistMatrix <- as.matrix(sampleDists)

pdf("DESeq2_results/sample_distance_heatmap.pdf", width = 8, height = 6)
pheatmap(sampleDistMatrix,
         clustering_distance_rows = sampleDists,
         clustering_distance_cols = sampleDists,
         col = colorRampPalette(c("blue", "white", "red"))(100),
         annotation_col = sample_info[, "condition", drop = FALSE])
dev.off()

# 12.2 PCA图
pdf("DESeq2_results/PCA_plot.pdf", width = 8, height = 6)
plotPCA(rld, intgroup = "condition", returnData = FALSE) +
  geom_text_repel(aes(label = name), size = 3) +
  theme_bw() +
  ggtitle("PCA Plot of Samples")
dev.off()

# 12.3 MA图
res <- results(dds, contrast = c("condition", "nucleus", "cytoplasm"))
pdf("DESeq2_results/MA_plot.pdf", width = 8, height = 6)
plotMA(res, ylim = c(-5, 5), main = "MA Plot: Nucleus vs Cytoplasm")
dev.off()

###------------------------------------ valcona map 火山图 (13) ----------------------------------------------------------------------- ###
# DESeq2_valcano.R
###------------------------------------ heatmap 热图 (14) ----------------------------------------------------------------------------- ###
# heatmap_label_star.R
# heatmap_label_specific_genes.R
# heatmap_label_specific_genes_star.R
# heatmap_label_significant_genes.R

# 从GENECODE 的注释文件gencode.v48.chr_patch_hapl_scaff.annotation.gtf.gz。
# 怎样提取最长转录本对应的基因长度(exon only), 用于FPKM,TPM计算？

# 设置工作目录
setwd("/Users/lixia/Data/database/ref_RNA/GENCODE")
# 创建结果保存目录
dir.create("Gene_length", showWarnings = FALSE)

# 加载必要包 - 增加GenomicFeatures用于高效计算
library(rtracklayer)
library(GenomicRanges)
library(GenomicFeatures)  
library(dplyr)

# 1. 导入GTF并创建TxDb对象（一次性处理）
txdb <- makeTxDbFromGFF("gencode.v48.chr_patch_hapl_scaff.annotation.gtf.gz")
exons_by_tx <- exonsBy(txdb, by = "tx", use.names = TRUE)  # 按转录本分组exons
exons_by_tx2 <- as.data.frame(exons_by_tx) 

# 2. 计算转录本长度（向量化操作）
reduced_exons <- reduce(exons_by_tx)  # 合并重叠exon
transcript_lengths <- sum(width(reduced_exons))  # 向量化计算长度

# 3. 构建结果数据框
transcript_df <- data.frame(
  transcript_id = names(transcript_lengths),
  transcript_length = transcript_lengths
)

# 4. 获取基因ID与基因名的对应关系（包含版本号）
# 从GTF文件直接解析基因信息, 耗时1min
gtf <- import("gencode.v48.chr_patch_hapl_scaff.annotation.gtf.gz")
transcript_info <- gtf[gtf$type == "transcript"]
tx_to_gene <- data.frame(
  gene_id = transcript_info$gene_id,
  gene_name = transcript_info$gene_name,
  transcript_id = transcript_info$transcript_id,
  stringsAsFactors = FALSE
)

# 合并长度数据
tx_gene_length <- merge(transcript_df, tx_to_gene, by = "transcript_id")

# 5. 找到每个基因的最长转录本
gene_lengths <- tx_gene_length %>%
  group_by(gene_id) %>% 
  summarise(
    longest_tx = transcript_id[which.max(transcript_length)],
    gene_name = gene_name[which.max(transcript_length)],
    gene_length = max(transcript_length),
    .groups = 'drop'
  )

gene_lengths <-  as.data.frame(gene_lengths[, c("gene_id", "gene_name", "gene_length")])

# 6. 查看统计信息
print(paste("总共处理了", nrow(gene_lengths), "个基因"))
print(paste("基因长度范围：", min(gene_lengths$gene_length), "-", max(gene_lengths$gene_length)))
print(paste("平均基因长度：", round(mean(gene_lengths$gene_length), 2)))

# 7. 保存结果
write.csv(gene_lengths, "Gene_length/gencode.v48.gene_longest_transcript_lengths_GenomicFeatures_rtracklayer.csv", row.names = FALSE)
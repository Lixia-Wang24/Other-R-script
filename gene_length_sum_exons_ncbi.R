# 从NCBI 的注释文件GCF_000001405.40_GRCh38.p14_genomic.gtf.gz。
# 怎样提取最长转录本对应的基因长度(exon only), 用于FPKM,TPM计算？

# 设置工作目录
setwd("/Users/lixia/Data/database/ref_RNA/NCBI")
# 创建结果保存目录
dir.create("Gene_length", showWarnings = FALSE)

# 加载必要包 - 增加GenomicFeatures用于高效计算
library(rtracklayer)
library(GenomicRanges)
library(GenomicFeatures)  
library(dplyr)

# 1. 从GTF文件创建TxDb对象
gtf_file <- "GCF_000001405.40_GRCh38.p14_genomic.gtf.gz"
txdb <- makeTxDbFromGFF(gtf_file, format="gtf")

# 计算基因有效长度
exons_by_genes <- exonsBy(txdb, by = "gene")
gene_lengths <- sum(width(reduce(exons_by_genes)))


# 2. 获取基因的外显子信息
exons_by_gene <- exonsBy(txdb, by="gene")
# 3. 计算每个基因的非冗余外显子长度, 耗时5min
gene_effective_lengths <- sapply(exons_by_gene, function(gene_exons) { 
  reduced_exons <- reduce(gene_exons)  # 合并重叠的外显子区间
  sum(width(reduced_exons))  # 计算总长度
})




# 从GENECODE 的注释文件gencode.v48.chr_patch_hapl_scaff.annotation.gtf.gz。
# 怎样获取基因有效长度（非冗余外显子长度之和）用于fpkm，tpm计算？

# 设置工作目录
setwd("/Users/lixia/Data/database/ref_RNA/GENCODE")
# 创建结果保存目录
dir.create("Gene_length", showWarnings = FALSE)

library(rtracklayer)
library(GenomicRanges)
library(dplyr)

# 读取GTF文件
gtf_file <- "gencode.v48.chr_patch_hapl_scaff.annotation.gtf.gz"
gtf <- import(gtf_file)

# 提取外显子信息
exons <- gtf[gtf$type == "exon"]

# 设置“计算每个基因的有效长度（非冗余外显子长度之和)"的函数
calculate_gene_effective_length <- function(exons) {
  exons_by_gene <- split(exons, exons$gene_id)                  # 按基因ID分组
  gene_lengths <- sapply(exons_by_gene, function(gene_exons) {  # 对每个基因计算非冗余外显子长度
  merged_exons <- reduce(gene_exons)                            # 将同一基因的所有外显子合并，去除重叠
  sum(width(merged_exons))                                      # 计算总长度
  })
  return(gene_lengths)
}

# 计算基因有效长度
gene_effective_lengths <- calculate_gene_effective_length(exons)

# 转换为数据框，包含基因ID（带版本号）和有效长度
gene_length_df <- data.frame(
  gene_id = names(gene_effective_lengths),
  effective_length = as.numeric(gene_effective_lengths),
  stringsAsFactors = FALSE
)

# 查看结果
head(gene_length_df)

# 保存结果
write.table(gene_length_df, 
            file = "gene_effective_lengths.txt", 
            sep = "\t", 
            row.names = FALSE, 
            quote = FALSE)

# 输出统计信息
cat("总基因数:", nrow(gene_length_df), "\n")
cat("平均有效长度:", round(mean(gene_length_df$effective_length), 2), "bp\n")
cat("有效长度范围:", min(gene_length_df$effective_length), "-", 
    max(gene_length_df$effective_length), "bp\n")

# 可选：创建一个命名向量用于FPKM计算
gene_lengths_vector <- setNames(gene_length_df$effective_length, 
                                gene_length_df$gene_id)

# 示例：如何在FPKM计算中使用
# fpkm_value <- (raw_count * 1e9) / (gene_lengths_vector[gene_id] * total_reads)
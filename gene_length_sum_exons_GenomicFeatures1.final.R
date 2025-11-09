# 从GENECODE 的注释文件gencode.v48.chr_patch_hapl_scaff.annotation.gtf.gz。
# 怎样获取基因有效长度（非冗余外显子长度之和）用于fpkm，tpm计算？

# 设置工作目录
setwd("/Users/lixia/Data/database/ref_RNA/GENCODE")
# 创建结果保存目录
dir.create("Gene_length", showWarnings = FALSE)

# 计算基因有效长度（非冗余外显子长度之和）
# 用于FPKM计算

library(GenomicFeatures)
library(rtracklayer)

# 1. 从GTF文件创建TxDb对象
gtf_file <- "gencode.v48.chr_patch_hapl_scaff.annotation.gtf.gz"
txdb <- makeTxDbFromGFF(gtf_file, format="gtf")

# 2. 获取基因的外显子信息
exons_by_gene <- exonsBy(txdb, by="gene")

exons_by_gene2 <- as.data.frame(exons_by_gene)
head(exons_by_gene2)
# group         group_name seqnames     start       end width strand exon_id         exon_name
# 1     1 ENSG00000000003.16     chrX 100627108 100629986  2879      - 1017954 ENSE00001459322.5

# 3. 计算每个基因的非冗余外显子长度, 耗时5min
gene_effective_lengths <- sapply(exons_by_gene, function(gene_exons) { 
  reduced_exons <- reduce(gene_exons)  # 合并重叠的外显子区间
  sum(width(reduced_exons))  # 计算总长度
})

# 4. 整合结果
gene_lengths <- data.frame(
  gene_id = names(gene_effective_lengths),
  gene_length = gene_effective_lengths,
  stringsAsFactors = FALSE
)

# 5. 获取基因ID与基因名的对应关系（包含版本号）
# 从GTF文件直接解析基因信息, 耗时1min
gtf <- import(gtf_file)
gene_info <- gtf[gtf$type == "gene"]
gene_mapping <- data.frame(
  gene_id = gene_info$gene_id,
  gene_name = gene_info$gene_name,
  stringsAsFactors = FALSE
)

# 添加基因名信息
gene_lengths <- merge(gene_lengths, gene_mapping, 
                            by = "gene_id", all.x = TRUE)

# 6. 查看统计信息
print(paste("总共处理了", nrow(gene_lengths), "个基因"))
print(paste("基因长度范围：", min(gene_lengths$gene_length), "-", max(gene_lengths$gene_length)))
print(paste("平均基因长度：", round(mean(gene_lengths$gene_length), 2)))

# 7. 保存结果
write.csv(gene_lengths, "Gene_length/gencode.v48.gene_length_sum_exons_GenomicFeatures_1.csv", row.names = FALSE)

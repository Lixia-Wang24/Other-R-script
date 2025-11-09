# 从GENECODE 的注释文件gencode.v48.chr_patch_hapl_scaff.annotation.gtf.gz。
# 怎样提取最长转录本对应的基因长度(exon only), 用于FPKM,TPM计算？

# 设置工作目录
setwd("/Users/lixia/Data/database/ref_RNA/GENCODE")
# 创建结果保存目录
dir.create("Gene_length", showWarnings = FALSE)

# 加载必要的包
library(GenomicFeatures)
library(rtracklayer)
library(dplyr)

# 1. 读取GTF文件并创建TxDb对象
gtf_file <- "gencode.v48.chr_patch_hapl_scaff.annotation.gtf.gz"
txdb <- makeTxDbFromGFF(gtf_file, format="gtf")
#txdb <- makeTxDbFromGFF("gencode.v48.chr_patch_hapl_scaff.annotation.gtf.gz", format="gtf")
# 查看具体字段内容
columns(txdb)
#  "CDSCHROM"   "CDSEND"     "CDSID"      "CDSNAME"    "CDSPHASE"   "CDSSTART"   "CDSSTRAND"  
#  "EXONCHROM"  "EXONEND"    "EXONID"     "EXONNAME"   "EXONRANK"  "EXONSTART"  "EXONSTRAND"
#  "TXCHROM"    "TXEND"      "TXID"   "TXNAME"     "TXSTART"    "TXSTRAND"   "TXTYPE"    
#  "GENEID"     

# 2. 计算每个转录本的长度（包括所有外显子）
# 获取外显子按转录本分组
exons_by_tx <- exonsBy(txdb, by="tx", use.names=TRUE)
exons_by_tx2 <- as.data.frame(exons_by_tx)

# 计算每个转录本的总长度
tx_lengths <- sum(width(exons_by_tx))
tx_lengths <- data.frame(
  transcript_id = names(tx_lengths),
  transcript_length = as.numeric(tx_lengths)
)

# 3. 获取转录本到基因的映射关系（包含版本号）
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
tx_gene_length <- merge(tx_lengths, tx_to_gene, by = "transcript_id")

# 4. 找到每个基因的最长转录本
gene_lengths <- tx_gene_length %>%
  group_by(gene_id) %>% 
  summarise(
    longest_tx = transcript_id[which.max(transcript_length)],
    gene_name = gene_name[which.max(transcript_length)],
    gene_length = max(transcript_length),
    .groups = 'drop'
  )

gene_lengths <-  as.data.frame(gene_lengths[, c("gene_id", "gene_name", "gene_length")])

# 5. 查看统计信息
print(paste("总共处理了", nrow(gene_lengths), "个基因"))
print(paste("基因长度范围：", min(gene_lengths$gene_length), "-", max(gene_lengths$gene_length)))
print(paste("平均基因长度：", round(mean(gene_lengths$gene_length), 2)))

# 6. 保存结果
write.csv(gene_lengths, "Gene_length/gencode.v48.gene_longest_transcript_lengths_GenomicFeatures.csv", row.names = FALSE)
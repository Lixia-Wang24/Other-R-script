# 从GENECODE 的注释文件gencode.v48.chr_patch_hapl_scaff.annotation.gtf.gz。
# 怎样提取最长转录本对应的基因长度(exon only), 用于FPKM,TPM计算？

# 设置工作目录
setwd("/Users/lixia/Data/database/ref_RNA/GENCODE")
# 创建结果保存目录
dir.create("Gene_length", showWarnings = FALSE)

# 加载必要包
library(rtracklayer)  
library(GenomicRanges)
library(dplyr)

# 导入GTF文件
#import(): 来自rtracklayer包，自动解压并导入GTF/GFF文件为GRanges对象（基因组坐标的容器)
gtf <- import("gencode.v48.chr_patch_hapl_scaff.annotation.gtf.gz")
# 查看GTF文件结构
print(head(gtf))

# 2. 筛选出exon区域（用于计算转录本长度）
exons <- gtf[gtf$type == "exon"] # [ ]子集化: 基于type列筛选

# 查看exon数据结构
print("Exon data structure:")
print(head(mcols(exons))) # mcols(): 提取GRanges对象的元数据列（metadata columns）

# 3. 按转录本分组计算每个转录本的总长度 耗时30-40min
# 使用GenomicRanges包计算
transcript_lengths <- exons %>%                # 管道(%>%): 来自dplyr包，传递对象
  split(., .$transcript_id) %>%                # split(): 按transcript_id分组创建列表
  lapply(function(x) {
    reduced_exons <- reduce(x)                 # reduce(): 合并重叠/相邻的exon区域
    total_length <- sum(width(reduced_exons))  # sum(width()): 计算合并后的总长度
    return(data.frame(
      transcript_id = unique(x$transcript_id),
      gene_id = unique(x$gene_id),
      gene_name = unique(x$gene_name),
      transcript_length = total_length
    ))
  }) %>%
  do.call(rbind, .)  # 合并所有数据框为单个数据框

# 4. 为每个基因找到最长的转录本
longest_transcripts <- transcript_lengths %>%
  group_by(gene_id) %>%                                        # group_by(): 按基因ID分组
  slice_max(transcript_length, n = 1, with_ties = FALSE) %>%   # slice_max(): 取每组中长度最大值（n=1取一条，with_ties=FALSE随机处理并列情况）
  ungroup() # `dplyr`包中的一个函数，用于移除数据框的分组结构，停止分组操作，将数据框恢复到未分组状态

# 5. 创建最终的基因长度表，用于FPKM计算
gene_lengths <- data.frame(
  gene_id = longest_transcripts$gene_id,
  gene_name = longest_transcripts$gene_name,
  gene_length = longest_transcripts$transcript_length,
  stringsAsFactors = FALSE
)

# 查看统计信息
print(paste("总共处理了", nrow(gene_lengths), "个基因"))
print(paste("基因长度范围：", min(gene_lengths$gene_length), "-", max(gene_lengths$gene_length)))
print(paste("平均基因长度：", round(mean(gene_lengths$gene_length), 2)))

# 6. 保存结果
write.csv(gene_lengths, "Gene_length/gencode.v48.gene_longest_transcript_lengths_rtracklayer_timecosted.csv", row.names = FALSE)

# 如果有并列最长的转录本，上面的代码会随机选择一个
# 如果你想要更明确的选择策略，可以使用以下代码：
longest_transcripts_explicit <- transcript_lengths %>%
  group_by(gene_id) %>%
  filter(transcript_length == max(transcript_length)) %>%
  slice(1) %>%  # 如果有并列，选择第一个
  ungroup()

# 8. 可选：如果你需要处理基因的不同版本号
# 有时基因ID可能包含版本号（如ENSG00000290825.2），你可能需要去掉版本号
gene_lengths$gene_id_no_version <- gsub("\\.\\d+$", "", gene_lengths$gene_id) # gsub(): 正则表达式移除版本号（如ENSG000001.2 → ENSG000001）


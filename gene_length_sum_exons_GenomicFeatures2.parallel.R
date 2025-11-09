
# 设置工作目录
setwd("/Users/lixia/Data/database/ref_RNA/GENCODE")
# 创建结果保存目录
dir.create("Gene_length", showWarnings = FALSE)

library(parallel) #并行计算  parApply parLapply parSaplly 
cl <- makeCluster(0.75*detectCores())  #设计启用计算机3/4的核
# 在工作节点加载必要的包
clusterEvalQ(cl, {
  library(GenomicRanges)  # 提供 width() 和 reduce() 函数
  library(GenomicFeatures)
})

## 利用GenomicFeatures包导入gtf处理
library(GenomicFeatures)
txdb <- makeTxDbFromGFF("gencode.v48.chr_patch_hapl_scaff.annotation.gtf.gz",
                        format="gtf") 
exons_gene <- exonsBy(txdb, by = "gene") ###提取基因外显子
head(exons_gene)

##计算总外显子长度：用reduce去除掉重叠冗余的部分，,width统计长度，最后计算总长度
exons_gene_lens <- parLapply(cl,exons_gene,function(x){sum(width(reduce(x)))}) 
exons_gene_lens[1:10]

# 转换为数据框
gene_lengths <- data.frame(
  gene_id = names(exons_gene_lens),
  gene_length = as.numeric(exons_gene_lens)
)

# 停止集群
stopCluster(cl)

# 查看统计信息
print(paste("总共处理了", nrow(gene_lengths), "个基因"))
print(paste("基因长度范围：", min(gene_lengths$gene_length), "-", max(gene_lengths$gene_length)))
print(paste("平均基因长度：", round(mean(gene_lengths$gene_length), 2)))

# 6. 保存结果
write.csv(gene_lengths, "Gene_length/gencode.v48.gene_length_sum_exons_GenomicFeatures_2.csv", row.names = FALSE)
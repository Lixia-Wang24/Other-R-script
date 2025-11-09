# 载入必要的包
library(limma)
library(edgeR)
library(readr)
library(dplyr)

# 你已经完成的步骤：
# 设置工作目录
setwd("/Users/lixia/Data/data/Fractionation_seq/HTSeqcount_labE")

# 创建结果保存目录
dir.create("limma_results", showWarnings = FALSE)

# 读取count数据
count_file <- "prepare_counts/HTSeqcounts_counts_filter1.csv"
counts_filtered <- read_csv(count_file, locale = locale(encoding = "UTF-8"))
counts_matrax <- as.data.frame(counts_filtered[, 2:5]) 
rownames(counts_matrax) <- counts_filtered$gene_id

# 读取样本信息
sample_info <- as.data.frame(read_csv("prepare_counts/sample_info.csv", locale = locale(encoding = "UTF-8")))
rownames(sample_info) <- sample_info$sample 
sample_info$condition <- factor(sample_info$condition, levels = c("cytoplasm", "nucleus"))

# 对齐样本顺序
counts_matrax <- counts_matrax[, rownames(sample_info)]

# ==================== Limma分析步骤 ====================

# 第4步：创建DGEList对象（用于标准化）
dge <- DGEList(counts = counts_matrax)

# 第5步：过滤低表达基因
# 保留在至少2个样本中表达量 > 1的基因
keep <- filterByExpr(dge, min.count = 1, min.total.count = 2)
dge <- dge[keep, , keep.lib.sizes = FALSE]

# 第6步：标准化（TMM标准化）
dge <- calcNormFactors(dge, method = "TMM")

# 第7步：转换为log2-CPM值用于limma分析
logCPM <- cpm(dge, log = TRUE)

# 第8步：创建设计矩阵
design <- model.matrix(~ condition, data = sample_info)
colnames(design) <- c("Intercept", "nucleus_vs_cytoplasm")
print("设计矩阵:")
print(design)

# 第9步：拟合线性模型
fit <- lmFit(logCPM, design)

# 第10步：经验贝叶斯统计
fit <- eBayes(fit)

# 第11步：获取差异表达结果
# 获取nucleus vs cytoplasm的比较结果
results <- topTable(fit, coef = "nucleus_vs_cytoplasm", 
                    number = Inf, 
                    sort.by = "P")

# 添加基因ID列
results$gene_id <- rownames(results)
results <- results[, c("gene_id", "logFC", "AveExpr", "t", "P.Value", "adj.P.Val", "B")]

# 第12步：筛选差异表达基因
# 设置阈值：|logFC| > 1 且 adj.P.Val < 0.05
deg_threshold <- list(
  logFC_cutoff = 1,
  pvalue_cutoff = 0.05
)

# 标记差异表达基因
results$significant <- ifelse(abs(results$logFC) > deg_threshold$logFC_cutoff & 
                                results$adj.P.Val < deg_threshold$pvalue_cutoff, 
                              "YES", "NO")

results$regulation <- ifelse(results$logFC > deg_threshold$logFC_cutoff & 
                               results$adj.P.Val < deg_threshold$pvalue_cutoff, "Up",
                             ifelse(results$logFC < -deg_threshold$logFC_cutoff & 
                                      results$adj.P.Val < deg_threshold$pvalue_cutoff, "Down", "Not_significant"))

# 第13步：统计结果
cat("差异表达基因统计:\n")
cat("上调基因数量:", sum(results$regulation == "Up"), "\n")
cat("下调基因数量:", sum(results$regulation == "Down"), "\n")
cat("总差异基因数量:", sum(results$significant == "YES"), "\n")
cat("总检测基因数量:", nrow(results), "\n")

# 第14步：保存结果
# 保存完整结果
write.csv(results, "limma_results/limma_all_results.csv", row.names = FALSE)

# 保存显著差异基因
deg_results <- results[results$significant == "YES", ]
write.csv(deg_results, "limma_results/limma_DEG_results.csv", row.names = FALSE)

# 分别保存上调和下调基因
up_genes <- results[results$regulation == "Up", ]
down_genes <- results[results$regulation == "Down", ]
write.csv(up_genes, "limma_results/limma_up_regulated_genes.csv", row.names = FALSE)
write.csv(down_genes, "limma_results/limma_down_regulated_genes.csv", row.names = FALSE)

# 第15步：生成基本统计图
# 创建图形保存目录
dir.create("limma_results/plots", showWarnings = FALSE)

# MA图
pdf("limma_results/plots/MA_plot.pdf", width = 8, height = 6)
limma::plotMA(fit, coef = "nucleus_vs_cytoplasm", 
              main = "MA Plot: Nucleus vs Cytoplasm")
abline(h = c(-deg_threshold$logFC_cutoff, deg_threshold$logFC_cutoff), 
       col = "red", lty = 2)
dev.off()

# 火山图
pdf("limma_results/plots/volcano_plot.pdf", width = 8, height = 6)
plot(results$logFC, -log10(results$P.Value), 
     pch = 20, cex = 0.5, 
     xlab = "log2 Fold Change", 
     ylab = "-log10(P-value)",
     main = "Volcano Plot: Nucleus vs Cytoplasm")
# 标记显著差异基因
points(results$logFC[results$regulation == "Up"], 
       -log10(results$P.Value[results$regulation == "Up"]), 
       col = "red", pch = 20, cex = 0.5)
points(results$logFC[results$regulation == "Down"], 
       -log10(results$P.Value[results$regulation == "Down"]), 
       col = "blue", pch = 20, cex = 0.5)
# 添加阈值线
abline(v = c(-deg_threshold$logFC_cutoff, deg_threshold$logFC_cutoff), 
       col = "gray", lty = 2)
abline(h = -log10(deg_threshold$pvalue_cutoff), col = "gray", lty = 2)
legend("topright", legend = c("Up", "Down", "Not significant"), 
       col = c("red", "blue", "black"), pch = 20)
dev.off()

# 第16步：显示前20个最显著的差异基因
cat("\n前20个最显著的差异基因:\n")
print(head(results[results$significant == "YES", ], 20))

# 第17步：样本相关性分析
pdf("limma_results/plots/sample_correlation.pdf", width = 8, height = 6)
cor_matrix <- cor(logCPM)
heatmap(cor_matrix, main = "Sample Correlation Heatmap")
dev.off()

# 输出会话信息
cat("\n分析完成！结果已保存在 limma_results/ 目录中\n")
cat("主要输出文件:\n")
cat("- limma_all_results.csv: 所有基因的分析结果\n")
cat("- limma_DEG_results.csv: 显著差异基因\n")
cat("- limma_up_regulated_genes.csv: 上调基因\n")
cat("- limma_down_regulated_genes.csv: 下调基因\n")
cat("- plots/: 包含MA图、火山图和样本相关性热图\n")
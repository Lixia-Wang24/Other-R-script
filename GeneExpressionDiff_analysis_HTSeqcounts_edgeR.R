# edgeR差异基因表达分析完整流程
# 加载必要的包
library(limma)
library(edgeR)  # 需要同时加载limma
library(dplyr)
library(readr)
library(ggplot2)

# 设置工作目录
setwd("/Users/lixia/Data/data/Fractionation_seq/HTSeqcount_labE")
# 创建结果保存目录
dir.create("edgeR_results", showWarnings = FALSE) # 避免目录已存在时的警告

### ---------------------------- Step1: 准备count数据 ---------------------------------------------------------------------- ###

# (1) 读取count数据，gene_id不能有重复值
count_file <- "prepare_counts/HTSeqcounts_counts_filter1.csv"
counts_filtered <- read_csv(count_file,locale = locale(encoding = "UTF-8"))
counts_matrax <- as.data.frame(counts_filtered[, 2:5]) 
rownames(counts_matrax) <- counts_filtered$gene_id

# (2) 读取样本信息， 对照组在前，处理组在后。 前面的位置默认是control的位置
sample_info <- as.data.frame(read_csv("prepare_counts/sample_info.csv", locale = locale(encoding = "UTF-8")))
rownames(sample_info) <- sample_info$sample 

# 设置参考水平（细胞质作为参考）
sample_info$condition <- factor(sample_info$condition, levels = c("cytoplasm", "nucleus"))

#设置group信息
# group <- rep(c('control', 'treat'), each = 2)

# (3) 再次对齐样本顺序, 确保样本顺序一致
counts_matrax <- counts_matrax[, rownames(sample_info)]

# 查看数据维度和结构
print("计数矩阵维度:")
print(dim(counts_matrax))
print("样本信息:")
print(sample_info)

### ---------------------------- Step2: “创建DGEList对象” 和 “前期数据准备” ------------------------------------------------ ###

# (1). 创建DGEList对象（edgeR的核心数据结构）
dge <- DGEList(counts = counts_matrax, 
               samples = sample_info,
               group = sample_info$condition)
# 查看
print("DGEList对象信息:")
print(dge)
summary(dge)

# (2). 基因过滤 - 去除低表达基因, 避免影响分析准确性，edgeR无内置函数过滤基因，需手动；DESeq2内置函数过滤
# 保留（在≥2个样本中raw count >1）（count-per-million）
# 标准差异分析，自动适配实验设计，建议保留至少 10,000-15,000 个基因进行后续分析
keep <- filterByExpr(dge, min.count = 1, min.total.count = 2) #前者一个样品(1样品>=1)，后者所有样品总数(sum>=2)
dge <- dge[keep, , keep.lib.sizes = FALSE] # keep.lib.sizes = FALSE 重新计算库大小

# 上述过滤的实际应用脚本是
# 等价逻辑
# cpm = count / (lib.size/1e6)
# cpm_threshold <- min.count / (lib.size/1e6)  # 转换为样本特异的 CPM 阈值
# keep <- rowSums(cpm(dge) >= cpm_threshold) >= min.prop.samples # min.prop.samples默认为最小分组大小的样本数

print(paste("过滤前基因数:", nrow(counts_matrax)))
print(paste("过滤后基因数:", nrow(dge), ",剩余:", round((nrow(dge)/nrow(counts_matrax)*100)), "%"))

# 画图查看：可视化基因表达分布并标记过滤阈值
plotDensities(log2(cpm(dge)+1), legend = "topright") # X轴: log2(CPM + 1)
abline(v = log2(1+1), col = "red") # 阈值线 (abline)，垂直线位置: 对应 CPM = 1 的表达水平

# filterByExpr() 比固定 CPM 阈值更智能：
# 自动适应测序深度：对深度不同的样本使用不同阈值
# 考虑实验设计：根据分组大小确定需要表达的样本数
# 避免过度过滤：对低深度样本更宽容

# 手动过滤方法----------------------------------------------------------------
# CPM 阈值过滤
# cpm_matrix <- cpm(dge)   # cpm(dge) 输出input gene的cpm
# keep <- rowSums(cpm_matrix > 1) >= 2  # 在≥2个样本中CPM>1
# dge <- dge[keep, , keep.lib.sizes = FALSE]

# 表达量分位数过滤 # prior.count 默认为 2/L，L是平均库大小的百万分之一
# gene_means <- rowMeans(cpm(dge, log = TRUE)) # log表示对计算得到的 CPM 值进行对数转换log2(cpm_raw + prior.count)
# keep <- gene_means > quantile(gene_means, probs = 0.25) # 保留上75%表达的基因
# dge <- dge[keep, , keep.lib.sizes = FALSE]

# 方差稳定性过滤，寻找高变异基因。关注动态变化基因
# cpm_log <- cpm(dge, log = TRUE, prior.count = 1)
# gene_vars <- apply(cpm_log, 1, var)
# keep <- gene_vars > quantile(gene_vars, probs = 0.5) # 保留方差前50%的基因
# dge <- dge[keep, , keep.lib.sizes = FALSE]

# 基因长度过滤
# gene_lengths <- ... # 从GTF等文件获取基因长度
# keep <- gene_lengths > 200  # 保留长度>200bp的基因
# dge <- dge[keep, , keep.lib.sizes = FALSE]

# 组合过滤策略，单细胞数据
# keep1 <- filterByExpr(dge, design = design)
# keep2 <- rowSums(cpm(dge) > 0.5) >= 3
# keep3 <- rowMeans(cpm(dge)) > quantile(rowMeans(cpm(dge)), 0.2)

# keep <- keep1 & keep2 & keep3
# dge <- dge[keep, , keep.lib.sizes = FALSE]
# ----------------------------------------------------------------------------

# (3). 标准化 - 计算标准化因子（使用TMM方法计算标准化因子，校正样本间RNA组成差异）
dge <- calcNormFactors(dge, method = "TMM")

print("标准化因子:")
print(dge$samples)

# (4). 探索性数据分析
# A:  MDS图 - 查看样本间关系 (Multidimensional Scaling Plot - 多维尺度分析图)
plotMDS(dge, labels = sample_info$sample, 
        col = as.numeric(as.factor(sample_info$condition)))
legend("topright", 
       legend = levels(as.factor(sample_info$condition)), 
       col = 1:4, pch = 1)
title("MDS Plot - Sample Clustering")

# B: 计算CPM值用于可视化
cpm_values <- cpm(dge, log = TRUE) # 对数转换log2

# 箱线图查看样本分布
boxplot(cpm_values, 
        names = colnames(cpm_values),
        main = "Log2 CPM Distribution",
        las = 1, cex.axis = 0.5) # label size

### ---------------------------- Step3: 差异表达分析 ----------------------------------------------------------------------- ###

# （1） 设计矩阵
design <- model.matrix(~ condition, data = sample_info)
colnames(design) <- c("Intercept", "nucleus_vs_cytoplasm")
# Intercept 代表了细胞质样本的基准表达水平，而 nucleus_vs_cytoplasm 系数则量化了基因在细胞核中的表达变化
# 基因表达 = β₀ × Intercept + β₁ × (nucleus_vs_cytoplasm)
print("设计矩阵:")
print(design)

# design <- model.matrix(~ replicate + condition, data = sample_info) # 未考虑批次效应（replicate）


# （2）估计基因表达的离散度（dispersion）它量化了基因表达变异的程度（包括技术变异和生物变异
dge <- estimateDisp(dge, design)

print("离散度估计:")
print(paste("Common dispersion:", round(dge$common.dispersion, 4))) # 所有基因共享的离散度基准值
print(paste("Trended dispersion range:",                            # 随基因表达水平变化的离散度值范围
            round(min(dge$trended.dispersion), 4), "-", 
            round(max(dge$trended.dispersion), 4)))

# 查看离散度分布： 绘制BCV图
plotBCV(dge, main = "Biological Coefficient of Variation")

# 结果解读
# 1:数据质量评估：
# 公共离散度 0.0114 → √0.0114 ≈ 0.107 BCV
# 表示平均生物变异系数为 10.7%，属于优秀水平 (BCV < 20%：高质量数据; BCV > 40%：可能有严重批次效应)

# 2:技术噪音评估：
# 低表达基因离散度 0.0334 → BCV=18.3%
# 高表达基因离散度 0.0097 → BCV=9.8%
# 符合预期：低表达基因噪音更高
# 若最小值>0.05 可能表示重复不一致
# 若最大值<0.01 可能过度过滤

# （3）拟合负二项式广义线性模型
fit <- glmQLFit(dge, design) # 拟合准似然负二项模型

# 对每个基因的表达数据拟合一个广义线性模型（GLM）,计算基因水平的离散度（biological variation）和准似然离散度（quasi-dispersion）
# 生成拟合对象：返回一个包含模型参数、离散度估计等信息的对象（DGEGLM类），供后续检验使用
# fit <- glmQLFit(y, design, prior.count = 1)  # 直接修改伪计数， 默认使用 log2(n + 0.125),可改为log2(n + 1)，这和DEseq2一样

# （4） 进行似然比检验
qlf <- glmQLFTest(fit, coef = 2)  # 检验组间差异（coef=2对应nucleus vs cytoplasm）

# 基于glmQLFit的拟合结果，检验设计矩阵中特定系数（coefficient）是否显著不为0（即该因素是否引起差异表达）
# 准似然F检验：使用改良的F检验（结合准似然离散度），比传统检验更适应RNA-seq数据的特性（过度离散、小样本）
# 生成结果对象：返回DGELRT对象，包含基因的logFC、p值、FDR等统计量。

### ---------------------------- Step4: 提取差异分析结果 ------------------------------------------------------------------- ###

# （1）获取差异表达基因结果
results <- topTags(qlf, n = Inf, sort.by = "p.value") # 提取全部差异表达结果，FDR= "p.value"，从小到大排序
de_results <- results$table
de_results$gene_id <- rownames(de_results)

# 查看
print("差异表达分析结果摘要:")
print(summary(decideTests(qlf, p.value = 0.05, lfc = log2(1.5)))) 

# 保存结果
write.csv(de_results, "edgeR_results/differential_expression_results.csv", row.names = F)

# （2）添加显著性标记
de_results$significant <- "No"
de_results$significant[de_results$FDR < 0.05 & abs(de_results$logFC) > log2(1.5)] <- "Yes" # 标记显著基因（FDR<0.05且|logFC|>1）

# 查看
print("显著差异表达基因数量:")
print(table(de_results$significant))

# 保存结果
de_results_significant <- de_results[de_results$significant == "Yes", ]
write.csv(de_results_significant, "edgeR_results/differential_expression_results_significant.csv", row.names = F)


# (3) 上调基因和下调基因
# 上调基因
de_results_up <- de_results_significant[de_results_significant$logFC > 1, ]
de_results_down <- de_results_significant[de_results_significant$logFC < -1, ]

# 查看
print(paste("上调基因数量:", nrow(de_results_up)))
print(paste("下调基因数量:", nrow(de_results_down)))

# 保存结果
write.csv(de_results_up, "edgeR_results/differential_expression_results_up.csv", row.names = F)
write.csv(de_results_down, "edgeR_results/differential_expression_results_down.csv", row.names = F)

### ---------------------------- Step5: 差异分析结果-可视化 ---------------------------------------------------------------- ###

# (1) 火山图，展示显著性vs变化倍数
#保存路径
png("edgeR_results/differential_expression_volcanno.png", 
    width = 8, height = 6, units = "in", res = 300)
#开始画图
plot(de_results$logFC, -log10(de_results$FDR),
     pch = 16, cex = 0.6,
     col = ifelse(de_results$significant == "Yes", "red", "gray"),
     xlab = "log2 Fold Change (Nucleus vs Cytoplasm)",
     ylab = "-log10(FDR)",
     main = "Volcano Plot")
# 添加阈值线
abline(h = -log10(0.05), col = "blue", lty = 2)
abline(v = c(-log2(1.5), log2(1.5)), col = "blue", lty = 2)
# 添加图例
legend("topright", 
       legend = c("Significant", "Not significant"),
       col = c("red", "gray"), pch = 16)

dev.off()

# (2) MA图 
#保存路径
png("edgeR_results/differential_expression_MAplot.png", 
    width = 8, height = 6, units = "in", res = 300)
#开始画图
avg_log_cpm <- qlf$table$logCPM
log_fc <- qlf$table$logFC
plot(avg_log_cpm, log_fc,
     pch = 16, cex = 0.6,
     col = ifelse(abs(log_fc) > log2(1.5) & qlf$table$PValue < 0.05, "red", "gray"),
     xlab = "Average log2 CPM",
     ylab = "log2 Fold Change",
     main = "MA Plot")
abline(h = c(-log2(1.5), 0, log2(1.5)), col = c("blue", "black", "blue"), lty = c(2, 1, 2))

dev.off()

### ---------------------------- Step6: 差异分析完成，输出分析工具信息 ---------------------------------------------------- ###
print("分析完成!")
sessionInfo()
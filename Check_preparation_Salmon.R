# 加载必要的包
library(dplyr)
library(readr)
library(tximport)

# 设置工作目录
setwd("/Users/lixia/Data/data/Fractionation_seq/Salmon/Ensembl")
# 创建结果保存目录
dir.create("prepare_counts", showWarnings = FALSE)

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

# 保存样本信息列表
write.csv(sample_info, file = "prepare_counts/sample_info.csv", row.names = FALSE)

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
                countsFromAbundance = "no")  # 选择没有经过矫正的counts， 此处设置是为了用DESeq2的内置frpkm()函数计算fpkm，（此处TPM计算用不到）

### ---------------------------- Step3: gene-level TPM ，gene_name annotation--------------------------------------------------------------- ###

# (5) 提取gene_TPM; tix 矩阵中含有"gene-level TPM", "gene length","length&TPM 计算出的gene count";
#txi$abundances是gene-level TPM：同一个gene不同transcript isoform TPM 的“加和”
gene_TPM <- as.data.frame(txi$abundance)
gene_TPM$gene_id <- rownames(gene_TPM) # 一个transcript_id对应一个gene_id， 需要去除重复的gene_id
gene_TPM_unique <- gene_TPM %>%
  distinct(gene_id, .keep_all = TRUE)

# (6) 通过gene_id 注释gene_name；若无对应gene_name,则用gene_id代替gene_name
tx2gene_name <- tx2gene_ano[,c("gene_id","gene_name")]

gene_name_TPM <- gene_TPM_unique %>%
  left_join(
    tx2gene_name %>% distinct(gene_id, .keep_all = TRUE), 
    by = "gene_id"
  ) %>%
  mutate(gene_name = coalesce(gene_name, gene_id)) %>% # 若gene_name有值，保留原值；若 gene_name缺失，则使用gene_id
  relocate(gene_id, gene_name) 

# (7). 检查
print(paste("总基因数:", nrow(gene_TPM_unique)))
print(paste("所有样本中tpm > 0的基因数:", sum(rowSums(gene_TPM_unique[, 1:4]) > 0)))

# (8). 保存结果
write.csv(gene_name_TPM, file = "prepare_counts/Salmon_gene_name_TPM.csv", row.names = FALSE)


### ---------------------------- Step4: 过滤低表达基因 --------------------------------------------------------------------------- ###

####（筛选标准不唯一、依情况而定）
# (1) 筛选“至少在1个列中的tpm大于0” # gene_name_TPM[, 3:6]>1 创建一个逻辑矩阵，判断每个样本中每个基因的计数是否大于0
keep_feature <- rowSums(gene_name_TPM[, 3:6]>0) >= 1 # rowSums()对每行（每个基因）的TRUE值求和
table(keep_feature)  #查看筛选情况，FALSE为低表达基因数（行数），TURE为要保留基因数

# 筛选所有列的tpm之和大于0；现有cutoff counts_filter1和2结果相同
counts_filter2  <- gene_name_TPM %>%
  filter(rowSums(select(., 3:6), na.rm = TRUE) > 0)

# (2) 过滤
gene_name_TPM_filter1 <- gene_name_TPM[keep_feature, ] #替换tpm为筛选后的基因矩阵（保留较高表达量的基因)

# (3)检查
print(paste("保留基因数量数量为:", nrow(gene_name_TPM_filter1) )) 
print(paste("保留基因数量百分比为:", round((nrow(gene_name_TPM_filter1)/nrow(gene_name_TPM))*100), "%")) # round四舍五入取整

# (4) 保存数据
write.csv(gene_name_TPM_filter1, file = "prepare_counts/Salmon_gene_name_TPM_filter1.csv", row.names = F)
# 用于组合图形

### -------------差异分析前的准备——数据检查----------------------------------------------------------------------------------------------- ###
#预先判断样品间的差异： hclust 图、距离热图、PCA图(主要特征)、前500差异性大的基因热图、相关性热图

# BiocManager::install("x")
library(FactoMineR)
library(factoextra)  
library(tidyverse) # ggplot2 stringer dplyr tidyr readr purrr  tibble forcats
library(pheatmap)
library(RColorBrewer)
# library(DESeq2)
# library(edgeR)
library(ggpointdensity)  # 用于密度散点图
library(patchwork)   # 用于合并图片
library(ggpubr)  # `viridis` 色标
# library(hexbin)  # 用于创建六边形分箱图
library(gridExtra) # 用于组合图形

### ------------ Step5: 样本件归一化 （列举了3种归一化方法A-B-C）——————————----------———————----------------------------------------—--———— ###

# A: dat <- as.data.frame(log2(edgeR::cpm(counts)+1))  #简单归一化 CPM:Counts per million
# B: DESeq2_normalize： rld <- rlogTransformation(dds, blind = FALSE)  
# C: log2(TPM+1) 此处使用
dat <- log2((gene_name_TPM_filter1[, 3:6])+1)

### ------------ Step5: boxplot 查看样本的基因整体表达情况——————————----------———————---------------------------------------------------------------------———--———— ###
# R基础命令画图
num_samples <- ncol(dat) - 1
color <- rainbow(num_samples)  # 方法1：彩虹色系（默认）
# color <- heat.colors(num_samples)  # 方法2：热力图色系
# color <- terrain.colors(num_samples)  # 方法3：地形色系
# color <- c(rep("#1f77b4", 2), rep("#ff7f0e", 2)) #手动指定颜色，细胞质重复1-2: 蓝色系细胞核重复1-2: 橙色系 

boxplot(dat, col=color, ylab="log2(TPM+1)", main=" normalized data ",
        outline = F, notch = F)
dev.copy(png, "prepare_counts/TMP_boxplot.png", width = 2000, height = 1500, res = 300)
dev.off()


### ------------ Step5: boxplot 查看样本的基因整体表达情况——————————----------———————---------------------------------------------------------------------———--———— ###

# ggplot画图  

# 1. 准备数据
# 添加基因ID（如果dat中没有行名）
dat_with_id <- as.data.frame(dat)
dat_with_id$gene_id <- rownames(dat)  # 行名是编号

# 转换为长格式
dat_long <- pivot_longer(
  dat_with_id,
  cols = -gene_id,
  names_to = "sample",
  values_to = "tpm"
)

# 2. 从样本名称提取分组信息
dat_long$group <- ifelse(grepl("cyt", dat_long$sample), "Cytoplasm", "Nucleus")

# 3. 设置颜色映射
group_colors <- c("Cytoplasm" = "#1f78b4", "Nucleus" = "#e31a1c")  # 蓝/红配色

# 4. 创建ggplot箱线图
p <- ggplot(dat_long, aes(x = sample, y = tpm, fill = group)) +
  geom_boxplot(outlier.shape = NA, notch = FALSE) +  # 不显示异常点
  scale_fill_manual(values = group_colors) +        # 设置分组颜色
  labs(title = "Normalized TPM Expression",
       y = "TPM",
       x = "Samples",
       fill = "Fraction") +
  theme_minimal() +
  theme(
    # 自定义x轴标签
    axis.text.x = element_text(
      angle = 45,          # 45度倾斜
      hjust = 1,           # 右对齐
      size = 12,           # 字体大小
      face = "italic",     # 斜体
      color = "black"      # 字体颜色
    ),
    panel.border = element_rect(color = "black", fill = NA, size = 1),
    # 自定义y轴标签
    axis.text.y = element_text(
      size = 12,
      face = "italic"
    ),
    # 标题和轴标题
    plot.title = element_text(
      size = 16,
      face = "bold",
      hjust = 0.5
    ),
    axis.title = element_text(
      size = 14,
      face = "bold"
    ),
    # 图例
    legend.position = "top",
    legend.title = element_text(face = "bold")
  )

# 5. 显示图形
print(p)

# 6. 保存图形
ggsave("prepare_counts/tpm_boxplot_ggplot.pdf", plot = p, width = 10, height = 8, dpi = 300)

### ------------ step6. 基因表达相关性： A: 组内pearson correlation——————————-------------------------------------------------———--———— ###

# (1). 绘制细胞质重复样本的相关性图 (K562.LabE.polyA.cyt.rep1 vs rep2)
cyt_plot <- ggplot(as.data.frame(dat), 
                   aes(x = K562.LabE.polyA.cyt.rep1, 
                       y = K562.LabE.polyA.cyt.rep2)) +
  geom_pointdensity(size = 0.6) +  # 密度点图
  scale_color_viridis_c(option = "C") +  # 使用viridis颜色方案
  geom_smooth(method = "lm", color = "red", se = FALSE) +  # 添加回归线
  labs(title = "Cytoplasmic Replicates Correlation",
       x = "rep1 (log2(TPM+1))", 
       y = "rep2 (log2(TPM+1))") +
  theme_bw() +
  stat_cor(method = "pearson",  # 自动添加Pearson相关系数
           label.x.npc = "left",
           label.y.npc = "top")

# (2). 绘制细胞核重复样本的相关性图 (K562.LabE.polyA.nuc.rep1 vs rep2)
nuc_plot <- ggplot(as.data.frame(dat), 
                   aes(x = K562.LabE.polyA.nuc.rep1, 
                       y = K562.LabE.polyA.nuc.rep2)) +
  geom_pointdensity(size = 0.6) +
  scale_color_viridis_c(option = "D") +
  geom_smooth(method = "lm", color = "red", se = FALSE) +
  labs(title = "Nuclear Replicates Correlation",
       x = "rep1 (log2(TPM+1))", 
       y = "rep2 (log2(TPM+1))") +
  theme_bw() +
  stat_cor(method = "pearson",
           label.x.npc = "left",
           label.y.npc = "top")

# 显示图形
cyt_plot+nuc_plot

# 保存图形（可选）
ggsave("prepare_counts/cyt_replicates_corr.pdf", plot = cyt_plot, width = 6, height = 5)
ggsave("prepare_counts/nuc_replicates_corr.pdf", plot = nuc_plot, width = 6, height = 5)

# `viridis` 色标有五种选项（A到E），分别代表不同的颜色渐变方案：
# option = "A"：viridis（默认，绿色到紫色）
# option = "B"：magma（深紫色到亮黄色）
# option = "C"：plasma（深紫色到亮橙色）
# option = "D"：inferno（黑色到亮黄色）
# option = "E"：cividis（深蓝色到亮黄色）

### ------------ step6. 基因表达相关性： B: 组间pearson correlation——————————----------———————--------------------------------———--———— ###

# 创建绘图函数
create_density_plot <- function(df, x_col, y_col) {
  # 计算斯皮尔曼相关系数
  corr <- cor(df[[x_col]], df[[y_col]], method = "spearman")
  corr_label <- paste0("ρ = ", round(corr, 3))
  
  ggplot(df, aes(x = x_col, y = y_col)) +
    geom_density_2d(aes(color = ..level..)) +
    scale_color_viridis_c(option = "plasma")
  
  # 创建密度散点图
  ggplot(df, aes(x = .data[[x_col]], y = .data[[y_col]])) +
    #geom_hex(bins = 50) +  # 使用六边形分箱
    geom_pointdensity(size = 0.8, alpha = 0.5) +
    scale_fill_viridis_c(option = "plasma", trans = "log10") + # 使用对数颜色标度
    geom_smooth(method = "lm", color = "red", se = FALSE) + # 添加回归线
    # 将标签移动到左上角
    annotate("text", x = -Inf, y = Inf, label = corr_label,
             hjust = -0.5, vjust = 4, size = 5, color = "black") +
    # 将标签移动到右下角
    #annotate("text", x = Inf, y = -Inf, label = corr_label,
    #         hjust = 1.1, vjust = -1, size = 5, color = "black") +
    labs(title = paste(y_col, "vs", x_col),
         x = paste("Expression:", x_col),
         y = paste("Expression:", y_col)) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 5),
      panel.border = element_rect(color = "black", fill = NA, size = 1),  # 添加黑色边框
      panel.grid.major = element_line(color = "gray90"),
      panel.grid.minor = element_blank(),
      legend.position = "right"
    )
}

# 创建第一组图（rep1）
p1 <- create_density_plot(dat, 
                          "K562.LabE.polyA.cyt.rep1", 
                          "K562.LabE.polyA.nuc.rep1")

# 创建第二组图（rep2）
p2 <- create_density_plot(dat, 
                          "K562.LabE.polyA.cyt.rep2", 
                          "K562.LabE.polyA.nuc.rep2")
# 查看图形
p1+p2

# 保存图形（可选）
ggsave("prepare_counts/samples_rep1_corr.png", plot = p1, width = 6, height = 5)
ggsave("prepare_counts/samples_rep2_corr.png", plot = p2, width = 6, height = 5) #pdf 不能正确输出p = 0.824
# 保存图形（可选）
#ggsave("prepare_counts/correlation_plots.png", arrangeGrob(p1, p2, ncol=2), 
#       width = 12, height = 6, dpi = 300)


# 并排显示
correlation_scaterplot <- (cyt_plot+nuc_plot) / (p1+p2) + 
  plot_annotation(title = "Comprehensive PCA Analysis",
                  theme = theme(plot.title = element_text(size = 20, hjust = 0.5)))
print(correlation_scaterplot)
ggsave(correlation_scaterplot, filename = 'prepare_counts/correlation_scaterplot.png', width = 12, height = 8)
dev.off()

### ------------ step6. 基因表达相关性： C:  correlation heatmap (取500高表达基因) -——————————--------------------------------———--———— ###

# 取500高表达基因做相关性检测，避免低表达基因的干扰
dat_500 <- dat[names(sort(apply(dat,1,mad),decreasing = T)[1:500]),] 
M <- cor(dat_500) #默认method = "pearson"， 组内用pearson
# M <- cor(dat_500, method = "spearman")  # 组间用spearman， 但从结果看pearson更好

# 准备列注释 (确保与表达矩阵列名匹配)
annotation_col <- sample_info[, c("condition", "replicate")]

correlation_heatmap <-pheatmap::pheatmap(M,
                                         show_rownames = T,
                                         angle_col=45,
                                         fontsize=7,
                                         annotation_col=annotation_col) 

ggsave(correlation_heatmap,filename = 'prepare_counts/check_cor_top500.pdf',width = 7.5,height =6)

### ------------ step7. 样本距离：hclust and Heatmap of the sample-to-sample distances——————————------------------------------———--———— ###

sampleDists <- dist(t(dat))   #dist默认计算矩阵行与行的距离， 因此需要转置
sampleDistMatrix <- as.matrix(sampleDists)  
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)  #选取热图的颜色

# 画 distance heatmap
p0 <- pheatmap::pheatmap(sampleDistMatrix,
                         fontsize=7,
                         clustering_distance_rows=sampleDists,
                         clustering_distance_cols=sampleDists,
                         angle_col=45,
                         col=colors)
ggsave(p0,filename = 'prepare_counts/check_dist.pdf',width = 7.5,height =6)
dev.off()

# 画cluster
pdf("prepare_counts/check_hclust.pdf")
plot(hclust(sampleDists))
dev.off()

### ------------ step8. 主成分分析：PCA (DESeq2自带/PCA函数) 查看样品分组； A: PC1-2-3分别画——————————---------------———-----------———— ###

# （1）计算PCA
dat.pca <- PCA(t(dat), graph = FALSE)

# （2）获取方差百分比
eig.val <- get_eigenvalue(dat.pca)
percentVar <- eig.val[, "variance.percent"]

# （3）定义分组
group_exp <- c("cyt", "cyt", "nuc", "nuc")
group_batch <- c("batch1", "batch2", "batch1", "batch2")

#（4） 设置动态坐标范围计算函数-------------------------------------------------
get_dynamic_limits <- function(pca_result, pc_index, buffer_percent = 0.2) {
  # 获取指定PC的坐标值
  pc_values <- get_pca_ind(pca_result)$coord[, pc_index]
  # 计算范围
  pc_range <- range(pc_values)
  range_diff <- diff(pc_range)
  # 添加缓冲区 (默认10%)
  buffer <- range_diff * buffer_percent
  limits <- c(pc_range[1] - buffer, pc_range[2] + buffer)
  return(limits)
}
# 设置动态坐标范围计算函数------------------------------------------------------

# （5）设置画图函数-------------------------------------------------------------
plot_custom_pca_ggplot <- function(pc_x, pc_y, title) {
  # 获取PCA坐标
  pca_data <- as.data.frame(get_pca_ind(dat.pca)$coord)
  colnames(pca_data) <- paste0("PC", 1:ncol(pca_data))
  pca_data$Group <- group_exp
  pca_data$Batch <- group_batch
  
  # 动态计算坐标范围
  xlim <- get_dynamic_limits(dat.pca, pc_x)
  ylim <- get_dynamic_limits(dat.pca, pc_y)
  
  # 创建图形
  ggplot(pca_data, aes(x = .data[[paste0("PC", pc_x)]], 
                       y = .data[[paste0("PC", pc_y)]], 
                       color = Group,
                       shape = Batch)) +
    geom_point(size = 4) +
    geom_text(aes(label = rownames(pca_data)), 
              hjust = +0.5, vjust = 2, size = 2, show.legend = FALSE) +
    ggtitle(title) +
    coord_fixed(ratio = 1) + 
    xlab(paste0("PC", pc_x, " (", round(percentVar[pc_x], 1), "%)")) +
    ylab(paste0("PC", pc_y, " (", round(percentVar[pc_y], 1), "%)")) +
    coord_cartesian(xlim = xlim, ylim = ylim) +
    theme_bw() +
    theme(panel.border = element_rect(color = "black", size = 1),
          legend.position = "right") +
    scale_color_manual(values = c("cyt" = "blue", "nuc" = "red")) +
    scale_shape_manual(values = c("batch1" = 16, "batch2" = 17)) +
    guides(color = guide_legend(title = "Cell Location"),
           shape = guide_legend(title = "Batch"))
}
# 设置画图函数------------------------------------------------------------------

# （6）开始画图
p1_2 <- plot_custom_pca_ggplot(1, 2, "PC1 vs PC2: Biological Variation")
p1_3 <- plot_custom_pca_ggplot(1, 3, "PC1 vs PC3: Batch Effect Check")
p2_3 <- plot_custom_pca_ggplot(2, 3, "PC2 vs PC3: Secondary Patterns")

# （7）保存图形 
ggsave(p1_2, filename = 'prepare_counts/p1_2_PCA.pdf', width = 8, height = 7)
ggsave(p1_3, filename = 'prepare_counts/p1_3_PCA.pdf', width = 8, height = 7)
ggsave(p2_3, filename = 'prepare_counts/p2_3_PCA.pdf', width = 8, height = 7)
dev.off()

# （8）并排显示
combined_1_2_3_PCA <- (p1_2 / p1_3 / p2_3)
print(combined_1_2_3_PCA)
ggsave(combined_1_2_3_PCA, filename = 'prepare_counts/combined_1_2_3_PCA.pdf', width = 8, height = 8)
dev.off()

### ------------ step8. 主成分分析：PCA (DESeq2自带/PCA函数) 查看样品分组； B: 经典PC1-2——————————-------------------———-----------———— ###

# PCA函数画图
dat.pca <- PCA(t(dat) , graph = F)  # t()转置，graph = F 不自动生成图形（后续手动定制）

# 计算主成分的方差贡献百分比（新增关键步骤）
eig.val <- get_eigenvalue(dat.pca)  # 获取PCA结果的特征值（每个主成分的方差）
percentVar <- eig.val[, "variance.percent"]  # 特征值包含多个列，这里提取了主成分解释的方差百分比
# group_list <- colnames(Salmon_tpm3) # 将原始数据的列名（即样本名）作为分组信息。
group_list <- ifelse(grepl("cyt", colnames(dat)), "cyt", "nuc")
# group_list <- c("cyt", "cyt", "nuc", "nuc")

# 打印可用的主成分数
n_pc <- ncol(dat.pca$ind$coord)
# 最大主成分数量 = min(样本数-1, 基因数) = min(3, 22378) = 3
# 因此，有效的主成分编号只能是1,2,3
# 但如果数据只有2个有效PC（某些特殊情况），PC3就不存在

# 获取样本坐标数据
# `get_pca_ind`（来自`factoextra`）获取样本（个体）在主成分上的坐标。
# `$coord`提取坐标矩阵，每一行是一个样本，每一列是一个主成分。
pca_data <- get_pca_ind(dat.pca)$coord
pc1_range <- range(pca_data[, 1])  # 第一列是PC1
pc2_range <- range(pca_data[, 2])  # 第二列是PC2

# 计算扩展比例 (在图形上，增加10%的边界空间)
# 为了图形美观，我们希望在坐标轴两端留一些空白。这里计算了PC1和PC2范围宽度的10%作为扩展量。
expand_factor <- 0.20
x_expand <- diff(pc1_range) * expand_factor 
y_expand <- diff(pc2_range) * expand_factor

# 绘制PCA图并动态调整范围
pca <- fviz_pca_ind(dat.pca,
                    title = "Principal Component Analysis",
                    legend.title = "Groups",
                    geom.ind = c("point", "text"), #指定样本的几何形状，这里同时显示点（point）和文本标签（text）
                    pointsize = 1.5, #点的大小
                    labelsize = 2, # 文本标签的大小
                    col.ind = group_list,  # 分组上色
                    axes.linetype = NA,    # 移除坐标轴线
                    mean.point = FALSE     # 去除分组中心点
) + 
  coord_fixed(ratio = 1) +  # # 固定坐标轴比例，确保x轴和y轴单位长度相等
  xlab(paste0("PC1 (", round(percentVar[1], 1), "%)")) + # 设置x轴标签，显示PC1的方差贡献百分比
  ylab(paste0("PC2 (", round(percentVar[2], 1), "%)")) + # 设置y轴标签，显示PC2的方差贡献百分比
  # 动态设置坐标轴范围
  xlim(pc1_range[1] - x_expand, pc1_range[2] + x_expand) + # 设置x轴范围，扩展10%
  ylim(pc2_range[1] - y_expand, pc2_range[2] + y_expand) + # 设置y轴范围，扩展10%
  # 添加黑色边框
  theme(
    panel.border = element_rect(colour = "black", fill = NA, size = 0.5), # 添加黑色边框
    # plot.background = element_rect(colour = "black", size = 0.5)
  )

pca
# 保存图形 - 增加尺寸
ggsave(pca, filename = 'prepare_counts/check_PCA.pdf', width = 8, height = 7)
dev.off()

### ------------ 结束语 --------------------------------------------------------——————————-------------------———-----------———— ###
print("基因表达counts预处理已完成")

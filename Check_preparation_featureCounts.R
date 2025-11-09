# Load necessary packages
library(dplyr)
library(readr)

# Set working directory
setwd("/Users/lixia/Data/data/Fractionation_seq/featureCounts_labE")
# Create results directory
dir.create("prepare_counts", showWarnings = FALSE)

### ---------------------------- Step1: Read and Check Count Data --------------------------------------------------------------- ###

# (1) Read counts data (CSV, tab-separated)
counts_data <- as.data.frame(read_delim("K562.LabE.polyA.CytvsNuc.featureCounts.count.csv", delim = "\t", col_names = TRUE))

# (2) Check gene count distribution
print(paste("Total number of genes:", nrow(counts_data)))
print(paste("Number of genes with counts > 0 across all samples:", sum(rowSums(counts_data[, 3:6]) > 0)))
str(counts_data) # Check the data type of each column

# (3) Create sample metadata
sample_info <- data.frame(
  sample = colnames(counts_data[, 3:6]),
  condition = factor(c("cytoplasm", "cytoplasm", "nucleus", "nucleus")),
  replicate = c("rep1", "rep2", "rep1", "rep2"),
  row.names = colnames(counts_data[, 3:6])
)

# (4) Save sample information table
write.csv(sample_info, file = "prepare_counts/sample_info.csv", row.names = FALSE)

### ---------------------------- Step2: Calculate TPM --------------------------------------------------------------------------- ###

# (1) Pre-process counts data. Allows: counts = 0 (unless all are 0). Prohibits: effLen = 0 (leads to NaN/Inf), NULL (causes error), NA (results in all NA). Ensure effLen > 0 and no missing values.
# Check if gene_length is <=0 or NA
invalid_lengths <- counts_data$gene_length <= 0 | is.na(counts_data$gene_length)
# View number of invalid values
sum(invalid_lengths)
# Handle invalid values (Example: set to NA and exclude in subsequent calculations)
counts_data$gene_length[invalid_lengths] <- NA

# (2) Define TPM function
# Input: "raw count vector of genes" and "effective length vector of genes"
countToTpm <- function(counts, effLen)
{
  rate <- log(counts) - log(effLen) # Calculate log(count/effective length) for each gene, i.e., log(RPK)
  denom <- log(sum(exp(rate)))     # Take exponent of RPK, sum, and take log, i.e., log(total RPK)
  exp(rate - denom + log(1e6))     # Calculate TPM: (RPK/total RPK) * 10^6
}

# (3) Calculate TPM for each sample column separately

# a. Create results container
TPM_results <- data.frame(gene_id = counts_data$gene_id)

# b. Sample list
samples <- c("K562.LabE.polyA.cyt.rep1", "K562.LabE.polyA.cyt.rep2",
             "K562.LabE.polyA.nuc.rep1", "K562.LabE.polyA.nuc.rep2")

# c. Calculate TPM for each sample
for (sample in samples) {
  tpm_values <- countToTpm(
    counts = counts_data[[sample]],
    effLen = counts_data$gene_length
  )
  TPM_results[[paste0("TPM_", sample)]] <- tpm_values
}

# d. Check
print(paste("Total number of genes:", nrow(TPM_results)))
print(paste("Number of genes with tpm > 0 across all samples:", sum(rowSums(TPM_results[, 2:5]) > 0)))

# e. Save results
write.csv(TPM_results, file = "prepare_counts/featureCounts_gene_TPM.csv", row.names = FALSE)

### ---------------------------- Step3: Annotate gene_name --------------------------------------------------------------------------- ###

# (1) Read gene_id_name mapping table
map_file <- "/Users/lixia/Data/database/ref_RNA/Ensembl/GeneID_GeneName_Transcript_mapping_table/gene_id_name_map_Ensembl_rtracklayer.csv"
map <- as.data.frame(read_csv(file = map_file, locale = locale(encoding = "UTF-8")))

# (2) Annotate gene_name using gene_id; if no corresponding gene_name, use gene_id as substitute

TPM_geneNmae <- TPM_results %>%
  left_join(
    map %>% distinct(gene_id, .keep_all = TRUE),
    by = "gene_id"
  ) %>%
  mutate(gene_name = coalesce(gene_name, gene_id)) %>% # If gene_name exists, keep the original value; if gene_name is missing, use gene_id
  select(gene_id, gene_name, everything())

# (3) Check how many gene_ids were annotated with gene_name
success_anno <- sum(TPM_geneNmae$gene_name != TPM_geneNmae$gene_id, na.rm = TRUE)
print(paste(nrow(TPM_results),"genes were annotated with", success_anno, "gene_names"))

# (4). Save results
write.csv(TPM_geneNmae, file = "prepare_counts/featureCounts_gene_TPM_geneName.csv", row.names = FALSE)

### ---------------------------- Step4: Filter low-expression genes --------------------------------------------------------------------------- ###

#### (Filtering criteria are not unique and depend on the situation)
# (1) Filter: "count >= 1 in at least 1 column" # counts_data[, 3:6]>=1 creates a logical matrix to check if the count for each gene in each sample is >= 1
keep_feature <- rowSums(counts_data[, 3:6]>=1) >= 1 # rowSums() sums the TRUE values for each row (each gene)
table(keep_feature)  # Check filtering status. FALSE: number of low-expression genes (rows), TRUE: number of genes to keep

# Filter: sum of counts across all columns is greater than 1; current cutoff gives the same results for counts_filter1 and 2
counts_filter2  <- counts_data %>%
  filter(rowSums(select(., 3:6), na.rm = TRUE) >= 1)

# (2) Filter
counts_filter1 <- counts_data[keep_feature, ] # Replace counts with the filtered gene matrix (keeping higher expressed genes)
TPM_geneNmae_filter1 <- TPM_geneNmae[keep_feature, ]

# (3) Check
print(paste("Number of retained genes:", nrow(counts_filter1) ))
# round to the nearest integer
print(paste("Percentage of retained genes:", round((nrow(counts_filter1)/nrow(counts_data))*100), "%"))

# (4) Save data
write.csv(counts_filter1, file = "prepare_counts/featureCounts_counts_filter1.csv", row.names = FALSE)
write.csv(TPM_geneNmae_filter1, file = "prepare_counts/featureCounts_gene_TPM_geneName_filter1.csv", row.names = FALSE)

### ------------- Preparation before Differential Analysis - Data Checks ----------------------------------------------------------------------------------------------- ###
# Pre-assess differences between samples: hclust plot, distance heatmap, PCA plot (major features), heatmap of top 500 highly variable genes, correlation heatmap

# BiocManager::install("x")
library(FactoMineR)
library(factoextra)
library(tidyverse) # ggplot2 stringer dplyr tidyr readr purrr  tibble forcats
library(pheatmap)
library(RColorBrewer)
# library(DESeq2)
# library(edgeR)
library(ggpointdensity)  # Used for density scatter plots
library(patchwork)     # Used for combining plots
library(ggpubr)  # 'viridis' color scale
# library(hexbin)  # Used for creating hexagonal binning plots
library(gridExtra) # Used for combining plots


# Used for combining plots
### ------------ Step5: Sample Normalization (Listing 3 normalization methods A-B-C) ——————————----------———————----------------------------------------—--———— ###

# A: dat <- as.data.frame(log2(edgeR::cpm(counts)+1)) # Simple normalization: CPM: Counts per million
# B: DESeq2_normalize: rld <- rlogTransformation(dds, blind = FALSE)
# C: log2(TPM+1) Used here
dat <- log2((TPM_geneNmae_filter1[, 3:6])+1)

rawCounts <- log2((counts_filter1[, 3:6])+1)

### ------------ Step5: Boxplot to check overall gene expression of samples ——————————----------———————---------------------------------------------------------------------———--———— ###
# Plotting with R base commands
num_samples <- ncol(dat) - 1
color <- rainbow(num_samples)  # Method 1: Rainbow color palette (default)
# color <- heat.colors(num_samples)  # Method 2: Heat color palette
# color <- terrain.colors(num_samples)  # Method 3: Terrain color palette
# color <- c(rep("#1f77b4", 2), rep("#ff7f0e", 2)) # Manually specify colors, cytoplasm rep1-2: blue palette, nucleus rep1-2: orange palette

boxplot(dat, col=color, ylab="log2(TPM+1)", main="Normalized Data",
        outline = F, notch = F)
dev.copy(png, "prepare_counts/TMP_boxplot.png", width = 2000, height = 1500, res = 300)
dev.off()

boxplot(rawCounts, col=color, ylab="rawCounts", main="Raw Data",
        outline = F, notch = F)
dev.copy(png, "prepare_counts/rawCounts_boxplot.png", width = 2000, height = 1500, res = 300)
dev.off()

### ------------ Step5: Boxplot to check overall gene expression of samples ——————————----------———————---------------------------------------------------------------------———--———— ###

# Plotting with ggplot

# 1. Prepare data
# Add gene ID (if no row names in dat)
dat_with_id <- as.data.frame(dat)
dat_with_id$gene_id <- rownames(dat)  # If row names are gene IDs

# Convert to long format
dat_long <- pivot_longer(
  dat_with_id,
  cols = -gene_id,
  names_to = "sample",
  values_to = "tpm"
)

# 2. Extract grouping information from sample names
dat_long$group <- ifelse(grepl("cyt", dat_long$sample), "Cytoplasm", "Nucleus")

# 3. Set color mapping
group_colors <- c("Cytoplasm" = "#1f78b4", "Nucleus" = "#e31a1c")  # Blue/Red color scheme

# 4. Create ggplot boxplot
p <- ggplot(dat_long, aes(x = sample, y = tpm, fill = group)) +
  geom_boxplot(outlier.shape = NA, notch = FALSE) +  # Do not show outliers
  scale_fill_manual(values = group_colors) +          # Set group colors
  labs(title = "Normalized TPM Expression",
       y = "log2(TPM+1)",
       x = "Samples",
       fill = "Fraction") +
  theme_minimal() +
  theme(
    # Custom x-axis labels
    axis.text.x = element_text(
      angle = 45,          # 45-degree angle
      hjust = 1,           # Right-justified
      size = 12,           # Font size
      face = "italic",     # Italic font face
      color = "black"      # Font color
    ),
    panel.border = element_rect(color = "black", fill = NA, size = 1),
    # Custom y-axis labels
    axis.text.y = element_text(
      size = 12,
      face = "italic"
    ),
    # Title and axis titles
    plot.title = element_text(
      size = 16,
      face = "bold",
      hjust = 0.5
    ),
    axis.title = element_text(
      size = 14,
      face = "bold"
    ),
    # Legend
    legend.position = "top",
    legend.title = element_text(face = "bold")
  )

# 5. Display plot
print(p)

# 6. Save plot
ggsave("prepare_counts/tpm_boxplot_ggplot.pdf", plot = p, width = 10, height = 8, dpi = 300)

### ------------ step6. Gene Expression Correlation: A: Intra-group Pearson Correlation ——————————-------------------------------------------------———--———— ###

# (1). Plot correlation for cytoplasmic replicate samples (K562.LabE.polyA.cyt.rep1 vs rep2)
cyt_plot <- ggplot(as.data.frame(dat),
                   aes(x = TPM_K562.LabE.polyA.cyt.rep1,
                       y = TPM_K562.LabE.polyA.cyt.rep2)) +
  geom_pointdensity(size = 0.6) + # Density point plot
  scale_color_viridis_c(option = "C") + # Use viridis color scheme
  geom_smooth(method = "lm", color = "red", se = FALSE) +  # Add regression line
  labs(title = "Cytoplasmic Replicates Correlation",
       x = "rep1 (log2(TPM+1))",
       y = "rep2 (log2(TPM+1))") +
  theme_bw() +
  stat_cor(method = "pearson",  # Automatically add Pearson correlation coefficient
           label.x.npc = "left",
           label.y.npc = "top")

# (2). Plot correlation for nuclear replicate samples (K562.LabE.polyA.nuc.rep1 vs rep2)
nuc_plot <- ggplot(as.data.frame(dat),
                   aes(x = TPM_K562.LabE.polyA.nuc.rep1,
                       y = TPM_K562.LabE.polyA.nuc.rep2)) +
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

# Display plot
print(cyt_plot)
print(nuc_plot)
cyt_plot+nuc_plot
# Save plot (optional)
ggsave("prepare_counts/cyt_replicates_corr.pdf", plot = cyt_plot, width = 6, height = 5)
ggsave("prepare_counts/nuc_replicates_corr.pdf", plot = nuc_plot, width = 6, height = 5)

# The 'viridis' color scale has five options (A to E), each representing a different color gradient scheme:
# option = "A": viridis (default, green to purple)
# option = "B": magma (deep purple to bright yellow)
# option = "C": plasma (deep purple to bright orange)
# option = "D": inferno (black to bright yellow)
# option = "E": cividis (deep blue to bright yellow)

### ------------ step6. Gene Expression Correlation: B: Inter-group Pearson Correlation ——————————----------———————--------------------------------———--———— ###

# Create plotting function
create_density_plot <- function(df, x_col, y_col) {
  # Calculate Spearman correlation coefficient
  corr <- cor(df[[x_col]], df[[y_col]], method = "spearman")
  corr_label <- paste0("ρ = ", round(corr, 3))

  ggplot(df, aes(x = .data[[x_col]], y = .data[[y_col]])) +
    #geom_hex(bins = 50) +  # Use hexagonal binning
    geom_pointdensity(size = 0.8, alpha = 0.5) +
    scale_fill_viridis_c(option = "plasma", trans = "log10") + # Use logarithmic color scale
    geom_smooth(method = "lm", color = "red", se = FALSE) + # Add regression line
    # Move label to top-left corner
    annotate("text", x = -Inf, y = Inf, label = corr_label,
             hjust = -0.5, vjust = 4, size = 5, color = "black") +
    # Move label to bottom-right corner
    #annotate("text", x = Inf, y = -Inf, label = corr_label,
    #         hjust = 1.1, vjust = -1, size = 5, color = "black") +
    labs(title = paste(y_col, "vs", x_col),
         x = paste("Expression:", x_col),
         y = paste("Expression:", y_col)) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 5),
      panel.border = element_rect(color = "black", fill = NA, size = 1), # Add black border
      panel.grid.major = element_line(color = "gray90"),
      panel.grid.minor = element_blank(),
      legend.position = "right"
    )
}

# Create first plot set (rep1)
p1 <- create_density_plot(dat,
                          "TPM_K562.LabE.polyA.cyt.rep1",
                          "TPM_K562.LabE.polyA.nuc.rep1")

# Create second plot set (rep2)
p2 <- create_density_plot(dat,
                          "TPM_K562.LabE.polyA.cyt.rep2",
                          "TPM_K562.LabE.polyA.nuc.rep2")

# Save plot (optional)
ggsave("prepare_counts/samples_rep1_corr.png", plot = p1, width = 6, height = 5)
ggsave("prepare_counts/samples_rep2_corr.png", plot = p2, width = 6, height = 5) #pdf cannot correctly output p = 0.824
# Save plot (optional)
#ggsave("prepare_counts/correlation_plots.png", arrangeGrob(p1, p2, ncol=2),
#       width = 12, height = 6, dpi = 300)
p1+p2

# Display side-by-side
correlation_scaterplot <- (cyt_plot+nuc_plot) / (p1+p2) +
  plot_annotation(title = "Comprehensive Correlation Analysis",
                  theme = theme(plot.title = element_text(size = 20, hjust = 0.5)))
print(correlation_scaterplot)
ggsave(correlation_scaterplot, filename = 'prepare_counts/correlation_scaterplot.png', width = 12, height = 8)
dev.off()

### ------------ step6. Gene Expression Correlation: C: correlation heatmap (Top 500 highly variable genes) -——————————--------------------------------———--———— ###

# Select top 500 highly variable genes for correlation check to avoid interference from low-expression genes
colnames(dat) <- c("K562.LabE.polyA.cyt.rep1", "K562.LabE.polyA.cyt.rep2",
                   "K562.LabE.polyA.nuc.rep1", "K562.LabE.polyA.nuc.rep2")

# mad (Median Absolute Deviation) is used as a robust measure of variability
dat_500 <- dat[names(sort(apply(dat,1,mad),decreasing = T)[1:500]),]
M <- cor(dat_500) # Default method = "pearson", use pearson for intra-group
# M <- cor(dat_500, method = "spearman")  # Use spearman for inter-group, but pearson looks better based on results

# Prepare column annotation (ensure it matches expression matrix column names)
annotation_col <- sample_info[, c("condition", "replicate")]

correlation_heatmap <-pheatmap::pheatmap(M,
                                         show_rownames = T,
                                         angle_col=45,
                                         fontsize=7,
                                         annotation_col=annotation_col)

ggsave(correlation_heatmap,filename = 'prepare_counts/check_cor_top500.pdf',width = 7.5,height =6)

### ------------ step7. Sample Distance: hclust and Heatmap of the sample-to-sample distances ——————————------------------------------———--———— ###

# dist calculates distance between rows by default, so transposition is required
sampleDists <- dist(t(dat))
sampleDistMatrix <- as.matrix(sampleDists)
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)  # Select heatmap colors

# Plot distance heatmap
p0 <- pheatmap::pheatmap(sampleDistMatrix,
                         fontsize=7,
                         clustering_distance_rows=sampleDists,
                         clustering_distance_cols=sampleDists,
                         angle_col=45,
                         col=colors)
ggsave(p0,filename = 'prepare_counts/check_dist.pdf',width = 7.5,height =6)
dev.off()

# Plot cluster
pdf("prepare_counts/check_hclust.pdf")
plot(hclust(sampleDists))
dev.off()

### ------------ step8. Principal Component Analysis: PCA (DESeq2 native/PCA function) for sample grouping; A: Plot PC1-2-3 individually ——————————---------------———-----------———— ###

# (1) Calculate PCA
dat.pca <- PCA(t(dat), graph = FALSE)

# (2) Get variance percentage
eig.val <- get_eigenvalue(dat.pca)
percentVar <- eig.val[, "variance.percent"]

# (3) Define grouping
group_exp <- c("cyt", "cyt", "nuc", "nuc")
group_batch <- c("batch1", "batch2", "batch1", "batch2")

# (4) Define function to calculate dynamic coordinate limits -------------------------------------------------
get_dynamic_limits <- function(pca_result, pc_index, buffer_percent = 0.2) {
  # Get coordinate values for the specified PC
  pc_values <- get_pca_ind(pca_result)$coord[, pc_index]
  # Calculate range
  pc_range <- range(pc_values)
  range_diff <- diff(pc_range)
  # Add buffer (default 10%)
  buffer <- range_diff * buffer_percent
  limits <- c(pc_range[1] - buffer, pc_range[2] + buffer)
  return(limits)
}
# Define function to calculate dynamic coordinate limits -------------------------------------------------

# (5) Define plotting function -------------------------------------------------------------
plot_custom_pca_ggplot <- function(pc_x, pc_y, title) {
  # Get PCA coordinates
  pca_data <- as.data.frame(get_pca_ind(dat.pca)$coord)
  colnames(pca_data) <- paste0("PC", 1:ncol(pca_data))
  pca_data$Group <- group_exp
  pca_data$Batch <- group_batch

  # Dynamically calculate coordinate limits
  xlim <- get_dynamic_limits(dat.pca, pc_x)
  ylim <- get_dynamic_limits(dat.pca, pc_y)

  # Create plot
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
# Define plotting function ------------------------------------------------------------------

# (6) Start plotting
p1_2 <- plot_custom_pca_ggplot(1, 2, "PC1 vs PC2: Biological Variation")
p1_3 <- plot_custom_pca_ggplot(1, 3, "PC1 vs PC3: Batch Effect Check")
p2_3 <- plot_custom_pca_ggplot(2, 3, "PC2 vs PC3: Secondary Patterns")

# (7) Save plot
ggsave(p1_2, filename = 'prepare_counts/p1_2_PCA.pdf', width = 8, height = 7)
ggsave(p1_3, filename = 'prepare_counts/p1_3_PCA.pdf', width = 8, height = 7)
ggsave(p2_3, filename = 'prepare_counts/p2_3_PCA.pdf', width = 8, height = 7)


# (8) Display side-by-side
combined_1_2_3_PCA <- (p1_2 / p1_3 / p2_3)
ggsave(combined_1_2_3_PCA, filename = 'prepare_counts/combined_1_2_3_PCA.pdf', width = 8, height = 8)


### ------------ step8. Principal Component Analysis: PCA (DESeq2 native/PCA function) for sample grouping; B: Classic PC1-2 ——————————-------------------———-----------———— ###

# t() transpose, graph = F prevents automatic plot generation (manual customization follows)
dat.pca <- PCA(t(dat) , graph = F)

# Calculate variance contribution percentage of principal components (New key step)
eig.val <- get_eigenvalue(dat.pca)  # Get eigenvalues of the PCA result (variance of each principal component)
# Eigenvalues include multiple columns; here, we extract the percentage of variance explained by the principal components
percentVar <- eig.val[, "variance.percent"]
# Use the original data column names (i.e., sample names) as grouping information.
group_list <- ifelse(grepl("cyt", colnames(dat)), "cyt", "nuc")
# group_list <- c("cyt", "cyt", "nuc", "nuc")

# Print available number of principal components
n_pc <- ncol(dat.pca$ind$coord)
# Maximum number of principal components = min(number of samples - 1, number of genes) = min(3, 22378) = 3
# Therefore, the valid principal component indices can only be 1, 2, 3
# But if the data only has 2 effective PCs (certain special cases), PC3 will not exist

# Get sample coordinate data
# `get_pca_ind` (from `factoextra`) retrieves the coordinates of samples (individuals) on the principal components.
# `$coord` extracts the coordinate matrix; each row is a sample, each column is a principal component.
pca_data <- get_pca_ind(dat.pca)$coord
pc1_range <- range(pca_data[, 1])  # First column is PC1
pc2_range <- range(pca_data[, 2])  # Second column is PC2

# Calculate expansion ratio (add 20% boundary space on the plot)
# For aesthetic reasons, we want some whitespace at both ends of the axes. Here, 20% of the PC1 and PC2 range width is calculated as the expansion amount.
expand_factor <- 0.20
x_expand <- diff(pc1_range) * expand_factor
y_expand <- diff(pc2_range) * expand_factor

# Plot PCA and dynamically adjust range
pca <- fviz_pca_ind(dat.pca,
                    title = "Principal Component Analysis (PC1 vs PC2)",
                    legend.title = "Groups",
                    geom.ind = c("point", "text"), # Specify the geometric shape of the samples, displaying both points and text labels here
                    pointsize = 1.5, # Size of the points
                    labelsize = 2, # Size of the text labels
                    col.ind = group_list,  # Color by group
                    axes.linetype = NA,    # Remove axis lines
                    mean.point = FALSE     # Remove group center point
) +
  coord_fixed(ratio = 1) +  # Fix coordinate ratio to ensure equal unit length on x and y axes
  # Set x-axis label, showing PC1's variance contribution percentage
  xlab(paste0("PC1 (", round(percentVar[1], 1), "%)")) +
  # Set y-axis label, showing PC2's variance contribution percentage
  ylab(paste0("PC2 (", round(percentVar[2], 1), "%)")) +
  # Dynamically set coordinate axis range
  xlim(pc1_range[1] - x_expand, pc1_range[2] + x_expand) + # Set x-axis range, expanded by 20%
  ylim(pc2_range[1] - y_expand, pc2_range[2] + y_expand) + # Set y-axis range, expanded by 20%
  # Add black border
  theme(
    panel.border = element_rect(colour = "black", fill = NA, size = 0.5), # Add black border
    # plot.background = element_rect(colour = "black", size = 0.5)
  )
# Save plot - increase size
ggsave(pca, filename = 'prepare_counts/check_PCA.pdf', width = 8, height = 7)
dev.off()

### ------------ Conclusion --------------------------------------------------------——————————-------------------———-----------———— ###
print("Gene expression counts pre-processing is complete")

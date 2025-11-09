# Load necessary packages
library(dplyr)
library(readr)
library(tximport)
library(FactoMineR)
library(factoextra) 
library(tidyverse)
library(pheatmap)
library(RColorBrewer)
library(ggpointdensity)
library(patchwork) 
library(ggpubr)
library(gridExtra)

# Set working directory (Placeholder path)
setwd("~/RNA_Seq_Analysis/Salmon_Ensembl") 

# Create results directory
dir.create("prepare_counts", showWarnings = FALSE)

### ---------------------------- Step1: Read and inspect Salmon transcript-level count data ------------------------------------------------------- ###

#(1) Prepare sample information
samples <- c("K562.LabE.polyA.cyt.rep1", "K562.LabE.polyA.cyt.rep2", 
             "K562.LabE.polyA.nuc.rep1", "K562.LabE.polyA.nuc.rep2")
# Create sample information table
sample_info <- data.frame(
  sample = samples,
  condition = c("cytoplasm", "cytoplasm", "nucleus", "nucleus"),
  replicate = c("rep1", "rep2", "rep1", "rep2"),
  row.names = samples,
  stringsAsFactors = FALSE)

# Save sample information table
write.csv(sample_info, file = "prepare_counts/sample_info.csv", row.names = FALSE)

# (2) Prepare Salmon file paths
salmon_files <- file.path(samples, "quant.sf")
names(salmon_files) <- samples # Required

# Check if files exist
if(!all(file.exists(salmon_files))) {
  stop("Some quant.sf files are missing. Please check the file paths.")
}

### ---------------------------- Step2: Read Salmon gene-level data -------------------------------------------------------------------- ###

# (3) Import tx2gene mapping table
tx2gene_file <- "~/Data/annotation/gene_id_name_tx_id_map_Ensembl_rtracklayer.csv" # Placeholder path
if(!file.exists(tx2gene_file)) {
  stop("The tx2gene file does not exist. Please check the path.")
}

# gene_id, gene_name, transcript_id
tx2gene_ano <- as.data.frame(read_csv(tx2gene_file, col_names = TRUE, locale = locale(encoding = "UTF-8")))
# Select transcript_id and gene_id columns
tx2gene <- tx2gene_ano[,c("transcript_id", "gene_id")]

# (4) Import Salmon quantification data and summarize to gene level
txi <- tximport(salmon_files, 
                type = "salmon", 
                tx2gene = tx2gene,
                countsFromAbundance = "no")  # Choose uncorrected counts. This is set to calculate FPKM using DESeq2's internal frpkm() function later (not used for TPM calculation here).

### ---------------------------- Step3: Gene-level TPM and gene_name annotation --------------------------------------------------------------- ###

# (5) Extract gene_TPM; txi matrix contains "gene-level TPM", "gene length", and "gene count calculated from length & TPM";
# txi$abundances is gene-level TPM: the sum of TPMs of different transcript isoforms for the same gene
gene_TPM <- as.data.frame(txi$abundance)
gene_TPM$gene_id <- rownames(gene_TPM) 
# Ensure unique gene_id rows by removing potential duplicates (though tximport should handle this)
gene_TPM_unique <- gene_TPM %>%
  distinct(gene_id, .keep_all = TRUE)

# (6) Annotate gene_name using gene_id. If gene_name is missing, use gene_id instead.
tx2gene_name <- tx2gene_ano[,c("gene_id","gene_name")]

gene_name_TPM <- gene_TPM_unique %>%
  left_join(
    tx2gene_name %>% distinct(gene_id, .keep_all = TRUE), 
    by = "gene_id"
  ) %>%
  mutate(gene_name = coalesce(gene_name, gene_id)) %>% # If gene_name is present, keep it; if missing, use gene_id
  relocate(gene_id, gene_name) 

# (7). Check summary statistics
print(paste("Total number of genes:", nrow(gene_TPM_unique)))
print(paste("Number of genes with TPM > 0 across all samples:", sum(rowSums(gene_TPM_unique[, 1:4]) > 0)))

# (8). Save results
write.csv(gene_name_TPM, file = "prepare_counts/Salmon_gene_name_TPM.csv", row.names = FALSE)


### ---------------------------- Step4: Filter out low-expressed genes --------------------------------------------------------------------------- ###

#### (Filtering criteria is not unique and depends on the situation)
# (1) Filter: "at least 1 column has tpm > 0" 
# gene_name_TPM[, 3:6] > 0 creates a logical matrix checking if tpm is > 0 in each sample
keep_feature <- rowSums(gene_name_TPM[, 3:6] > 0) >= 1 # rowSums() sums the TRUE values (genes with TPM > 0) per row
table(keep_feature) # Check filtering results: FALSE is number of low-expressed genes (rows), TRUE is number of genes to keep

# Another cutoff: filter where the sum of TPM across all columns is > 0
counts_filter2  <- gene_name_TPM %>%
  filter(rowSums(select(., 3:6), na.rm = TRUE) > 0) # This should give the same result as keep_feature for this cutoff

# (2) Apply filter
gene_name_TPM_filter1 <- gene_name_TPM[keep_feature, ] # Replace TPM matrix with filtered genes (retaining higher expressed genes)

# (3) Check filtered count
print(paste("Number of retained genes:", nrow(gene_name_TPM_filter1) )) 
print(paste("Percentage of retained genes:", round((nrow(gene_name_TPM_filter1)/nrow(gene_name_TPM))*100), "%")) # Round to nearest integer

# (4) Save filtered data
write.csv(gene_name_TPM_filter1, file = "prepare_counts/Salmon_gene_name_TPM_filter1.csv", row.names = F)
# Used for combining plots

### ------------- Pre-differential analysis check - Data QC ----------------------------------------------------------------------------------------------- ###
# Pre-assess differences between samples: hclust plot, distance heatmap, PCA plot (major features), heatmap of top 500 variable genes, correlation heatmap

# Libraries are already loaded at the top

### ------------ Step5: Sample Normalization (Listing 3 normalization methods A-B-C) ------------------------------------------------------ ###

# A: dat <- as.data.frame(log2(edgeR::cpm(counts)+1)) # Simple normalization CPM: Counts per million
# B: DESeq2_normalize: rld <- rlogTransformation(dds, blind = FALSE) 
# C: log2(TPM+1) - Using this method here
dat <- log2((gene_name_TPM_filter1[, 3:6])+1)

### ------------ Step5: Boxplot to check overall gene expression distribution ---------------------------------------------------------- ###
# Plot using base R commands
num_samples <- ncol(dat)
color <- rainbow(num_samples) # Method 1: Rainbow color scheme (default)
# color <- heat.colors(num_samples) # Method 2: Heat color scheme
# color <- terrain.colors(num_samples) # Method 3: Terrain color scheme
# color <- c(rep("#1f77b4", 2), rep("#ff7f0e", 2)) # Manually specify colors: Cytoplasm reps 1-2: blue; Nucleus reps 1-2: orange

boxplot(dat, col=color, ylab="log2(TPM+1)", main="Normalized Data - Boxplot",
        outline = F, notch = F)
dev.copy(png, "prepare_counts/TMP_boxplot.png", width = 2000, height = 1500, res = 300)
dev.off()


### ------------ Step5: Boxplot to check overall gene expression distribution (ggplot) ------------------------------------------------- ###

# Plot using ggplot 

# 1. Prepare data
# Add gene ID (if dat does not have row names)
dat_with_id <- as.data.frame(dat)
dat_with_id$gene_id <- rownames(dat) 

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
group_colors <- c("Cytoplasm" = "#1f78b4", "Nucleus" = "#e31a1c") # Blue/Red scheme

# 4. Create ggplot boxplot
p <- ggplot(dat_long, aes(x = sample, y = tpm, fill = group)) +
  geom_boxplot(outlier.shape = NA, notch = FALSE) + # Do not show outliers
  scale_fill_manual(values = group_colors) +        # Set group colors
  labs(title = "Normalized TPM Expression Distribution",
       y = "log2(TPM+1)",
       x = "Samples",
       fill = "Fraction") +
  theme_minimal() +
  theme(
    # Custom x-axis labels
    axis.text.x = element_text(
      angle = 45,            # 45 degree tilt
      hjust = 1,             # Right alignment
      size = 12,             # Font size
      face = "italic",       # Italic font
      color = "black"        # Font color
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

### ------------ Step6. Gene Expression Correlation: A: Intra-group Pearson Correlation ------------------------------------------------- ###

# (1). Plot correlation for Cytoplasmic replicates (K562.LabE.polyA.cyt.rep1 vs rep2)
cyt_plot <- ggplot(as.data.frame(dat), 
                   aes(x = K562.LabE.polyA.cyt.rep1, 
                       y = K562.LabE.polyA.cyt.rep2)) +
  geom_pointdensity(size = 0.6) + # Density scatter plot
  scale_color_viridis_c(option = "C") + # Use viridis color scheme
  geom_smooth(method = "lm", color = "red", se = FALSE) + # Add regression line
  labs(title = "Cytoplasmic Replicates Correlation",
       x = "Rep 1 (log2(TPM+1))", 
       y = "Rep 2 (log2(TPM+1))") +
  theme_bw() +
  stat_cor(method = "pearson", # Automatically add Pearson correlation coefficient
           label.x.npc = "left",
           label.y.npc = "top")

# (2). Plot correlation for Nuclear replicates (K562.LabE.polyA.nuc.rep1 vs rep2)
nuc_plot <- ggplot(as.data.frame(dat), 
                   aes(x = K562.LabE.polyA.nuc.rep1, 
                       y = K562.LabE.polyA.nuc.rep2)) +
  geom_pointdensity(size = 0.6) +
  scale_color_viridis_c(option = "D") +
  geom_smooth(method = "lm", color = "red", se = FALSE) +
  labs(title = "Nuclear Replicates Correlation",
       x = "Rep 1 (log2(TPM+1))", 
       y = "Rep 2 (log2(TPM+1))") +
  theme_bw() +
  stat_cor(method = "pearson",
           label.x.npc = "left",
           label.y.npc = "top")

# Display plots side-by-side
cyt_plot + nuc_plot

# Save plots (optional)
ggsave("prepare_counts/cyt_replicates_corr.pdf", plot = cyt_plot, width = 6, height = 5)
ggsave("prepare_counts/nuc_replicates_corr.pdf", plot = nuc_plot, width = 6, height = 5)

# viridis color scales have five options (A to E), representing different color gradient schemes:
# option = "A": viridis (default, green to purple)
# option = "B": magma (deep purple to bright yellow)
# option = "C": plasma (deep purple to bright orange)
# option = "D": inferno (black to bright yellow)
# option = "E": cividis (deep blue to bright yellow)

### ------------ Step6. Gene Expression Correlation: B: Inter-group Pearson Correlation ------------------------------------------------- ###

# Create plotting function
create_density_plot <- function(df, x_col, y_col) {
  # Calculate Spearman correlation coefficient
  corr <- cor(df[[x_col]], df[[y_col]], method = "spearman")
  corr_label <- paste0("ρ = ", round(corr, 3))
  
  # The original code has two ggplot structures here, I will retain the density scatter plot structure which is used for saving
  
  # Create density scatter plot
  ggplot(df, aes(x = .data[[x_col]], y = .data[[y_col]])) +
    #geom_hex(bins = 50) + # Use hexagonal binning
    geom_pointdensity(size = 0.8, alpha = 0.5) +
    scale_fill_viridis_c(option = "plasma", trans = "log10") + # Use log color scale
    geom_smooth(method = "lm", color = "red", se = FALSE) + # Add regression line
    # Move label to top left corner
    annotate("text", x = -Inf, y = Inf, label = corr_label,
             hjust = -0.5, vjust = 4, size = 5, color = "black") +
    # Move label to bottom right corner (Original second option commented out)
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

# Create first set of plots (rep1)
p1 <- create_density_plot(dat, 
                          "K562.LabE.polyA.cyt.rep1", 
                          "K562.LabE.polyA.nuc.rep1")

# Create second set of plots (rep2)
p2 <- create_density_plot(dat, 
                          "K562.LabE.polyA.cyt.rep2", 
                          "K562.LabE.polyA.nuc.rep2")
# View plots
p1 + p2

# Save plots (optional)
ggsave("prepare_counts/samples_rep1_corr.png", plot = p1, width = 6, height = 5)
ggsave("prepare_counts/samples_rep2_corr.png", plot = p2, width = 6, height = 5) 

# Combine all scatterplots
correlation_scaterplot <- (cyt_plot + nuc_plot) / (p1 + p2) + 
  plot_annotation(title = "Comprehensive Correlation Analysis",
                  theme = theme(plot.title = element_text(size = 20, hjust = 0.5)))
print(correlation_scaterplot)
ggsave(correlation_scaterplot, filename = 'prepare_counts/correlation_scaterplot.png', width = 12, height = 8)
dev.off()

### ------------ Step6. Gene Expression Correlation: C: Correlation Heatmap (Top 500 High Expressed Genes) --------------------------------- ###

# Use the top 500 most variable genes to check correlation, avoiding noise from low-expressed genes
dat_500 <- dat[names(sort(apply(dat,1,mad),decreasing = T)[1:500]),] 
M <- cor(dat_500) # Default method = "pearson" (use Pearson for intra-group check)
# M <- cor(dat_500, method = "spearman") # Use Spearman for inter-group check, but Pearson often looks better here

# Prepare column annotation (ensure it matches the expression matrix column names)
annotation_col <- sample_info[, c("condition", "replicate")]

correlation_heatmap <-pheatmap::pheatmap(M,
                                         show_rownames = T,
                                         angle_col=45,
                                         fontsize=7,
                                         annotation_col=annotation_col) 

ggsave(correlation_heatmap,filename = 'prepare_counts/check_cor_top500.pdf',width = 7.5,height =6)

### ------------ Step7. Sample Distance: Hclust and Heatmap of the sample-to-sample distances ------------------------------------------------ ###

sampleDists <- dist(t(dat))    # dist defaults to calculating distance between rows, so transpose is needed
sampleDistMatrix <- as.matrix(sampleDists)  
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255) # Select colors for the heatmap

# Plot distance heatmap
p0 <- pheatmap::pheatmap(sampleDistMatrix,
                         fontsize=7,
                         clustering_distance_rows=sampleDists,
                         clustering_distance_cols=sampleDists,
                         angle_col=45,
                         col=colors)
ggsave(p0,filename = 'prepare_counts/check_dist.pdf',width = 7.5,height =6)
dev.off()

# Plot cluster dendrogram (hclust)
pdf("prepare_counts/check_hclust.pdf")
plot(hclust(sampleDists))
dev.off()

### ------------ Step8. Principal Component Analysis: PCA (DESeq2 internal / PCA function) to check sample grouping; A: PC1-2-3 plots ---------------- ###

# (1) Calculate PCA
dat.pca <- PCA(t(dat), graph = FALSE)

# (2) Get variance percentage
eig.val <- get_eigenvalue(dat.pca)
percentVar <- eig.val[, "variance.percent"]

# (3) Define groupings
group_exp <- c("cyt", "cyt", "nuc", "nuc")
group_batch <- c("batch1", "batch2", "batch1", "batch2")

# (4) Function to calculate dynamic coordinate limits -------------------------------------------------
get_dynamic_limits <- function(pca_result, pc_index, buffer_percent = 0.2) {
  # Get coordinate values for the specified PC
  pc_values <- get_pca_ind(pca_result)$coord[, pc_index]
  # Calculate range
  pc_range <- range(pc_values)
  range_diff <- diff(pc_range)
  # Add buffer (default 20%)
  buffer <- range_diff * buffer_percent
  limits <- c(pc_range[1] - buffer, pc_range[2] + buffer)
  return(limits)
}
# End of dynamic coordinate limits function ------------------------------------------------------

# (5) Plotting function -------------------------------------------------------------
plot_custom_pca_ggplot <- function(pc_x, pc_y, title) {
  # Get PCA coordinates
  pca_data <- as.data.frame(get_pca_ind(dat.pca)$coord)
  colnames(pca_data) <- paste0("PC", 1:ncol(pca_data))
  pca_data$Group <- group_exp
  pca_data$Batch <- group_batch
  
  # Calculate dynamic coordinate ranges
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
# End of plotting function ------------------------------------------------------------------

# (6) Start plotting
p1_2 <- plot_custom_pca_ggplot(1, 2, "PC1 vs PC2: Biological Variation")
p1_3 <- plot_custom_pca_ggplot(1, 3, "PC1 vs PC3: Batch Effect Check")
p2_3 <- plot_custom_pca_ggplot(2, 3, "PC2 vs PC3: Secondary Patterns")

# (7) Save plots 
ggsave(p1_2, filename = 'prepare_counts/p1_2_PCA.pdf', width = 8, height = 7)
ggsave(p1_3, filename = 'prepare_counts/p1_3_PCA.pdf', width = 8, height = 7)
ggsave(p2_3, filename = 'prepare_counts/p2_3_PCA.pdf', width = 8, height = 7)
dev.off()

# (8) Display combined plots
combined_1_2_3_PCA <- (p1_2 / p1_3 / p2_3)
print(combined_1_2_3_PCA)
ggsave(combined_1_2_3_PCA, filename = 'prepare_counts/combined_1_2_3_PCA.pdf', width = 8, height = 8)
dev.off()

### ------------ Step8. Principal Component Analysis: PCA (DESeq2 internal / PCA function) to check sample grouping; B: Classic PC1-2 ---------------- ###

# PCA plot using the PCA function
dat.pca <- PCA(t(dat) , graph = F) # t() transpose, graph = F prevents automatic plotting (for manual customization later)

# Calculate the percentage of variance explained by principal components (New critical step)
eig.val <- get_eigenvalue(dat.pca) # Get eigenvalues of the PCA result (variance for each PC)
percentVar <- eig.val[, "variance.percent"] # Eigenvalues contain multiple columns; extract the variance percentage explained by PCs
# group_list <- colnames(Salmon_tpm3) # Use the original data column names (sample names) as grouping information.
group_list <- ifelse(grepl("cyt", colnames(dat)), "cyt", "nuc")
# group_list <- c("cyt", "cyt", "nuc", "nuc")

# Print the number of available principal components
n_pc <- ncol(dat.pca$ind$coord)
# Max number of principal components = min(num_samples-1, num_genes) = min(3, 22378) = 3
# Therefore, valid PC indices can only be 1, 2, 3.
# However, if the data only has 2 effective PCs (in certain cases), PC3 will not exist.

# Get sample coordinate data
# `get_pca_ind` (from `factoextra`) retrieves the coordinates of samples (individuals) on the principal components.
# `$coord` extracts the coordinate matrix, where each row is a sample and each column is a principal component.
pca_data <- get_pca_ind(dat.pca)$coord
pc1_range <- range(pca_data[, 1]) # First column is PC1
pc2_range <- range(pca_data[, 2]) # Second column is PC2

# Calculate expansion factor (add 20% border space on the plot)
# For aesthetic purposes, we want to leave some white space at both ends of the axes.
# Calculate 20% of the range width for PC1 and PC2 as the expansion amount.
expand_factor <- 0.20
x_expand <- diff(pc1_range) * expand_factor  
y_expand <- diff(pc2_range) * expand_factor

# Plot PCA and dynamically adjust range
pca <- fviz_pca_ind(dat.pca,
                    title = "Principal Component Analysis (PC1 vs PC2)",
                    legend.title = "Groups",
                    geom.ind = c("point", "text"), # Specify geometric shapes for samples: show both points and text labels
                    pointsize = 1.5, # Size of points
                    labelsize = 2, # Size of text labels
                    col.ind = group_list, # Color by group
                    axes.linetype = NA,  # Remove axis lines
                    mean.point = FALSE   # Remove group center points
) + 
  coord_fixed(ratio = 1) + # Fix coordinate aspect ratio to ensure equal unit lengths for x and y axes
  xlab(paste0("PC1 (", round(percentVar[1], 1), "%)")) + # Set x-axis label, showing PC1's variance contribution percentage
  ylab(paste0("PC2 (", round(percentVar[2], 1), "%)")) + # Set y-axis label, showing PC2's variance contribution percentage
  # Dynamically set axis ranges
  xlim(pc1_range[1] - x_expand, pc1_range[2] + x_expand) + # Set x-axis range, expanded by 20%
  ylim(pc2_range[1] - y_expand, pc2_range[2] + y_expand) + # Set y-axis range, expanded by 20%
  # Add black border
  theme(
    panel.border = element_rect(colour = "black", fill = NA, size = 0.5) # Add black border
  )

pca
# Save plot - increase size
ggsave(pca, filename = 'prepare_counts/check_PCA.pdf', width = 8, height = 7)
dev.off()

### ------------ Conclusion ---------------------------------------------------- ###
print("Gene expression counts pre-processing is complete.")

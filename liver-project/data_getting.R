# load the required libraries
library(gplots)
library(ggplot2)
library(curatedTCGAData)
library(MultiAssayExperiment)
library(TCGAutils)
library(pheatmap)
library(SummarizedExperiment)
library(TCGAbiolinks)
library(DESeq2)
library(data.table)
library(dplyr)


######## 1 make an TCGA object #######
# build a query to retrieve gene expression data ------------
query_TCGA <- GDCquery(project = 'TCGA-LIHC',data.category = c('Transcriptome Profiling') ,experimental.strategy = 'RNA-Seq',
                       workflow.type = 'STAR - Counts', access = 'open')



# download data - GDCdownload
####### 2 download the object ########
GDCdownload(query_TCGA, directory = "/home/balbio/Hackbio-Internship")


####### 3 prepare data #########
liver.tcga.data <- GDCprepare(query_TCGA, summarizedExperiment = TRUE,directory = "/home/balbio/Hackbio-Internship")

meta.data <- as.data.frame(liver.tcga.data@colData@listData)
write.csv(meta.data, "LIHC_Clinical_data.csv")
table(meta.data$tumor_descriptor)
######## 4 get the assay ########
liver.raw <- assay(liver.tcga.data, 'stranded_first')

liver.raw <- as.data.frame(liver.raw)

# Assuming you have your data in a dataframe called df and the column names are the sample IDs
sample_ids <- colnames(liver.raw)

# Function to classify tumor or normal based on sample ID
classify_sample <- function(sample_id) {
  sample_type_code <- substr(sample_id, 14, 15)  # Extract 14th and 15th characters
  if (sample_type_code %in% c("01", "02")) {
    return("Tumor")
  } else if (sample_type_code %in% c("10", "11")) {
    return("Normal")
  } else {
    return("Other")
  }
}

# Apply the classify function to each sample ID
sample_info <- data.frame(
  Sample_ID = sample_ids,
  Sample_Type = sapply(sample_ids, classify_sample)
)

# Count number of normal and tumor samples
sample_counts <- table(sample_info$Sample_Type)



######## 5 convert Ensemble Ids to gene symbols ########

# Install and load org.Hs.eg.db package
library(org.Hs.eg.db)

# Convert Ensembl IDs to gene symbols
ensembl_ids <- rownames(liver.raw)

# Remove version numbers for compatibility
ensembl_ids <- sub("\\..*", "", ensembl_ids)

# Map IDs
gene_symbols <- mapIds(org.Hs.eg.db, keys = ensembl_ids, column = "SYMBOL", keytype = "ENSEMBL", multiVals = "first")

# Remove NA values and keep only mapped genes
mapped_genes <- !is.na(gene_symbols)
liver_matrix <- liver.raw[mapped_genes, ]
gene_symbols <- gene_symbols[mapped_genes]

unique_gene_symbols <- make.unique(gene_symbols, sep = "_")
length(unique_gene_symbols)

# Assign gene symbols to prostate_matrix
rownames(liver_matrix) <- unique_gene_symbols
#==================================================================================
# DESeq2 Analysis 
#==================================================================================
library(DESeq2)
# making the rownames and column names identical
all(rownames(sample_info) %in% colnames(liver_matrix))
all(rownames(sample_info) == colnames(liver_matrix))

expr_matrix <- round(liver_matrix)


# Create DESeqDataSet object
dds <- DESeqDataSetFromMatrix(countData = expr_matrix,
                              colData = sample_info,
                              design = ~ Sample_Type)
# Filter out low count genes
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

# Set the reference level for the factor of interest
dds$Sample_Type <- relevel(dds$Sample_Type, ref = "Normal")

# Run DESeq2 analysis
dds <- DESeq(dds)

# Extract results
res <- results(dds)

# Explore results
summary(res)

# Filter significant genes
sig_genes <- subset(res, padj < 0.05 & abs(log2FoldChange) > 2)

summary(sig_genes)
# Visualize results
plotMA(res, ylim=c(-3,3))


plotCounts(dds, gene=which.min(res$padj), intgroup="Sample_Type")

d <- plotCounts(dds, gene=which.min(res$padj), intgroup="Sample_Type", 
                returnData=TRUE)
library("ggplot2")
ggplot(d, aes(x=Sample_Type, y=count)) + 
  geom_point(position=position_jitter(w=0.1,h=0)) + 
  scale_y_log10(breaks=c(25,100,400))

write.csv(as.data.frame(sig_genes), 
          file="sample_type_differential_results.csv")

# Load the necessary libraries
library("pheatmap")
library("DESeq2")

# Assuming you have already created the DESeq2 object 'dds'

# Perform variance stabilizing transformation (VST) on the dataset
ntd <- vst(dds, blind=FALSE)

# Alternatively, you could use the rlog transformation, depending on your preference
# ntd <- rlog(dds, blind=FALSE)

# Select the top 20 most highly expressed genes
select <- order(rowMeans(assay(ntd)), decreasing=TRUE)[1:20]

# Create a dataframe for sample annotations
df <- as.data.frame(colData(dds)[, c("Sample_Type")])
rownames(df) <- rownames(colData(dds))
# Plot the heatmap using the transformed counts
pheatmap(assay(ntd)[select, ], 
         cluster_rows=FALSE, 
         show_rownames=FALSE, 
         cluster_cols=FALSE, 
         annotation_col=df)

vsd <- vst(dds, blind=FALSE)
plotPCA(vsd, intgroup=c("Sample_Type"))

#==============================================================================
#-----------01# BoxPLot for the 20 upregulated genes in tumor vs normal ----------
#------------------------------------------------------------------------------
#-=============================================================================
install.packages("tidyverse")
sample_info$Sample_ID
library(DESeq2)
library(tidyverse)
library(ggplot2)
library(tidyr)
# Assuming you've already run DESeq2 and have the results in 'res'

# Filter for upregulated genes (log2FoldChange > 1 and adjusted p-value < 0.05)
res_upregulated <- res[which(res$log2FoldChange > 2 & res$padj < 0.05),]

# Order by log2FoldChange
res_upregulated <- res_upregulated[order(res_upregulated$log2FoldChange, decreasing = TRUE),]

# Select top 18 upregulated genes
top_upregulated <- head(res_upregulated, 20)

# Extract normalized counts for these genes
normalized_counts <- counts(dds, normalized=TRUE)[rownames(top_upregulated),]

# Prepare data for plotting
plot_data <- normalized_counts %>%
  as.data.frame() %>%
  rownames_to_column("gene") %>%
  pivot_longer(cols = -gene, names_to = "sample", values_to = "count") %>%
  left_join(sample_info, by = c("sample" = "Sample_ID"))

# Add log2FoldChange and padj values
plot_data <- plot_data %>%
  left_join(as.data.frame(res_upregulated) %>% 
              rownames_to_column("gene") %>%
              dplyr::select(gene, log2FoldChange, padj),
            by = "gene")

# Create the plot
ggplot(plot_data, aes(x = Sample_Type, y = log2(count + 1), fill = Sample_Type)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.2, size = 0.5, alpha = 0.5) +
  facet_wrap(~ gene, scales = "free_y", ncol = 6) +
  scale_fill_manual(values = c("Normal" = "lightgray", "Tumor" = "darkgreen")) +
  labs(y = "log2CPM", x = NULL) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        strip.background = element_blank(),
        strip.text = element_text(face = "bold"),
        legend.position = "none") +
  geom_text(data = plot_data %>% group_by(gene) %>% slice(1),
            aes(x = Sample_Type[1], y = Inf, label = sprintf("q = %.2e", padj)),
            vjust = 1.5, size = 3, inherit.aes = FALSE)


# Save the plot
ggsave("top_20_upregulated_genes.png", width = 15, height = 10, dpi = 300)

write.csv(res_upregulated, "top_20_upregulated_genes.csv")


#==============================================================================
#----- #02 BoxPLot for the top 20  downregulated genes in tumor vs normal -----
#------------------------------------------------------------------------------
#-=============================================================================
library(DESeq2)
library(tidyverse)
library(ggplot2)

# Assuming you've already run DESeq2 and have the results in 'res'

# Filter for upregulated genes (log2FoldChange > 1 and adjusted p-value < 0.05)
res_downregulated <- res[which(res$log2FoldChange < 2 & res$padj < 0.05),]

# Order by log2FoldChange
res_downregulated <- res_downregulated[order(res_downregulated$log2FoldChange, decreasing = FALSE),]

# Select top 18 upregulated genes
top_downregulated <- head(res_downregulated, 20)

# Extract normalized counts for these genes
normalized_counts <- counts(dds, normalized=TRUE)[rownames(top_downregulated),]

# Prepare data for plotting
plot_data <- normalized_counts %>%
  as.data.frame() %>%
  rownames_to_column("gene") %>%
  pivot_longer(cols = -gene, names_to = "sample", values_to = "count") %>%
  left_join(sample_info, by = c("sample" = "Sample_ID"))

# Add log2FoldChange and padj values
plot_data <- plot_data %>%
  left_join(as.data.frame(top_downregulated) %>% 
              rownames_to_column("gene") %>%
              dplyr::select(gene, log2FoldChange, padj),
            by = "gene")

# Create the plot
ggplot(plot_data, aes(x = Sample_Type, y = log2(count + 1), fill = Sample_Type)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.2, size = 0.5, alpha = 0.5) +
  facet_wrap(~ gene, scales = "free_y", ncol = 6) +
  scale_fill_manual(values = c("Normal" = "lightpink", "Tumor" = "purple")) +
  labs(y = "log2CPM", x = NULL) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        strip.background = element_blank(),
        strip.text = element_text(face = "bold"),
        legend.position = "none") +
  geom_text(data = plot_data %>% group_by(gene) %>% slice(1),
            aes(x = Sample_Type[1], y = Inf, label = sprintf("q = %.2e", padj)),
            vjust = 1.5, size = 3, inherit.aes = FALSE)


# Save the plot
ggsave("top_20_downregulated_genes.png", width = 15, height = 10, dpi = 300)
write.csv(res_downregulated,"Downregulated_20genes.csv")


#==============================================================================
#----- #02 BoxPLot for ARNT2 gene in downregulated deseq2 result -----
#------------------------------------------------------------------------------
#-=============================================================================
library(DESeq2)
library(tidyverse)
library(ggplot2)

# Assuming you've already run DESeq2 and have the results in 'res'

# Filter for downregulated genes (log2FoldChange < -3 and adjusted p-value < 0.05)
res_downregulated <- res[which(res$log2FoldChange < -1 & res$padj < 0.05),]


normalized_counts <- counts(dds, normalized=TRUE)["ARNT2",]

# Prepare data for plotting
plot_data <- normalized_counts %>%
  as.data.frame() %>%
  rownames_to_column("gene") %>%
  pivot_longer(cols = -gene, names_to = "sample", values_to = "count") %>%
  left_join(coldata150, by = c("sample" = "sample"))

# Add log2FoldChange and padj values
plot_data <- plot_data %>%
  left_join(as.data.frame(res_filtered) %>% 
              rownames_to_column("gene") %>%
              dplyr::select(gene, log2FoldChange, padj),
            by = "gene")

ggplot(plot_data, aes(x = condition, y = log2(count + 1), fill = condition)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.2, size = 0.5, alpha = 0.5) +
  facet_wrap(~ gene, scales = "free_y", ncol = 6) +
  scale_fill_manual(values = c("normal" = "lightpink", "aml" = "purple")) +
  labs(y = "log2CPM", x = NULL) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        strip.background = element_blank(),
        strip.text = element_text(face = "bold"),
        legend.position = "none") +
  geom_text(data = plot_data %>% group_by(gene) %>% slice(1),
            aes(x = condition[1], y = Inf, label = sprintf("q = %.2e", padj)),
            vjust = 1.5, size = 3, inherit.aes = FALSE)

# Assuming you've already run DESeq2 and have the results in 'res'

# Filter for downregulated genes (log2FoldChange < -1 and adjusted p-value < 0.05)
res_downregulated <- res_filtered[which(res_filtered$log2FoldChange < -1 & res_filtered$padj < 0.05),]

# Find ARNT2 in the downregulated results
arnt2_result <- res_downregulated["ARNT2",]

# Extract normalized counts for ARNT2
normalized_counts <- counts(dds, normalized=TRUE)["ARNT2",]

# Prepare data for plotting
plot_data <- normalized_counts %>%
  as.data.frame() %>%
  rownames_to_column("gene") %>%
  pivot_longer(cols = -gene, names_to = "sample", values_to = "count") %>%
  left_join(coldata150, by = c("sample" = "sample"))

# Add log2FoldChange and padj values
plot_data <- plot_data %>%
  left_join(as.data.frame(res_filtered) %>% 
              rownames_to_column("gene") %>%
              dplyr::select(gene, log2FoldChange, padj),
            by = "gene")

# Create the plot
ggplot(plot_data, aes(x = condition, y = log2(count + 1), fill = condition)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.2, size = 1, alpha = 0.5) +
  scale_fill_manual(values = c("normal" = "lightpink", "aml" = "purple")) +
  labs(y = "log2CPM", x = NULL, title = "ARNT2 Expression") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        strip.background = element_blank(),
        strip.text = element_text(face = "bold"),
        legend.position = "none") +
  geom_text(data = plot_data %>% slice(1),
            aes(x = condition[1], y = Inf, 
                label = sprintf("log2FC = %.2f\nq = %.2e", log2FoldChange, padj)),
            vjust = 1.5, size = 3, inherit.aes = FALSE)

# Save the plot











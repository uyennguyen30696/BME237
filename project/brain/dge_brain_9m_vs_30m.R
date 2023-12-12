# Follow this tutorial https://www.reneshbedre.com/blog/edger-tutorial.html

library(edgeR)

# Read the count data from your text file
# 9 months
b9m_rep1 <- read.table('Documents/ucsc_bme/BME237_RNA_bioinformatics/project/brain/9m/rep1/rep1_counts.txt', header = FALSE, col.names = c("GeneID", "Counts1_9"))
b9m_rep2 <- read.table('Documents/ucsc_bme/BME237_RNA_bioinformatics/project/brain/9m/rep2/rep2_counts.txt', header = FALSE, col.names = c("GeneID", "Counts2_9"))
b9m_rep3 <- read.table('Documents/ucsc_bme/BME237_RNA_bioinformatics/project/brain/9m/rep3/rep3_counts.txt', header = FALSE, col.names = c("GeneID", "Counts3_9"))
# Remove the last 5 rows
b9m_rep1 <- b9m_rep1[1:(nrow(b9m_rep1) - 5), ]
b9m_rep2 <- b9m_rep2[1:(nrow(b9m_rep2) - 5), ]
b9m_rep3 <- b9m_rep3[1:(nrow(b9m_rep3) - 5), ]
# Merge tables based on the GeneID column
b9m_counts <- merge(b9m_rep1, b9m_rep2, by = "GeneID", all = TRUE)
b9m_counts <- merge(b9m_counts, b9m_rep3, by = "GeneID", all = TRUE)

# 30 months
b30m_rep1 <- read.table('Documents/ucsc_bme/BME237_RNA_bioinformatics/project/brain/30m/rep1/rep1_counts.txt', header = FALSE, col.names = c("GeneID", "Counts1_30"))
b30m_rep2 <- read.table('Documents/ucsc_bme/BME237_RNA_bioinformatics/project/brain/30m/rep2/rep2_counts.txt', header = FALSE, col.names = c("GeneID", "Counts2_30"))
b30m_rep3 <- read.table('Documents/ucsc_bme/BME237_RNA_bioinformatics/project/brain/30m/rep3/rep3_counts.txt', header = FALSE, col.names = c("GeneID", "Counts3_30"))
# Remove the last 5 rows
b30m_rep1 <- b30m_rep1[1:(nrow(b30m_rep1) - 5), ]
b30m_rep2 <- b30m_rep2[1:(nrow(b30m_rep2) - 5), ]
b30m_rep3 <- b30m_rep3[1:(nrow(b30m_rep3) - 5), ]
# Merge tables based on the GeneID column
b30m_counts <- merge(b30m_rep1, b30m_rep2, by = "GeneID", all = TRUE)
b30m_counts <- merge(b30m_counts, b30m_rep3, by = "GeneID", all = TRUE)

# Combine
count_data = merge(b9m_counts, b30m_counts, by = "GeneID", all = TRUE)
# Set 'GeneID' as row names
rownames(count_data) <- count_data$GeneID

# Create count matrix
counts_matrix <- count_data[, c("Counts1_9", "Counts2_9", "Counts3_9", "Counts1_30", "Counts2_30", "Counts3_30")]
counts_matrix

# Create a sample information for the count data
sample_info <- c("9m", "9m", "9m", "30m", "30m", "30m")

# Create DGEList data class for count and sample information
dge <- DGEList(counts = counts_matrix, group = factor(sample_info))

# Filter out genes with low counts
keep <- filterByExpr(y = dge)
dge <- dge[keep, , keep.lib.sizes=FALSE]

# Normalize
dge <- calcNormFactors(object = dge)

# Estimate dispersion
dge <- estimateDisp(y = dge)

# Find DGE
et <- exactTest(object = dge)

# Add ajusted p values
top_degs = topTags(object = et, n = "Inf")

# Print the reference of the comparision
topTags(et)$comparison

# Get a summary DGE table (returns significant genes with absolute log fold change at least 1 and adjusted p value < 0.05)
summary(decideTests(object = et, lfc = 1))

# Export to CSV file 
write.csv(as.data.frame(top_degs), file="Documents/ucsc_bme/BME237_RNA_bioinformatics/project/brain/dge_brain_9m_vs_30m.csv")



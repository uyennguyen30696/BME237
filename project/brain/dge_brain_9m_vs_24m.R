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

# 24 months
b24m_rep1 <- read.table('Documents/ucsc_bme/BME237_RNA_bioinformatics/project/brain/24m/rep1/rep1_counts.txt', header = FALSE, col.names = c("GeneID", "Counts1_24"))
b24m_rep2 <- read.table('Documents/ucsc_bme/BME237_RNA_bioinformatics/project/brain/24m/rep2/rep2_counts.txt', header = FALSE, col.names = c("GeneID", "Counts2_24"))
b24m_rep3 <- read.table('Documents/ucsc_bme/BME237_RNA_bioinformatics/project/brain/24m/rep3/rep3_counts.txt', header = FALSE, col.names = c("GeneID", "Counts3_24"))
# Remove the last 5 rows
b24m_rep1 <- b24m_rep1[1:(nrow(b24m_rep1) - 5), ]
b24m_rep2 <- b24m_rep2[1:(nrow(b24m_rep2) - 5), ]
b24m_rep3 <- b24m_rep3[1:(nrow(b24m_rep3) - 5), ]
# Merge tables based on the GeneID column
b24m_counts <- merge(b24m_rep1, b24m_rep2, by = "GeneID", all = TRUE)
b24m_counts <- merge(b24m_counts, b24m_rep3, by = "GeneID", all = TRUE)

# Combine
count_data = merge(b9m_counts, b24m_counts, by = "GeneID", all = TRUE)
# Set 'GeneID' as row names
rownames(count_data) <- count_data$GeneID

# Create count matrix
counts_matrix <- count_data[, c("Counts1_9", "Counts2_9", "Counts3_9", "Counts1_24", "Counts2_24", "Counts3_24")]
counts_matrix

# Create a sample information for the count data
sample_info <- c("9m", "9m", "9m", "24m", "24m", "24m")

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
write.csv(as.data.frame(top_degs), file="Documents/ucsc_bme/BME237_RNA_bioinformatics/project/brain/dge_brain_9m_vs_24m.csv")




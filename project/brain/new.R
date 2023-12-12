# Follow this tutorial https://www.reneshbedre.com/blog/edger-tutorial.html

library(edgeR)

# Read the count data from your text file
# 9 months
l9m_rep1 <- read.table('Documents/ucsc_bme/BME237_RNA_bioinformatics/project/liver/9m/rep1/rep1_counts.txt', header = FALSE, col.names = c("GeneID", "Counts1_9"))
l9m_rep2 <- read.table('Documents/ucsc_bme/BME237_RNA_bioinformatics/project/liver/9m/rep2/rep2_counts.txt', header = FALSE, col.names = c("GeneID", "Counts2_9"))
l9m_rep3 <- read.table('Documents/ucsc_bme/BME237_RNA_bioinformatics/project/liver/9m/rep3/rep3_counts.txt', header = FALSE, col.names = c("GeneID", "Counts3_9"))
# Remove the last 5 rows
l9m_rep1 <- l9m_rep1[1:(nrow(l9m_rep1) - 5), ]
l9m_rep2 <- l9m_rep2[1:(nrow(l9m_rep2) - 5), ]
l9m_rep3 <- l9m_rep3[1:(nrow(l9m_rep3) - 5), ]
# Merge tables based on the GeneID column
l9m_counts <- merge(l9m_rep1, l9m_rep2, by = "GeneID", all = TRUE)
l9m_counts <- merge(l9m_counts, l9m_rep3, by = "GeneID", all = TRUE)

# 24 months
l24m_rep1 <- read.table('Documents/ucsc_bme/BME237_RNA_bioinformatics/project/liver/24m/rep1/rep1_counts.txt', header = FALSE, col.names = c("GeneID", "Counts1_24"))
l24m_rep2 <- read.table('Documents/ucsc_bme/BME237_RNA_bioinformatics/project/liver/24m/rep2/rep2_counts.txt', header = FALSE, col.names = c("GeneID", "Counts2_24"))
l24m_rep3 <- read.table('Documents/ucsc_bme/BME237_RNA_bioinformatics/project/liver/24m/rep3/rep3_counts.txt', header = FALSE, col.names = c("GeneID", "Counts3_24"))
# Remove the last 5 rows
l24m_rep1 <- l24m_rep1[1:(nrow(l24m_rep1) - 5), ]
l24m_rep2 <- l24m_rep2[1:(nrow(l24m_rep2) - 5), ]
l24m_rep3 <- l24m_rep3[1:(nrow(l24m_rep3) - 5), ]
# Merge tables based on the GeneID column
l24m_counts <- merge(l24m_rep1, l24m_rep2, by = "GeneID", all = TRUE)
l24m_counts <- merge(l24m_counts, l24m_rep3, by = "GeneID", all = TRUE)

# Combine
count_data = merge(l9m_counts, l24m_counts, by = "GeneID", all = TRUE)
nrow(count_data)
nrow(l9m_counts)
nrow(l24m_counts)


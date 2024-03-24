setwd("/Volumes/USB/learn_transcriptomics")
library("dplyr")
library("knitr")
library(ggplot2)
library("magrittr")
library(tibble)
library(edgeR)

# Read in the count matrix
counts <- as.matrix(read.csv("dge_counts.csv", row.names = 1))
metadata_1 = read.csv("metadata.csv",row.names = 1)
# Identify the row name to remove
row_name_to_remove <- "LG45"  # Specify the row name you want to remove

# Check if the row name exists in the row names of the data frame
if (row_name_to_remove %in% rownames(metadata_1)) {
  # Remove the corresponding row from metadata
  metadata_1 <- metadata_1[rownames(metadata_1) != row_name_to_remove, , drop = FALSE]
} else {
  # If the row name does not exist, print a message
  message(paste("Row name", row_name_to_remove, "not found. No rows removed."))
}
dge <- DGEList(counts=counts, group=metadata_1$genotype, samples = metadata_1)
## filterByExpr function in edgeR is designed to filter out genes that are not expressed at a significant level across the samples in a dataset.
keep <- filterByExpr(dge)
dge <- dge[keep, , keep.lib.sizes=FALSE]
# Normalizing the Data
#Normalization is crucial for removing differences in library sizes and other technical variations.
dge <- calcNormFactors(dge)


# Design Matrix and Model Fitting
#Creating a design matrix reflecting the experimental setup and fit a model to the normalized counts.
#By defining the levels and their order, we are implicitly setting the first level as the baseline category against which comparisons will be made in our analysis.
dge$samples$genotype <- factor(dge$samples$genotype, 
                               levels = c("tt8fae1", "da1tt8fae1", 
                                          "dar1tt8fae1","da1dar1tt8fae1" ,"upl3tt8fae1"))

design <- model.matrix(~ genotype, data = dge$samples)

dge <- estimateDisp(dge, design)
fit <- glmQLFit(dge, design)

# Performing the Differential Expression Analysis
# glmQLFTest to test for differential expression.
# da1tt8fae1 compared to the baseline tt8fae1, we would use coef=2
qlf_1 <- glmQLFTest(fit, coef=2)
# next one  
#"dar1tt8fae1" to control 
qlf_2 <- glmQLFTest(fit, coef=3)

# next one "da1dar1tt8fae1" 
qlf_3 <- glmQLFTest(fit, coef=4)

# next one 
# "upl3tt8fae1"
qlf_4 <- glmQLFTest(fit, coef=5)


## saving important files 
save(fit, file="edgeR_fit_object.RData")
## qlf file 
save(qlf, file="edgeR_qlf_object.RData")


# Function to save qlf results to CSV
## for a single file 
results1 <- topTags(qlf_1, n=Inf)$table
write.csv(results1, "DGE_results1.csv")

##. for all file at once make a function 
save_qlf_results <- function(qlf_list, prefix="DGE_results") {
  for (i in seq_along(qlf_list)) {
    qlf_obj <- qlf_list[[i]]
    results <- topTags(qlf_obj, n=Inf)$table
    file_name <- paste(prefix, i, ".csv", sep="")
    write.csv(results, file_name)
  }
}

# Assuming qlf_1, qlf_2, qlf_3, qlf_4 exist and are your qlf objects
# Generate the names of these qlf objects
qlf_names <- paste("qlf_", 1:4, sep="")

# Retrieve the qlf objects by their names and create a list
qlf_list <- mget(qlf_names)

# Save the results of each qlf analysis to a CSV file
save_qlf_results(qlf_list)



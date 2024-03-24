library(ggplot2)
library("magrittr")
library(tibble)
library(edgeR)

library(ggplot2)

plot_volcano_ggplot <- function(file_path) {
  results <- read.csv(file_path)
  # Assume we have a data frame called 'results' with columns 'logFC', 'PValue', and 'FDR'
  threshold_logFC <- 2 # Threshold for log fold change
  threshold_pval <- 0.05 # Threshold for P-value
  threshold_fdr <- 0.05
  # Create a 'Significant' column based on both logFC and PValue/FDR
  # Creating the 'Significant' column based on logFC threshold and P-value
  # Create a vector to store the classification
  results$Significant <- ifelse(results$logFC <= -2 & results$PValue <=0.05, "Down", 
                           ifelse(results$logFC >= 2 & results$PValue <=0.05, "Up", "Normal"))
  # Assume 'Significant' is a column in your results that is TRUE for significant genes
  # and 'Up' is TRUE for upregulated genes and FALSE for downregulated genes.
  # If not, you would need to create this column based on your criteria for significance and direction of change.

  ggplot(results, aes(x = logFC, y = -log10(PValue), color = Significant)) +
    geom_point(size = 2) +
    scale_color_manual(values = c("Down" = "blue", "Up" = "red", "Normal" = "grey")) +
    labs(x="Log2 Fold Change", y="-Log10 P-value", title=paste("Volcano Plot for", basename(file_path))) +
    geom_hline(yintercept=-log10(0.05), linetype="dashed", color="black") +
    geom_vline(xintercept=2,color="red",linetype="dashed") +  # Threshold for log2 fold change
    geom_vline(xintercept=-2,color="blue",linetype="dashed") +  # Threshold for log2 fold change
    
    theme_minimal() +
    guides(color=guide_legend(override.aes=list(alpha=1)))

}

# Usage
plot_volcano_ggplot("DGE_results1.csv")
output_file <-plot_volcano_ggplot("DGE_results1.csv")



# Save the plot as PDF
ggsave(filename = paste0("volcano_plot_", "DGE_1", ".pdf"),
       plot = output_file,
       device = "pdf")


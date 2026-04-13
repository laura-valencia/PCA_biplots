
# This script is used to produce PCA biplots, in which differential gene expression (the variables) are shown as vectors.
# These are useful to capture variance and observe which genes add most.

library(tidyverse)

# Read raw-data (counts obtained from feature counts. These will be used to obtain a dds object in DESeq2, only to create a PCA biplot later.
# Used in step 1: Create dds object
data_cc <- read.csv(file= "data/gene_matrix_counts.csv")

# make a vector of column names per sample using colnames
#colnames(data_cc)

CC_colnames <-  c("C1N1B",  "C1N5B",  "C1N4B",  "C2A1B",  "C1A4B",  "C1A3B",
                       "C2N5B",  "C2N1B",  "C2N3B",  "C2A5B",  "C1A1B",  "C2A3B" )


# Read samples
# These are the DEG obtained from DESeq2 in R (after running a model)
# Used in step 2: plotting 
DEG_CC <- read.csv(file="DEG_filtered/DESeq2_output_CC_species") 
# this was obtained with my previous model with more details.



### PCA biplot ###
# Step 1: Data normalization
# need dds objects

# Design: Condition only. We run 
get_dds <- function(col.names, data){ 
  #libraries
  library("apeglm")
  library("DESeq2")
  
  
  # Create metadata file
  # add to a vector
  column1 <- col.names
  
  # Extract the last character of each element in the first column
  column2 <- substring(column1, 3, 3)
  
  # Extract the second character of each element in the first column
  column3 <- substring(column1, 2, 2)
  
  # Create a dataframe using the three vectors
  metaData_df <- data.frame(sample = column1, condition = column2, days = column3)
  
  
  #Defining dds object from matrix
  
  dds <- DESeqDataSetFromMatrix(countData=data, 
                                colData=metaData_df, 
                                design=~condition, tidy = TRUE)
  
  
  # pre filtering
  
  # pre-filtering: Here we perform pre-filtering to keep only rows that have a count 
  # of at least 10 for a minimal number of samples. 
  # The count of 10 is a reasonable choice for bulk RNA-seq.
  
  
  smallestGroupSize <- 3
  keep <- rowSums(counts(dds) >= 10) >= smallestGroupSize
  dds <- dds[keep,]
  
  
  # Estimate Normalization factors
  
  dds <- estimateSizeFactors(dds)
  
  # First, check the current levels of the factor
  levels(dds$condition)
  
  ####  Make and format DESeqObject and results object. 
  ## set Normoxia as Reference level
  dds$condition <- relevel(dds$condition, ref = "N") 
  
  
  # After, check the current levels of the factor
  levels(dds$condition)
  
  # running DESeq
  # DESes2 does default normaixation by median of ratios method
  
  dds <- DESeq(dds)
  
  # For each species, identify DEGs associated with the treatment by using contrast
  # select cutoff for pvalue here with alpha
  
  
  return(dds)
}


# get dds.

dds_CC <- get_dds(CC_colnames, data_cc)

# get counts, This uses DESeq2 dunction counts()
# Extract normalized counts
CC_counts <- counts(dds_CC, normalized=TRUE)

# Restart the R session before next step because DESeq2 creates conflicts.

## Step 2: Prepare top-genes list, those genes that contribute most to variation, as identified by PCA biplot
# This function also creates the PCA biplot.
# Notice that the PCA biplot needs a transposition because it will represent individual genes as vectors.
# This function will also produce a table with the top 5 % genes that contribute to variation, saved to folder.

make_biplot_PCA <- function(dds, normalized_counts, results_dataframe, col.names, label){
  

  #test
  #normalized_counts <- CC_counts
  #results_dataframe <- DEG_CC
  
  # Obtain Metadata df, the same used for dds
  
  #add to a vector
  column1 <- col.names
  
  # Extract the last character of each element in the first column
  column2 <- substring(column1, 3, 3)
  
  # Extract the second character of each element in the first column
  column3 <- substring(column1, 2, 2)
  
  # Create a dataframe using the three vectors
  metaData_df <- data.frame(sample = column1, condition = column2, days = column3)
  
  
  # Filter normalized counts to keep only significant genes
  significant_genes <- results_dataframe$Symbol
  normalized_counts_filtered <- normalized_counts[significant_genes, ]
  
  
  # Scale and center the data
  pca_input <- t(scale(normalized_counts_filtered))  # Transpose because PCA works on samples
  pca_result <- prcomp(pca_input, center=TRUE, scale=TRUE)
  
  # Get scores for the first two principal components
  scores <- as.data.frame(pca_result$x[, 1:2])
  # Add sample identifiers if needed
  scores$sample <- rownames(scores)
  
  # Merge scores with metadata
  scores_with_metadata <- merge(scores, metaData_df, by = "sample")
  
  # Get loadings (genes' contribution to principal components)
  loadings <- pca_result$rotation[, 1:2]  # First two PCs
  
  # Convert loadings to a format that can be visualized
  loadings_df <- as.data.frame(loadings)
  loadings_df$Gene <- rownames(loadings_df)
  
  # Filter loadings for significant genes
  loadings_df_filtered <- loadings_df[rownames(loadings_df) %in% significant_genes, ]
  

  # Calculate absolute loadings
  loadings_df_filtered$AbsLoading <- abs(loadings_df_filtered$PC1)  # Use PC1 or PC2 as needed
  
  # Determine the number of top genes to extract (5%)
  number_of_genes <- nrow(loadings_df_filtered)
  top_gene_count <- ceiling(number_of_genes * 0.05)  # Calculate 5% of total genes
  
  # Get the top 5% based on absolute loadings
  top_genes <- loadings_df_filtered[order(loadings_df_filtered$AbsLoading, decreasing=TRUE), ][1:top_gene_count, ]
  
  # Store only the gene names or IDs for easy reference
  top_gene_symbols <- top_genes$Gene
  
  # Display the top genes
  print(top_genes)
  
  write.csv(top_genes, file=paste("biplot_PCA_top_genes_",label, sep=""))
  
  
  ## PCA plotting
  
  # Using PCA and loadings computed above
  
  # Load required libraries
  library(ggplot2)
  
  # Get absolute loadings for PC1
  loadings_df_filtered$AbsLoading <- abs(loadings_df_filtered$PC1)
  
  # Calculate 5% of significant genes
  number_of_genes <- nrow(loadings_df_filtered)
  top_gene_count <- ceiling(number_of_genes * 0.05)
  
  # Select top 5% genes
  top_genes <- loadings_df_filtered[order(loadings_df_filtered$AbsLoading, decreasing=TRUE), ][1:top_gene_count, ]
  
  # Store gene symbols for plotting
  top_gene_symbols <- top_genes$Gene
  
 
  ## Adding noise between vector lines for improved readability
  
  # Setting a random seed for reproducibility
  set.seed(123) 
  
  # Create a jittered version of the loadings for non-top genes
  loadings_jittered <- loadings_df_filtered
  loadings_jittered$PC1 <- loadings_jittered$PC1 + runif(nrow(loadings_jittered), -0.3, 0.3)  # Adjust the range as needed
  loadings_jittered$PC2 <- loadings_jittered$PC2 + runif(nrow(loadings_jittered), -0.3, 0.3)  # Adjust the range as needed
  
  # Create a jittered version of the top genes for additional visual separation
  top_genes_jittered <- top_genes
  top_genes_jittered$PC1 <- top_genes_jittered$PC1 + runif(nrow(top_genes_jittered), -0.3, 0.3)
  top_genes_jittered$PC2 <- top_genes_jittered$PC2 + runif(nrow(top_genes_jittered), -0.3, 0.3)
  
  # Create the PCA biplot
  ggplot() +
    # Plot all gene vectors in gray with jitter
    geom_segment(data=loadings_jittered, 
                 aes(x=0, y=0, xend=PC1, yend=PC2), 
                 arrow=arrow(type="closed", length=unit(0.1, "inches")), 
                 color="gray", size=0.5) +
    # Highlight top genes in black with jitter
    geom_segment(data=top_genes_jittered, 
                 aes(x=0, y=0, xend=PC1, yend=PC2), 
                 arrow=arrow(type="closed", length=unit(0.1, "inches")), 
                 color="black", size=0.8) +  # Adjust the size if needed
    # Add labels for top genes
    geom_text(data=top_genes_jittered, aes(x=PC1, y=PC2, label=Gene), 
              vjust=1, hjust=1, size=3, color="black") +
    theme_minimal() +
    labs(title=paste("PCA Biplot - Top Genes Highlighted with Jitter: ", label), 
         x="Principal Component 1", 
         y="Principal Component 2")
  
  
}


biplot_CC_B_2 <- make_biplot_PCA(dds_CC, CC_counts, DEG_CC, CC_colnames, "plotting_label_title")



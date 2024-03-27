## Draw genes from N=5.5K samples using gggenes

# Load necessary libraries
library(gggenes)
library(readr)
library(ggplot2)
library(RColorBrewer)

display.brewer.all(colorblindFriendly = TRUE)

# Read the gene details from the file
gene_data1 <- read_delim("output.tsv", delim = "\t", show_col_types = FALSE)
gene_data1$gene_label[gene_data1$gene_label == "Beta-lactamase CTX-M-1"] <- "Beta-lactamase CTX-M-15"
gene_data1$gene_label[gene_data1$gene_label == "Periplasmic murein peptide-binding protein"] <- "Periplasmic murein peptide-binding protein mppA"
#gene_data1$gene_label[gene_data1$gene_label == "hypothetical protein"] <- "Spacer sequence"

# Aggregate the gene data by cluster name
aggregated_gene_data <- aggregate(cbind(gene_start, gene_end) ~ cluster_name + gene_label, data = gene_data1, FUN = function(x) c(min(x), max(x)))

# Plot the collapsed genes for each cluster name
gene_plot <- ggplot(aggregated_gene_data, aes(xmin = gene_start[,1], xmax = gene_end[,1], y = cluster_name, fill = gene_label)) +
  geom_gene_arrow(arrowhead_width = unit(0.1, "inches")) +
  scale_y_discrete(expand = c(0.02, 0.02)) +
  scale_fill_brewer(palette="Set3") +
  theme_genes() +
  labs(y = "Sample Name", x = "Position", fill = "Gene/Gene product")


## Anchor the alignments on blaCTX-M-15

# Set larger plot dimensions
#plot_width <- 8  # Width in inches
#plot_height <- 6  # Height in inches

dummies <- make_alignment_dummies(
  gene_data1,
  aes(xmin = gene_start, xmax = gene_end, y = cluster_name, id = gene_label),
  on = "Beta-lactamase CTX-M-15"
)

aligned_genes <- ggplot(gene_data1, aes(xmin = gene_start, xmax = gene_end, y = cluster_name, fill = gene_label)) +
  geom_gene_arrow() +
  geom_blank(data = dummies) +
  facet_wrap(~ cluster_name, scales = "free", ncol = 1) +
  #scale_fill_brewer(palette = "RdYlBu") +
  scale_color_gradient(low = "blue", high = "red") +
  theme_genes()

enlarge_output_pane <- function(height. = 700, width. = 1300){
  
  # detect current output pane dimensions
  dim_px <- grDevices::dev.size("px")
  
  # if detected width is less than 'width.', widen output pane with RStudio's layoutZoomRightColumn command
  if(dim_px[1] < width.){ rstudioapi::executeCommand("layoutZoomRightColumn") }
  
  # if detected height is less than 'height.', switch to Viewer windom and set height
  if(dim_px[2] < height.){
    
    viewer <- getOption("viewer")
    viewer("http://localhost:8100", height = height.)
    
  }
  
}

enlarge_output_pane

ggsave("aligned_genes_plot.pdf", aligned_genes, limitsize = FALSE)

dev.new(width = x, height = y)

#dev.off()


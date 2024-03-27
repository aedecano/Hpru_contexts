## Plot root-to-tip distances in the final timetree vs Year of isolation

rt <- read.csv("dist_root-to-tip.tsv", sep = "\t", header = TRUE)
colnames(rt) <- c("ID","Root_to_tip_distance")
metad <- read.csv("HPRU_metadata-date_context_n5516.tsv", sep = "\t", header = FALSE)
colnames_metad <- c("ID", "Year", "Context", "MLST",	"Host",	"Specimen",	"Country", "Sample")
colnames(metad) <- colnames_metad

rtWithmetad <- data.frame(inner_join(metad, rt, by = "ID"))
#View(rtWithmetad)

#write.table(rtWithmetad, file = "Joined_roottotip_year.tsv", sep = "\t", row.names = FALSE)

#Drop entries with unknown year

#rtWithmetad <- rtWithmetad[rtWithmetad$Year != "Unknown", ]

# Ensure Year is numeric
rtWithmetad$Year <- as.numeric(as.character(rtWithmetad$Year))

# Sort data by Year
rtWithmetad <- rtWithmetad[order(rtWithmetad$Year), ]

# Plot
cbPalette <- c("#999999", "#E69F00")

ggplot(rtWithmetad, aes(y = factor(Year), x = Root_to_tip_distance, color = Context)) +
  geom_point() +
  scale_y_discrete(labels = function(y) as.character(y)) +
  theme_light() +
  scale_color_manual(values=cbPalette) +
  labs(y = "Year", x = "Root-to-tip distance", color = "Context") +
  ggtitle("Root-to-tip distance vs. Year")

##
# Sample data representing root-to-tip distances
root_to_tip_distances <- rt$Root_to_tip_distance

# Compute pairwise distances (optional, if not already computed)
dist_matrix <- dist(root_to_tip_distances)

# Perform hierarchical clustering
hc <- hclust(dist_matrix)

# Plot the dendrogram
plot(hc, main = "Dendrogram of Root-to-Tip Distances")

# Cut the dendrogram to obtain clusters
# You can specify the desired number of clusters or a height threshold
 clusters <- cutree(hc, k = 5)
 
 # Now, let's visualize the clusters in a phylogenetic tree
 # Replace this with your actual phylogenetic tree visualization code
 # Here, we use a simple example of plotting a dendrogram
 # Plot the dendrogram
 plot(hc, main = "Dendrogram with Clustered Tips")
 
 # Color the tips according to their clusters
 # Here, we assume clusters is a vector containing cluster assignments
 tiplabels <- rep(1, length(clusters))
 tiplabels[clusters == 2] <- 2
 tiplabels[clusters == 3] <- 3
 # You can continue this pattern for more clusters if needed
 

 
 
 
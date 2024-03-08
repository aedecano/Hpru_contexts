#install.packages("data.table")
# Load ggplot2 library
library(ggplot2)
library(stringr)
library(dplyr)
library(tidyverse)
library(data.table)
library(vtable)
library(ggpubr)
library(gridExtra)
library(RColorBrewer)

## Colour-blind pallettes
# The palette with grey:
cbPalette <- c("#999999", "#E69F00")
#cbPalette1 <- c("#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

# The palette with black:
#cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

display.brewer.all(colorblindFriendly = TRUE)


# To use for fills, add
#scale_fill_manual(values=cbPalette)

# To use for line and point colors, add
#scale_colour_manual(values=cbPalette)

##Unpub with context
contig_data <- read.table("contigLength_n117plus1_Unpub.tsv", header = TRUE, sep = "\t")
context <- read.table("context_unpub.tsv", header = TRUE, sep = "\t")
merged_data <- dplyr::full_join(contig_data, context, by = "Fasta_File_Name")
#View(merged_data)
write.table(merged_data, file = "context_contigLength_n117plus1_Unpub_ctxm15positive.tsv", sep = "\t", row.names = FALSE)
sorted_contig_data <- merged_data[order(contig_data$Contig_length), ]
sorted_contig_data$Context[sorted_contig_data$Context == "chrom"] <- "chromosome"

# Plotting a bar plot of contig lengths

ggplot(sorted_contig_data , aes(x = reorder(Fasta_File_Name, -Contig_length), y = Contig_length, fill = Context)) + 
  geom_bar(stat = "identity") +
  #geom_vline(xintercept = nrow(sorted_contig_data) - 90, linetype = "dashed", color = "blue", size = 0.5) +
  labs(title = "Distribution of Contig Lengths (Unpublished, N=117)",
       y = "Contig Length",
       x = "Sample") +
  theme(axis.text.x  = element_blank(),  # Remove x-axis labels
        axis.ticks.x = element_blank())


### Enterobase with context
contig_data <- read.table("context_contigLength_enterobase_n404_ctxm15positive.tsv", header = TRUE, sep = "\t")

sorted_contig_data <- contig_data[order(contig_data$Contig_length), ]
sorted_contig_data$Context[sorted_contig_data$Context == "chrom"] <- "chromosome"

# Plotting a bar plot of contig lengths

ggplot(sorted_contig_data , aes(x = reorder(Fasta_File_Name, -Contig_length), y = Contig_length, fill = Context)) + 
  geom_bar(stat = "identity") +
  #geom_vline(xintercept = nrow(sorted_contig_data) - 90, linetype = "dashed", color = "blue", size = 0.5) +
  labs(title = "Distribution of Contig Lengths (Enterobase, N=404)",
       y = "Contig Length",
       x = "Sample") +
  theme(axis.text.x  = element_blank(),  # Remove x-axis labels
        axis.ticks.x = element_blank())

### Kev - with context
contig_data <- read.table("contigLength_n16_Kev_ctxm15positive.tsv", header = TRUE, sep = "\t")

sorted_contig_data <- contig_data[order(contig_data$Contig_length), ]
sorted_contig_data$Context[sorted_contig_data$Context == "chrom"] <- "chromosome"

# Plotting a bar plot of contig lengths

ggplot(sorted_contig_data , aes(x = reorder(Fasta_File_Name, -Contig_length), y = Contig_length, fill = Context)) + 
  geom_bar(stat = "identity") +
  #geom_vline(xintercept = nrow(sorted_contig_data) - 90, linetype = "dashed", color = "blue", size = 0.5) +
  labs(title = "Distribution of Contig Lengths (KChau, N=16)",
       y = "Contig Length",
       x = "Sample") +
  theme(axis.text.x  = element_blank(),  # Remove x-axis labels
        axis.ticks.x = element_blank())

### REHAB - with context

contig_data <- read.table("context_contigLength_n23_REHAB_ctxm15positive.tsv", header = TRUE, sep = "\t")

sorted_contig_data <- contig_data[order(contig_data$Contig_length), ]
sorted_contig_data$Context[sorted_contig_data$Context == "chrom"] <- "chromosome"

# Plotting a bar plot of contig lengths

ggplot(sorted_contig_data , aes(x = reorder(Fasta_File_Name, -Contig_length), y = Contig_length, fill = Context)) + 
  geom_bar(stat = "identity") +
  #geom_vline(xintercept = nrow(sorted_contig_data) - 90, linetype = "dashed", color = "blue", size = 0.5) +
  labs(title = "Distribution of Contig Lengths (REHAB, N=23)",
       y = "Contig Length",
       x = "Sample") +
  theme(axis.text.x  = element_blank(),  # Remove x-axis labels
        axis.ticks.x = element_blank())

##Blackwell
contig_data <- read.table("blackwell_fastaName_contigLength.tsv", header = TRUE, sep = "\t")
context <- read.table("contigName_context_blackwell_n5005plus3.tsv", header = TRUE, sep = "\t")
merged_data <- dplyr::full_join(contig_data, context, by = "Fasta_File_Name")
View(merged_data)
df <- distinct(merged_data[c(1,2,3)])
#View(df)
write.table(df, file = "context_contigLength_n5005plus3_Blackwell_ctxm15positive.tsv", row.names = FALSE)

sorted_contig_data <- merged_data[order(contig_data$Contig_length), ]
sorted_contig_data$Context[sorted_contig_data$Context == "chrom"] <- "chromosome"

# Plotting a bar plot of contig lengths

ggplot(sorted_contig_data , aes(x = reorder(Fasta_File_Name, -Contig_length), y = Contig_length, fill = Context)) + 
  geom_bar(stat = "identity") +
  #geom_vline(xintercept = nrow(sorted_contig_data) - 90, linetype = "dashed", color = "blue", size = 0.5) +
  labs(title = "Distribution of Contig Lengths (Blackwell, N=5005)",
       y = "Contig Length",
       x = "Sample") +
  theme(axis.text.x  = element_blank(),  # Remove x-axis labels
        axis.ticks.x = element_blank())

###Concatenate all datasets
entero <- read.table("context_contigLength_enterobase_n404_ctxm15positive.tsv", header = TRUE, sep = "\t")
kchau <- read.table("context_contigLength_n16_Kev_ctxm15positive.tsv", header = TRUE, sep = "\t")
rehab <- read.table("context_contigLength_n23_REHAB_ctxm15positive.tsv", header = TRUE, sep = "\t")
unpub <- read.table("context_contigLength_n117plus1_Unpub_ctxm15positive.tsv", header = TRUE, sep = "\t")
blackwell <- read.table("context_contigLength_n5005plus3_Blackwell_ctxm15positive.tsv", header = TRUE, sep = "\t")

# Merge tables
concat <- do.call(rbind, list(entero, kchau, rehab, unpub, blackwell))
concat$Context[concat$Context == "chrom"] <- "chromosome"
#View(concat)
uniq <- distinct(concat)
#View(uniq)
write.table(uniq, file ="N5565_fastaName_context_contigLength_withStudy.tsv", sep = "\t", row.names = FALSE)

# Summary stats

sumtable(uniq)

#Plot

ggplot(uniq , aes(x = reorder(Fasta_File_Name, -Contig_length), y = Contig_length, fill = Context)) + 
  geom_bar(stat = "identity") +
  #geom_vline(xintercept = nrow(sorted_contig_data) - 90, linetype = "dashed", color = "blue", size = 0.5) +
  scale_y_continuous(breaks = round(seq(min(0), max(uniq$Contig_length), by = 250000))) +
  labs(title = "Distribution of Contig Lengths (ALL, N=5565)",
       y = "Contig Length",
       x = "Sample") +
  theme(axis.text.x  = element_blank(),  # Remove x-axis labels
        axis.ticks.x = element_blank()) +
  facet_wrap(~Study, nrow = 2)

hist(uniq$Contig_length, col="red")

###
meta1 <- read.table("SUMMARY_HPRU_CTXM15postive_metadata_geoloc_year_host_source_mlst.tsv", header = TRUE, sep = "\t")
uniq1 <- read.table("N5565_fastaName_context_contigLength_withStudy.tsv",header = TRUE, sep = "\t")
summary <- dplyr::full_join(meta1, uniq1, by = "ID")
colnames(summary) <- c("ID","Year", "Country","Host","Specimen","MLST", "Context","Contig length (bp)","Study" )

###Further re-codes

#Country
summary$Country[summary$Country == "UK"] <- "United Kigndom"
summary$Country[summary$Country == "United Kigndom"] <- "United Kingdom"

#Year

#Country
summary$Year[summary$Year == "2012/2014"] <- "2012"

#MLST

summary$MLST <- ifelse(summary$MLST != "131" & summary$MLST != "410" & summary$MLST != "648" & summary$MLST != "405" & summary$MLST != "10" &
                         summary$MLST != "38" & summary$MLST != "156" & summary$MLST != "617" & summary$MLST != "69" & summary$MLST != "448" & 
                         summary$MLST != "1193" & summary$MLST != "58" & summary$MLST != "44", "Other", summary$MLST)


update_summary <- function(summary_df) {
  # Update Specimen
  specimen_lookup <- list(
    "unknown" = "Unknown",
    "canine" = "Unknown",
    "bovine intestine" = "Intra-abdominal",
    "bacter" = "Unknown",
    "animal feces" = "Faeces",
    "abscess" = "Skin/wound",
    "pleural cavity" = "Respiratory",
    "JP drainage" = "Unknown",
    "gall bladder" = "Intra-abdominal",
    "calf" = "Unknown",
    "respiratory tract" = "Respiratory",
    "pus" = "Skin/wound",
    "bile" = "Intra-abdominal",
    "abdominal" = "fluid Intra-abdominal",
    "sputum" = "Respiratory",
    "raw milk cheese" = "Food",
    "clinical isolate" = "Unknown",
    "animal" = "Unknown",
    "de-identified patient sample" = "Unknown",
    "wastewater" = "Wastewater",
    "stool" = "Faeces",
    "food" = "Food",
    "rectal swab" = "Faeces",
    "human" = "Unknown",
    "urine" = "Urine",
    "blood" = "Blood",
    "whale" = "Unknown",
    "cattle" = "Unknown",
    "cow" = "Unknown",
    "feces" = "Faeces",
    "clinical setting" = "Unknown",
    "dog lung" = "Respiratory",
    "drainage" = "Wastewater",
    "fluid" = "Unknown",
    "induced sputum" = "Respiratory",
    "poultry" = "Unknown",
    "wound abscess" = "Skin/wound",
    "wound secretion" = "Skin/wound",
    "skin" = "Skin/wound",
    "Skin" = "Skin/wound",
    "surgical wound" = "Skin/wound",
    "ulcer" = "Skin/wound",
    "tracheal aspirate" = "Respiratory",
    "abdominal fluid" = "Intra-abdominal",
    "lung" = "Lung",
    "liver" = "Liver",
    "kidney" = "Kidney",
    "intestine" = "Intestine",
    "ear swab" = "Ear swab"
  )
  
  for (old_value in names(specimen_lookup)) {
    summary_df$Specimen[summary_df$Specimen == old_value] <- specimen_lookup[[old_value]]
  }
  
  # Update Host
  host_lookup <- list(
    "HUMAN/FOOD" = "Human/Food",
    "Sus" = "Sus spp.",
    "Sus scrofa" = "Sus spp.",
    "Sus scrofa domesticus" = "Sus spp.",
    "Animalia" = "Animal - unknown",
    "Anas platyrhynchos" = "Duck",
    "Corvus frugilegus" = "Rook",
    "Equus caballus" = "Horse",
    "Environmental" = "Environment - unknown",
    "Canis lupus familaris" = "Dog",
    "Chroicocephalus" = "Gull",
    "Homo sapiens" = "Human",
    "Bos taurus" = "Boar",
    "Sus spp." = "Pig"
  )
  
  for (old_value in names(host_lookup)) {
    summary_df$Host[summary_df$Host == old_value] <- host_lookup[[old_value]]
  }
  
  return(summary_df)
}
summary <- update_summary(summary)
#View(summary)
sumtable(summary)
write.table(summary, file="SUMMARY_HPRU_CTXM15postive_metadata_geoloc_year_host_source_mlst_context_contigLength_study1.tsv", sep = "\t", row.names = FALSE)

### By Context ###

#ALL
ggplot(summary, aes(y=fct_infreq(Host))) +
  geom_histogram(aes(fill = Context),stat = "count", binwidth = 500) +
  theme_light() +
  scale_fill_manual(values=cbPalette) +
  labs( title = "Per Host Distribution",
      x = "Count",
       y  = "Host") 


ggplot(summary, aes(y=fct_infreq(Country))) +
  geom_histogram(aes(fill = Context),stat = "count", binwidth = 500) +
  theme_light() +
  scale_fill_manual(values=cbPalette) +
  labs( title = "Per Country Distribution",
        x = "Count",
        y  = "Country") 


ggplot(summary, aes(y=fct_infreq(Specimen))) +
  geom_histogram(aes(fill = Context),stat = "count", binwidth = 500) +
  theme_light() +
  scale_fill_manual(values=cbPalette) +
  labs(title = "Per Specimen Distribution", x = "Count",
       y  = "Specimen")


ggplot(summary, aes(y=fct_infreq(MLST))) +
  geom_histogram(aes(fill = Context),stat = "count", binwidth = 500) +
  theme_light() +
  scale_fill_manual(values=cbPalette) +
  labs(title = "Per MLST Distribution", x = "Count",
       y  = "MLST")

# contig length distribution

ggplot(summary, aes(x = reorder(ID, -Contig_length), y = Contig_length, fill = Context)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values=cbPalette) +
  labs(x = "Sample",
       y  = "Contig Length") +
  ggtitle("Contig Length Distribution") +
  theme(plot.title = element_text(vjust = 1, size = 18, face = "bold"),
        axis.text.x  = element_blank(),  # Remove x-axis labels
        axis.ticks.x = element_blank()) 


t <- tibble(summary$ID,summary$Contig_length)
colnames(t) <- c("ID", "Contig_length")
t
sumtable(t)


g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}

mylegend<-g_legend(a)

grid.arrange(arrangeGrob(a + theme(legend.position="none"),
                         b + theme(legend.position="none"),
                         c + theme(legend.position="none"),
                         d + theme(legend.position="none"),
                         nrow=1),
             mylegend, nrow=1,heights=c(10))


#grid.arrange(a,b,c,d, heights=c(3,3,3,4))


#Per group

cl <- ggplot(summary, aes(x = reorder(ID, -Contig_length), y = Contig_length, fill = Context)) +
  geom_bar(stat = "identity") +
  facet_wrap(~Study, nrow = 1, scales = "free_x") +
  theme_light() +
  scale_fill_manual(values=cbPalette) +
  labs(x = "Sample",
       y  = "Contig Length") +
  ggtitle("Contig Length Distribution per Study") +
  theme(plot.title = element_text(vjust = 1, size = 18, face = "bold"),
        axis.text.x  = element_blank(),  # Remove x-axis labels
        axis.ticks.x = element_blank())

cl

p1 <- ggplot(summary, aes(y=fct_infreq(Year))) +
  geom_histogram(aes(fill = Context),stat = "count", binwidth = 500) +
  #facet_wrap(~Study, nrow = 1, scales = "free_x") +
  theme_light() +
  scale_fill_manual(values=cbPalette) +
  labs(x = "Count",
       y  = "Year") +
  ggtitle("A") +
  theme(plot.title = element_text(vjust = 1, size = 18, face = "bold"))
p1

p2 <- ggplot(summary, aes(y=fct_infreq(Country))) +
  geom_histogram(aes(fill = Context),stat = "count", binwidth = 500) +
  #facet_wrap(~Study, nrow = 1, scales = "free_x") +
  theme_light() +
  scale_fill_manual(values=cbPalette) +
  labs(x = "Count",
       y  = "Country") +
  ggtitle("B") +
  theme(plot.title = element_text(vjust = 1, size = 18, face = "bold"), legend.position="none") 
  
p2

p3 <- ggplot(summary, aes(y=fct_infreq(Host))) +
  geom_histogram(aes(fill = Context),stat = "count", binwidth = 500) +
  facet_wrap(~Study, nrow = 1, scales = "free_x") +
  theme_light() +
  scale_fill_manual(values=cbPalette) +
  labs(x = "Count",
      y  = "Host") +
  ggtitle("C") +
  theme(plot.title = element_text(vjust = 1, size = 18, face = "bold"), legend.position="none") 
p3

ggplot(summary, aes(y=fct_infreq(Host))) +
  geom_histogram(aes(fill = Context), stat = "count", binwidth = 500) +
  facet_wrap(~Study, nrow = 1, scales = "free_x") +
  theme_light() +
  scale_fill_manual(values=cbPalette) +
  labs(x = "Count",
       y  = "Host") +
  ggtitle("C") +
  theme(plot.title = element_text(vjust = 1, size = 18, face = "bold"), legend.position="none") 

## Unknown Host
unknown_host <- data.frame(subset(summary, summary$Host == "Unknown"))
colnames(unknown_host) <- c("ID", "Year", "Country", "Host", "Specimen", "MLST", "Context", "Contig_length", "Study")
#View(unknown_host)

ggplot(unknown_host, aes(y=fct_infreq(Specimen))) +
  geom_histogram(aes(fill = Context), stat = "count", binwidth = 500) +
  #facet_wrap(~Study, nrow = 1, scales = "free_x") +
  theme_light() +
  scale_fill_manual(values=cbPalette) +
  labs(x = "Count",
       y  = "Specimen") +
  ggtitle("Unknown Host") +
  theme(plot.title = element_text(vjust = 1, size = 18, face = "bold"), legend.position="none") 

ggplot(unknown_host, aes(y=fct_infreq(Year))) +
  geom_histogram(aes(fill = Context), stat = "count", binwidth = 500) +
  #facet_wrap(~Study, nrow = 1, scales = "free_x") +
  theme_light() +
  scale_fill_manual(values=cbPalette) +
  labs(x = "Count",
       y  = "Year") +
  ggtitle("Unknown Host") +
  theme(plot.title = element_text(vjust = 1, size = 18, face = "bold"), legend.position="none") 

ggplot(unknown_host, aes(y=fct_infreq(Country))) +
  geom_histogram(aes(fill = Context), stat = "count", binwidth = 500) +
  #facet_wrap(~Study, nrow = 1, scales = "free_x") +
  theme_light() +
  scale_fill_manual(values=cbPalette) +
  labs(x = "Count",
       y  = "Country") +
  ggtitle("Unknown Host") +
  theme(plot.title = element_text(vjust = 1, size = 18, face = "bold"), legend.position="none") 

ggplot(unknown_host, aes(y=fct_infreq(MLST))) +
  geom_histogram(aes(fill = Context), stat = "count", binwidth = 500) +
  #facet_wrap(~Study, nrow = 1, scales = "free_x") +
  theme_light() +
  scale_fill_manual(values=cbPalette) +
  labs(x = "Count",
       y  = "MLST") +
  ggtitle("Unknown Host") +
  theme(plot.title = element_text(vjust = 1, size = 18, face = "bold"), legend.position="none") 

##

p4 <- ggplot(summary, aes(y=fct_infreq(Specimen))) +
  geom_histogram(aes(fill = Context),stat = "count", binwidth = 500) +
  facet_wrap(~Study, nrow = 1, scales = "free_x") +
  theme_light() +
  scale_fill_manual(values=cbPalette) +
  labs(x = "Count",
       y  = "Specimen") +
  ggtitle("D") +
  theme(plot.title = element_text(vjust = 1, size = 18, face = "bold"), legend.position="none") 
p4

## Unknown Specimen
unknown_specimen <- data.frame(subset(summary, summary$Specimen == "Unknown"))
colnames(unknown_specimen) <- c("ID", "Year", "Country", "Host", "Specimen", "MLST", "Context", "Contig_length", "Study")
#View(unknown_specimen)

ggplot(unknown_specimen, aes(y=fct_infreq(Host))) +
  geom_histogram(aes(fill = Context),stat = "count", binwidth = 500) +
  #facet_wrap(~Study, nrow = 1, scales = "free_x") +
  theme_light() +
  scale_fill_manual(values=cbPalette) +
  labs(x = "Count",
       y  = "Host") +
  ggtitle("Unknown Specimen") +
  theme(plot.title = element_text(vjust = 1, size = 18, face = "bold"), legend.position="none") 

ggplot(unknown_specimen, aes(y=fct_infreq(Year))) +
  geom_histogram(aes(fill = Context),stat = "count", binwidth = 500) +
  #facet_wrap(~Study, nrow = 1, scales = "free_x") +
  theme_light() +
  scale_fill_manual(values=cbPalette) +
  labs(x = "Count",
       y  = "Year") +
  ggtitle("Unknown Specimen") +
  theme(plot.title = element_text(vjust = 1, size = 18, face = "bold"), legend.position="none") 

ggplot(unknown_specimen, aes(y=fct_infreq(Country))) +
  geom_histogram(aes(fill = Context),stat = "count", binwidth = 500) +
  #facet_wrap(~Study, nrow = 1, scales = "free_x") +
  theme_light() +
  scale_fill_manual(values=cbPalette) +
  labs(x = "Count",
       y  = "Country") +
  ggtitle("Unknown Specimen") +
  theme(plot.title = element_text(vjust = 1, size = 18, face = "bold"), legend.position="none") 

ggplot(unknown_specimen, aes(y=fct_infreq(MLST))) +
  geom_histogram(aes(fill = Context),stat = "count", binwidth = 500) +
  #facet_wrap(~Study, nrow = 1, scales = "free_x") +
  theme_light() +
  scale_fill_manual(values=cbPalette) +
  labs(x = "Count",
       y  = "MLST") +
  ggtitle("Unknown Specimen") +
  theme(plot.title = element_text(vjust = 1, size = 18, face = "bold"), legend.position="none") 

##

p5 <- ggplot(summary, aes(y=fct_infreq(MLST))) +
  geom_histogram(aes(fill = Context),stat = "count", binwidth = 500) +
  facet_wrap(~Study, nrow = 1, scales = "free_y") +
  theme_light() +
  scale_fill_manual(values=cbPalette) +
  labs(y = "MLST",
       x  = "Count") +
  ggtitle("E") +
  theme(plot.title = element_text(vjust = 1, size = 18, face = "bold"), legend.position="none") 
p5


## Function for plotting target ST

# subset data frame and plot histograms
plot_st <- function(summary, mlst) {
  cbPalette <- c("#999999", "#E69F00")
  # Subset the data frame based on MLST value
  st <- subset(summary, MLST == mlst)
  # Rename columns
  colnames(st) <- c("ID", "Year", "Country", "Host", "Specimen", "MLST", "Context", "Contig_length", "Study")
  
  # Plot histogram for Host
  plot_host <- ggplot(st, aes(y = fct_infreq(Host))) +
    geom_histogram(aes(fill = Context), stat = "count", binwidth = 500) +
    facet_wrap(~Study, nrow = 1, scales = "free_y") +
    theme_light() +
    scale_fill_manual(values = cbPalette) +
    labs(y = "Host", x = "Count") +
    ggtitle(paste("ST", mlst, "Host")) +
    theme(plot.title = element_text(vjust = 1, size = 18, face = "bold"))
  
  # Plot histogram for Specimen
  plot_specimen <- ggplot(st, aes(y = fct_infreq(Specimen))) +
    geom_histogram(aes(fill = Context), stat = "count", binwidth = 500) +
    facet_wrap(~Study, nrow = 1, scales = "free_y") +
    theme_light() +
    scale_fill_manual(values = cbPalette) +
    labs(y = "Specimen", x = "Count") +
    ggtitle(paste("ST", mlst, "Specimen")) +
    theme(plot.title = element_text(vjust = 1, size = 18, face = "bold"))
  
  # Plot histogram for Country
  
  plot_country <- ggplot(st, aes(y = fct_infreq(Country))) +
    geom_histogram(aes(fill = Context), stat = "count", binwidth = 500) +
    facet_wrap(~Study, nrow = 1, scales = "free_y") +
    theme_light() +
    scale_fill_manual(values = cbPalette) +
    labs(y = "Country", x = "Count") +
    ggtitle(paste("ST", mlst, "Country")) +
    theme(plot.title = element_text(vjust = 1, size = 18, face = "bold"))
  
  # Plot histogram for Year
  
  plot_year <- ggplot(st, aes(y = fct_infreq(Year))) +
    geom_histogram(aes(fill = Context), stat = "count", binwidth = 500) +
    facet_wrap(~Study, nrow = 1, scales = "free_y") +
    theme_light() +
    scale_fill_manual(values = cbPalette) +
    labs(y = "Year", x = "Count") +
    ggtitle(paste("ST", mlst, "Year")) +
    theme(plot.title = element_text(vjust = 1, size = 18, face = "bold"))
  
  list(host_plot = plot_host, specimen_plot = plot_specimen, country_plot = plot_country, year_plot = plot_year)
}

# Specify ST
plots <- plot_st(summary, "405")
plots$host_plot
plots$specimen_plot
plots$country_plot
plots$year_plot

## Context
plot_context <- function(summary, context) {
  # Subset the data frame based on MLST value
  context <- subset(summary, Context == context)
  # Rename columns
  colnames(context) <- c("ID", "Year", "Country", "Host", "Specimen", "MLST", "Context", "Contig_length", "Study")
  
  # Plot histogram for Host
  plot_host <- ggplot(context, aes(y = fct_infreq(Host))) +
    geom_histogram(stat = "count", binwidth = 500) +
    #facet_wrap(~Study, nrow = 1, scales = "free_y") +
    theme_light() +
    labs(y = "Host", x = "Count") +
    theme(plot.title = element_text(vjust = 1, size = 18, face = "bold"))
  
  # Plot histogram for Specimen
  plot_specimen <- ggplot(context, aes(y = fct_infreq(Specimen))) +
    geom_histogram(stat = "count", binwidth = 500) +
    #facet_wrap(~Study, nrow = 1, scales = "free_y") +
    theme_light() +
    labs(y = "Specimen", x = "Count") +
    theme(plot.title = element_text(vjust = 1, size = 18, face = "bold"))
  
  # Plot histogram for Country
  
  plot_country <- ggplot(context, aes(y = fct_infreq(Country))) +
    geom_histogram(stat = "count", binwidth = 500) +
    #facet_wrap(~Study, nrow = 1, scales = "free_y") +
    theme_light() +
    labs(y = "Country", x = "Count") +
    theme(plot.title = element_text(vjust = 1, size = 18, face = "bold"))
  
  # Plot histogram for Year
  
  plot_year <- ggplot(context, aes(x = Year)) +
    geom_histogram(stat = "count", binwidth = 500) +
    #facet_wrap(~Study, nrow = 1, scales = "free_y") +
    theme_light() +
    labs(x = "Year", y = "Count") +
    theme(plot.title = element_text(vjust = 1, size = 18, face = "bold"))
  
  list(host_plot = plot_host, specimen_plot = plot_specimen, country_plot = plot_country, year_plot = plot_year)
}

#chromosome
chr <- plot_context(summary, "chromosome")
chr$host_plot
chr$specimen_plot
chr$country_plot
chr$year_plot

#plasmid
plm <- plot_context(summary, "plasmid")
plm$host_plot
plm$specimen_plot
plm$country_plot
plm$year_plot

#mylegend<-g_legend(p1)

#figure <- ggarrange(p1,p2,p3,p4,p5,
#          ncol = 1, nrow = 1,
#          common.legend = TRUE, legend = "bottom"
#          )
#figure

library(ggplot2)
library(gggenes)

# Create a dataframe with gene information
gene_data <- read.table("output.tsv", sep = "\t", header = TRUE)
genes <- data.frame(
  gene = gene_data$gene_label,
  start = gene_data$gene_start,
  end = gene_data$gene_end,
  strand = gene_data$gene_strand
)

genes$gene[genes$gene == "Beta-lactamase CTX-M-1"] <- "Beta-lactamase CTX-M-15"
genes$gene[genes$gene == "Periplasmic murein peptide-binding protein"] <- "Periplasmic murein peptide-binding protein mppA"
# Create the plot
ggplot(genes, aes(xmin = start, xmax = end, y = gene, strand = strand)) +
  geom_gene_arrow() +
  theme_genes()

# Count the frequency of each gene
gene_counts <- table(genes$gene)

# Convert the table to a data frame
gene_counts_df <- data.frame(gene = names(gene_counts), count = as.vector(gene_counts))

# Plot the gene frequencies
ggplot(gene_counts_df, aes(x = reorder(gene, count), y = count)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  labs(title = "Gene Frequencies", x = "Gene", y = "Count") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 60, hjust = 1))

## Draw genes using gggenes

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
# Aggregate the gene data by cluster name
aggregated_gene_data <- aggregate(cbind(gene_start, gene_end) ~ cluster_name + gene_label, data = gene_data1, FUN = function(x) c(min(x), max(x)))

# Plot the collapsed genes for each cluster name
gene_plot <- ggplot(aggregated_gene_data, aes(xmin = gene_start[,1], xmax = gene_end[,1], y = cluster_name, fill = gene_label)) +
  geom_gene_arrow(arrowhead_width = unit(0.1, "inches")) +
  scale_y_discrete(expand = c(0.02, 0.02)) +
  scale_fill_brewer(palette="RdYlBu") +
  theme_genes() +
  labs(y = "Cluster Name", fill = "Gene/Gene product")

# Display the plot
gene_plot



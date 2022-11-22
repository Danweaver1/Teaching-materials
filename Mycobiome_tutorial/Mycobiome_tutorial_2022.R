### Setup ###
#load packages
library(devtools)
library("stringr")
library("phyloseq"); packageVersion("phyloseq")
library("ggplot2"); packageVersion("ggplot2")
library("Biostrings", quietly = T)
library("tidyr")

#set plot theme
theme_set(theme_light())

##set working directory (the folder we want to work in)
setwd("C:/Users/Public/Mycobiome_tutorial/")
### Read in data ####
## metadata - this holds information about the samples
samples <- read.table(file="metadata.csv",
                      header=F,sep=",",row.names = 1, fill=TRUE)
#set the column names
colnames(samples) <- c("sample_ID","extraction_method","sample_type")

## ASV table - this is the counts for all the samples
otu_mat <- as.matrix(read.table("asv.csv", sep =",", header=TRUE, row.names = 1))

## Taxa table - this is all the fungal taxa found in the dataset
tax_mat <- as.matrix(read.table("taxa.csv", sep =",", header=TRUE, row.names = 1))

# Construct phyloseq object - this brings all the above information together: sample info, counts (otus) and taxa
physeq <- phyloseq(otu_table(otu_mat, taxa_are_rows=FALSE), 
                   sample_data(samples), 
                   tax_table(tax_mat))


#remove samples with no data 
ps.raw1 <- prune_samples(sample_sums(physeq) > 0, physeq)
#report empty samples
empty_samples <- which(sample_sums(physeq) == 0 )

if (length(empty_samples) > 0){
  cat("Following samples have no taxa identified:", sample_names(physeq)[empty_samples], "\n")
}

### Perform sample count filter - removes any low yield samples with < 5000 reads ####
ps.raw <- phyloseq::subset_samples(ps.raw1, phyloseq::sample_sums(ps.raw1) > 5000)

low_yield<- which(sample_sums(ps.raw1) < 5000 )
if (length(low_yield) > 0){
  cat("Following samples have fewer than 5K counts (raw data):", sample_names(ps.raw1)[low_yield], "\n")
}

#visualise read counts per sample - this plot will appear in the bottom right window 
plot_bar(ps.raw) + 
    geom_bar(stat="identity", position="stack") +
    facet_grid(extraction_method~sample_type,scales="free", space="free") +
    #theme(legend.position="none") +
    geom_hline(yintercept = 5000) +
    ylab("Abundance") + #set colours using colour table, set labels using formatted species names
    theme(panel.background = element_rect(fill = "white"), #change background to white
          axis.line = element_line(size = 1, colour = "black"),
          axis.ticks = element_line(size = 2),
          axis.text.x = element_text(color = "grey20", size = 10, angle =90, hjust = .5, vjust = .5, face = "plain"),
          axis.text.y = element_text(color = "grey20", size = 10, angle = 0, hjust = 1, vjust = .5, face = "plain"),
          axis.title.x = element_text(face = "bold", color = "grey20", size = 10, angle = 0, hjust = .5, vjust = 0),
          axis.title.y = element_text(face = "bold", color = "grey20", size = 10, angle = 90, hjust = .5, vjust = .5),
          plot.title = element_text(size = 29, face = "bold",hjust=0.5, vjust=-.5),
          legend.title=element_text(size=29),
          legend.text=element_text(size=27),
          legend.direction="vertical",
          legend.position ="bottom") +
    guides(fill = guide_legend(ncol = 4))

### Normalise counts - standardise abundances to the median sequencing depth ####
# As all sample read counts are variable, we are going to normalise them all 
#   to the median read counts
total = median(sample_sums(ps.raw))
standf = function(x, t=total) round(t * (x / sum(x)))
physeqn = transform_sample_counts(ps.raw, standf)

#### filter low level species - this is filter for species occurring above 0.2% in ANY sample ####
# this step removes species which are present at a very low abundance
physeq_abund <- filter_taxa(physeqn, function(x) sum(x > total*0.002) > 0, TRUE)

#remove any taxa that are no longer present after the filter 
taxa_to_keep <- names(which(taxa_sums(physeq_abund)>0))
physeq_abund <- prune_taxa(taxa_to_keep, physeq_abund)

#### Explore raw and filtered data ####
#send to file - this will create a text file with a summary of the data
sink("data_summary.txt")
##ASV level 
cat("Number of taxa in raw data: ", ntaxa(physeq), "\n")
cat("Number of taxa in filtered data: ", ntaxa(physeq_abund), "\n")

cat("Number of samples in raw data: ", nsamples(physeq), "\n")
cat("Number of samples in filtered data: ", nsamples(physeq_abund), "\n")

#ntaxa per sample
taxa_per_sample <- rowSums(otu_table(physeq_abund) > 0 )
cat("Median taxa per sample in filtered data: ",round(median(taxa_per_sample)), "\n")
cat("Range of number of taxa per sample in filtered data: ",range(taxa_per_sample), "\n")

sink()

#### Visualise fungal taxa identified in filtered dataset ####

#send the plot to file 
tiff("sample_heatmap.tiff",width = 900, height = 500)
print(
#plot species level taxa for each sample as heatmap, but group samples based on extraction method
phyloseq::plot_heatmap(physeq_abund, taxa.label = "Species", method = "NMDS",
                       distance = "bray", taxa.order = "Genus",  low="beige", high="red", na.value="white") +
  facet_grid(~extraction_method,scales="free", space="free") +
  theme(panel.background = element_rect(fill = "white"), #change background to white
        axis.line = element_line(size = 1, colour = "black"),
        axis.ticks = element_line(size = 2),
        axis.text.x = element_text(color = "grey20", size = 15, angle =90, hjust = .5, vjust = .5, face = "plain"),
        axis.text.y = element_text(color = "grey20", size = 15, angle = 0, hjust = 1, vjust = .5, face = "plain"),
        axis.title.x = element_text(face = "bold", color = "grey20", size = 15, angle = 0, hjust = .5, vjust = 0),
        axis.title.y = element_text(face = "bold", color = "grey20", size = 15, angle = 90, hjust = .5, vjust = .5),
        plot.title = element_text(size = 29, face = "bold",hjust=0.5, vjust=-.5),
        legend.title=element_text(size=27),
        legend.direction="vertical")
)
dev.off()
#you should now see that a plot file has appeared in the folder we are working in

#### Identify the top 10 most abundant species ####
#sort the total counts of all taxa and take the top 10
top10 <- names(sort(taxa_sums(physeq_abund),TRUE)[1:10])
#extract the sample data for the top 10 taxa found above
t10_data <- phyloseq::psmelt(physeq_abund) %>% dplyr::filter(OTU %in% top10 )

### plot top 10 species abundances - comparing the species abundances for each extraction method ###
tiff("top_10_species_abundances.tiff",width = 900, height = 700)
print(
  phyloseq::psmelt(physeq_abund) %>% #this reformats data fthe way ggplot likes it
    dplyr::filter(OTU %in% top10 ) %>% #filters data to top 10 taxa
    ggplot(data = ., aes(x = Species, y = Abundance, fill=extraction_method)) + #plots data using ggplot
    geom_boxplot(width=0.6) + #makes plot into a boxplot
    ggtitle(paste("Top 10 taxa abundance",sep=" ")) + #this and lines below are all changing how the plot looks/labels etc.
      ggtitle(paste("Top 10 taxa abundance",sep=" ")) +
      labs(y = "Abundance\n") +
      scale_fill_manual(values=c("#10639B","#808080", "#EFEF93")) +
      theme(panel.background = element_blank(),#change background to white
            panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
            axis.line = element_line(size = 1, colour = "black"),
            axis.ticks = element_line(size = 2),
            axis.text.x = element_text(color = "grey20", size = 20, angle =45, hjust = 1, vjust = 1, face = "plain"),
            axis.text.y = element_text(color = "grey20", size = 20, angle = 0, hjust = 1, vjust = .5, face = "plain"),
            axis.title.x = element_text(face = "bold", color = "grey20", size = 20, angle = 0, hjust = .5, vjust = 0),
            axis.title.y = element_text(face = "bold", color = "grey20", size = 20, angle = 90, hjust = .5, vjust = .5),
            plot.title = element_text(size = 15, face = "bold",hjust=0.5, vjust=-.5),
            legend.title=element_text(size=10),
            legend.text=element_text(size=10),
            legend.direction="vertical",
            legend.position ="bottom",
            plot.margin = ggplot2::margin(1,1,1,1, "cm")) 
)
dev.off()
  
### Diversity tests & stats  ####

#create richness estimate data to readout
Obs_abund <- estimate_richness(physeq_abund,measures="Observed")
Chao_abund <- estimate_richness(physeq_abund,measures="Chao1")
Shann_abund <- estimate_richness(physeq_abund,measures="Shannon")

#compare diversity between sample extraction methods
pairwise.wilcox.test(Obs_abund$Observed, sample_data(physeq_abund)$extraction_method, exact=F)
pairwise.wilcox.test(Shann_abund$Shannon, sample_data(physeq_abund)$extraction_method)
pairwise.wilcox.test(Chao_abund$Chao1, sample_data(physeq_abund)$extraction_method, exact=F)
#what are the p values for each comparison?
# are the diversity measures significantly different between any of the extraction methods?

#plot the diversity measures (to file)
tiff("richness_by_extraction_method.tiff",width = 800, height = 600)
print(plot_richness(physeq_abund,measures=c("Chao1", "Shannon","Observed"), x="extraction_method")+ geom_boxplot() +
        theme(panel.background = element_rect(fill = "white"), #change background to white
              axis.line = element_line(size = 1, colour = "black"),
              axis.ticks = element_line(size = 2),
              axis.text.x = element_text(color = "grey20", size = 30, angle =90, hjust = .5, vjust = .5, face = "plain"),
              axis.text.y = element_text(color = "grey20", size = 30, angle = 0, hjust = 1, vjust = .5, face = "plain"),
              axis.title.x = element_text(face = "bold", color = "grey20", size = 32, angle = 0, hjust = .5, vjust = 0),
              axis.title.y = element_text(face = "bold", color = "grey20", size = 32, angle = 90, hjust = .5, vjust = .5),
              plot.title = element_text(size = 29, face = "bold",hjust=0.5, vjust=-.5),
              legend.title=element_text(size=27),
              legend.text=element_text(size=25, face="italic"),
              legend.direction="vertical",
              legend.position="bottom",
              strip.text.x = element_text(size = 32, colour = "black"))
      
)
dev.off()
#which extraction method gave the highest diversity?

### Plot the taxa found at the genus level - as a barplot ####
###set colours 
genera_names <- tax_table(physeq_abund)[,6]
n <- length(genera_names)
##use Rcolorbrewer set and just sample from it
library(tidyr)
library(RColorBrewer)
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
genus_colour_table <- tibble(species=genera_names, Colour=sample(col_vector, n, replace=TRUE))


#plot taxa found (at genus level) as a stacked barplot
#plot to file
tiff("genera_stacked_barplot.tiff",width = 1000, height = 600)
print(
plot_bar(physeq_abund, fill = "Genus") + 
  geom_bar(aes(fill=Genus), stat="identity", position="stack") +
  facet_grid(sample_type~extraction_method,scales="free_x", space="free") +
  scale_fill_manual(values= genus_colour_table$Colour,breaks = genus_colour_table$species) +
  theme(panel.background = element_rect(fill = "white"), #change background to white
      axis.line = element_line(size = 1, colour = "black"),
      axis.ticks = element_line(size = 2),
      axis.text.x = element_text(color = "grey20", size = 15, angle =90, hjust = .5, vjust = .5, face = "plain"),
      axis.text.y = element_text(color = "grey20", size = 15, angle = 0, hjust = 1, vjust = .5, face = "plain"),
      axis.title.x = element_text(face = "bold", color = "grey20", size = 15, angle = 0, hjust = .5, vjust = 0),
      axis.title.y = element_text(face = "bold", color = "grey20", size = 15, angle = 90, hjust = .5, vjust = .5),
      legend.title=element_text(size=15),
      legend.text=element_text(size=15, face="italic"),
      legend.direction="vertical",
      legend.position="right",
      strip.text.x = element_text(size = 15, colour = "black"))
)
dev.off()
#which extraction methods had contamination? (ie. fungal taxa identified in the 
#   negative extraction controls?)

## Ordination ####
# calculate beta diversity - bray curtis
physeq.ord <- ordinate(physeq_abund, "NMDS", "bray")

#visualise the beta diversity clustering
#plot to file
tiff("sample_ordination.tiff",width = 600, height = 600)
print(
plot_ordination(physeq_abund, physeq.ord, type="samples", color="extraction_method", 
                shape="sample_type", title="Samples") + theme(aspect.ratio=1) + geom_point(size=3) +
  theme(panel.background = element_blank(), #change background to white
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.line = element_line(size = 1, colour = "black"),
        axis.ticks = element_line(size = 2),
        axis.text.x = element_text(color = "grey20", size = 30, angle =0, hjust = .5, vjust = .5, face = "plain"),
        axis.text.y = element_text(color = "grey20", size = 30, angle = 0, hjust = 1, vjust = .5, face = "plain"),
        axis.title.x = element_text(face = "bold", color = "grey20", size = 32, angle = 0, hjust = .5, vjust = 0),
        axis.title.y = element_text(face = "bold", color = "grey20", size = 32, angle = 90, hjust = .5, vjust = 0),
        legend.title=element_text(size=27),
        legend.text=element_text(size=25, face="italic"),
        legend.direction="vertical",
        legend.position = "bottom",
        strip.text.x = element_text(size = 32, colour = "black"))
)
dev.off()


#### Run permanova test ####
# this is testing if the there is sig. difference in beta diversity ###
#load library
library(vegan)
##remove NA samples for plotting
var_values <- sample_data(physeq_abund)[["extraction_method"]]
clean <- prune_samples(!is.na(var_values), physeq_abund)
#calculate distances
data <- phyloseq::distance(clean, method = "bray")

metadata <- as(sample_data(clean), "data.frame")
#run permanova test
res.permanova <- adonis2(data~ extraction_method,
                         data = metadata)
#report the results
if (res.permanova$`Pr(>F)`[1] < 0.05 ) {
  print(paste("Permanova pval is significant, pval = ",res.permanova$`Pr(>F)`[1], sep=""))
} else {
  print("no significance")
}
#################################################################################################


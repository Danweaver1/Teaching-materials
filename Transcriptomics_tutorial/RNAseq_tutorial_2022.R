#Code by Paul Bowyer & Danielle Weaver
### Initial setup ####
# Load DESeq2 library
library("DESeq2")
library("ggplot2")

theme_set(theme_bw())
# Set the working directory
#first define the directory
directory <- "C:/Users/Public/Transcriptomics_tutorial/"

#then tell R what directory to use
setwd(directory)

# Set the prefix for each output file name
outputPrefix <- "DESeq2"

#now tell R which sequence files you will use
sampleFiles<- c("Untreated_1.counts",
                "Untreated_2.counts",
                "Untreated_3.counts",
                "Treated_1.counts",
                "Treated_2.counts",
                "Treated_3.counts")

#and give the actual names
sampleNames <- c("Af 1","Af 2","Af 3","Af 4","Af 5","Af 6")

#tell R what the conditions are
sampleCondition <- c("no_drug","no_drug","no_drug","treated","treated","treated")

#and now get R to make the names and conditions into a table for DESeq2
sampleTable <- data.frame(sampleName = sampleNames, fileName = sampleFiles, condition = sampleCondition)

#And tell R which treatments you will be comparing
treatments = c("no_drug","treated")


### Run DESeq2 #######################################################################

# now define the experiment in DESeq2 - here we are adding the table, 
#   directory, the data from the alignment and the conditions we will test to ddsHTSeq - a dataset in R
ddsHTSeq <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable,
                                       directory = directory,
                                       design = ~ condition)

colData(ddsHTSeq)$condition <- factor(colData(ddsHTSeq)$condition,
                                      levels = treatments)

# run DESeq2 - this puts the results into a dataset (dds) and extracts more human readable results to res
dds <- DESeq(ddsHTSeq)
res <- results(dds)

#adapted from: https://benchtobioinformatics.wordpress.com/category/dexseq/


### Extract results from DESeq2 ##################################################

# filter data to include significant hits only
res = subset(res, padj<0.05)
# order results by padj value (most significant to least)
res <- res[order(res$padj),]

# view ordered results table
res
#should see DataFrame of baseMean, log2Foldchange, stat, pval, padj

#then send the ordered results (plus normalised count data) to a csv file that Excel can read
resdata <- merge(as.data.frame(res), as.data.frame(counts(dds,normalized =TRUE)), by = 'row.names', sort = FALSE)
names(resdata)[1] <- 'gene'

write.csv(resdata, file = paste0(outputPrefix, "_results_with_normalized.csv"), row.names = F)

#or you can send normalized counts to tab delimited file for GSEA, etc.
write.table(as.data.frame(counts(dds),normalized=T), file = paste0(outputPrefix, "_normalized_counts.txt"), sep = '\t', row.names = F)


### Exploratory data analysis of RNAseq data with DESeq2 #########################
#
# the next sections of code are for a variety of visualization, QC and other plots to
# get a sense of what the RNAseq data looks like based on DESeq2 analysis
#
# 1) MA plot
# 2) PCA plot
# 3) Volcano plot
# 4) Heatmap of most highly expressed genes
#

### MA plot ###############################################################################

# MA plot of RNAseq data for entire dataset
# http://en.wikipedia.org/wiki/MA_plot
# Log2 fold change (LFC) is plotted against mean normalised counts
# genes with padj < 0.05 are colored blue

png(paste0(outputPrefix, "-MAplot.png"), width = 500, height = 500, bg = "white")
plotMA(dds, ylim=c(-10,10),main = "RNAseq experiment", alpha=0.05)
dev.off()


### Heatmap of sample clustering #########################################
# excerpts from http://dwheelerau.com/2014/02/17/how-to-use-deseq2-to-analyse-rnaseq-data/

# perform variance stabilization (vsd) transformation to make it easier to visualise
vsd <- varianceStabilizingTransformation(dds, blind=T)

#load required packages for plotting
library("RColorBrewer")
library("gplots")

#create distance matrix
distsRL <- dist(t(assay(vsd)))
mat <- as.matrix(distsRL)

#add conditions to sample names to improve plot labels
rownames(mat) <- colnames(mat) <- with(colData(dds), paste(condition, sampleNames, sep=" : "))

#setup colour palette
hmcol <- colorRampPalette(brewer.pal(9, "GnBu"))(100)

#print plot to file 
png(paste0(outputPrefix, "-clustering.png"), width = 700, height = 600, bg = "white")
heatmap.2(mat, lhei = c(1,4), trace = "none", main="Sample to Sample Correlation (vst)", col = rev(hmcol), margin = c(13,13))
dev.off()


### PCA plot ################################################################################

#Principal components plot shows additional but rough clustering of samples

#print plot to file
png(paste0(outputPrefix, "-PCAplot.png"), width = 600, height = 500, bg = "white")
{
  #get PCA data
  pcaData <- plotPCA(vsd, intgroup=c("condition"), returnData=TRUE)
  #calculate percent variance per priniciple component
  percentVar <- round(100 * attr(pcaData, "percentVar"), digits=1)
  #make pretty PCA plot
  ggplot(pcaData, aes(PC1, PC2, color=condition)) +
    theme(panel.background = element_blank(),
        text=element_text(size=20)) +
    geom_point(size=5) +
    xlab(paste0("PC1: ",percentVar[1],"% variance")) +
    ylab(paste0("PC2: ",percentVar[2],"% variance")) 
}
dev.off()

### Volcano plot #############################################################################
#code from: https://support.bioconductor.org/p/78962/

tab = data.frame(logFC = res$log2FoldChange, negLogPval = -log10(res$pvalue))
head(tab)
png(paste0(outputPrefix, "-volcanoplot.png"), width = 550, height = 600, bg = "white")
{
  par(mar = c(5, 4, 4, 4))
  plot(tab, pch = 16, cex = 0.5, xlab = expression(log[2]~fold~change), ylab = expression(-log[10]~pvalue))
  lfc = 2
  pval = 0.05
  signGenes = (abs(tab$logFC) > lfc & tab$negLogPval > -log10(pval))
  points(tab[signGenes, ], pch = 16, cex = 0.5, col = "red") 
  abline(h = -log10(pval), col = "green3", lty = 2) 
  abline(v = c(-lfc, lfc), col = "blue", lty = 2) 
  mtext(paste("pval =", pval), side = 4, at = -log10(pval), cex = 0.8, line = 0.5, las = 1) 
  mtext(c(paste("-", lfc, "fold"), paste("+", lfc, "fold")), side = 3, at = c(-lfc, lfc), cex = 0.8, line = 0.5)
}
dev.off()


### Heatmap of most highly expressed genes ##################################################
library("RColorBrewer")
library("gplots")
# select the top 1000 expressed genes with heatmap.2
select <- order(rowMeans(counts(dds,normalized=T)),decreasing=T)[1:1000]
top_1000 <- assay(vsd)[select,]

#set colour palette
my_palette <- colorRampPalette(c("blue",'white','red'))(n=100)

png(paste0(outputPrefix, "-heatmap.png"), width = 600, height = 800, bg = "white")
heatmap.2(top_1000, col=my_palette, lhei = c(0.5,4),
          scale="row", key=T, keysize=1, symkey=T,
          density.info="none", trace="none",
          cexCol=1, labRow=F,
          main="Top 1000 Expressed Genes Heatmap")
dev.off()


### Final filtered significant hit list ###########################################################
#results have already been filtered to include only those with a significant adjusted p value 
# also want to filter using a LFC threshold of 2:

res.lfc = subset(res, abs(log2FoldChange)>2)
#then send the ordered results (plus normalised count data) to a csv file that Excel can read
resdata.lfc <- as.data.frame(res.lfc)

write.csv(resdata.lfc, file = paste0(outputPrefix, "_results_LFC_filtered.csv"), row.names = T)

# You now have your list of DEGs!
#   now go back to the tutorial sheet and follow the instructions 
#     to perform some pathway analysis on the DEGs using FungiDB
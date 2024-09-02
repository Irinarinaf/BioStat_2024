###################################################################################################
#
#     Prelude: You only need to run this code the first time
#
###################################################################################################

#source("https://bioconductor.org/biocLite.R")
#biocLite("tximport")
#biocLite("readr")
#biocLite("tximportData")
#biocLite("biomaRt")
#biocLite("GenomicFeatures")
#biocLite("DESeq2")
#install.packages("ggplot2")
#install.packages("gplots")
#install.packages("ggrepel")


###################################################################################################
#
#     Set up working environment and import data
#
###################################################################################################

## Load libraries (install if necessary)
library(tximport)

library(readr)
library(GenomicFeatures)
library(biomaRt)
library(DESeq2)
library(ggplot2)
library(gplots)
library(ggrepel)
library(dplyr)
library(Seurat)

## Change to location of results directory and samples.csv file
#setwd("/Users/Matthias/Dropbox/praveen/sequencing_projects/RNAseq_tutorial/")
base_dir = getwd()

# Create directory for output results
dir.create(file.path(getwd(), 'DESeq_output'), showWarnings = FALSE)

## Import sample and condition file
samples = read.csv(file.path(base_dir, "samples.csv"), header = TRUE, stringsAsFactors=FALSE)
samples$condition <- factor(samples$condition)
samples$condition <- relevel()
samples        # Prints the sample / condition list

## Make a TxDb object from transcript annotations
## - available as a GFF3 or GTF file
## - can be downloaded from gencode for mouse / human, ensembl for other species
## - this dataset was mapped against mm10

gtf="Homo_sapiens.GRCh38.101.gtf"   # Download this file from here; ftp://cbsuftp.tc.cornell.edu/pr46ftp/tutorial_files/DESeq/
# then move to the working directory (where the script is)
txdb=makeTxDbFromGFF(gtf,
                     format="gtf",
                     organism="Homo sapiens",
                     taxonomyId=9606)

k <- keys(txdb, keytype = "GENEID")
df <- select(txdb, keys = k, keytype = "GENEID", columns = "TXNAME")
tx2gene <- df[, 2:1]
head(tx2gene)

## Import Salmon quant files and create counts table
files <- file.path(base_dir, "salmon", samples$sample, "quant.sf")
all(file.exists(files))        # Verify names of files match names in samples.csv, should return True
names(files)=samples$sample
txi <- tximport(files, 
                type = "salmon", 
                tx2gene = tx2gene, ignoreTxVersion = TRUE)
head(txi$counts)               # This is the counts table for all of our samples

## Now to import the data into a DESeq Data Set (dds)
## Verify that sample names and colnames are the same
identical(samples$sample,colnames(txi$counts))

## Create a DEseqDataSet from txi count table
dds <- DESeqDataSetFromTximport(txi, samples, ~condition)

###################################################################################################
#
#    EXPLORATORY DATA ANALYSIS
#
###################################################################################################

library(tidyverse)
library("RColorBrewer")
library(pheatmap)

## Set color palette for figures
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)

## Perform a rlog transformation on count data (essentially a puts on a log2 scale)
## This helps our data assume a normal distribution and is good to do before these analyses
rld <- rlog(dds, blind=TRUE)

## Setup annotation file to show the conditions on the figures
treat_ann <- samples
rownames(treat_ann) <- treat_ann$sample
treat_ann$sample <- NULL
treat_ann

## SAMPLE TO SAMPLE DISTANCE & CORRELATION HEATMAPS

## Sample correlation heatmap
corr_samps <- cor(as.matrix(assay(rld)))      # Computes pairwise correlations between samples based on gene expression
png(filename="DESeq_output/DESeq_sampleCorr_HM.png", units = 'in', width = 8, height = 8, res = 250)
pheatmap(corr_samps,
         annotation = treat_ann,
         col=colors,
         main="Sample Correlations")
dev.off()

# Sample distance heatmap
sampleDists <- dist(t(assay(rld)))            # Computes Euclidean distance between samples based on gene expression
sampleDistMatrix <- as.matrix(sampleDists)

png(filename="DESeq_output/DESeq_sampleDist_HM.png", units = 'in', width = 8, height = 8, res = 250)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         annotation = treat_ann,
         col=colors,
         main="Sample to Sample Distances")
dev.off()

## Principal Component Analysis
## Separates samples based on variation between sample's gene expression
## Greater variation will affect separation to a greater degree

data <- plotPCA(rld, intgroup=c("condition"), returnData=TRUE)
percentVar <- round(100 * attr(data, "percentVar"))

png('DESeq_output/DESeq_PCA.png', units='in', width=8, height=6, res=250)
ggplot(data, aes(PC1, PC2, color=condition)) +
  geom_point(size=3.5) +
  geom_text_repel(aes(label=name)) +
  scale_colour_manual(values = c("orange", "steelblue", 'red')) +
  theme_bw() +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  ggtitle("")
dev.off()

###################################################################################################
#
#   Get table to convert names
#
###################################################################################################

## Convert the ensembl gene ID to gene name
## This will require an active internet connection
## Need to change according to species of interest

mart<- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))

ensembl_2_geneName <- getBM(
  attributes = c("ensembl_gene_id", "external_gene_name", "description"),
  mart = mart)
head(ensembl_2_geneName)
names(ensembl_2_geneName) <- c("GENEID", "geneName", "description")


###################################################################################################
#
#     DIFFERENTIAL EXPRESSSION ANALYSIS 
#
###################################################################################################
## DESeq = fx to calculate DE
## Combines multiple steps from DESeq
dds <- DESeq(dds)

## Write normalized counts to file
normalized.counts <- as.data.frame(counts(dds, normalized=TRUE ))
rownames(normalized.counts)<-sub("\\.[0-9]*", "", rownames(normalized.counts))
head(normalized.counts)

# Add gene name column
idx <- match( rownames(normalized.counts), ensembl_2_geneName$GENEID )
normalized.counts$geneName <- ensembl_2_geneName$geneName[ idx ]
head(normalized.counts)
write.table(normalized.counts, file = 'DESeq_output/DESeq_normalized_counts.csv', qmethod = NULL, sep = ',', quote = F)


###################################################################################################
#
#  Getting fold changes from direct comparisons with control
#
###################################################################################################

## Volcano plot function
## Will be used to visualize the differential expression of genes in the next section
volcanoPlot <- function(df, line_val, fc_cut = 0, pv_cut = .05, padj_cut = .1) {
  log2_lim <- 10
  pval_lim <- 20
  
  df_plt <- as.data.frame(df) %>%
    mutate(threshold = ifelse(pvalue <= pv_cut & log2FoldChange < fc_cut, -1,
                              ifelse(pvalue <= pv_cut & log2FoldChange > fc_cut, 1, 0))) %>%
    mutate(threshold = as.factor(threshold)) %>%
    mutate(log_pval = -log10(pvalue)) %>%
    mutate(shape = ifelse(log_pval > pval_lim | abs(log2FoldChange) > log2_lim, 17, 16)) %>%
    mutate(log_pval = ifelse(log_pval > pval_lim, pval_lim * .99, log_pval),
           log2FoldChange = ifelse(abs(log2FoldChange) > log2_lim, log2_lim * .99 * sign(log2FoldChange), log2FoldChange))
  
  
  
  ##Construct the plot object
  g = ggplot(df_plt, aes(log2FoldChange, y=log_pval, color=threshold)) +
    geom_point(alpha=0.4, size=1.75, shape = df_plt$shape) +
    scale_colour_manual(values = c("blue", "gray", "red")) +
    geom_hline(yintercept = -log10(pv_cut), linetype = 'dashed') +
    annotate('text', label = 'p-value', x = -log2_lim * .98, y = -log10(pv_cut) + .15, vjust = 0, hjust = 0) +
    ggtitle('') +
    scale_x_continuous("log2 fold change",
                       limit = c(-log2_lim, log2_lim),
                       expand = c(0,0)) +
    scale_y_continuous("-log10 p-value",
                       limit = c(0, pval_lim),
                       expand = c(0, 0)) +
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background = element_rect(colour = "black", fill = "white"), 
          text = element_text(size = 24),
          legend.position = "none")
  return(g)
}

###################################################################################################
#
#  Getting fold changes from direct comparisons with control
#
###################################################################################################

## Make comparison between conditions
samples
res <- results( dds, contrast = c("condition", "Flow", "Static") )
head(res)

## Add gene name and gene description to res file
rownames(res)<-sub("\\.[0-9]*", "", rownames(res))
idx <- match( rownames(res), ensembl_2_geneName$GENEID )
res$geneName <- ensembl_2_geneName$geneName[ idx ]
res$description <- ensembl_2_geneName$description[idx]

## Filter to remove genes with a baseMean of 5 or less
## baseMean is the average expression for that gene across all samples
res.5<-res[res$baseMean>5, ]

## Adjust p-value according to Benjamini & Hochberg method (need to do this since we filtered out genes by base mean 5 above)
res$padj <- p.adjust(res$pvalue, method="BH")

## Remove lines where pvalue is NA
res <- res[!is.na(res$pvalue),]
resdf = merge(as.data.frame(res), normalized.counts, by = 0 )
resdf = merge(resdf, markers, by.x='Row.names', by.y='V1', all.x = T)
resdf <- mutate(resdf, staticMean = rowMeans(select(resdf, starts_with("static")), na.rm = TRUE))
resdf = merge(resdf, ensembl_2_geneName, by.x = 'Row.names', by.y = 'GENEID', all.x = T)
resdf_tabula = merge(tabulagenes, resdf, by.x = 'row')

resdf<-resdf%>%mutate(ecde = case_when(
  EC == 1 & ((log2FoldChange > 0 & log2FoldChange < 1)|(padj > 0.05 & log2FoldChange > 0))  ~ 'Up',
  EC == 1 & padj <=0.05  & log2FoldChange > 1 ~ 'Up DE',
  EC == 1 & ((log2FoldChange < 0 & log2FoldChange > -1)|(padj > 0.05 & log2FoldChange < 0))  ~ 'Down',
  EC == 1 & padj <=0.05  & log2FoldChange < -1 ~ 'Down DE'
  
))

resdf <- mutate(resdf, flowMean = rowMeans(select(resdf, starts_with("flow")), na.rm = TRUE))

resdf<-resdf%>%mutate(DE = case_when(
  padj <=0.05  & log2FoldChange > 1 ~ 1,
  padj <=0.05  & log2FoldChange < -1 ~ 1
))

brainec = read.table('brainEC_TMP.txt', header = T)
resdf = merge(resdf, brainec, by.x = 'geneName', by.y = 'Bgene', all.x = T)
resdf <- mutate(resdf, BECMean = rowMeans(select(resdf, starts_with("BEC")), na.rm = TRUE))

tabula = readRDS('tabulaEC.rds')
kidney <- colnames(tabula)[grepl("*Kidney*", colnames(tabula))]
liver <- colnames(tabula)[grepl("*Liver*", colnames(tabula))]
salivary <- colnames(tabula)[grepl("*Salivary*", colnames(tabula))]
thymus <- colnames(tabula)[grepl("*Thymus*", colnames(tabula))]
spleen <- colnames(tabula)[grepl("*Spleen*", colnames(tabula))]
trachea <- colnames(tabula)[grepl("*Trachea*", colnames(tabula))]
tongue <- colnames(tabula)[grepl("*Tongue*", colnames(tabula))]
uterus <- colnames(tabula)[grepl("*Uterus*", colnames(tabula))]
mammary <- colnames(tabula)[grepl("*Mammary*", colnames(tabula))]
adipose <- c(colnames(tabula)[grepl("*FAT*", colnames(tabula))],
                 colnames(tabula)[grepl("*Fat*", colnames(tabula))])
eye <- colnames(tabula)[grepl("*Eye*", colnames(tabula))]
skin <- colnames(tabula)[grepl("*Skin*", colnames(tabula))]
heart <- colnames(tabula)[grepl("*Heart*", colnames(tabula))]
prostate <- colnames(tabula)[grepl("*Prostate*", colnames(tabula))]
pancreas <- colnames(tabula)[grepl("*ancreas*", colnames(tabula))]
muscle <- sample(colnames(tabula)[grepl("*uscle*", colnames(tabula))],2000)
lung <- c(colnames(tabula)[grepl("*Lung*", colnames(tabula))], 
              colnames(tabula)[grepl("*lung*", colnames(tabula))])
large_int <- colnames(tabula)[grepl("*_LI_*", colnames(tabula))]
small_int <- colnames(tabula)[grepl("*_SI_*", colnames(tabula))]
vasc <- colnames(tabula)[grepl("*Vasculature*", colnames(tabula))]
bladder <- colnames(tabula)[grepl("*ladder*", colnames(tabula))]
lymph <- colnames(tabula)[grepl("*LymphNodes*", colnames(tabula))]

mean_by_tissue <- function(tissue, tissue_name) {
  tabula_tissue  <- as.data.frame(as.matrix(tabula@assays$RNA@data[,tissue]))
  tabula_tissue <- mutate(tabula_tissue, tissue = rowMeans(tabula_tissue), na.rm = TRUE)
  tabula_means[,tissue_name] = tabula_tissue$tissue
  return(tabula_means)
}
tabula_means<-mean_by_tissue(eye,'eye')
tabula_means<-mean_by_tissue(kidney,'kidney')
tabula_means<-mean_by_tissue(liver,'liver')
tabula_means<-mean_by_tissue(skin,'skin')
tabula_means<-mean_by_tissue(muscle,'muscle')
tabula_means<-mean_by_tissue(bladder,'bladder')
tabula_means<-mean_by_tissue(vasc,'vasc')
tabula_means<-mean_by_tissue(large_int,'large_int')
tabula_means<-mean_by_tissue(small_int,'small_int')
tabula_means<-mean_by_tissue(lung,'lung')
tabula_means<-mean_by_tissue(heart,'heart')
tabula_means<-mean_by_tissue(prostate,'prostate')
tabula_means<-mean_by_tissue(pancreas,'pancreas')
tabula_means<-mean_by_tissue(lymph,'lymph')
tabula_means<-mean_by_tissue(adipose,'adipose')
tabula_means<-mean_by_tissue(spleen,'spleen')
tabula_means<-mean_by_tissue(mammary,'mammary')
tabula_means<-mean_by_tissue(thymus,'thymus')
tabula_means<-mean_by_tissue(tongue,'tongue')
tabula_means<-mean_by_tissue(salivary,'salivary')
tabula_means<-mean_by_tissue(uterus,'uterus')
tabula_means<-mean_by_tissue(trachea,'trachea')

## Check number of differentially expressed genes
as.data.frame(res) %>% filter(pvalue < .05) %>% dim()                       # Number of sig genes by p-value
as.data.frame(res) %>% filter(pvalue < .05 & log2FoldChange < 0) %>% dim()  # Number of sig down genes by p-value
as.data.frame(res) %>% filter(pvalue < .05 & log2FoldChange > 0) %>% dim()  # Number of sig up genes by p-value

## Write res.cont DESeq data to output file
write.csv(res, file="DESeq_output/DESeq_CONDAvsCONTROL.csv", quote=F, row.names = T)

## Construct gene list for miRhub
as.data.frame(res.5) %>%
  filter(pvalue < .05) %>%                  # p-value must be less than .05
  filter(log2FoldChange > 0) %>%            # log2FoldChange must be positive (aka upregulated genes)
  mutate(Name = toupper(geneName)) %>%      # Capitalize gene names
  pull(Name) %>%                            # Take only gene names column
  as.vector() %>%                           # Set as vector
  write.table(., file = 'DESeq_output/DEG_CONDAvsCONTROL_up_PVAL.05.txt', quote = F, row.name = F, col.names = F)     # Write to file

## Volcano plot
png('DESeq_output/DEG_CONDAvsCONTROL_VolcanoPlot.png', units = 'in', width = 6, height = 6, res = 250)
g <- volcanoPlot(res.5)
g
dev.off()

rownames(tabula_res) <- tabula_res$GeneID
cor_names <- c(colnames(tabula_means), c('staticMean', 'flowMean'))
cor_tabula <- tabula_res[,cor_names]
#cor_tabula <- tabula_res[which(tabula_res$DE == 1), cor_names]
corr_samps <- cor(as.matrix(cor_tabula), method = 'spearman', use = 'pairwise.complete.obs')  
corr_samps <- corr_samps[1:22,23:24]
pheatmap(corr_samps,
         col=colors,
         main="Tabula Sapiens Correlations")

tabula_res <- merge(tabula_res, brain, by.x = 'GeneID', by.y = 'V1', all.x = T)
cor_names <- c('brain', cor_names)
cor_tabula <- tabula_res[,cor_names]
corr_samps <- cor(as.matrix(cor_tabula), method = 'spearman', use = 'pairwise.complete.obs')
corr_samps <- corr_samps[1:23,24:25]
###########################

boxplot(tabula_res[which(tabula_res$EC == 1), 'liver'], tabula_res[which(is.na(tabula_res$EC)), 'liver'],
        tabula_res[which(tabula_res$EC == 1), 'kidney'], tabula_res[which(is.na(tabula_res$EC)), 'kidney']
        , tabula_res[which(tabula_res$EC == 1), 'brain'], tabula_res[which(is.na(tabula_res$EC)), 'brain'],
        na.rm = TRUE, outline = F, col=(c("white","grey", "white","grey","white","grey")),
        names = c('Liver EC', 'Liver', 'Kidney EC', 'Kidney', 'Brain EC', 'Brain'),
        ylab= 'Mean expression')

#Очитска от мертвых клеток
#надо заменить енсембл идентификаторы на имена генов, посчитать процент митохондриальных
#https://www.biostars.org/p/445418/
#https://satijalab.org/seurat/articles/pbmc3k_tutorial.html

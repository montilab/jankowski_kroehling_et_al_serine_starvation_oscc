library(DESeq2)
library(Biobase)
library(AUCell)
library(tximport)
library(hypeR)
setwd('/restricted/projectnb/montilab-p/projects/oralcancer/serine_starvation/')
setwd('/restricted/projectnb/montilab-p/CBMrepositoryData/otherStudies/RNAseq/2023_03_14_KukuruzinskaM_JankowskiS/scripts/fastq_human/')
stacy <- readRDS('./Output/Expression/JS_Gene_Expression.rds')

###PART 1: DE with DESEQ2 on HSC3 cell line

#read in files
files <- file.path(c("./JS_KM_1_2HC1/JS_KM_1_2HC1/RSEM/JS_KM_1_2HC1.genes.results", "./JS_KM_2_2HS1/JS_KM_2_2HS1/RSEM/JS_KM_2_2HS1.genes.results", 
                     "./JS_KM_3_3HC2/JS_KM_3_3HC2/RSEM/JS_KM_3_3HC2.genes.results", "./JS_KM_4_4HS2/JS_KM_4_4HS2/RSEM/JS_KM_4_4HS2.genes.results", 
                     "./JS_KM_5_5HC3/JS_KM_5_5HC3/RSEM/JS_KM_5_5HC3.genes.results", "./JS_KM_6_6HS3/JS_KM_6_6HS3/RSEM/JS_KM_6_6HS3.genes.results"))

names(files) <- c("JS_KM_1_2HC1", "JS_KM_2_2HS1", "JS_KM_3_3HC2", "JS_KM_4_4HS2", "JS_KM_5_5HC3", "JS_KM_6_6HS3")

txi.rsemhh <- tximport(files, type = "rsem", txIn = FALSE, txOut = FALSE)

geneWiseCountshh <- apply(txi.rsemhh$counts, 1, sum)
head(geneWiseCountshh)# What this looks like

# Find genes with at least 25 counts across all samples
CountGreater0hh <- geneWiseCountshh > 25
# How many genes will be removed
sum(!CountGreater0hh)

# Subset for CountGreater0 == TRUE
txi.rsemhh$counts <- txi.rsemhh$counts[ CountGreater0hh, ]
txi.rsemhh$length <- txi.rsemhh$length[ CountGreater0hh, ]
txi.rsemhh$abundance <- txi.rsemhh$abundance[ CountGreater0hh, ]

head(txi.rsemhh$counts)
colDatahh <- data.frame(condition=as.character(c("HC", "HS", "HC", "HS", "HC", "HS")))
ddshh <- DESeqDataSetFromTximport(txi.rsemhh, colDatahh, ~condition)
#dds$condition <- factor(dds$condition, levels=c("HC","HS"))

#run deseq2
#3 steps:
#1. estimate size factors
#estimate dispersion
#3. negative binomial GLM fitting and wald test
dds_reshh <- DESeq(ddshh)
dds_reshh <- estimateSizeFactors(dds_reshh)
reshh <- results(dds_reshh)
reshh$dispersion <- dispersions(dds_reshh)
#saveRDS(ddshh, '/restricted/projectnb/montilab-p/projects/oralcancer/serine_starvation/results/DESeq2/20230501_H_dds.rds')
#saveRDS(reshh, '/restricted/projectnb/montilab-p/projects/oralcancer/serine_starvation/results/DESeq2/20230501_H_DEout.rds')


reshh <- readRDS('/restricted/projectnb/montilab-p/projects/oralcancer/serine_starvation/results/DESeq2/20230501_H_DEout.rds')
plotMA(reshh, ylim=c(-7,7))
#shrunken log2fc
#plotMA(resLFC, ylim=c(-2,2))

resfilthh <- na.omit(reshh)

conversion <- as.data.frame(convert_ENSG_to_genename(rownames(resfilthh)))
rownames(conversion) <- conversion$Gene.stable.ID
resfilthh$Gene.stable.ID <- rownames(resfilthh)
table <- data.frame(resfilthh)
finaltable <- merge(x = table, y = conversion, by = 'Gene.stable.ID', all = T)
finaltable1 <- finaltable[finaltable$Gene.name!='',]
finaltable2 <- na.omit(finaltable1)

sigreshhall <- finaltable2[finaltable2$padj <=0.05,]
sigreshh <- finaltable2[(finaltable2$log2FoldChange >=1 | finaltable2$log2FoldChange<=-1) & finaltable2$padj <=0.05,]

#write.csv(finaltable2, '/restricted/projectnb/montilab-p/projects/oralcancer/serine_starvation/results/DESeq2/20230501_H_DEout.csv')
finaltable2 <- read_csv('/restricted/projectnb/montilab-p/projects/oralcancer/serine_starvation/results/DESeq2/20230501_H_DEout.csv')

convert_ENSG_to_genename <- function( humanex ) 
{
  ## input checks
  stopifnot( is.character(humanex) )
  stopifnot( require("biomaRt") )
  
  #human <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  #mouse <- useEnsembl("ensembl", dataset = "mmusculus_gene_ensembl")
  ## the following is a temporary workaround as the above doesn't currently work
  human <- useEnsembl(biomart="genes", dataset = "hsapiens_gene_ensembl", 
                      host = "https://dec2021.archive.ensembl.org/", mirror="useast")
  #mouse <- useEnsembl(biomart="genes", dataset = "mmusculus_gene_ensembl", 
  #                    host = "https://dec2021.archive.ensembl.org/", mirror="useast")
  genesV2 <- getLDS(attributes = c("ensembl_gene_id"), 
                    filters = "ensembl_gene_id", 
                    values = humanex, 
                    mart = human, 
                    attributesL = c("external_gene_name"), 
                    martL = human, 
                    uniqueRows=TRUE)
  humanx <- unique(genesV2[, 2])
  #return(humanx)
  return(genesV2)
}

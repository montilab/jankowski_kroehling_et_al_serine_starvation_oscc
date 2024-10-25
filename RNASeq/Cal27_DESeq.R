library(DESeq2)
library(Biobase)
library(AUCell)
library(tximport)
library(hypeR)
setwd('/restricted/projectnb/montilab-p/projects/oralcancer/serine_starvation/')
setwd('/restricted/projectnb/montilab-p/CBMrepositoryData/otherStudies/RNAseq/2023_03_14_KukuruzinskaM_JankowskiS/scripts/fastq_human/')
stacy <- readRDS('./Output/Expression/JS_Gene_Expression.rds')

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

#CS v CC
files <- file.path(c("./JS_KM_10_4CS2/JS_KM_10_4CS2/RSEM/JS_KM_10_4CS2.genes.results", "./JS_KM_11_5CC3/JS_KM_11_5CC3/RSEM/JS_KM_11_5CC3.genes.results", 
                     "./JS_KM_12_6CS3/JS_KM_12_6CS3/RSEM/JS_KM_12_6CS3.genes.results", "./JS_KM_7_1CC1/JS_KM_7_1CC1/RSEM/JS_KM_7_1CC1.genes.results",
                     "./JS_KM_8_2CS1/JS_KM_8_2CS1/RSEM/JS_KM_8_2CS1.genes.results", "./JS_KM_9_3CC2/JS_KM_9_3CC2/RSEM/JS_KM_9_3CC2.genes.results"))

names(files) <- c("JS_KM_10_4CS2", "JS_KM_11_5CC3", "JS_KM_12_6CS3",
                  "JS_KM_7_1CC1", "JS_KM_8_2CS1", "JS_KM_9_3CC2")

#names(files) <- paste0("sample", 1:6)
txi.rsemcc <- tximport(files, type = "rsem", txIn = FALSE, txOut = FALSE)

geneWiseCountscc <- apply(txi.rsemcc$counts, 1, sum)
head(geneWiseCountscc)# What this looks like

# Find genes with at least 1 counts across all samples
CountGreater0cc <- geneWiseCountscc > 25
# How many genes will be removed (Notice the !)
sum(!CountGreater0cc)

# Subset for CountGreater0 == TRUE
txi.rsemcc$counts <- txi.rsemcc$counts[ CountGreater0cc, ]
txi.rsemcc$length <- txi.rsemcc$length[ CountGreater0cc, ]
txi.rsemcc$abundance <- txi.rsemcc$abundance[ CountGreater0cc, ]


#continue on
head(txi.rsemcc$counts)
colDatacc <- data.frame(condition=as.character(c("CS", "CC", "CS", "CC", "CS", "CC")))
ddscc <- DESeqDataSetFromTximport(txi.rsemcc, colDatacc, ~condition)


#run deseq2
#3 steps:
#1. estimate size factors
#estimate dispersion
#3. negative binomial GLM fitting and wald test
dds_rescc <- DESeq(ddscc)
dds_rescc <- DESeq2::estimateSizeFactors(dds_rescc) 
rescc <- results(dds_rescc)
rescc$dispersion <- dispersions(dds_rescc)
#saveRDS(ddscc, '/restricted/projectnb/montilab-p/projects/oralcancer/serine_starvation/results/DESeq2/20230501_C_dds.rds')
#saveRDS(rescc, '/restricted/projectnb/montilab-p/projects/oralcancer/serine_starvation/results/DESeq2/20230501_C_DEout.rds')
#rescc <- readRDS('/restricted/projectnb/montilab-p/projects/oralcancer/serine_starvation/results/DESeq2/20230501_C_DEout.rds')


#make eset
eSet <- ExpressionSet(assayData=as.matrix(txi.rsemcc$counts))
colDatacc <- data.frame(condition=as.character(c("CS", "CC", "CS", "CC", "CS", "CC")), samples = c(colnames(eSet)))
pData(eSet) <- as.data.frame(colDatacc)
colnames(eSet) <- colnames(txi.rsemcc$counts)
sexp <- SummarizedExperiment::makeSummarizedExperimentFromExpressionSet(eSet[,order(sampleNames(eSet))])
names(assays(sexp)) <- "rawcounts"
assays(sexp)$log2norm <- log2(DESeq2::counts(dds_rescc, normalized=TRUE) + 1)
#saveRDS(sexp, '/restricted/projectnb/montilab-p/projects/oralcancer/serine_starvation/data/RNASeq/Cal27_sexp.rds')



res_summarycc <- data.frame(deseq2_padj=rescc$padj,
                            deseq2_logfc=rescc$log2FoldChange,
                            mean_exprs=rowMeans(log10(txi.rsemcc$counts+1)))

plotMA(rescc, ylim=c(-7,7))
#shrunken log2fc
#plotMA(resLFC, ylim=c(-2,2))

resfiltcc <- na.omit(rescc)

conversionc <- as.data.frame(convert_ENSG_to_genename(rownames(resfiltcc)))
rownames(conversionc) <- conversionc$Gene.stable.ID
resfiltcc$Gene.stable.ID <- rownames(resfiltcc)
tablec <- data.frame(resfiltcc)
finaltablecc <- merge(x = tablec, y = conversionc, by = 'Gene.stable.ID', all = T)
finaltable1cc <- finaltablecc[finaltablecc$Gene.name!='',]
finaltable2cc <- na.omit(finaltable1cc)

sigresccall <- finaltable2cc[finaltable2cc$padj <=0.05,]
sigrescc <- finaltable2cc[(finaltable2cc$log2FoldChange >=1 | finaltable2cc$log2FoldChange<=-1) & finaltable2cc$padj <=0.05,]

write.csv(finaltable2cc, '/restricted/projectnb/montilab-p/projects/oralcancer/serine_starvation/results/DESeq2/20230501_C_DEout.csv')

finaltable2cc <- read.csv('/restricted/projectnb/montilab-p/projects/oralcancer/serine_starvation/results/DESeq2/20230501_C_DEout.csv')
#upreg
uprescc <- finaltable2cc[finaltable2cc$log2FoldChange >=1 & finaltable2cc$padj <=0.05,]$Gene.name

#downreg
downrescc <- finaltable2cc[finaltable2cc$log2FoldChange <=-1 & finaltable2cc$padj <=0.05,]$Gene.name

signaturescc <- c(list(uprescc), list(downrescc))
names(signaturescc) <- c("C-SS", "C-C")

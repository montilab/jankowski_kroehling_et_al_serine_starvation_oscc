library(DiffBind)
library(parallel)
library(ggplot2)
library(GenomicRanges)
library(rtracklayer)
setwd('/restricted/projectnb/montilab-p/projects/oralcancer/serine_starvation/')

# Diffbind analysis to compare differential peak abundance between SS and C in H3K27me3.

setwd('/restricted/projectnb/montilab-p/projects/oralcancer/serine_starvation/results/cutnrun/diffbind/')

# We have already found the set of consensus peaks for SS and C separately, we want to use the union of these peaks 

# load consensus peaks
C_consensus <- import("/restricted/projectnb/montilab-p/projects/oralcancer/serine_starvation/results/cutnrun/peaksnoSMPR_broad/consensuspeaks/Cmplt_H3K27stitch_consensusPeaks.bed")   # GRanges
SS_consensus <- import("/restricted/projectnb/montilab-p/projects/oralcancer/serine_starvation/results/cutnrun/peaksnoSMPR_broad/consensuspeaks/SS_H3K27_consensusPeaks.bed")

union_peaks <- reduce(c(C_consensus, SS_consensus))          # union of both peak sets - reduce() function merges any any overlap (even 1bp) into one item
length(C_consensus)
length(SS_consensus)      # e.g. 20k + 18k = 38k
length(union_peaks)   
export(union_peaks, "H3K27me3_union_consensus.bed")

#read in sample sheet
samples <- read.csv('./samplesheet27_only1batch_20260209.csv')

db <- dba(sampleSheet = samples)
plot(db)

# count reads over !!Self-specified union peakset
# Even though you supplied union_peaks coordinates, DiffBind's dba.count(peaks = union_peaks) by default:
# Counts reads across your union regions
# Finds the highest-read summit (point of max pileup) within each union region
# Re-centers and shrinks the region around that summit to a fixed width (default 250bp)
# to Disable: db_union <- dba.count(db, peaks = union_peaks, summits = FALSE)
db_union <- dba.count(db, peaks = union_peaks, summits = FALSE) #change filter = 0 to not filter

#validate number of peaks in union peaks is the same that were counted for:
# Before counting
length(union_peaks)  # Your expected regions
# After counting  
dba.show(db_union)   # #peaks should match exactly

# use default normalization (don't need to define w separate step)
# Define contrast with batch blocking
db_union <- dba.contrast(db_union, categories = DBA_CONDITION, minMembers = 2) #same batch so don't need to include covariate
db_union <- dba.analyze(db_union)
dba.show(db_union, bContrasts=TRUE)
dba.plotMA(db_union)
dba.plotVolcano(db_union)
plot(db_union, contrast=1)
pdf('./202604_figures/H3K27_VlnPlot.pdf', height = 2.5, width = 4)
dba.plotVolcano(db_union, bFlip = T)
dev.off()

#get diff bound sites
db_union.DB.sig <- dba.report(db_union)
db_union.DB.sig[db_union.DB.sig$FDR < 0.05 & db_union.DB.sig$Fold > 0,] #up in C
db_union.DB.sig[db_union.DB.sig$FDR < 0.05 & db_union.DB.sig$Fold < (-2),]

db_union.DB.all <- dba.report(db_union, th=1, bCalled=TRUE)
head(db_union.DB.sig)
head(db_union.DB.all)
dba.plotVenn(db_union, contrast=1, bDB=TRUE,
             bGain=TRUE, bLoss=TRUE, bAll=FALSE)
dba.plotPCA(db_union,DBA_CONDITION,label=DBA_TREATMENT)
dba.plotPCA(db_union,DBA_TREATMENT,label=DBA_CONDITION)
dba.plotPCA(db_union, contrast=1, label=DBA_CONDITION)

dba.plotMA(db_union)
dba.plotVolcano(db_union)

sum(db_union.DB.all$Fold<0)
sum(db_union.DB.all$Fold>0)
sum(db_union.DB.sig$Fold<0)
sum(db_union.DB.sig$Fold>0)
pvals <- dba.plotBox(db_union)

# peak to gene overlap annotation
library(ChIPseeker)
library(org.Hs.eg.db)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)

# Annotation with promoter = -2kb:+1kb from TSS
peak_anno <- annotatePeak(db_union.DB.all, 
                          TxDb = TxDb.Hsapiens.UCSC.hg38.knownGene, #TxDb.Hs.eg.db
                          annoDb = "org.Hs.eg.db",
                          tssRegion = c(-2000, 100))

# Summary plot of genomic features (promoter/distal/intergenic)
plotAnnoPie(peak_anno)

anno_df <- as.data.frame(peak_anno@anno)
anno_df[anno_df$Fold > 0 & anno_df$FDR < 0.05,]
anno_df[anno_df$Fold < (-2) & anno_df$FDR < 0.05,]

# save
write.table(anno_df, "./202604_figures/H3K27_CUTRUN_differential_peaks_annotated_Cpos.txt", 
            sep = "\t", quote = FALSE, row.names = FALSE)


# run gsea (hypergeometric - KS is below)
sigs <- list("SS_up" = unique(anno_df[anno_df$Fold < 0 & anno_df$FDR < 0.05,]$SYMBOL),
             "C_up" = unique(anno_df[anno_df$Fold > 0 & anno_df$FDR < 0.05,]$SYMBOL))

lengths(sigs)

library(msigdbr)
library(hypeR)

background <- 20000
dogsea <- function(signatures = signatures) {
  enrichments1 <- list()
  genesets <- msigdb_gsets("Homo sapiens", "C2", "CP:KEGG_LEGACY", clean=TRUE)
  mhyp <- hypeR(signatures, genesets, test="hypergeometric", background = background) 
  #hyp_dots(mhyp, merge = T, top = 100, fdr = 0.05)
  enrichments1[[1]] <- mhyp
  names(enrichments1)[[1]] <- "kegg"
  genesets <- msigdb_gsets("Homo sapiens", "H", clean=TRUE)
  mhyp <- hypeR(signatures, genesets, test="hypergeometric", background = background) 
  #hyp_dots(mhyp, merge = T, top = 100, fdr = 0.05)
  enrichments1[[2]] <- mhyp
  names(enrichments1)[[2]] <- "hallmarks"
  REACTOME <- msigdb_gsets(species="Homo sapiens", "C2", "CP:REACTOME")
  mhyp <- hypeR(signatures, REACTOME, test="hypergeometric", background = background) 
  #hyp_dots(mhyp, merge=TRUE, fdr=0.01, title = "Reactome")
  enrichments1[[3]] <- mhyp
  names(enrichments1)[[3]] <- "reactome"
  biocarta <- msigdb_gsets(species="Homo sapiens", "C2", "CP:BIOCARTA")
  mhyp <- hypeR(signatures, biocarta, test="hypergeometric", background = background) 
  #hyp_dots(mhyp, merge=TRUE, fdr=0.05, title = "biocarta")
  enrichments1[[4]] <- mhyp
  names(enrichments1)[[4]] <- "biocarta"
  wiki <- msigdb_gsets(species="Homo sapiens", "C2", "CP:WIKIPATHWAYS")
  mhyp <- hypeR(signatures, wiki, test="hypergeometric", background = background) 
  #hyp_dots(mhyp, merge=TRUE, fdr=0.05, title = "wiki")
  enrichments1[[5]] <- mhyp
  names(enrichments1)[[5]] <- "wiki"
  bp <- msigdb_gsets(species="Homo sapiens", "C5", "BP")
  mhyp <- hypeR(signatures, bp, test="hypergeometric", background= background) 
  #hyp_dots(mhyp, merge=TRUE, fdr=0.05, abrv = 100, title = "bp")
  enrichments1[[6]] <- mhyp
  names(enrichments1)[[6]] <- "biologicalpathway"
  mf <- msigdb_gsets(species="Homo sapiens", "C5", "MF")
  mhyp <- hypeR(signatures, mf, test="hypergeometric", background= background)
  #hyp_dots(mhyp, merge=TRUE, fdr=0.05, title = "KEGG pathways enriched in each celltype")
  enrichments1[[7]] <- mhyp
  names(enrichments1)[[7]] <- "molecularfunction"
  cc <- msigdb_gsets(species="Homo sapiens", "C5", "CC")
  mhyp <- hypeR(signatures, cc, test="hypergeometric", background= background) 
  #hyp_dots(mhyp, merge=TRUE, fdr=0.05, title = " pathways enriched in each celltype")
  enrichments1[[8]] <- mhyp
  names(enrichments1)[[8]] <- "cellularcompartment"
  return(enrichments1)
}


gseas <- dogsea(signatures = sigs)

hyp_dots(gseas[[1]], fdr = 0.05, merge = T)
hyp_dots(gseas[[2]], fdr = 0.05, merge = T)
hyp_dots(gseas[[3]], fdr = 0.05, merge = T)
hyp_dots(gseas[[4]], fdr = 0.05, merge = T)
hyp_dots(gseas[[5]], fdr = 0.05, merge = T)
hyp_dots(gseas[[6]], fdr = 0.05, merge = T)
hyp_dots(gseas[[7]], fdr = 0.05, merge = T)
hyp_dots(gseas[[8]], fdr = 0.05, merge = T)

pdf('./202604_figures/H3K27_enrichment.pdf', height = 4, width = 6)
hyp_dots(gseas[[1]], fdr = 0.05, merge = T)
hyp_dots(gseas[[2]], fdr = 0.05, merge = T)
hyp_dots(gseas[[5]], fdr = 0.05, merge = T)
dev.off()

library(openxlsx)

# Create workbook
wb <- createWorkbook()
# Add each data frame as a separate sheet
addWorksheet(wb, "KEGG_SS")
addWorksheet(wb, "KEGG_C")
addWorksheet(wb, "Hallmarks_SS")
addWorksheet(wb, "Hallmarks_C")
addWorksheet(wb, "WP_SS")
addWorksheet(wb, "WP_C")
writeData(wb, "KEGG_SS", gseas[[1]]$data$SS_up$data)
writeData(wb, "KEGG_C", gseas[[1]]$data$C_up$data)
writeData(wb, "Hallmarks_SS", gseas[[2]]$data$SS_up$data)
writeData(wb, "Hallmarks_C", gseas[[2]]$data$C_up$data)
writeData(wb, "WP_SS", gseas[[5]]$data$SS_up$data)
writeData(wb, "WP_C", gseas[[5]]$data$C_up$data)
# Save
saveWorkbook(wb, "./202604_figures/H3K27_ORA.xlsx", overwrite=TRUE)


# gsea with only peaks in promoter regions (<=1kb)

sigs <- list("SS_up" = unique(anno_df[anno_df$Fold < 0 & anno_df$FDR < 0.05 & anno_df$annotation == 'Promoter (<=1kb)',]$SYMBOL),
             "C_up" = unique(anno_df[anno_df$Fold > 0 & anno_df$FDR < 0.05 & anno_df$annotation == 'Promoter (<=1kb)',]$SYMBOL))

lengths(sigs)

library(msigdbr)
library(hypeR)

background <- 20000
dogsea <- function(signatures = signatures) {
  enrichments1 <- list()
  genesets <- msigdb_gsets("Homo sapiens", "C2", "CP:KEGG_LEGACY", clean=TRUE)
  mhyp <- hypeR(signatures, genesets, test="hypergeometric", background = background) 
  #hyp_dots(mhyp, merge = T, top = 100, fdr = 0.05)
  enrichments1[[1]] <- mhyp
  names(enrichments1)[[1]] <- "kegg"
  genesets <- msigdb_gsets("Homo sapiens", "H", clean=TRUE)
  mhyp <- hypeR(signatures, genesets, test="hypergeometric", background = background) 
  #hyp_dots(mhyp, merge = T, top = 100, fdr = 0.05)
  enrichments1[[2]] <- mhyp
  names(enrichments1)[[2]] <- "hallmarks"
  REACTOME <- msigdb_gsets(species="Homo sapiens", "C2", "CP:REACTOME")
  mhyp <- hypeR(signatures, REACTOME, test="hypergeometric", background = background) 
  #hyp_dots(mhyp, merge=TRUE, fdr=0.01, title = "Reactome")
  enrichments1[[3]] <- mhyp
  names(enrichments1)[[3]] <- "reactome"
  biocarta <- msigdb_gsets(species="Homo sapiens", "C2", "CP:BIOCARTA")
  mhyp <- hypeR(signatures, biocarta, test="hypergeometric", background = background) 
  #hyp_dots(mhyp, merge=TRUE, fdr=0.05, title = "biocarta")
  enrichments1[[4]] <- mhyp
  names(enrichments1)[[4]] <- "biocarta"
  wiki <- msigdb_gsets(species="Homo sapiens", "C2", "CP:WIKIPATHWAYS")
  mhyp <- hypeR(signatures, wiki, test="hypergeometric", background = background) 
  #hyp_dots(mhyp, merge=TRUE, fdr=0.05, title = "wiki")
  enrichments1[[5]] <- mhyp
  names(enrichments1)[[5]] <- "wiki"
  bp <- msigdb_gsets(species="Homo sapiens", "C5", "BP")
  mhyp <- hypeR(signatures, bp, test="hypergeometric", background= background) 
  #hyp_dots(mhyp, merge=TRUE, fdr=0.05, abrv = 100, title = "bp")
  enrichments1[[6]] <- mhyp
  names(enrichments1)[[6]] <- "biologicalpathway"
  mf <- msigdb_gsets(species="Homo sapiens", "C5", "MF")
  mhyp <- hypeR(signatures, mf, test="hypergeometric", background= background)
  #hyp_dots(mhyp, merge=TRUE, fdr=0.05, title = "KEGG pathways enriched in each celltype")
  enrichments1[[7]] <- mhyp
  names(enrichments1)[[7]] <- "molecularfunction"
  cc <- msigdb_gsets(species="Homo sapiens", "C5", "CC")
  mhyp <- hypeR(signatures, cc, test="hypergeometric", background= background) 
  #hyp_dots(mhyp, merge=TRUE, fdr=0.05, title = " pathways enriched in each celltype")
  enrichments1[[8]] <- mhyp
  names(enrichments1)[[8]] <- "cellularcompartment"
  return(enrichments1)
}


gseas <- dogsea(signatures = sigs)

hyp_dots(gseas[[1]], fdr = 0.05, merge = T)
hyp_dots(gseas[[2]], fdr = 0.05, merge = T)
hyp_dots(gseas[[3]], fdr = 0.05, merge = T)
hyp_dots(gseas[[4]], fdr = 0.05, merge = T)
hyp_dots(gseas[[5]], fdr = 0.05, merge = T)
hyp_dots(gseas[[6]], fdr = 0.05, merge = T)
hyp_dots(gseas[[7]], fdr = 0.05, merge = T)
hyp_dots(gseas[[8]], fdr = 0.05, merge = T)

pdf('./202604_figures/H3K27_enrichment_PromoterOnly.pdf', height = 4, width = 6)
hyp_dots(gseas[[1]], fdr = 0.05, merge = T)
hyp_dots(gseas[[2]], fdr = 0.05, merge = T)
hyp_dots(gseas[[5]], fdr = 0.05, merge = T)
dev.off()

library(openxlsx)

# Create workbook
wb <- createWorkbook()
# Add each data frame as a separate sheet
addWorksheet(wb, "KEGG_SS")
addWorksheet(wb, "KEGG_C")
addWorksheet(wb, "Hallmarks_SS")
addWorksheet(wb, "Hallmarks_C")
addWorksheet(wb, "WP_SS")
addWorksheet(wb, "WP_C")
writeData(wb, "KEGG_SS", gseas[[1]]$data$SS_up$data)
writeData(wb, "KEGG_C", gseas[[1]]$data$C_up$data)
writeData(wb, "Hallmarks_SS", gseas[[2]]$data$SS_up$data)
writeData(wb, "Hallmarks_C", gseas[[2]]$data$C_up$data)
writeData(wb, "WP_SS", gseas[[5]]$data$SS_up$data)
writeData(wb, "WP_C", gseas[[5]]$data$C_up$data)
# Save
saveWorkbook(wb, "./202604_figures/H3K27_ORA_PromoterOnly.xlsx", overwrite=TRUE)


# KS gsea
anno_df1 <- anno_df[anno_df$annotation == 'Promoter (<=1kb)',]
# order peaks by fold change
# positive is up in SS, so peak anno order will start w C on top
peak_anno_order <- anno_df1[order(anno_df1$Fold, decreasing = TRUE),]
# or order with weighted pval (this doesn't change much)
# anno_df$weightedFold <- anno_df$Fold * (-log10(anno_df$p.value))
# peak_anno_order <- anno_df[order(anno_df$weightedFold, decreasing = TRUE),]


sigs <- list("C_up" = peak_anno_order$SYMBOL,
             "SS_up" = rev(peak_anno_order$SYMBOL))

lengths(sigs)

# pull all possible genes from genome ref
dogsea <- function(signatures = signatures) {
  enrichments1 <- list()
  genesets <- msigdb_gsets("Homo sapiens", "C2", "CP:KEGG_LEGACY", clean=TRUE)
  mhyp <- hypeR(signatures, genesets, test="kstest", fdr=0.05) 
  #hyp_dots(mhyp, merge = T, top = 100, fdr = 0.05)
  enrichments1[[1]] <- mhyp
  names(enrichments1)[[1]] <- "kegg"
  genesets <- msigdb_gsets("Homo sapiens", "H", clean=TRUE)
  mhyp <- hypeR(signatures, genesets, test="kstest", fdr=0.05) 
  #hyp_dots(mhyp, merge = T, top = 100, fdr = 0.05)
  enrichments1[[2]] <- mhyp
  names(enrichments1)[[2]] <- "hallmarks"
  REACTOME <- msigdb_gsets(species="Homo sapiens", "C2", "CP:REACTOME")
  mhyp <- hypeR(signatures, REACTOME, test="kstest", fdr=0.05) 
  #hyp_dots(mhyp, merge=TRUE, fdr=0.01, title = "Reactome")
  enrichments1[[3]] <- mhyp
  names(enrichments1)[[3]] <- "reactome"
  biocarta <- msigdb_gsets(species="Homo sapiens", "C2", "CP:BIOCARTA")
  mhyp <- hypeR(signatures, biocarta, test="kstest", fdr=0.05) 
  #hyp_dots(mhyp, merge=TRUE, fdr=0.05, title = "biocarta")
  enrichments1[[4]] <- mhyp
  names(enrichments1)[[4]] <- "biocarta"
  wiki <- msigdb_gsets(species="Homo sapiens", "C2", "CP:WIKIPATHWAYS")
  mhyp <- hypeR(signatures, wiki, test="kstest", fdr=0.05) 
  #hyp_dots(mhyp, merge=TRUE, fdr=0.05, title = "wiki")
  enrichments1[[5]] <- mhyp
  names(enrichments1)[[5]] <- "wiki"
  bp <- msigdb_gsets(species="Homo sapiens", "C5", "BP")
  mhyp <- hypeR(signatures, bp, test="kstest", fdr=0.05) 
  #hyp_dots(mhyp, merge=TRUE, fdr=0.05, abrv = 100, title = "bp")
  enrichments1[[6]] <- mhyp
  names(enrichments1)[[6]] <- "biologicalpathway"
  mf <- msigdb_gsets(species="Homo sapiens", "C5", "MF")
  mhyp <- hypeR(signatures, mf, test="kstest", fdr=0.05)
  #hyp_dots(mhyp, merge=TRUE, fdr=0.05, title = "KEGG pathways enriched in each celltype")
  enrichments1[[7]] <- mhyp
  names(enrichments1)[[7]] <- "molecularfunction"
  cc <- msigdb_gsets(species="Homo sapiens", "C5", "CC")
  mhyp <- hypeR(signatures, cc, test="kstest", fdr=0.05) 
  #hyp_dots(mhyp, merge=TRUE, fdr=0.05, title = " pathways enriched in each celltype")
  enrichments1[[8]] <- mhyp
  names(enrichments1)[[8]] <- "cellularcompartment"
  return(enrichments1)
}

gseas <- dogsea(signatures = sigs)

hyp_dots(gseas[[1]], fdr = 0.05, merge = T)
hyp_dots(gseas[[2]], fdr = 0.05, merge = T)
hyp_dots(gseas[[3]], fdr = 0.05, merge = T)
hyp_dots(gseas[[4]], fdr = 0.05, merge = T)
hyp_dots(gseas[[5]], fdr = 0.05, merge = T)
hyp_dots(gseas[[6]], fdr = 0.05, merge = T)
hyp_dots(gseas[[7]], fdr = 0.05, merge = T)
hyp_dots(gseas[[8]], fdr = 0.05, merge = T)


pdf('./202604_figures/H3K27_enrichment_KStest_promoter.pdf', height = 4, width = 6)
hyp_dots(gseas[[1]], fdr = 0.05, merge = T)
hyp_dots(gseas[[2]], fdr = 0.05, merge = T)
hyp_dots(gseas[[5]], fdr = 0.05, merge = T)
dev.off()


library(openxlsx)

# Create workbook
wb <- createWorkbook()
# Add each data frame as a separate sheet
addWorksheet(wb, "KEGG_SS")
addWorksheet(wb, "KEGG_C")
addWorksheet(wb, "Hallmarks_SS")
addWorksheet(wb, "Hallmarks_C")
addWorksheet(wb, "Wiki_SS")
addWorksheet(wb, "Wiki_C")
writeData(wb, "KEGG_SS", gseas[[1]]$data$SS_up$data)
writeData(wb, "KEGG_C", gseas[[1]]$data$C_up$data)
writeData(wb, "Hallmarks_SS", gseas[[2]]$data$SS_up$data)
writeData(wb, "Hallmarks_C", gseas[[2]]$data$C_up$data)
writeData(wb, "Wiki_SS", gseas[[5]]$data$SS_up$data)
writeData(wb, "Wiki_C", gseas[[5]]$data$C_up$data)
# Save
saveWorkbook(wb, "./202604_figures/H3K27_ORA_KS_promoter.xlsx", overwrite=TRUE)


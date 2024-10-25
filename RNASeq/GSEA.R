library(Biobase)
library(hypeR)

#read in DESeq2 outputs
hh <- read.csv('./results/DESeq2/20230501_H_DEout.csv')
cc <- read.csv('./results/DESeq2/20230501_C_DEout.csv')

#HSC3
finaltable2 <- read_csv('/restricted/projectnb/montilab-p/projects/oralcancer/serine_starvation/results/DESeq2/20230501_H_DEout.csv')

#get signatures
#upreg
upreghh <- finaltable2[finaltable2$log2FoldChange >=1 & finaltable2$padj <= 0.05,]$Gene.name
#downreg
downreghh <- finaltable2[finaltable2$log2FoldChange <=-1 & finaltable2$padj <= 0.05,]$Gene.name

#Cal27
##
finaltable2cc <- read.csv('/restricted/projectnb/montilab-p/projects/oralcancer/serine_starvation/results/DESeq2/20230501_C_DEout.csv')
#upreg
uprescc <- finaltable2cc[finaltable2cc$log2FoldChange >=1 & finaltable2cc$padj <=0.05,]$Gene.name

#downreg
downrescc <- finaltable2cc[finaltable2cc$log2FoldChange <=-1 & finaltable2cc$padj <=0.05,]$Gene.name

signaturescc <- c(list(uprescc), list(downrescc))
names(signaturescc) <- c("C-SS", "C-C")

#both
## 
signatures <- c(list(uprescc), list(downrescc), list(upreghh), list(downreghh))
names(signatures) <- c("C-SS", "C-C", "H-SS", "H-C")
#saveRDS(signatures, './results/20230501_signatureslfc1.rds')



pdf('./results/20240524_hallmarksenrichment.pdf', height = 5.5, width = 7)
fctable <- readRDS('./results/DESeq2/20230501_signatureslfc1.rds')
sigs <- c(list("C-C" = fctable$`C-C`),list("C-SS" = fctable$`C-SS`))
hallmarks <- msigdb_gsets("Homo sapiens", "H", clean=TRUE)
mhyph <- hypeR(sigs, hallmarks, test="hypergeometric", background = cc$Gene.name)
hyp_dots(mhyph, merge = T, title = 'hallmarks', top = 5, fdr=0.05) 
hyp_dots(mhyph, merge = T, title = 'hallmarks', top = 5, fdr=0.07) 
hyp_dots(mhyph, merge = T, title = 'hallmarks', top = 5, fdr=0.07) + ggstyle()

sigs <- c(list("H-C" = fctable$`H-C`),list("H-SS" = fctable$`H-SS`))
hallmarks <- msigdb_gsets("Homo sapiens", "H", clean=TRUE)
mhyph <- hypeR(sigs, hallmarks, test="hypergeometric", background = hh$Gene.name)
hyp_dots(mhyph, merge = T, title = 'hallmarks', top = 10, fdr=0.05) + ggstyle()
hyp_dots(mhyph, merge = T, title = 'hallmarks', top = 10, fdr=0.05) 
dev.off()

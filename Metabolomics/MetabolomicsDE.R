library(DESeq2)
library(Biobase)
library(AUCell)
library(tximport)
library(hypeR)
library(readr)
library(limma)
library(edgeR)
library(ggplot2)

setwd('/restricted/projectnb/montilab-p/projects/oralcancer/serine_starvation/')

dat <- read_delim("data/083122QsjankowskiMets1_6.txt", 
                  delim = "\t", escape_double = FALSE, 
                  col_types = cols(`083022QSJ1` = col_number(), 
                                   `083022QSJ2` = col_number(), `083022QSJ3` = col_number(), 
                                   `083022QSJ4` = col_number(), `083022QSJ5` = col_number(), 
                                   `083022QSJ6` = col_number()), trim_ws = TRUE)
dat1 <- dat[-c(1,2),]
names(dat1) <- c("metabolite", "s1", "s2", "s3", "s4", "s5", "s6")
dat2 <- data.frame(dat1)
rownames(dat2) <- dat1$metabolite
dat2 <- dat2[,-1]

omitnas <- na.omit(dat2)
#eliminates 65 metabolites

# Plots
par(mfrow=c(1,2))
#look as distribution of counts
hist(dat2$s1, breaks = 100000, xlim = c(0,100000), main = "s1 counts")
hist(log(dat2$s1,2), breaks = 10, main = "log2 s1 counts")
par(mfrow=c(1,2))
hist(dat2$s1, breaks = 100000, xlim = c(0,100000), main = "s2 counts")
hist(log(dat2$s1,2), breaks = 10, main = "log2 s2 counts")
par(mfrow=c(1,2))
hist(dat2$s1, breaks = 100000, xlim = c(0,100000), main = "s3 counts")
hist(log(dat2$s1,2), breaks = 10, main = "log2 s3 counts")
par(mfrow=c(1,2))
hist(dat2$s1, breaks = 100000, xlim = c(0,100000), main = "s4 counts")
hist(log(dat2$s1,2), breaks = 10, main = "log2 s4 counts")
par(mfrow=c(1,2))
hist(dat2$s1, breaks = 100000, xlim = c(0,100000), main = "s5 counts")
hist(log(dat2$s1,2), breaks = 10, main = "log2 s5 counts")
par(mfrow=c(1,2))
hist(dat2$s1, breaks = 100000, xlim = c(0,100000), main = "s6 counts")
hist(log(dat2$s1,2), breaks = 10, main = "log2 s6 counts")

#boxplots
labels <- c("s1_C", "s2_M", "s3_C", "s4_M", "s5_C", "s6_M")
boxplot(dat2$s1, dat2$s2, dat2$s3, dat2$s4, dat2$s5, dat2$s6, col="orange",
        border="brown", ylim = c(0,2000000), names = labels)
boxplot(log(dat2$s1,2), log(dat2$s2,2), log(dat2$s3,2), log(dat2$s4,2), log(dat2$s5,2), log(dat2$s6,2), col="orange",
        border="brown", ylim = c(0,30), names = labels)


colData <- data.frame(condition=as.character(c("COMPLETE", "MINSG", "COMPLETE", "MINSG", "COMPLETE", "MINSG")))
rownames(colData) <- colnames(omitnas)


#limma
## Organize data and design
cnt <- DGEList(omitnas)
cnt$samples$group <- group <- colData$condition

design <- model.matrix(~0+group)


## Normalization options
# 1) log2 and quantile normalization
data.norm <- normalizeQuantiles(log2(omitnas))

#make eset
#make eset
eSet <- ExpressionSet(assayData=as.matrix(omitnas))
colData <- data.frame(condition=as.character(c("COMPLETE", "MINSG", "COMPLETE", "MINSG", "COMPLETE", "MINSG")),
                      labels = c("s1_C", "s2_M", "s3_C", "s4_M", "s5_C", "s6_M"))
pData(eSet) <- as.data.frame(colData)
colnames(eSet) <- colnames(omitnas)
sexp <- SummarizedExperiment::makeSummarizedExperimentFromExpressionSet(eSet[,order(sampleNames(eSet))])
names(assays(sexp)) <- "rawcounts"
assays(sexp)$log2norm <- data.norm

#saveRDS(sexp, '/restricted/projectnb/montilab-p/projects/oralcancer/serine_starvation/data/metabolomics/Metab_sexp.rds')


boxplot(data.norm$s1, data.norm$s2, data.norm$s3, data.norm$s4, data.norm$s5, data.norm$s6, col="orange",
        border="brown", ylim = c(0,30), names = labels)


## Define contrasts and fit model
contr.matrix <- makeContrasts(CUP=groupCOMPLETE-groupMINSG, MUP=groupMINSG-groupCOMPLETE,  levels = colnames(design))

tfit <- lmFit(data.norm, design)
head(coef(tfit))
tfit <- contrasts.fit(tfit, contr.matrix)
efit <- eBayes(tfit, trend=TRUE, robust=TRUE)
top.table <- topTable(efit, sort.by = "F", n = Inf)
head(top.table, 20)

length(which(top.table$adj.P.Val < 0.05))

library(EnhancedVolcano)

keyvals <-  ifelse(
  top.table$MUP <= -1 & top.table$adj.P.Val < 0.05, 'blue', 
  ifelse(
    top.table$MUP >= 1 & top.table$adj.P.Val < 0.05, 'red3',
    'grey'))

names(keyvals)[keyvals == 'grey'] <- 'Not significant'
names(keyvals)[keyvals == 'blue'] <- 'Log2 FC <= -1 & padj < 0.05'
names(keyvals)[keyvals == 'red3'] <- 'Log2 FC >= 1 & padj < 0.05'

metablabels <- c("serine", "glycine","3-phospho-serine", "3-phosphoglycerate","purine","Cystine", "glutathione disulfide-nega")
pdf('./results/metabolomics/20241025volcano1_quantilenorm.pdf', height = 5, width = 10)
EnhancedVolcano(top.table,
                lab = rownames(top.table),
                x = 'MUP',
                y = 'adj.P.Val',
                pCutoff = 0.05,
                colCustom = keyvals,selectLab = metablabels,
                pointSize = 3.0,
                labSize = 5.0, 
                drawConnectors = TRUE,
                widthConnectors = 0.5, 
                max.overlaps = Inf) + ggstyle()
dev.off()

#bargraph
library(dplyr)
pdf('./metabolomics_aa_barplot.pdf', height = 5, width = 10)
top.table[aa,] %>%
  mutate(sig = case_when(adj.P.Val < 0.01 ~ "adj.p.val < 0.01",
                         adj.P.Val >= 0.01 ~ "adj.p.val >= 0.01")) %>%
  ggplot() + 
  geom_col(aes(x=reorder(metabolite, -adj.P.Val), y=MUP, fill=sig)) +
  scale_fill_manual(values=c("blue","grey")) + 
  xlab("Log2 fold change") +
  ylab("amino acid") +
  theme(legend.position="none") +
  ggstyle() + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1.05, hjust=1))
dev.off()

library(EnhancedVolcano)



keyvals <-  ifelse(
      cc$log2FoldChange <= -1 & cc$padj < 0.05, 'blue', 
      ifelse(
        cc$log2FoldChange >= 1 & cc$padj < 0.05, 'red3',
        'grey'))
        
names(keyvals)[keyvals == 'grey'] <- 'Not significant'
names(keyvals)[keyvals == 'red3'] <- 'Log2 FC => 1 & padj < 0.05'
names(keyvals)[keyvals == 'blue'] <- 'Log2 FC <= 1 & padj < 0.05'

clab <- c("SERPINE1",
          "TGM2",
          "IGFBP3",
          "LGALS1",
          "MKI67",
          "SPARC",
          "VIM",
          "MATN2",
          "CDH2",
          "PHGDH",
          "PSAT1",
          "PSPH",
          "CDH1",
          "KRT14",
          "KRT5",
          "ALDH7A1",
          "BMI1",
          "ALDH3A1",
          "CD24",
          "ASNS",
          "PPP1R15A",
          "TRIB3",
          "LCN2",
          "SLC7A11",
          "SQSTM1",
          "DDIT4",
          "GAS5",
          "KRT80")

pdf('./results/DESeq2/20240429figures/20240525_ccVolvano.pdf', width = 8, height = 8)
EnhancedVolcano(cc,
                lab = cc$Gene.name,
                x = 'log2FoldChange',
                y = 'padj',
                title = 'Cal27',
                pCutoffCol = 'padj',
                pCutoff = 0.05,
                FCcutoff = 1, 
                #xlim = c(-3,4),
                colCustom = keyvals,
                selectLab = clab,
                #col=c('grey', 'grey', 'grey', 'red3'),
                pointSize = 3.0,
                labSize = 4.0,
                drawConnectors = TRUE,
                widthConnectors = 0.5, 
                max.overlaps = Inf) #selectLab
dev.off()

keyvals <-  ifelse(
  hh$log2FoldChange <= -1 & hh$padj < 0.05, 'blue', 
  ifelse(
    hh$log2FoldChange >= 1 & hh$padj < 0.05, 'red3',
    'grey'))
    
names(keyvals)[keyvals == 'grey'] <- 'Not significant'
names(keyvals)[keyvals == 'red3'] <- 'Log2 FC => 1 & padj < 0.05'
names(keyvals)[keyvals == 'blue'] <- 'Log2 FC <= 1 & padj < 0.05'

hlab <- c("SERPINE1",
          "COL12A1",
          "TGM2",
          "RHOB",
          "IGFBP3",
          "LGALS1",
          "MKI67",
          "SPARC",
          "VIM",
          "MATN2",
          "CDH2",
          "PHGDH",
          "PSAT1",
          "PSPH",
          "CDH1",
          "KRT10",
          "KRT14",
          "KRT5",
          "ALDH7A1",
          "BMI1",
          "ALDH3A1",
          "CD24",
          "ASNS",
          "PPP1R15A",
          "TRIB3",
          "LCN2",
          "SLC7A11",
          "SQSTM1",
          "DDIT4",
          "GAS5",
          "KRT80")

pdf('./results/DESeq2/20240429figures/20240525_HscVolvano.pdf', width = 8, height = 8)
EnhancedVolcano(hh,
                lab = hh$Gene.name,
                x = 'log2FoldChange',
                y = 'padj',
                title = 'HSC-3',
                pCutoffCol = 'padj',
                pCutoff = 0.05,
                colCustom = keyvals,
                selectLab = hlab,
                FCcutoff = 1,
                pointSize = 3.0,
                labSize = 4.0,
                drawConnectors = TRUE,
                widthConnectors = 0.5, 
                max.overlaps = Inf) #selectLab
dev.off()

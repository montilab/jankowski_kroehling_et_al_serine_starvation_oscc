
### comparison to keratinization sig
ksig <- c("A2ML1", "A4GALT", "AADACL3", "ABCC5", "ABCG2", "ABHD2", "ABO", "ACER2", "ACKR2", "ACOT11", 
          "ACOXL", "ACP6", "ACP7", "ACSF2", "ADAM8", "ADAMTSL4", "ADGRF1", "ADGRF4", "ADH7", "ADRB2", 
          "AGTRAP", "AHNAK2", "AIG1", "AJUBA", "AKTIP", "ALDH1A3", "ALDH1L1", "ALOX12B", "ALOXE3", 
          "ALPP", "ALS2CL", "AMDHD1", "AMOTL1", "AMOTL2", "AMTN", "AMY1B", "AMY1C", "ANKK1", "ANKUB1", 
          "ANO1", "ANO6", "ANXA11", "ANXA3", "ANXA8", "ANXA8L1", "AOPEP", "AP1M2", "AP1S3", "APLP2", 
          "APOBEC3A", "AQP3", "AQP6", "AREG", "ARHGAP12", "ARHGAP23", "ARHGAP26", "ARHGAP27", "ARHGAP8", 
          "ARHGEF35", "ARHGEF37", "ARHGEF4", "ARL14.00", "ARRDC4", "ARSI", "ARTN", "ASAP2", "ASPG", 
          "ASPRV1", "ATG16L1", "ATOSB", "ATP13A4", "ATP1B3", "ATP2C2", "ATP6V1B1", "B3GALT4", "B3GALT5",
          "B4GALNT3", "B4GALT4", "BAG3", "BAK1", "BCL2L1", "BICD2", "BICDL2", "BLCAP", "BNC1", "BNIPL",
          "BPGM", "BSPRY", "BTBD16", "BTC", "BTG4", "C10orf67", "C12orf54", "C15orf48", "C15orf62", 
          "C16orf74", "C17orf50", "C19orf33", "C19orf84", "C1QTNF12", "C1QTNF6", "C1orf116", "C1orf68", 
          "C1orf74", "C22orf23", "C3", "C3orf52", "C4BPB", "CAB39", "CAGE1", "CALB2", "CALHM3", "CALML3",
          "CAMK2G", "CAPN1", "CAPN6", "CAPNS2", "CARD10", "CASP14", "CAST", "CAV2", "CBLC", "CBY2", "CCDC107",
          "CCDC120", "CCDC17", "CCDC68", "CCDC87", "CCDC9B", "CCL20", "CCL21", "CCL27", "CCN6", "CD200R1L", 
          "CD300E", "CD55", "CDA", "CDC42BPG", "CDC42EP3", "CDCP1", "CDH1", "CDH3", "CDHR2", "CDKN2B", "CDS1", 
          "CEACAM6", "CERS3", "CFAP251", "CFB", "CGB7", "CGN", "CHIA", "CHST4", "CLCA2", "CLCA4", "CLDN12",
          "CLDN17", "CLDN4", "CLDN7", "CLDN8", "CLIC3", "CLTA", "CLTB", "CNFN", "CNGA1", "CNKSR1", "COL17A1", 
          "COL4A5", "COL4A6", "CORO2A", "CPAMD8", "CPD", "CRB3", "CRCT1", "CREB3", "CRNN", "CRYBB3", "CRYBG2",
          "CSN3", "CST2", "CST5", "CST6", "CSTA", "CTAGE4", "CTDSPL", "CTNNBIP1", "CTNND1", "CTSH", "CXCL1", 
          "CXCL11", "CXCL14", "CXCL16", "CXCL17", "CXCL2", "CXCL3", "CXCL5", "CXCL6", "CXCL8", "CXorf58", 
          "CYB5R1", "CYP1A1", "CYP1A2", "CYP2F1", "CYP8B1", "CYSRT1", "DAB2IP", "DDR1", "DEFB103A", "DEFB103B", 
          "DEFB110", "DEFB112", "DEFB123", "DEFB4A", "DEFB4B", "DGKH", "DHDH", "DHRS1", "DHRS9", "DLK2", "DNAH2",
          "DNAH5", "DNTTIP1", "DOCK5", "DPPA5", "DPRX", "DRD1", "DSC2", "DSC3", "DSE", "DSG2", "DSG3", "DSP", 
          "DTX2", "DTX4", "DUOX1", "DUOX2", "DUOXA1", "DUSP11", "DUSP5", "DUSP7", "EDAR", "EFHD2", "EFNA1", 
          "EFNB1", "EFNB2", "EGFR", "EHBP1L1", "ELF3", "ELMO3", "ENC1", "ENSG00000258691", "ENSG00000266202", 
          "ENSG00000268674", "ENSG00000271810", "ENSG00000280778", "ENSG00000285130", "ENSG00000285446",
          "ENSG00000289700", "ENSG00000290317", "ENTPD3", "EPAS1", "EPB41L1", "EPHA1", "EPHA2", "EPHX3", 
          "EPHX4", "EPPK1", "EPS8L2", "EQTN", "EREG", "ERMP1", "ERO1A", "ERP27", "ERRFI1", "ESRP1", "ESRP2",
          "EVPL", "EXOC6B", "EXPH5", "EZR", "F11R", "F3", "FABP5", "FABP9", "FAM110C", "FAM166B", "FAM181B",
          "FAM183A", "FAM228A", "FAM237B", "FAM25A", "FAM25C", "FAM25G", "FAM3C", "FAM83A", "FAM83B", "FAM83C",
          "FAM83F", "FAM83H", "FAM90A27P", "FAT2", "FBXO47", "FERMT1", "FETUB", "FGD6", "FGFBP1", "FHIP1A", 
          "FLNB", "FOLR3", "FOS", "FOXE1", "FOXL1", "FOXL2", "FOXN1", "FRRS1", "FSTL4", "FZD6", "GADD45A", 
          "GALNT12", "GALNT18", "GALNT3", "GARIN2", "GARIN4", "GAST", "GBA2", "GBP6", "GCNT2", "GCNT4", "GCOM1",
          "GDPD2", "GIPC1", "GJA5", "GJB2", "GJB3", "GJB4", "GJB5", "GJB6", "GKN1", "GLTP", "GNE", "GPD2", 
          "GPR87", "GPRC5A", "GRAMD2A", "GRHL2", "GRHL3", "GSDMC", "GSDME", "GSTA4", "GSTP1", "H1-0", "HACD2",
          "HAS3", "HBEGF", "HCAR1", "HCAR2", "HCAR3", "HEPHL1", "HES2", "HLA-G", "HNRNPCL2", "HOOK2", "HOPX", 
          "HORMAD1", "HOXB1", "HOXC13", "HPSE", "HRAS", "HS3ST1", "HSD17B3", "HSD3B1", "HSF2BP", "HSFY1", "HSPA2",
          "ICAM1", "ID1", "IDS", "IER3", "IGFL1", "IL1A", "IL1B", "IL1F10", "IL1RL2", "IL1RN", "IL20RA", "IL20RB", 
          "IL23A", "IL36G", "IL36RN", "IL37", "IL6", "INAVA", "INPP4B", "INSIG2", "IQANK1", "IRAK2", "IRF6", 
          "IRX4", "ITGA2", "ITGA6", "ITGB4", "ITGB6", "ITPKC", "ITPRID2", "ITPRIPL2", "IVL", "JAG1", "JAG2", 
          "JDP2", "JMJD7", "JMJD7-PLA2G4B", "JPH2", "JUP", "KCNJ15", "KCNK7", "KCNQ3", "KCNS3", "KCTD1", "KCTD14",
          "KDF1", "KHDC1L", "KIAA1522", "KIF13A", "KIF13B", "KLF5", "KLF8", "KLHL10", "KLHL5", "KLK10", "KLK14",
          "KLK5", "KLK6", "KLK7", "KLK8", "KLK9", "KPNA7", "KRT1", "KRT13", "KRT14", "KRT15", "KRT16", "KRT17", 
          "KRT18", "KRT19", "KRT23", "KRT31", "KRT36", "KRT4", "KRT5", "KRT6A", "KRT6B", "KRT6C", "KRT7", "KRT75",
          "KRT76", "KRT77", "KRT8", "KRT80", "KRT84", "KRT9", "KRTAP20-3", "KRTAP21-2", "KRTAP3-1", "KRTAP3-3",
          "KRTAP4-6", "KRTCAP3", "KRTDAP", "LAD1", "LAMA3", "LAMB3", "LAMB4", "LAMC2", "LBX2", "LCE1A", "LCE1B",
          "LCE1C", "LCE1F", "LCE2D", "LCE3A", "LCE3C", "LCE3D", "LCE3E", "LCE6A", "LCE7A", "LCN2", "LDAF1", 
          "LDLR", "LGALS2", "LGALS7", "LGALS7B", "LGALS9B", "LHX5", "LIPH", "LIPK", "LITAFD", "LONRF3", "LPCAT4",
          "LRIG3", "LRRC1", "LRRC20", "LRRC8A", "LSR", "LTBP3", "LTBR", "LTO1", "LY6D", "LY6G6C", "LY6K", 
          "LYPD2", "LYPD3", "MACC1", "MAL2", "MALL", "MAML2", "MAP3K9", "MAPK13", "MAPKBP1", "MARCOL", "MARK1",
          "MARVELD2", "MARVELD3", "MAST4", "MBOAT1", "MBOAT2", "MBOAT7", "MCIDAS", "MEAK7", "MET", "METTL27", 
          "MFSD4A", "MIB2", "MICALL1", "MKRN2OS", "MMP13", "MMP28", "MMP9", "MOB3B", "MOB3C", "MPP7", "MPRIP",
          "MPZL2", "MPZL3", "MROH6", "MSANTD5", "MSMO1", "MST1R", "MT-ATP6", "MT-ATP8", "MT-CO1", "MT-CO2", 
          "MT-CO3", "MT-CYB", "MT-ND1", "MT-ND2", "MT-ND3", "MT-ND4", "MT-ND4L", "MT-ND5", "MT-ND6", "MT1A", 
          "MUC1", "MUC21", "MUC22", "MUC4", "MVP", "MXD1", "MYH14", "MYO5C", "MYOF", "MYZAP", "N4BP1", "NECTIN1",
          "NECTIN4", "NEDD8-MDP1", "NEDD9", "NET1", "NETO2", "NEU2", "NFE2L3", "NGEF", "NIPAL1", "NIPAL4", "NKPD1",
          "NKX1-2", "NKX2-1", "NMUR2", "NOTCH4", "NPBWR1", "NPBWR2", "NPC2", "NR1D1", "NR2F2", "NRK", "NSUN7", 
          "NT5C2", "NT5DC4", "NTF4", "NTSR1", "NXN", "OBSL1", "OCIAD2", "ODAPH", "ODF1", "ODF3L1", "OLR1", 
          "OR10H1", "OR11H7", "OR1B1", "OR1L8", "OR1P1", "OR2A7", "OR4D6", "OR5B12", "OR5D13", "OR5H14", "OR5M11", 
          "OR6C2", "OR6C70", "OSMR", "OTOP3", "OTUB2", "OTUD1", "OVCH2", "OVOL1", "OVOL2", "PACSIN3", "PAK6", 
          "PAQR7", "PARD6G", "PATJ", "PCDH1", "PDLIM1", "PDLIM4", "PDXK", "PERP", "PFN3", "PGLYRP3", "PGLYRP4", 
          "PHLDB2", "PHLDB3", "PI3", "PIK3C2B", "PINK1", "PIP4K2C", "PITPNM3", "PKD1L2", "PKP2", "PKP3", "PLA2G2F", 
          "PLA2G4E", "PLA2R1", "PLAAT3", "PLAAT4", "PLAU", "PLBD1", "PLCB3", "PLCD3", "PLCH2", "PLEKHA6", "PLEKHA7", 
          "PLEKHF1", "PLEKHG3", "PLEKHN1", "PLK3", "PLXNA1", "PNLIPRP3", "PPARD", "PPARG", "PPL", "PPP1R13L", 
          "PPP1R14C", "PPP2R3A", "PPP4R4", "PRAG1", "PRB2", "PRDM11", "PRDM9", "PRKCG", "PRKCZ", "PRMT8", "PROM2",
          "PROSER2", "PRR20G", "PRR23D1", "PRR5-ARHGAP8", "PRRG2", "PRRG4", "PRSS16", "PRSS22", "PRSS27", "PRSS8", 
          "PSCA", "PSG8", "PSORS1C1", "PSORS1C2", "PTGES", "PTGFRN", "PTGS2", "PTHLH", "PTPRF", "PTPRK", "PTPRU", 
          "PWWP2B", "RAB19", "RAB25", "RAET1E", "RAET1G", "RAET1L", "RAP2B", "RARG", "RARRES1", "RASSF7", "RASSF9",
          "RBMS2", "RDH16", "REN", "RGP1", "RHBDF2", "RHBDL2", "RHCG", "RHOBTB3", "RHOD", "RHOF", "RIN2", "RIPK4",
          "RNASE10", "RNASE7", "RND1", "RNF103", "RNF141", "RNF152", "RNF222", "RNF223", "RNF225", "RNF227", "RNF39",
          "RPL3L", "RPTN", "RXRA", "S100A10", "S100A11", "S100A12", "S100A2", "S100A3", "S100A4", "S100A7", "S100A7A",
          "S100A8", "S100A9", "S1PR5", "SAA1", "SAA2", "SAA2-SAA4", "SAT1", "SBSN", "SCEL", "SCGB1A1", "SCPEP1", 
          "SDC4", "SDR42E1", "SDR9C7", "SERINC2", "SERPINB13", "SERPINB5", "SESN3", "SFN", "SFR1", "SFTA2", 
          "SFTPA1", "SFTPA2", "SH2D3A", "SH3BP1", "SH3RF2", "SH3YL1", "SHB", "SKAP2", "SLC10A6", "SLC15A1", 
          "SLC16A5", "SLC1A3", "SLC1A6", "SLC22A31", "SLC26A9", "SLC29A3", "SLC2A9", "SLC31A2", "SLC35G3", 
          "SLC37A2", "SLC39A2", "SLC3A2", "SLC44A2", "SLC44A3", "SLC52A1", "SLC6A11", "SLC6A14", "SLC7A5", 
          "SLC7A7", "SLC8B1", "SLC9A1", "SLFN5", "SLK", "SLURP1", "SLURP2", "SMAD3", "SMAGP", "SMCO2", "SMCP",
          "SMIM5", "SMPDL3B", "SNCG", "SNX31", "SOD2", "SOWAHB", "SOX15", "SOX9", "SPATS2L", "SPINK6", "SPINT1",
          "SPINT4", "SPRR1A", "SPRR1B", "SPRR2A", "SPRR2B", "SPRR2D", "SPRR2E", "SPRR2G", "SPRR3", "SPRR4", "SPRR5",
          "SPTBN1", "SPTLC2", "SPZ1", "SRC", "SRGAP3", "SSH3", "SSX4", "ST14", "ST3GAL1", "STAP2", "STK39", "STN1.00",
          "STX19", "STYK1", "SULF2", "SYT15", "SYT16", "SYT8", "TAB3", "TACSTD2", "TBC1D2", "TBC1D26", "TBC1D3F", 
          "TBC1D3G", "TBPL2", "TCF7L2", "TCHH", "TCHHL1", "TCIRG1", "TEAD3", "TEF", "TENT5B", "TESK1", "TGM1", "TGM2",
          "THBD", "TIGD2", "TINAGL1", "TINCR", "TIPARP", "TJP1", "TJP2", "TLN1", "TM4SF1", "TMBIM1", "TMC5", "TMC7", 
          "TMEM125", "TMEM156", "TMEM238", "TMEM265", "TMEM30B", "TMEM40", "TMEM51", "TMEM53", "TMEM79", "TMEM87B",
          "TMEM8B", "TMIGD1", "TMOD3", "TMPRSS11A", "TMPRSS11B", "TMPRSS11E", "TMPRSS12", "TMPRSS13", "TMPRSS4", 
          "TMPRSS7", "TNFAIP2", "TNFRSF10A", "TNFRSF10B", "TNFRSF21", "TNFSF13", "TNIP1", "TNK2", "TNKS1BP1", 
          "TNNT3", "TNS4", "TOM1L2", "TOR4A", "TP53AIP1", "TP63", "TPBG", "TPCN1", "TPRA1", "TREM2", "TRIM29", 
          "TRIM31", "TRIM6-TRIM34", "TSBP1", "TSLP", "TSPAN15", "TSPAN6", "TTC9", "TTLL11", "TUFT1", "TXNDC8", 
          "TYMP", "UBAP1", "UGT1A10", "UGT1A8", "ULBP2", "UNC13D", "UNC93B1", "UPK1A", "UPK1B", "UPK2", "UPP1",
          "USH1G", "USP17L23", "USP31", "VGLL1", "VSIG8", "VSIR", "VSTM2L", "VXN", "WFDC10B", "WFDC12", "WFDC13",
          "WFDC3", "WFDC5", "WNT10A", "WNT3A", "WNT6", "WNT7B", "WNT9A", "XAGE2", "XAGE3", "XDH", "ZBTB7C", "ZC3H12A",
          "ZDHHC3", "ZNF114", "ZNF165", "ZNF185", "ZNF311", "ZNF488", "ZNF750", "ZSCAN5C")

ksGenescore <- function(
    n.x,               # length of ranked list
    y,                 # positions of geneset items in ranked list (basically, ranks)
    do.pval=TRUE,      # compute asymptotic p-value
    alternative=c("two.sided","greater","less"),
    do.plot=F,         # draw the ES plot
    bare=FALSE,        # return score & p-value only (a 2-tuple)
    cls.lev=c(0,1),    # class labels to display
    absolute=FALSE,    # takes max - min score rather than the maximum deviation from null
    plot.labels=FALSE, # hits' labels
    ...                # additional plot arguments
)
{
  # efficient version of ks.score (should give same results as ks.test, when weight=NULL)
  #
  alternative <- match.arg(alternative)
  n.y <- length(y)
  if ( n.y < 1 )  stop("Not enough y data")
  if ( any(y>n.x) ) stop( "y must be <= n.x: ", max(y) )
  if ( any(y<1) ) stop( "y must be positive: ", min(y) )
  x.axis <- y.axis <- NULL
  
  # KS score
  #
  y <- sort(y)
  n <- n.x * n.y/(n.x + n.y)
  hit <- 1/n.y
  mis <- 1/n.x
  
  ## to compute score, only the y positions and their immediate preceding
  ## ..positions are needed
  ##
  Y <- sort(c(y-1,y)); Y <- Y[diff(Y)!=0]; y.match <- match(y,Y)
  D <- rep( 0, length(Y) ); D[y.match] <- (1:n.y)
  zero <- which(D==0)[-1]; D[zero] <- D[zero-1]
  
  z <- D*hit - Y*mis
  
  score <- if (absolute) max(z)-min(z) else z[which.max(abs(z))]
  names(score) <- "D"
  
  if (do.plot) {
    x.axis <- Y;
    y.axis <- z;
    if(Y[1]>0) {
      x.axis <- c(0,x.axis);
      y.axis <- c(0,y.axis);
    }
    if ( max(Y)<n.x ) {
      x.axis <- c(x.axis,n.x)
      y.axis <- c(y.axis,0)
    }
    plot( x.axis, y.axis, type="l",
          xlab=paste("up-regulated for class ", cls.lev[2], " (KS>0) vs ",
                     "up-regulated for class ", cls.lev[1], " (KS<0)", sep="" ),
          ylab="gene hits",...)
    abline(h=0)
    abline(v=n.x/2,lty=3)
    axis(1,at=y,labels=plot.labels,tcl=0.25,las=2)
    i.max <- which.max(abs(y.axis))
    points( x.axis[i.max], y.axis[i.max], pch=20, col="red")
    text(x.axis[i.max]+n.x/20,y.axis[i.max],round(y.axis[i.max],2))
  }
  if ( !do.pval ) {
    return(score)
  }
  ## ELSE compute p-value as in function ks.test but return signed statistic
  ##
  tmp <- suppressWarnings(ks.test(1:n.x,y,alternative=alternative))
  tmp$statistic <- score # use the signed statistic
  return( if (bare) c(tmp$statistic, tmp$p.value) else tmp )
}


hh
cc

hhks <- hh
ccks <- cc

#now k sig onto hh
hhkssort <- hhks[order(hhks$log2FoldChange, decreasing = TRUE),]
hhkssort$bulkorder <- seq.int(nrow(hhkssort))
hhkssort1 <- hhkssort[hhkssort$Gene.name %in% intersect(hhkssort$Gene.name, ksig),]
print("ksig projected on hh bulk. X-axis left is bulk upreg in hh, right X-axis is downregin hh")
ksGenescore(nrow(hhkssort), hhkssort1$bulkorder, do.plot = T)

#k sig onto cc
cckssort <- ccks[order(ccks$log2FoldChange, decreasing = TRUE),]
cckssort$bulkorder <- seq.int(nrow(cckssort))
cckssort1 <- cckssort[cckssort$Gene.name %in% intersect(cckssort$Gene.name, ksig),]
print("ksig projected on cc bulk. X-axis left is bulk upreg in hh, right X-axis is downregin cc")
ksGenescore(nrow(cckssort), cckssort1$bulkorder, do.plot = T)

          

library(DESeq2)
library(Biobase)
library(AUCell)
library(tximport)
library(hypeR)
library(PMCMRplus)
library(GSVA)
setwd('/restricted/projectnb/montilab-p/projects/oralcancer/serine_starvation/')

#gene name sigs
sigs <- readRDS('/restricted/projectnb/montilab-p/projects/oralcancer/serine_starvation/results/DESeq2/20230501_signatureslfc1.rds')
#ENS name sigs (same signatures just different name for grade)
sigs <- readRDS('/restricted/projectnb/montilab-p/projects/oralcancer/serine_starvation/results/DESeq2/20230901_signatureslfc1_ENSID.rds')
names(sigs) <- c("CSS", "CC", "HSS", "HC")

data <- readRDS("/restricted/projectnb/montilab-p/CBMrepositoryData/TCGA/TCGA-GDC/RNAseq/esets_postprocessed/TCGA-OSCC_2020-03-22_wLegacyPData_DESeq2_log_eSet.rds")# TCGA-GDC/RNAseq/esets_postprocessed/TCGA-HNSC_2020-03-22_wLegacyPData_DESeq2_log_eSet.rds")


neg <- data[, data$hpv_status == "negative" ]

#score with signatures
hnscgsva <- gsva(neg, sigs, verbose = FALSE)

#########SURVIVAL
#survival info
pData(hnscgsva) <- pData(hnscgsva) %>% mutate(time = case_when(days_to_death != "--" ~ hnscgsva$days_to_death, 
                                                               days_to_death == "--" ~ hnscgsva$days_to_last_follow_up)
)

hnscgsva$vital_status1[hnscgsva$vital_status == "Alive"] = 1
hnscgsva$vital_status1[hnscgsva$vital_status == "Dead"] = 2

pData(hnscgsva) <- pData(hnscgsva) %>% mutate(x5year = case_when(vital_status == "Alive" ~ as.character(hnscgsva$days_to_last_follow_up),
                                                                 (vital_status == "Alive" & as.numeric(days_to_last_follow_up) >= 1825.0) ~ as.character(1825),
                                                                 (vital_status == "Dead" & as.numeric(days_to_death) < 1825.0) ~ as.character(hnscgsva$time),
                                                                 (vital_status == "Dead" & as.numeric(days_to_death) >= 1825.0) ~ as.character(1825))
)

pData(hnscgsva) <- pData(hnscgsva) %>% mutate(vital_status5 = case_when(vital_status == "Alive" ~ 1,
                                                                        (vital_status == "Dead" & as.numeric(days_to_death) < 1825.0) ~ 2,
                                                                        (vital_status == "Dead" & as.numeric(days_to_death) >= 1825.0) ~ 1)
)


addgsvadata <- function(obj, sigsxx) {
  for (sig in names(sigsxx)) {
    print(sig)
    threshold <- median(exprs(obj[sig,]))
    print(threshold)
    obj[[paste0(sig,"val")]] <- t(exprs(obj[sig,]))
    obj[[paste0(sig,"stat")]] <- with(obj, ifelse(obj[[paste0(sig,"val")]] <= threshold, "low", "high"))
  }
  return(obj)
}

test <- addgsvadata(hnscgsva, sigs)
test$Cdiff <- test$CCval - test$CSSval
test$Hdiff <- test$HCval - test$HSSval

test$age <- as.numeric(test$age_at_diagnosis)/365
test$age1 <- scale(test$age)

threshold1 <- 66
theshold2 <- 40
test$agestat <- with(test, ifelse(test$age <= threshold1, "young", "old"))
test$agestat1 <- with(test, ifelse(test$age >= theshold2 & test$agestat == 'young', "mid", test$agestat))

test1 <- test[,!is.na(test$age)]

Cdiffcoxph <- coxph(formula = Surv(as.numeric(test1$time), test1$vital_status1) ~ test1$Cdiff + test1$age)
summary(Cdiffcoxph)
Cdiffcoxph <- coxph(formula = Surv(as.numeric(test1$time), test1$vital_status1) ~ test1$Cdiff + strata(test1$agestat1))
summary(Cdiffcoxph)
Hdiffcoxph <- coxph(formula = Surv(as.numeric(test1$time), test1$vital_status1) ~ test1$Hdiff + test1$age)
summary(Hdiffcoxph)
Hdiffcoxph <- coxph(formula = Surv(as.numeric(test1$time), test1$vital_status1) ~ test1$Hdiff + strata(test1$agestat1))
summary(Hdiffcoxph)

#######GRADE AND STAGE
mat <- pData(test1)
my_comparisons <- list( c("AN","stage i"), c("AN","stage ii"), c("AN","stage iii"), c("AN","stage iv"), c("stage i","stage ii"), c("stage i","stage iii"), c("stage i","stage iv"), c("stage ii","stage iii"), c("stage ii","stage iv") , c("stage iii","stage iv"))

mat$my_stage <- factor(mat$my_stage, levels = c("AN", "stage i", "stage ii", "stage iii", "stage iv"))

p <- ggplot(mat, aes(y=Hdiff, x=my_stage)) + 
  geom_boxplot() + scale_fill_brewer(palette="Blues") + theme_classic()
p + geom_jitter(shape=16, position=position_jitter(0.2)) + 
  ggpubr::stat_compare_means(label.y=1) +
  ggpubr::stat_compare_means(comparisons = my_comparisons) + geom_jitter(shape=16, position=position_jitter(0.2)) 

jonckheereTest(mat$my_stage, mat$Hdiff)

p <- ggplot(pData(test), aes(y=Cdiff, x=my_stage)) + 
  geom_boxplot() + theme_classic()
p + geom_jitter(shape=16, position=position_jitter(0.2)) + 
  ggpubr::stat_compare_means(label.y=1) +
  ggpubr::stat_compare_means(comparisons = my_comparisons) + geom_jitter(shape=16, position=position_jitter(0.2)) 

jonckheereTest(mat$my_stage, mat$Cdiff)


#grade
unique(mat$my_grade)
my_comparisons <- list( c("AN","g1"), c("AN","g2"), c("AN","g3"), c("g1","g2"), c("g1", "g3"), c("g2","g3"))

mat$my_grade <- factor(mat$my_grade, levels = c("AN", "g1", "g2", "g3", "gx"))

p <- ggplot(mat, aes(y=Hdiff, x=my_grade)) + 
  geom_boxplot() + scale_fill_brewer(palette="Blues") + theme_classic()
p + geom_jitter(shape=16, position=position_jitter(0.2)) + 
  ggpubr::stat_compare_means(label.y=1) +
  ggpubr::stat_compare_means(comparisons = my_comparisons) + geom_jitter(shape=16, position=position_jitter(0.2)) 

jonckheereTest(mat$my_grade, mat$Hdiff)

p <- ggplot(pData(test), aes(y=Cdiff, x=my_grade)) + 
  geom_boxplot() + theme_classic()
p + geom_jitter(shape=16, position=position_jitter(0.2)) + 
  ggpubr::stat_compare_means(label.y=1) +
  ggpubr::stat_compare_means(comparisons = my_comparisons) + geom_jitter(shape=16, position=position_jitter(0.2)) 

jonckheereTest(mat$my_grade, mat$Cdiff)

#
newdf <- pData(test1)[,88:99]
newdf <- pData(test1)[,88:100]
newdf <- pData(test1)[,c(88,96,98,100,102)]
library(reshape2)
melt_data <- melt(newdf, scores = c("CC","CSS","HCval", "HSSval")) 
melt_data <- melt_data %>% mutate(cellline =
                                    case_when(variable == 'CCval' ~ "Cal27",
                                              variable == 'CSSval' ~ "Cal27",
                                              variable == 'HCval' ~ "HSC3",
                                              variable == 'HSSval' ~ "HSC3"))
melt_data <- melt_data %>% mutate(condition =
                                    case_when(variable == 'CCval' ~ "Complete",
                                              variable == 'CSSval' ~ "SerineStarvation",
                                              variable == 'HCval' ~ "Complete",
                                              variable == 'HSSval' ~ "SerineStarvation"))
melt_data <- melt_data %>% mutate(variable1 =
                                    case_when(variable == 'CCval' ~ "Cal27 Complete",
                                              variable == 'CSSval' ~ "Cal27 SS",
                                              variable == 'HCval' ~ "HSC3 Complete",
                                              variable == 'HSSval' ~ "HSC3 SS"))
melt_data$variable <- factor(melt_data$variable, levels = c("CCval", "CSSval", "HCval", "HSSval"))
melt_data$variable1 <- factor(melt_data$variable1, levels = c("Cal27 Complete", "Cal27 SS", "HSC3 Complete", "HSC3 SS"))
#melt_data1 <- melt_data[melt_data$my_grade!="NA",]
melt_data1 <- melt_data[!is.na(melt_data$my_grade),]
melt_data2 <- melt_data1[melt_data1$my_grade!="gx",]

ggplot(melt_data2, aes(y=value, x=my_grade, fill = condition))+
  geom_boxplot()+ theme_classic()+
  facet_wrap(.~variable1, scales = "free", ncol = 2)+
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_manual(values=c("#666666", "#009900")) +
  ylim(c(-0.5,0.6)) + 
  xlab("Grade") + ylab("Signature Score")
mat1 <- mat[mat$my_grade!="NA",]
mat2 <- mat[mat$my_grade!="gx",]
jonckheereTest(mat$my_grade, mat$CCval)
jonckheereTest(mat$my_grade, mat$CSSval)
jonckheereTest(mat$my_grade, mat$HCval)
jonckheereTest(mat$my_grade, mat$HSSval)

ggplot(melt_data2, aes(y=value, x=my_grade, fill = condition))+
  geom_boxplot()+ theme_classic()+
  facet_wrap(.~cellline, scales = "free", ncol = 1)+
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_manual(values=c("#666666", "#009900"))

pdf('/restricted/projectnb/montilab-p/projects/oralcancer/serine_starvation/results/DESeq2/TCGAgrade.pdf', width = 6, height = 5)
dev.off()

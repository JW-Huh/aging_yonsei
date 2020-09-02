

# This script is written for actual display 
# of human gut microbiome data

##########      Aging      ##########

rm(list=ls())
options(java.parameters = "-Xmx64g", stringsAsFactors = F) 
#setwd(dir = "D:/1_Data/2020/0530_Aging_combined/")

library(vegan) ; library(Hmisc) ; library(MASS) ; library(dplyr)

memory.size(max=TRUE) ; memory.size(max=FALSE) ; memory.limit(size=NA) ; memory.limit(size=56000) ; memory.limit(size=NA) 


######  1) Firmicutes / Bacteroidetes Ratio   ######

# F/B ratio is one of common indicator of gut microbiome

# First of all, we needs a phylum level abundance data

p_AGP <- 
  read.csv("input/AGP_phylum.csv", row.names = 1) ; rownames(p_AGP) <- 
  gsub(rownames(p_AGP), pattern = "X", replacement = "")

p_China <- 
  read.csv("input/China_phylum.csv", row.names = 1) ; rownames(p_China) <- gsub(rownames(p_China), pattern = "X", replacement = "")

# m_agp <- meta.sub ; m_chi <- meta.ch
m_agp <- read.csv("input/meta_AGP_rarefied.csv", row.names = 1)
m_chi <- read.csv("input/meta_China_rarefied.csv", row.names = 1)

p_AGP <- p_AGP %>%
  mutate(Age=m_agp$AGE_CORRECTED,
         Sex=m_agp$SEX,
         Project="AGP")

p_China <- p_China %>%
  mutate(Age=m_chi$age,
         Sex=m_chi$sex,
         Project="China")

intersect(colnames(p_China), colnames(p_AGP))

p_all <- rbind(p_AGP, p_China) ;  head(p_all)

# Add F/B ratio
p_all <- p_all %>% 
    mutate(Age_3class=case_when(Age<=40 ~ "Young",
                                Age<=60 & Age>40 ~ "Middle",
                                Age>60 ~ "Old"),
           Age_2class=ifelse(Age<=45, "Young", "Old"),
           FB.Ratio=(Firmicutes+1)/(Bacteroidetes+1)) 

tiff("figures/Phylum_FB_ratio.tiff", res = 150, 
     width=600, height = 800, units = "px");ggplot(p_all, aes(Age, log10(FB.Ratio))) +
  xlab("Age") + ylab("Firmicutes/Bacteroidetes")+
  geom_point(alpha=0.01)+
  geom_smooth(method = "lm") +
  theme_classic() +
  stat_cor(method = "spearman", size=3.5)+
  facet_grid(Sex~Project) +
  theme(axis.text = element_text(size = rel(0.8), color="black")); dev.off()


######  2) PCoA   ######

# To visualize overall microbiome at the genus + unassigned_family level

# 2-1) Family + Genus level PCA
fg_all <- read.csv("input/Combined_family_N_genus.csv", row.names = 1)

# fg_dist <- vegdist(fg_all, method="bray")
# fg_pcoa <- cmdscale(fg_dist, k=10, eig=T)

# metadata modification_1) Alcohol consumption
m_chi <- m_chi %>%
  mutate(Alcohol=ifelse(low_alcohol_liquor==0 & high_alcohol_liquor==0,
                         "No", "Yes")) 

m_agp <- m_agp %>%
  mutate(Alcohol=ifelse(ALCOHOL_CONSUMPTION %in% c("FALSE", "No"),
                        "No", "Yes")) ; m_agp$Alcohol[m_agp$ALCOHOL_CONSUMPTION=="Unknown"] <- "Unknown"


# metadata modification_2) Sugar consumption
m_agp$Sugar <- 
  "Yes" ; m_agp$Sugar[m_agp$SUGAR_SWEETENED_DRINK_FREQUENCY=="Never" & m_agp$SUGARY_SWEETS_FREQUENCY=="Never"] <- 
  "No"  ; m_agp$Sugar[m_agp$SUGAR_SWEETENED_DRINK_FREQUENCY %in% 
                        c("Unknown", "Unspecified") | 
                        m_agp$SUGARY_SWEETS_FREQUENCY %in% 
                        c("Unspecified", "Unknown")] <- 
  "Unknown"

m_chi$Sugar <- "Yes" ; m_chi$Sugar[m_chi$sugar==0 &
                                     m_chi$soft_drinks==0] <- "No"

# metadata modification_3) Meat consumption (livestock meat)
m_agp$Meat <- "No" ; m_agp$Meat[m_agp$DIET_TYPE %in% c("Omnivore", "Omnivore but do not eat red meat")] <- "Yes"

m_chi$Meat <- "Yes" ; m_chi$Meat[m_chi$livestock_meat==0] <- "No"


# metadata modification_4) Smoking
table(m_chi$smoking, m_chi$sex) # Factor number 4 might be non-smoker
m_chi$Smoking <- ifelse(m_chi$smoking==4, "No", "Yes")

table(m_agp$SMOKING_FREQUENCY)
m_agp$Smoking <- 
  "Yes" ; m_agp$Smoking[m_agp$SMOKING_FREQUENCY=="Never"] <- 
  "No"  ; m_agp$Smoking[m_agp$SMOKING_FREQUENCY %in%
                          c("Unknown", "Unspecified")] <- "Unknown"


m_all <- rbind(data.frame(ID=m_agp$SampleID,
                          Age=m_agp$AGE_CORRECTED,
                          Sex=m_agp$SEX,
                          BMI=m_agp$BMI,
                          Meat=m_agp$Meat,
                          Alcohol=m_agp$Alcohol,
                          Sugar=m_agp$Sugar,
                          Smoking=m_agp$Smoking,
                          Project="AGP"),
               data.frame(ID=m_chi$sample_name,
                          Age=m_chi$age,
                          Sex=m_chi$sex,
                          BMI=m_chi$bmi,
                          Meat=m_chi$Meat,
                          Alcohol=m_chi$Alcohol,
                          Sugar=m_chi$Sugar,
                          Smoking=m_chi$Smoking,
                          Project="China"))

fg_summary <- data.frame(fg_pcoa$points, m_all)
colnames(fg_summary)[1:10] <- paste0("PC", 1:10)

fg_vars <- round((fg_pcoa$eig^2/sum(fg_pcoa$eig^2))[1:10]*100, 1)

# 1) Project
pcoa_short <- function(data, color) {
  
  attach(data)
  ggplot(data, aes(PC1, PC2, color=color)) +
    xlab(paste0("PC1 (", fg_vars[1], "%)"))+
    ylab(paste0("PC2 (", fg_vars[2], "%)"))+
    geom_point(alpha=0.5, na.rm = T) +
    theme_classic() +
    coord_fixed()+
    theme(axis.text = element_text(size=rel(0.9), 
                                   color="black"),
          legend.position = "top",
          legend.title = element_blank())
    
}

# 1) Project
tiff("figures/Both_PCoA_project.tiff", res=150, units = "px", width = 600, height = 650); pcoa_short(fg_summary, Project) ; dev.off()


# 2) Sex
tiff("figures/Both_PCoA_sex.tiff", res=150, units = "px", width = 600, height = 650); pcoa_short(fg_summary, Sex) ; dev.off()

# 3) Age
tiff("figures/Both_PCoA_Age.tiff", res=150, units = "px", width = 600, height = 650); pcoa_short(fg_summary, Age) +
  scale_color_viridis_c(direction = -1); dev.off()

# 4) Smoking
fg_summary_Smoking <- fg_summary %>% filter(Smoking %in% c("Yes", "No"))
tiff("figures/Both_PCoA_smoking.tiff", res=150, units = "px", width = 600, height = 650); pcoa_short(fg_summary_Smoking, Smoking) ; dev.off()

# 5) Shannon
tiff("figures/Both_PCoA_Shannon.tiff", res=150, units = "px", width = 600, height = 650); pcoa_short(fg_summary, Shannon) +
  scale_color_viridis_c(direction = -1); dev.off()

# 6) Richness
tiff("figures/Both_PCoA_Richness.tiff", res=150, units = "px", width = 600, height = 650); pcoa_short(fg_summary, log2(Richness)) +
  scale_color_viridis_c(direction = -1); dev.off()

# 7) Evenness
tiff("figures/Both_PCoA_Evenness.tiff", res=150, units = "px", width = 600, height = 650); pcoa_short(fg_summary, Evenness) +
  scale_color_viridis_c(direction = -1); dev.off()

# 8) Firmicutes/Bacteroidetes
tiff("figures/Both_PCoA_FB_ratio.tiff", res=150, units = "px", width = 600, height = 650); pcoa_short(fg_summary, FB.Ratio) +
  scale_color_viridis_c(direction = -1); dev.off()


tiff("figures/Both_PC1_Age_Sex.tiff", res=150, units = "px",
     width = 550, height = 650);ggplot(fg_summary, aes(Age, PC1)) +
  xlab("Age") + ylab("PC1 (51.1%)")+
  geom_jitter(alpha=0.01)+
  geom_smooth(method = "loess", aes(color=Sex)) +
  facet_wrap(~Sex) +
  theme_classic() +
  geom_vline(xintercept = c(40, 50), linetype="dotted", alpha=0.9)+
  theme(axis.text = element_text(size=rel(1.2), color="black"),
        axis.title.x = element_text(size=rel(1.6)),
        axis.title.y = element_text(size=rel(1.1)),
        strip.text = element_text(size=rel(1.6)),
        legend.position = "none"); dev.off()

  
# by Sex
tiff("figures/Both_Age_Sex_log-Abundance.tiff", res = 160,
     width = 2000, height = 650, units = "px"); ggplot(p_melt, aes(Age, log10(value+1e-5))) +
  geom_point(alpha=0.01)+
  xlab("Age") + ylab("Log10 (relative abundance+1e-5)") +
  geom_smooth(method="loess", aes(color=Sex)) +
  stat_cor(method = "spearman", size=2.6, label.y = 0.05) +
  theme_classic() +
  facet_grid(Sex~variable) +
  theme(axis.text = element_text(size=rel(1.2), color="black"),
        axis.title.x = element_text(size=rel(2.2)),
        axis.title.y = element_text(size=rel(1.1)),
        strip.text = element_text(size=rel(1)),
        legend.position = "none"); dev.off()





######  3) Diversity plolt   ######

# compare alpha-diversity 

#rm(list=ls())
fg_agp <- read.csv("input/AGP_family_N_genus.csv", row.names = 1)
fg_chi <- read.csv("input/China_family_N_genus.csv", row.names = 1)


div_AGP <- data.frame(SampleID=rownames(fg_agp),
                      Richness=rowSums(fg_agp>0), # Richness
                      Shannon=diversity(fg_agp), # "Shannon"
                      Age=m_agp$AGE_CORRECTED,
                      Sex=m_agp$SEX,
                      Project="AGP") ; div_AGP$Evenness <- 
  div_AGP$Shannon/log(div_AGP$Richness) ;rownames(div_AGP) <- NULL

div_CHI <- data.frame(SampleID=rownames(fg_chi),
                      Richness=rowSums(fg_chi>0), # Richness
                      Shannon=diversity(fg_chi), # "Shannon"
                      Age=m_chi$age,
                      Sex=m_chi$sex,
                      Project="China") ; div_CHI$Evenness <- 
  div_CHI$Shannon/log(div_CHI$Richness) ;rownames(div_CHI) <- NULL


div_all <- rbind(div_AGP, div_CHI) ; head(div_all)


# Change of microbial diversity by age or sex

# 1) Shannon
tiff("figures/F1_both_Shannon_diveristy_by_Age.tiff",
     res = 160, width = 800, height = 900,
     units = "px"); ggplot(div_all, aes(Age, Shannon, color=Sex)) +
  xlab("Age") + ylab("Shannon diversity")+
  geom_smooth(method = "loess", size=0.8)+
  facet_grid(Sex~Project, scales = "free") +
  stat_cor(method = "spearman", color="black", label.y = 2.6, size=4.4) +
  theme_classic() +
  theme(legend.position = "none",
        axis.text = element_text(color="black",
                                 size=rel(1.4)),
        axis.title = element_text(size=rel(2.4), vjust=1),
        strip.text = element_text(size=rel(1.4))); dev.off()

# 2) Richness
tiff("figures/F1_both_Richness_diveristy_by_Age.tiff",
     res = 160, width = 800, height = 900,
     units = "px"); ggplot(div_all, aes(Age, Richness, color=Sex)) +
  xlab("Age") + ylab("Richness")+
  geom_smooth(method = "loess", size=0.8)+
  facet_grid(Sex~Project, scales = "free") +
  stat_cor(method = "spearman", color="black", label.y = 115, size=4.4) +
  theme_classic() +
  theme(legend.position = "none",
        axis.text = element_text(color="black",
                                 size=rel(1.4)),
        axis.title = element_text(size=rel(2.4), vjust=1),
        strip.text = element_text(size=rel(1.4))); dev.off()

# 3) Evenness
tiff("figures/F1_both_Evenness_diveristy_by_Age.tiff",
     res = 160, width = 800, height = 900,
     units = "px"); ggplot(div_all, aes(Age, Evenness, color=Sex)) +
  xlab("Age") + ylab("Evenness")+
  geom_smooth(method = "loess", size=0.8)+
  facet_grid(Sex~Project, scales = "free") +
  stat_cor(method = "spearman", color="black", label.y = 0.6, size=4.4) +
  theme_classic() +
  theme(legend.position = "none",
        axis.text = element_text(color="black",
                                 size=rel(1.4)),
        axis.title = element_text(size=rel(2.4), vjust=1),
        strip.text = element_text(size=rel(1.4))); dev.off()


identical(rownames(fg_summary), div_all$SampleID) # TRUE
  fg_summary$Shannon <- div_all$Shannon
  fg_summary$Richness <- div_all$Richness
  fg_summary$Evenness <- div_all$Evenness


######  4) Association between metadata and PCs   ######
  
  # Spearman correlation test for continuous metadata (ex. age, BMI)
  # Spearman correlation test for ordinal metadata (ex. diet type)
  # Wilcoxon rank sum test for categorical metadata (ex. sex)
  # Detailed description is written in method section of manuscript
  
# BMI
fg_summary_BMI <- fg_summary[!fg_summary$BMI %in% 
                               c("no_data", "Unknown", "Unspecified"),]
fg_summary_BMI$BMI <- as.numeric(fg_summary_BMI$BMI) 
fg_summary_BMI <- subset(fg_summary_BMI, BMI<=60 & BMI>=10) # 12905

# Meat 
table(fg_summary$Meat) # Omnivore 12156 vs. No livestock meat 1030

# Sugar
table(fg_summary$Sugar) # 2538 Unknown ; yes 8512 ; no 2126
fg_summary_sugar <- subset(fg_summary, Sugar!="Unknown")

# Alcohol
table(fg_summary$Alcohol) # 18 Unknown 
fg_summary_EtOH <- subset(fg_summary, Alcohol!="Unknown")

# Smoking
table(fg_summary$Smoking) # 25 Unknown 


# making input filese to calculate p-value
PC_pro <- fg_summary[, c("Project", paste0("PC", 1:10))]
PC_age <- fg_summary[, c("Age", paste0("PC", 1:10))]
PC_bmi <- fg_summary_BMI[, c("BMI", paste0("PC", 1:10))]
PC_sex <- fg_summary[, c("Sex", paste0("PC", 1:10))]
PC_meat <- fg_summary[, c("Meat", paste0("PC", 1:10))]
PC_sugar <- fg_summary_sugar[, c("Sugar", paste0("PC", 1:10))]
PC_EtOH <- fg_summary_EtOH[, c("Alcohol", paste0("PC", 1:10))]
PC_smoke <- fg_summary_Smoking[, c("Smoking", paste0("PC", 1:10))]


# Age, BMI
PC_pvalue <- data.frame(Age=rcorr(as.matrix(PC_age), 
                                  type = "spearman")$P[,"Age"],
                        BMI=rcorr(as.matrix(PC_bmi), 
                                  type = "spearman")$P[,"BMI"])

# To calculate exact p-values that were described as 0.000000e+00
library("pspearman") 
test <- spearman.test(x = PC_age$Age, y = PC_age$PC1, approximation = "t-distribution") ; PC_pvalue$Age[2] <- test$p.value

test <- spearman.test(x = PC_age$Age, y = PC_age$PC4, approximation = "t-distribution") ; PC_pvalue$Age[5] <- test$p.value

test <- spearman.test(x = PC_age$Age, y = PC_age$PC5, approximation = "t-distribution") ; PC_pvalue$Age[6] <- test$p.value

test <- spearman.test(x = PC_bmi$BMI, y = PC_bmi$PC4, approximation = "t-distribution") ; PC_pvalue$BMI[5] <- test$p.value


# Sex
test <- apply(PC_sex[,2:11], 2, function(x) wilcox.test(x~PC_sex[,1])) ; PC_pvalue$Sex <- NA ; for (i in 1:10){ PC_pvalue$Sex[i+1] <- as.numeric(test[[i]][3])}

# Meat
test <- apply(PC_meat[,2:11], 2, function(x) wilcox.test(x~PC_meat[,1])) ; PC_pvalue$Meat <- NA ;for (i in 1:10){ PC_pvalue$Meat[i+1] <- as.numeric(test[[i]][3])}

# Sugar
test <- apply(PC_sugar[,2:11], 2, function(x) wilcox.test(x~PC_sugar[,1])) ; PC_pvalue$Sugar <- NA ;for (i in 1:10){ PC_pvalue$Sugar[i+1] <- as.numeric(test[[i]][3])}

# Alcohol
test <- apply(PC_EtOH[,2:11], 2, function(x) wilcox.test(x~PC_EtOH[,1])) ; PC_pvalue$Alcohol <- NA ;for (i in 1:10){ PC_pvalue$Alcohol[i+1] <- as.numeric(test[[i]][3])}

# Smoking
test <- apply(PC_smoke[,2:11], 2, function(x) wilcox.test(x~PC_smoke[,1])) ; PC_pvalue$Smoking <- NA ;for (i in 1:10){ PC_pvalue$Smoking[i+1] <- as.numeric(test[[i]][3])}

# # Project
# test <- apply(PC_pro[,2:11], 2, function(x) wilcox.test(x~PC_pro[,1])) ; PC_pvalue$Project <- NA ;for (i in 1:10){ PC_pvalue$Project[i+1] <- as.numeric(test[[i]][3])}


# FDR, p-value adjustment
PC_pvalue <- apply(PC_pvalue, 2, 
                   function(x) p.adjust(x, method = "fdr")) 
PC_pvalue <- PC_pvalue[-1,]
PC_pvalue$Variance <- 10^(-(fg_vars)) # This will be log transformed

#write.csv(PC_pvalue, "input/P.value_metadata_summary.csv")

# Melting
PC_pvalue$Alcohol <- 
  NULL ; PC_pvalue$Smoking <- 
  NULL ; PC_pvalue$PC <- 
  rownames(PC_pvalue)

PC_melt <- 
  reshape::melt(PC_pvalue, id="PC") ; PC_melt$PC <- 
  factor(PC_melt$PC, levels=c(paste0("PC", 10:1))) ; PC_melt$variable <-
  factor(PC_melt$variable, 
         levels=c("Variance", "Sex", "Age", "BMI", "Sugar", "Meat"))

PC_melt$FDR <- "<10" ; PC_melt$FDR[PC_melt$variable!="Variance" & (-log10(PC_melt$value)>10)] <- 
  ">10" ; PC_melt$FDR[PC_melt$variable!="Variance" & (-log10(PC_melt$value)>30)] <- 
  ">30" ; PC_melt$FDR <- 
  factor(PC_melt$FDR, levels=c("<10", ">10", ">30"))

# PCs 
tiff("figures/Both_PCs_Bray_All.tiff", units = "px",
     res = 180, width = 1800, height = 500);ggplot(PC_melt, aes(PC, -log10(value)))+
  geom_bar(stat = "identity", aes(fill=FDR))+
  xlab("") + ylab("")+
  scale_fill_manual(values=c("#4B4B4B", "#FF4300", "#FF0087"))+
  theme_classic()+
  coord_flip()+
  facet_wrap(.~variable, nrow = 1, scale="free_x") +
  theme(axis.text.x = element_text(size=rel(1.2), color="black"),
        axis.text.y = element_text(size=rel(1.3), color="black"),
        strip.text = element_text(size=rel(1.2)),
        legend.justification = c(0,0),
        legend.title = element_text(size=rel(1.1)))+
  labs(fill="-Log10 (FDR)") ; dev.off()





######  5) LEfSe results analysis   ######

# LEfSe analysis was performed in web using input data 
# that prepared in 01_data_prep.R script & excel
# In excel, a little modification such as transposition and conversion of csv files to tsv was manually applied. 

# Four different input data files for LEfSe was separately analyzed at
# https://huttenhower.sph.harvard.edu/galaxy/

# Among 4 results, all significant biomarkers were combined in one result file (LEfSe_results.csv).


lef <- read.csv("LEfSe/lefse_results_modified_name.csv") ; head(lef)

lef.both <- lef ; lef.both <- lef.both %>% 
  dplyr::filter(Multiple=="Yes") %>%
  select(microbe) %>% 
  duplicated() %>% 
  which() %>% 
  c(., .-1) %>%
  lef.both[.,] %>%
  arrange(microbe) ; head(lef.both)

lef.both$LDA[lef.both$Age=="Old"]<- 
  with(lef.both, LDA[Age=="Old"]*-1)

ggplot(lef.both, aes(microbe, LDA, fill=Sex)) +
  xlab("") + ylab("Sum of LDA effect size") +
  geom_bar(stat = "identity", position = "stack") +
  theme_bw()+
  theme(axis.text = element_text(size=rel(1.1), angle = 90, 
                                 vjust=0.3, hjust=1, color="black"),
        legend.position = "top")


lef.both$Condition <- paste(lef.both$Project, lef.both$Sex, sep = "_")

lef$Condition <- paste(lef$Project, lef$Sex, sep = "_")

lef$LDA[lef$Age=="Old"] <- lef$LDA[lef$Age=="Old"]*-1

lef.all <- reshape::melt(lef, 
                         id=c("microbe", "Condition", "Age", "Sex"),
                         measure="LDA")
lef.list <- lef.all %>%
  reshape::cast(data = ., ~microbe, sum) %>% t() 

lef.list <- data.frame(microbe=attributes(lef.list)$dimnames[[1]],
           LDA.sum=lef.list[,1]) %>% arrange(desc(LDA.sum))

lef.list <- subset(lef.list, abs(LDA.sum)>=3)
lef.all <- lef.all[lef.all$microbe %in% lef.list$microbe,]

lef.all$microbe <- factor(lef.all$microbe, levels = lef.list$microbe)

tiff("figures/LEfSe_6condition_summary.tiff", res = 160,
     width = 1800, height = 1300, units = "px"); ggplot(lef.all, aes(microbe, value, fill=Condition)) +
  xlab("") + ylab("Cumulative effect size") +
  geom_bar(stat = "identity", position = "stack") +
  theme_bw()+
  theme(axis.text.x = element_text(size=rel(1.6), 
                                 angle = 90, vjust=0.3, hjust=1, color="black"),
        legend.position = "top",
        legend.title = element_blank(),
        axis.text.y = element_text(color="black"),
        legend.text = element_text(size=rel(1.5)),
        axis.title.y = element_text(size = rel(1.5)));dev.off()

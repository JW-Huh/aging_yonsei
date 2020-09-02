

# This script is written for preparing input data file
# for visualization of human gut microbiome data
# from American Gut Project (AGP) 
# & Guangdong Gut Microbiome Project (GGMP)

# BIOM file for both AGP and GGMP was downloaded from QIITA database (https://qiita.ucsd.edu/)
  # 1) Accession number of AGP was 10317 (https://qiita.ucsd.edu/study/description/10317)
    # Download link for AGP BIOM files: https://qiita.ucsd.edu/download_study_bioms/10317
  # 2) Accession number of GGMP was 11757 (https://qiita.ucsd.edu/study/description/11757)
    # Download link for AGP BIOM files: https://qiita.ucsd.edu/download_study_bioms/11757



#  Requirement 
rm(list=ls())
options(java.parameters = "-Xmx64g", stringsAsFactors = F) 
# You need to set a proper working directory on your own.
# setwd(dir = "D:/1_Data/2020/0530_Aging_combined/")

library(dplyr) ; library(biomformat); library(vegan) #; library(readr)

memory.size(max=TRUE) ; memory.size(max=FALSE) ; memory.limit(size=NA) ; memory.limit(size=56000) ; memory.limit(size=NA) 



######## 1) American Gut Project  ###########

# Loadding Biom data ; Genus levels
# Coerce data into a matrix and then into a data.frame
# to Extract microbial abundance matrix
# 2082 microbes x 15158 samples

otu <- read_biom(biom_file = "D:/1_Data/2020/0511_Aging_AGP/03-otus/03-otus/100nt/gg-13_8-97-percent/otu_table.biom")
  # otu$generated_by # QIIME 1.9.1
  #otu$date # 2017-12-05
  #otu$shape # 36405 sequence variants x 19524 samples
  

# 1) Extract unique taxanomy data
taxa_agp <- otu$rows %>%
  unlist() %>%
  matrix(data = ., 
         nrow = length(otu$rows),
         ncol = length(otu$rows[[1]]$metadata$taxonomy)+1,
         byrow = TRUE) %>%
  as.data.frame() %>%
  set_colnames(c("ID", "Kingdom", "Phylum", "Class", 
                 "Order", "Family","Genus", "Species")) %>%
  mutate(paste=paste0(Phylum, ";", Class, ";", Order, ";",
                      Family, ";", Genus, ";", Species)) 


# Unique taxa
taxa_agp_2768 <- taxa_agp %>%
  distinct(paste, .keep_all = TRUE) %>%
  mutate(fg.name=paste0(Family, ";", Genus)) ; taxa_agp_unique <-
  taxa_agp_2768$paste

# 2) Metadata 
  # 15274 samples x 528 columns (fecal data)
  # readr::guess_encoding("D:/1_Data/2020/0511_Aging_AGP/04-meta/04-meta/ag-fecal-cleaned.csv")

meta <- read.csv("D:/1_Data/2020/0511_Aging_AGP/04-meta/04-meta/ag-fecal-cleaned.csv", header = T, fileEncoding = "euc-kr")


## Choose samples that contained information of age and sex
  # because not all samples have those
meta_agp <- meta %>%
  filter(!AGE_CORRECTED %in% c("Unknown", "Unspecified", "no_data") &
           SEX %in% c("female", "male")) %>%
  mutate(Age=as.numeric(AGE_CORRECTED)) %>%
  dplyr::select(-AGE_CORRECTED) %>% # 13839
  filter(!BMI_CORRECTED %in% 
           c("Unknown", "Unspecified", "no_data")) %>%
  mutate(BMI=as.numeric(BMI_CORRECTED)) %>%
  dplyr::select(-BMI_CORRECTED) # 7408

# Healthy, non-pregnant, no recent Abx, having info. of diet, Caucasian
# Race: Caucasian (88%)

table(meta_agp$IBS)
table(meta_agp$IBD)
table(meta_agp$CANCER)
table(meta_agp$DIABETES)
table(meta_agp$ANTIBIOTIC_HISTORY)

# Strict filtering
  # Criteria
    # Age, BMI, diseases (IBD, IBS, Cancer, diabetes), diet
    # Not pregnant, no antibiotic history within 6-months
    # Caucasian to ameliorate the effect of race

meta_agp <- meta_agp %>% 
  filter(Age >= 20 & Age <=70 &
           BMI>=18.5 & BMI<=30 & # 7159
           IBD=="I do not have this condition" & 
           IBS=="I do not have this condition" & 
           CANCER=="I do not have this condition" & 
           DIABETES=="I do not have this condition" & 
           DIET_TYPE %in% 
           c("Omnivore", "Omnivore but do not eat red meat", 
             "Vegan", "Vegetarian", "Vegetarian but eat seafood") &
           PREGNANT != "Yes" &
           ANTIBIOTIC_HISTORY %in% 
           c("Year", "6 months",
             "I have not taken antibiotics in the past year.") &
           RACE=="Caucasian") %>%
  mutate(Age_5class=case_when(Age<30 ~ "20s",
                              Age<40 & Age>=30 ~ "30s",
                              Age<50 & Age>=40 ~ "40s",
                              Age<60 & Age>=50 ~ "50s",
                              Age>=60 ~ "60+")) # 2654 survived

## Soft filtering
# meta_agp <- meta_agp %>% 
#   filter(Age >= 20 & Age <=70 &
#            BMI>=15 & BMI<=40 & # 7159
#            !IBD %like% "iagnosed" & # Remove Diagnosed or diagnosed 
#            !IBS %like% "iagnosed"  &
#            !CANCER %like% "iagnosed" &
#            !DIABETES %like% "iagnosed" &
#            DIET_TYPE %in% 
#            c("Omnivore", "Omnivore but do not eat red meat", 
#              "Vegan", "Vegetarian", "Vegetarian but eat seafood") &
#            PREGNANT != "Yes" &
#            ANTIBIOTIC_HISTORY != "Week" &
#            RACE=="Caucasian") # 3665 


# 3) OTU table
otu_tmp <- otu$data # length(otu_tmp) = 36405 ; number of OTUs

# Samples (row) x OTUs (column)
otu_agp <- 
  matrix(nrow = 19524, # otu$shape[2] ; number of samples
         ncol = 0) %>% 
  as.data.frame() ; for(i in 1:36405) { 

  tmp_agp <- as.data.frame(otu_tmp[[i]])
  otu_agp <- cbind(otu_agp, tmp_agp)

} ; colnames(otu_agp) <- taxa_agp$ID ; rm(i, tmp_agp, otu_tmp)


# Shared ID with metadata 
id_agp_intersect <- sort(intersect(meta_agp$SampleID,                  
                                   rownames(otu_agp))) # 2651 

# Filtering 
  # OTU matrix ; 2651 samples ordered x 36405 OTUs
  # metadata; 2651 sample sordered
  otu_agp <- otu_agp[id_agp_intersect,] 
  meta_agp <- meta_agp %>% 
    filter(SampleID %in% id_agp_intersect) %>%
    arrange(SampleID)



# Collapsing to 2768 unique species
otu_agp_2768 <- matrix(nrow = nrow(meta_agp),   # Number of samples 
                       ncol = 2768) %>%  # Number of unique OTUs
  as.data.frame() %>%
  set_rownames(meta_agp$SampleID) %>%
  set_colnames(taxa_agp_unique) ; for(i in 1:2768){
    
  # Loading unique bugs taxon, one-by-one (total 2768 taxa)
  tmp.taxa <- taxa_agp_unique[i]
    
  # Choose taxon ID(s) which share same unique taxon
  tmp.ID <- taxa_agp %>%
    filter(paste==tmp.taxa) %>%
    dplyr::select(ID) %>% 
    unlist()
    
  # Sum read counts of those taxon IDs 
  otu_agp_2768[,i] <- otu_agp[,tmp.ID] %>% 
    as.data.frame() %>%
    rowSums() 
      
  }

# # Saving point_1) After metadata filtering
# write.csv(otu_agp, "input/OTU_AGP_2651x36405_asv.csv")
# write.csv(otu_agp_2768, "input/OTU_AGP_2651x2768_taxa.csv")
# write.csv(meta_agp, "input/meta_AGP_2651.csv")
# write.csv(taxa_agp, "input/taxa_AGP_36405_asv.csv")
# write.csv(taxa_agp_2768, "input/taxa_AGP_2768_taxa.csv")
# rm(otu, meta)


# Sum of read = Sequencing depth
depth_agp <- data.frame(sum=rowSums(otu_agp_2768))
  
tiff("figures/AGP_depth.tiff", res = 150, units = "px",
     width = 400, height = 350);ggplot(depth_agp,
                                       aes(log10(sum))) +
  xlab("Log10 (Sequencing depth)") +
  geom_histogram(bins=100) +
  theme_bw(); dev.off()
  
identical(rownames(otu_agp_2768), rownames(otu_agp)) # TRUE
identical(rowSums(otu_agp_2768), rowSums(otu_agp)) # TRUE ; sample depth

rm(i, tmp.ID, tmp.taxa)



######## 2) Chinese dataset  ###########

#rm(llist=ls())
memory.size(max=TRUE) ; memory.size(max=FALSE) ; memory.limit(size=NA) ; memory.limit(size=56000) ; memory.limit(size=NA) 


# 1) Metadata: 7009 fecal samples 
  # Mannaully converted to csv file 
  # "82680_mapping_file.txt "

meta <- read.csv("D:/1_Data/2020/0511_Aging_China/meta_China.csv",
                     header = T) %>% dplyr::rename(Age=age, Sex=sex)

meta_chi <- meta %>%
  arrange(sample_name) %>%
  filter(colitis!="y" & 
           t1dm!="y" & t2dm!="y" & 
           antibiotics!="y" &  # yes 6392 or no 449
           ibs!="y" &
           Age>=20 & Age<=70 &
           bmi>=18.5 & bmi<=30) %>%
  mutate(Age_5class=case_when(Age<30 ~ "20s",
                              Age<40 & Age>=30 ~ "30s",
                              Age<50 & Age>=40 ~ "40s",
                              Age<60 & Age>=50 ~ "50s",
                              Age>=60 ~ "60+")) # 4548

## Bristol stool form scale 
meta_chi_stool <- 
  meta_chi %>%
  filter(bristol_stool_type!="NA") 
  
# meta_chi_stool$bristol_stool_type <- 
#   as.factor(meta_chi_stool$bristol_stool_type)

is.numeric(meta_chi_stool$bristol_stool_type) # TRUE


tiff("figures/China_stool_type_age.tiff", res=160, units = "px",
     width = 450, height = 600) ; ggplot(meta_chi_stool, 
                                         aes(Age, bristol_stool_type,
                                             color=Sex)) +
  xlab("") + ylab("")+
  geom_smooth(method="loess", span=0.8, alpha=0.15)+
  theme_classic() +
  theme(axis.text.x = element_text(size=rel(1.6), color="black"),
        axis.text.y = element_text(size=rel(1.2), color="black"),
        legend.position = "none");dev.off()

# tiff("figures/China_stool_type_age.tiff", res=160, units = "px",
#      width = 450, height = 600) ; ggplot(meta_chi_stool, 
#                                          aes(Age_5class))+
#   xlab("") + ylab("") +
#   geom_bar(stat = "count",
#            position = "fill",
#            aes(fill=bristol_stool_type)) +
#   theme_classic() +
#   theme(legend.position = "none",
#         axis.text.x = element_text(size=rel(1.6), color="black"),
#         axis.text.y = element_text(size=rel(1.5), color="black"));dev.off()

rm(meta_chi_stool)

# The most recent version was used (82680)
  # QIIME 1.9.1
  # 2019-11-21T10:32:32.032421

# 2) OTU table
# otu <- read_biom("D:/1_Data/2020/0511_Aging_China/study_11757_042620-183355/BIOM/82680/otu_table.biom")

# Extract OTU abundance matrix
otu_chi <- as.data.frame(as.matrix(biom_data(otu)))
otu_chi <-
  otu_chi[sort(rownames(otu_chi)), 
          sort(colnames(otu_chi))] %>%
  dplyr::select(meta_chi$sample_name)

identical(colnames(otu_chi), meta_chi$sample_name) # TRUE
otu_chi <- as.data.frame(t(otu_chi)) ; dim(otu_chi) # 4548 x 16070


# Shared ID with metadata 
# Every samples in metadata are also present in abundance matrix
  # TRUE ; 4548 == 4548
nrow(meta_chi)==
  length(intersect(meta_chi$sample_name, colnames(otu_chi))) 


# 3) Extract unique taxanomy data
taxa_chi <- otu$rows %>%
  unlist() %>%
  matrix(data = ., 
         nrow = length(otu$rows),
         ncol = length(otu$rows[[1]]$metadata$taxonomy)+1,
         byrow = TRUE) %>%
  as.data.frame() %>%
  set_colnames(c("ID", "Kingdom", "Phylum", "Class", 
                 "Order", "Family","Genus", "Species")) %>%
  mutate(paste=paste0(Phylum, ";", Class, ";", Order, ";",
                      Family, ";", Genus, ";", Species)) 

# Unique taxa
taxa_chi_1623 <- taxa_chi %>%
  distinct(paste, .keep_all = TRUE) %>%
  mutate(fg.name=paste0(Family, ";", Genus)) ; taxa_chi_unique <-
  taxa_chi_1623$paste

# Collapsing 16070 taxa into 1623 unique taxa
# 4548 samples x 1623 microbes
otu_chi_1623 <- matrix(nrow = nrow(meta_chi),    # Number of samples 
                       ncol=1623) %>%  # Number of unique OTUs
  as.data.frame() %>%
  set_rownames(meta_chi$sample_name) %>%
  set_colnames(taxa_chi_unique) ; for(i in 1:1623){
    
    # Loading unique bugs taxon, one-by-one (total 1623 taxa)    
    tmp.taxa <- taxa_chi_unique[i]
    
    # Choose taxon ID(s) which share same unique taxon
    tmp.ID <- taxa_chi %>%
      filter(paste==tmp.taxa) %>%
      dplyr::select(ID) %>% 
      unlist()
    
    # Sum read counts of those taxon IDs 
    otu_chi_1623[,i] <- otu_chi[,tmp.ID] %>% 
      as.data.frame() %>%
      rowSums() 
    
}; rm(i, tmp.ID, tmp.taxa) 

# # Saving point_1) After metadata filtering
# write.csv(otu_chi, "input/OTU_GGMP_4548x16070_asv.csv")
# write.csv(otu_chi_1623, "input/OTU_GGMP_4548x1623_taxa.csv")
# write.csv(meta_chi, "input/meta_GGMP_4548.csv")
# write.csv(taxa_chi, "input/taxa_GGMP_16070_asv.csv")
# write.csv(taxa_chi_1623, "input/taxa_GGMP_1623_taxa.csv")
# rm(otu, meta)

taxa_chi <- read.csv("input/taxa_GGMP_16070_asv.csv", 
                     row.names = 1)
taxa_agp <- read.csv("input/taxa_GGMP_16070_asv.csv", 
                     row.names = 1)

taxa_chi$ID <- as.character(taxa_chi$ID)

taxa_chi <- taxa_chi[order(taxa_chi$ID),]
taxa_agp <- taxa_agp[order(taxa_agp$ID),]


taxa_chi$ID[taxa_chi$paste %like% "Akkermansia"]
taxa_chi$paste[taxa_chi$ID=="174337"]
taxa_agp$paste[taxa_agp$ID=="174337"]


# Depth
depth_chi <- data.frame(sum=rowSums(otu_chi_1623))

tiff("figures/China_depth.tiff",
     res = 150, units = "px",
     width = 400, height = 350);ggplot(depth_chi,
                                       aes(log10(sum))) +
  xlab("Log10 (Sequencing depth)")+
  geom_histogram(bins=100) +
  theme_bw();dev.off()

identical(rownames(otu_chi_1623), rownames(otu_chi)) # TRUE  
identical(rowSums(otu_chi_1623), rowSums(otu_chi)) # TRUE 
# rm(otu, otu_chi, depth_chi, taxa_chi)




#############    Prevalence filtering    #############

# Remove undetected species
  # 1426 ; not detected species
pre_agp <- 
  data.frame(pre = colSums(otu_agp_2768>0)/
               nrow(otu_agp_2768)*100) ; sum(pre_agp==0) 

  # 107 ; not detected species
pre_chi <- 
  data.frame(pre = colSums(otu_chi_1623>0)/
               nrow(otu_chi_1623)*100) ; sum(pre_chi==0) 

# tiff("figures/prevalence_agp.tiff", units = "px", res=160,
#      width = 600, height = 450);ggplot(pre_agp, aes(pre)) +
#   xlab("Prevalence")+ geom_histogram(bins=100) +
#   scale_y_continuous(limits=c(0, 2200)) + theme_bw(); dev.off()
# 
# tiff("figures/prevalence_chi.tiff", units = "px", res=160,
#      width = 600, height = 450);ggplot(pre_chi, aes(pre)) +
#   xlab("Prevalence")+ geom_histogram(bins=100) +
#   scale_y_continuous(limits=c(0, 2200)) + theme_bw(); dev.off()

pre_agp <- rownames(pre_agp)[pre_agp>0] # 1342 
pre_chi <-rownames(pre_chi)[pre_chi>0] # 1516

# Bugs abundance filtering
otu_agp_prev <- otu_agp_2768[,pre_agp]
otu_chi_prev <- otu_chi_1623[,pre_chi]


# taxanomy filtering
taxa_agp_prev <- taxa_agp_2768 %>% filter(paste %in% pre_agp)
taxa_chi_prev <- taxa_chi_1623 %>% filter(paste %in% pre_chi)
  #rm(otu_agp_2768, otu_chi_1623, taxa_agp_2768, taxa_chi_1623)


# # Saving point_2) After collapsing & prevalence filtering 

# write.csv(otu_agp_prev, "input/OTU_AGP_2651x1342_taxa.csv")
# write.csv(taxa_agp_prev, "input/taxa_AGP_1342.csv")
# 
# write.csv(otu_chi_prev, "input/OTU_GGMP_4548x1516_taxa.csv")
# write.csv(taxa_chi_prev, "input/taxa_GGMP_1516_taxa.csv")




#############    Rarefying    #############
  # 90% of lower samples depth 
  # In this case, 2178 (shown at line number 501)

# # Import
# 
# # 1) AGP
# otu_agp_prev <- 
#   read.csv("input/OTU_AGP_2651x1342_taxa.csv", row.names = 1)
# taxa_agp_prev <- 
#   read.csv("input/taxa_AGP_1342.csv", row.names = 1)
# meta_agp <- read.csv("input/meta_AGP_2651.csv", row.names = 1)
# 
# # 2) China
# otu_chi_prev <- 
#   read.csv("input/OTU_GGMP_4548x1516_taxa.csv", row.names = 1)
# taxa_chi_prev <- 
#   read.csv("input/taxa_GGMP_1516_taxa.csv", row.names = 1)
# meta_chi <- read.csv("input/meta_GGMP_4548.csv", row.names = 1)

  
# 1) First, filtering by depth
# 12,750 sorted almost equal number of samples like below

  # depth_agp <- data.frame(sum=rowSums(otu_agp_prev))
  # depth_chi <- data.frame(sum=rowSums(otu_chi_prev))
  
  depth.thres = 5000
  sum(depth_agp>=depth.thres) # 2420
  sum(depth_chi>=depth.thres) # 4547
  
  
# 1-1) Samples with sequencing depth higher than 5000
depth_agp <- rownames(depth_agp)[depth_agp$sum>=depth.thres] 
depth_chi <- rownames(depth_chi)[depth_chi$sum>=depth.thres]


# 1-2) OTU matrix filtering
otu_agp_dep <- otu_agp_prev[depth_agp,]
otu_chi_dep <- otu_chi_prev[depth_chi,]
rm(otu_agp_prev, otu_chi_prev, otu_agp_2768, otu_chi_1623)


# 1-3) Metadata filtering
meta_agp <- meta_agp %>% filter(SampleID %in% depth_agp)
meta_chi <- meta_chi %>% filter(sample_name %in% depth_chi)


# tiff("figures/AGP_depth_afterFiltering.tiff",
#      res = 150, units = "px",
#      width = 400, height = 350); otu_agp_dep %>%
#   rowSums() %>% as.data.frame() %>%
#   set_colnames("Sum") %>%
#   ggplot() + geom_histogram(aes(log10(Sum))) + theme_bw() ; dev.off()
# 
# tiff("figures/China_depth_afterFiltering.tiff",
#      res = 150, units = "px",
#      width = 400, height = 350); otu_chi_dep %>%
#   rowSums() %>% as.data.frame() %>%
#   set_colnames("Sum") %>%
#   ggplot() + geom_histogram(aes(log10(Sum))) + theme_bw() ; dev.off()


# 1-4) Order check ; All TRUE
identical(rownames(otu_agp_dep), meta_agp$SampleID)   # TRUE
identical(rownames(otu_chi_dep), meta_chi$sample_name) # TRUE

# 1-5) Filtering 0 abundnace OTUs
otu_agp_dep <- otu_agp_dep[,colSums(otu_agp_dep)>0] ; dim(otu_agp_dep) 
otu_chi_dep <- otu_chi_dep[,colSums(otu_chi_dep)>0] ; dim(otu_chi_dep) 

# 2420 samples x 1326 OTUs ; AGP (Americans)
# 4547 samples x 1516 OTUs ; GGMP (Chinese)


# 2) Rarefying
rarefy.thres <- round(nrow(otu_agp_dep) * 0.9) # 2178

set.seed(123) ; otu_agp_rare <- 
  as.data.frame(rrarefy(otu_agp_dep, sample = rarefy.thres))

set.seed(123) ; otu_chi_rare <- 
  as.data.frame(rrarefy(otu_chi_dep, sample = rarefy.thres))

otu_agp_rare <- otu_agp_rare[,colSums(otu_agp_rare)>0] # 938 OTUs
otu_chi_rare <- otu_chi_rare[,colSums(otu_chi_rare)>0] # 1059 OTUs

  # Rarefying for OTU abundance matrix of ASV
  otu_agp_asv <- otu_agp[rownames(otu_agp_rare),] # 36405 OTUs
  otu_agp_asv <- otu_agp_asv[,colSums(otu_agp_asv)>0] # 12576 OTUs
  
  otu_chi_asv <- otu_chi[rownames(otu_chi_rare),] # 16070 OTUs
  otu_chi_asv <- otu_chi_asv[,colSums(otu_chi_asv)>0] # 14454 OTUs
  
  set.seed(123) ; otu_agp_asv <- 
    as.data.frame(rrarefy(otu_agp_asv, sample = rarefy.thres))
  otu_agp_asv <- otu_agp_asv[,colSums(otu_agp_asv)>0] # 8197 OTUs
  
  set.seed(123) ; otu_chi_asv <- 
    as.data.frame(rrarefy(otu_chi_asv, sample = rarefy.thres))
  otu_chi_asv <- otu_chi_asv[,colSums(otu_chi_asv)>0] # 8349 OTUs
  
  dim(otu_agp_asv) # 2420 samples x 8197 SVs
  dim(otu_chi_asv) # 4547 samples x 8349 SVs
  
# Taxon filtered
taxa_agp_rare <- taxa_agp_prev %>% 
  filter(paste %in% colnames(otu_agp_rare))

taxa_chi_rare <- taxa_chi_prev %>% 
  filter(paste %in% colnames(otu_chi_rare))


dim(otu_agp_rare)
dim(otu_chi_rare)

# # Saving point_2) After rarefying
# write.csv(otu_agp_rare, "input/OTU_AGP_rarefied_2420x938.csv")
# write.csv(taxa_agp_rare, "input/taxa_AGP_rarefied_938.csv")
# write.csv(meta_agp, "input/meta_AGP_2420.csv")
# 
# write.csv(otu_chi_rare, "input/OTU_GGMP_rarefied_4547x1059.csv")
# write.csv(taxa_chi_rare, "input/taxa_GGMP_rarefied_1059.csv")
# write.csv(meta_chi, "input/meta_GGMP_4547.csv")


# Import 
otu_agp_rare <- 
  read.csv("input/OTU_AGP_rarefied_2420x938.csv", row.names = 1)
otu_chi_rare <- 
  read.csv("input/OTU_GGMP_rarefied_4547x1059.csv", row.names = 1)

meta_agp <- read.csv("input/meta_AGP_2420.csv", row.names = 1)
meta_chi <- read.csv("input/meta_GGMP_4547.csv", row.names = 1)

taxa_agp_rare <- 
  read.csv("input/taxa_AGP_rarefied_938.csv", row.names = 1)
taxa_chi_rare <- 
  read.csv("input/taxa_GGMP_rarefied_1059.csv", row.names = 1)

# # Continuous age
# tiff("figures/AGP_distribution_sex_age_filtered.tiff",
#      res = 150, units = "px",
#      width = 600, height = 400); ggplot(meta_agp,
#                                         aes(Age, fill=SEX)) +
#   xlab("") + ylab("")+
#   geom_histogram(stat = "count") +
#   theme_classic() +
#   theme(legend.title = element_blank(),
#         legend.position = "none",
#         legend.key.size = unit(0.6, "line"),
#         axis.text = element_text(size=rel(1.2),color="black"),
#         axis.title = element_text(size=rel(1.5))); dev.off()
# 
# tiff("figures/GGMP_distribution_sex_age_filtered.tiff",
#      res = 150, units = "px",
#      width = 600, height = 400); ggplot(meta_chi,
#                                         aes(Age, fill=Sex)) +
#   xlab("") + ylab("")+
#   geom_histogram(stat = "count") +
#   theme_classic() +
#   theme(legend.title = element_blank(),
#         legend.position = "none",
#         legend.key.size = unit(0.6, "line"),
#         axis.text = element_text(size=rel(1.2),color="black"),
#         axis.title = element_text(size=rel(1.5))); dev.off()
# 
# # Categorized age
# tiff("figures/AGP_distribution_sex_categorized_age_filtered.tiff",
#      res = 150, units = "px",
#      width = 400, height = 400); ggplot(meta_agp,
#                                         aes(Age_5class, fill=SEX)) +
#   xlab("") + ylab("")+
#   geom_bar(stat = "count") +
#   theme_classic() +
#   theme(legend.title = element_blank(),
#         legend.position = "none",
#         legend.key.size = unit(0.6, "line"),
#         axis.text = element_text(size=rel(1.2),color="black"),
#         axis.title = element_text(size=rel(1.5))); dev.off()
# 
# 
# tiff("figures/GGMP_distribution_sex_categorized_age_filtered.tiff",
#      res = 150, units = "px",
#      width = 400, height = 400); ggplot(meta_chi,
#                                         aes(Age_5class, fill=Sex)) +
#   xlab("") + ylab("")+
#   geom_bar(stat = "count") +
#   theme_classic() +
#   theme(legend.title = element_blank(),
#         legend.position = "none",
#         legend.key.size = unit(0.6, "line"),
#         axis.text = element_text(size=rel(1.2),color="black"),
#         axis.title = element_text(size=rel(1.5))); dev.off()



#############    Collapsing    #############
  # To change the resolution of data
  # I collapsed microbial read count to specific taxonomic levels
  # 0) Extract the name of taxa
  # 1) Phylum: line number 649
  # 2) Family: line number 731
  # 3) Genus + genus-unassigned family: line nummber 791


# 0) Name of taxon----------

# 0-1) AGP
agp_name <- 
  strsplit(colnames(otu_agp_rare), ";") %>%
  unlist() %>%
  matrix(data=., ncol=6, byrow = T) %>%
  as.data.frame() 


# 0-2) China
chi_name <- 
  strsplit(colnames(otu_chi_rare), ";") %>%
  unlist() %>%
  matrix(data=., ncol=6, byrow = T) %>%
  as.data.frame() 


# 1) Phylum--------------

# 1-1) AGP

# Unique phylum
agp_p_all <- taxa_agp_rare %>% 
  distinct(Phylum) %>% filter(Phylum!="p__") %>% unlist()
  
# All phylum of 1485 
agp_p_name <- agp_name %>% dplyr::select("V1") %>% unlist()

# Collapsing into unique phylum
p_agp <- 
  as.data.frame(matrix(nrow = nrow(meta_agp), 
                       ncol = 0)) ; for(i in agp_p_all){

  p_tmp <- 
    otu_agp_rare[,agp_p_name==i] %>%
    as.data.frame() %>%
    rowSums() %>%
    as.data.frame()
    
  p_agp <- cbind(p_agp, p_tmp)
    
} ; colnames(p_agp) <- agp_p_all %>%
  gsub(pattern = "p__", replacement = "") ; head(p_agp)

identical(meta_agp$SampleID, rownames(p_agp)) # TRUE


# 1-2) GGMP

# Unique phylum
chi_p_all <- taxa_chi_rare %>% 
  distinct(Phylum) %>% filter(Phylum!="p__") %>% unlist()

# All phylum of 1231
chi_p_name <- chi_name %>% dplyr::select("V1") %>% unlist()

# Collapsing into unique phylum
p_chi <- 
  as.data.frame(matrix(nrow = nrow(meta_chi), 
                       ncol = 0)) ; for(i in chi_p_all){
                         
  p_tmp <- 
    otu_chi_rare[,chi_p_name==i] %>%
    as.data.frame() %>%
    rowSums() %>%
    as.data.frame()
                         
  p_chi <- cbind(p_chi, p_tmp)
                         
} ; colnames(p_chi) <- chi_p_all %>%
  gsub(pattern = "p__", replacement = "") ; head(p_chi)

identical(meta_chi$sample_name, rownames(p_chi)) # TRUE


# Top 9 phyla were retrieved in both dataset
p_mean1 <- apply(p_chi, 2, mean) %>% sort(decreasing = T)
p_mean2 <- apply(p_agp, 2, mean) %>% sort(decreasing = T)
  sort(names(p_mean1)[1:9])==sort(names(p_mean2)[1:9]) # All TRUE
  phylum.names <- names(p_mean1)[1:9]
  # Except Cyannobacteria
  phylum.names <- phylum.names[-which(phylum.names=="Cyanobacteria")]
  

# Top 8 + Other
p_agp <- data.frame(p_agp[,phylum.names],
                    Others=rowSums(p_agp[,!colnames(p_agp) %in% 
                                           phylum.names]))

p_chi <- data.frame(p_chi[,phylum.names],
                    Others=rowSums(p_chi[,!colnames(p_chi) %in% 
                                           phylum.names]))

# write.csv(p_agp, "input/phylum_AGP.csv")
# write.csv(p_chi, "input/phylum_GGMP.csv")
# rm(p_tmp, i, p_mean1, p_mean2, phylum.names)



# 2) Family--------------

# 2-1) AGP

# Unique family 
agp_fam <- taxa_agp_rare %>%
  distinct(Family) %>% filter(.!="f__") %>% unlist()


# All family in 2420 samples
agp_f_name <- agp_name %>%  dplyr::select("V4") %>%  unlist()

# Collapsing
f_agp <- 
  as.data.frame(matrix(nrow = nrow(meta_agp),
                       ncol = 0)) ; for(i in agp_fam){
    
                         
  f_tmp <- 
    otu_agp_rare[,agp_f_name==i] %>%
    as.data.frame() %>%
    rowSums() %>%
    as.data.frame()
                         
  f_agp <- cbind(f_agp, f_tmp)
                         
} ; colnames(f_agp) <- agp_fam ; f_agp$Others <- 4792-rowSums(f_agp)

# write.csv(f_agp, "input/Family_level_AGP.csv") ; rm(f_tmp, i)


# 2-2) China

# Unique family 
chi_fam <- taxa_chi_rare %>%
  distinct(Family) %>% filter(.!="f__") %>% unlist()


# All family in 4547 samples
chi_f_name <- chi_name %>% dplyr::select("V4") %>% unlist()

# Collapsing
f_chi <- 
  as.data.frame(matrix(nrow = nrow(meta_chi),
                       ncol = 0)) ; for(i in chi_fam){

  f_tmp <-
    otu_chi_rare[,chi_f_name==i] %>%
    as.data.frame() %>%
    rowSums() %>%
    as.data.frame()

  f_chi <- cbind(f_chi, f_tmp)

} ; colnames(f_chi) <- chi_fam ; f_chi$Others <- 4792-rowSums(f_chi)

# write.csv(f_chi, "input/Family_level_GGMP.csv") ; rm(f_tmp, i)



# 3) Genus + genus-unassigned Family--------------

# 3-1) AGP

# At least family levels  + Genus collapsed

agp_min_fam <- taxa_agp_rare %>% 
  filter(Family != "f__") %>% 
  dplyr::select("paste") %>%
  unlist() ; otu_agp_min <- otu_agp_rare[,agp_min_fam]  

length(agp_min_fam) # 874
dim(otu_agp_min) # 2420 samples x 874 OTUs 

# Rename OTU matrix
agp_min_fam_tag <- agp_min_fam %>%
  gsub(pattern = "\\[", replacement = ".") %>%
  gsub(pattern = "\\]", replacement = ".") ; colnames(otu_agp_min) <-
  agp_min_fam_tag

# paste Family + Genus name
agp_fg_name <- 
  paste0(taxa_agp_rare$Family, ";",
         taxa_agp_rare$Genus) %>% 
  gsub(pattern = "\\[", replacement = ".") %>%
  gsub(pattern = "\\]", replacement = ".") %>%
  as.data.frame() %>%
  unique() %>%
  filter(.!="f__;g__") %>% unlist() # 630 OTUs 

length(agp_fg_name) # 630

# Get genus collapsed + family, abundance matrix

otu_agp_fg  <- 
  as.data.frame(matrix(nrow=nrow(meta_agp), 
                       ncol=0)); for (i in 1:length(agp_fg_name)) {
                         
  tmp_genus <- agp_fg_name[[i]]
  
  if(tmp_genus %like% "g__$") { # Family, genus unassigned
             
  otu_tmp <- otu_agp_min %>%  
    dplyr::select(ends_with(paste0(tmp_genus, ";s__")))  %>%
    set_colnames(tmp_genus)
  
  } else { # Faminly + genus assigned
   
  otu_tmp <- otu_agp_min %>%  
    dplyr::select(contains(tmp_genus)) %>%
    rowSums() %>%
    as.data.frame() %>%
    set_colnames(tmp_genus)
    
  }
                         
  otu_agp_fg <- cbind(otu_agp_fg, otu_tmp)
                         
} ; otu_agp_fg <- otu_agp_fg %>%
  mutate(Others=(rarefy.thres-rowSums(.))) ; rownames(otu_agp_fg) <-
  rownames(otu_agp_min)

# Filtering by Family + Genus explanation levels
otu_agp_fg <- otu_agp_fg[otu_agp_fg$Others/rarefy.thres<0.1,] 
dim(otu_agp_fg) # 1652 samples x 631 OTUs (630 unique Fam+Gen + Others)

# write.csv(otu_agp_fg, "input/Family_Genus_AGP_1652x631_taxa.csv")  



# 3-2) China

chi_min_fam <- taxa_chi_rare %>% 
  filter(Family != "f__") %>% 
  dplyr::select("paste") %>%
  unlist() ; otu_chi_min <- otu_chi_rare[,chi_min_fam]  # 970 OTUs

dim(otu_chi_min) # 4547 samples x 970 OTUs

# Rename otu matrix
chi_min_fam_tag <- chi_min_fam %>%
  gsub(pattern = "\\[", replacement = ".") %>%
  gsub(pattern = "\\]", replacement = ".") ; colnames(otu_chi_min) <-
  chi_min_fam_tag


# paste Family + Genus name
chi_fg_name <- 
  paste0(taxa_chi_rare$Family, ";",
         taxa_chi_rare$Genus) %>% 
  gsub(pattern = "\\[", replacement = ".") %>%
  gsub(pattern = "\\]", replacement = ".") %>%
  as.data.frame() %>%
  unique() %>%
  filter(.!="f__;g__") %>% unlist() # 719 OTUs 


# Get genus collapsed + family, abundance matrix
  
otu_chi_fg  <- 
  as.data.frame(matrix(nrow=nrow(meta_chi), 
                       ncol=0)); for (i in 1:length(chi_fg_name)) {
                         
  tmp_genus <- chi_fg_name[[i]]
   
  if(tmp_genus %like% "g__$") { # Family, genus unassigned
                           
    otu_tmp <- otu_chi_min %>%  
      dplyr::select(ends_with(paste0(tmp_genus, ";s__")))  %>%
      set_colnames(tmp_genus)
                           
  } else { # Faminly + genus assigned
  
    otu_tmp <- otu_chi_min %>%  
      dplyr::select(contains(tmp_genus)) %>%
      rowSums() %>%
      as.data.frame() %>%
      set_colnames(tmp_genus)

    }
    
  otu_chi_fg <- cbind(otu_chi_fg, otu_tmp)
  
} ; otu_chi_fg <- otu_chi_fg %>%
  mutate(Others=(rarefy.thres-rowSums(.))) ; rownames(otu_chi_fg) <-
  rownames(otu_chi_min)


# Filtering by Family + Genus explanation levels
otu_chi_fg <- otu_chi_fg[otu_chi_fg$Others/rarefy.thres<=0.1,]
dim(otu_chi_fg) # 3743 samples x 720 OTUs (719 + Others)

# write.csv(otu_chi_fg, "input/Family_Genus_GGMP_3743x720_taxa")  




###### Project combine (Family + Genus) ######

otu_agp_fg <- 
  read.csv("input/Family_Genus_AGP_1652x631_taxa.csv", row.names = 1)
otu_chi_fg <- 
  read.csv("input/Family_Genus_GGMP_3743x720_taxa.csv", row.names = 1)


# Intersect ;
fg_sect <- sort(intersect(colnames(otu_agp_fg), 
                          colnames(otu_chi_fg))) # 483 OTUs (+ Others)

# 148 OTUs ; AGP-specific 
otu_agp_fg %>% dplyr::select(!fg_sect) %>% dim() 

# 237 OTUs ; GGMP-specific
otu_chi_fg %>% dplyr::select(!fg_sect) %>% dim() 


# 1) Intersect microbes
otu_all <- rbind(otu_agp_fg[,fg_sect], 
                 otu_chi_fg[,fg_sect]) 

dim(otu_all) # 5395 x 483

#---------------#

# 2) GGMP samples + AGP-specific OTUs
  # AGP-specific OTU ; columns
  # The number of GGMP samples : rows

chi_tmp <- matrix(data = 0, 
                  nrow = nrow(otu_chi_fg), 
                  ncol = dim(otu_agp_fg %>% dplyr::select(!fg_sect))[2]) %>% as.data.frame() %>%
  set_colnames(colnames(otu_agp_fg %>% dplyr::select(!fg_sect))) ; rownames(chi_tmp) <- rownames(otu_chi_fg)

# 3) AGP samples + Chinese-specific OTUs
  # GGMP-specific OTU ; columns
  # The number of AGP samples : rows
agp_tmp <- matrix(data = 0, 
                  nrow = nrow(otu_agp_fg), 
                  ncol = dim(otu_chi_fg %>% dplyr::select(!fg_sect))[2]) %>% as.data.frame() %>%
  set_colnames(colnames(otu_chi_fg %>% dplyr::select(!fg_sect))) ; rownames(agp_tmp) <- rownames(otu_agp_fg)

#---------------#

# Merge_1) GGMP-specific microbes
fg_all_agp <- rbind(agp_tmp,
                    otu_chi_fg[,-which(colnames(otu_chi_fg) %in% 
                                         fg_sect)])
dim(fg_all_agp) # 5395 x 237

# Merge_2) AGP-specific microbes
fg_all_chi <- rbind(otu_agp_fg[,-which(colnames(otu_agp_fg) %in%
                                         fg_sect)],
                    chi_tmp)
dim(fg_all_chi) # 5395 x 148


# Merge all 
fg_all <- cbind(otu_all, fg_all_chi)
fg_all <- cbind(fg_all, fg_all_agp)

rm(fg_all_agp, fg_all_chi, otu_all, fg_sect, chi_tmp, agp_tmp)

# write.csv(fg_all, "input/Family_Genus_Both_5395x868_taxa.csv")



#### LEfSe input ####

meta_agp <- read.csv("input/meta_AGP_2420.csv", row.names = 1)
otu_agp_rare <- 
  read.csv("input/OTU_AGP_rarefied_2420x938.csv", row.names = 1)
taxa_agp_rare <- 
  read.csv("input/taxa_AGP_rarefied_938.csv", row.names = 1)

meta_chi <- read.csv("input/meta_GGMP_4547.csv", row.names = 1)
otu_chi_rare <- 
  read.csv("input/OTU_GGMP_rarefied_4547x1059.csv", row.names = 1)
taxa_chi_rare <- 
  read.csv("input/taxa_GGMP_rarefied_1059.csv", row.names = 1)

# [, ], ; were imported as "."
colnames(otu_agp_rare) <- taxa_agp_rare$paste 
colnames(otu_chi_rare) <- taxa_chi_rare$paste 


# 0) Name of taxon: requisite
# 0-1) AGP
agp_name <- strsplit(colnames(otu_agp_rare), ";") %>%
  unlist() %>% matrix(data=., ncol=6, byrow = T) %>% as.data.frame() 

# 0-2) China
chi_name <- strsplit(colnames(otu_chi_rare), ";") %>% unlist() %>%
  matrix(data=., ncol=6, byrow = T) %>% as.data.frame() 


# 1) AGP
# Unique family excluding "f__"
agp_lef <- taxa_agp_rare %>%
  distinct(Family) %>% filter(Family!="f__") %>% unlist()


# All family in all samples
agp_f_name <- agp_name %>% dplyr::select("V4") %>% unlist()

# Taxa that are not resolved to family level
agp_above_fam <- otu_agp_rare[,agp_f_name=="f__"] 


# Collapsing
lefse_agp <- 
  as.data.frame(matrix(nrow = nrow(meta_agp),
                       ncol = 0)) ; for(i in agp_lef){
                        
f_tmp <- 
  otu_agp_rare[,agp_f_name==i] %>%
  as.data.frame() %>% rowSums() %>% as.data.frame() %>%
  set_colnames(i)

lefse_agp <- cbind(lefse_agp, f_tmp)
 
} ; lefse_agp <- cbind(lefse_agp, agp_above_fam)


# metadata
identical(meta_agp$SampleID, rownames(lefse_agp)) # TRUE

meta_agp$Age_2 <- ifelse(meta_agp$Age>=45, "Old", "Young")
lefse_agp %>%
  mutate(Sex=meta_agp$SEX,
         Age=meta_agp$Age_2,
         ID=meta_agp$SampleID) %>% write.csv("LEfSe/lefse_AGP.csv")


# 2) GGMP
# Unique family excluding "f__"
chi_lef <- taxa_chi_rare %>%
  distinct(Family) %>% filter(Family!="f__") %>% unlist()


# All family in 7162 samples
chi_f_name <- chi_name %>% dplyr::select("V4") %>% unlist()

# Taxa that are not resolved to family level
chi_above_fam <- otu_chi_rare[,chi_f_name=="f__"] 


# Collapsing
lefse_chi <- 
  as.data.frame(matrix(nrow = nrow(meta_chi),
                       ncol = 0)) ; for(i in chi_lef){
                         
  f_tmp <- 
    otu_chi_rare[,chi_f_name==i] %>%
    as.data.frame() %>% rowSums() %>% as.data.frame() %>%
    set_colnames(i)
                         
  lefse_chi <- cbind(lefse_chi, f_tmp)

} ; lefse_chi <- cbind(lefse_chi, chi_above_fam)


# metadata
identical(meta_chi$sample_name, rownames(lefse_chi)) # TRUE

meta_chi$Age_2 <- ifelse(meta_chi$Age>=45, "Old", "Young")
lefse_chi %>%
  mutate(Sex=meta_chi$Sex,
         Age=meta_chi$Age_2,
         ID=meta_chi$sample_name) %>% write.csv("LEfSe/lefse_GGMP.csv")




##############    For PICRUST2    ##############

# Input file

  # I am going to use the samples that remained after rarefying
  # Sequencing depth filtered (5000)

# OTU abundance table & metadata
rarefy.thres <- 2178

# 1) AGP
# asv_agp <- read.csv("input/OTU_AGP_2651x36405_asv.csv", 
#                     row.names = 1, check.names = F)
# meta_agp <- read.csv("input/meta_AGP_2420.csv", row.names = 1)
# 
# asv_agp <- asv_agp %>%
#   filter(rownames(.) %in% meta_agp$SampleID) %>%
#   dplyr::select(colnames(.)[colSums(.)>0])
#   rownames(asv_agp) <- meta_agp$SampleID
#   write.csv(asv_agp, "input/OTU_AGP_2420x12576_asv.csv")

asv_agp <- read.csv("input/OTU_AGP_2420x12576_asv.csv",
                    row.names = 1, check.names = F) ; colnames(asv_agp)

set.seed(123) ; asv_agp <- 
  as.data.frame(rrarefy(asv_agp, sample = rarefy.thres)) 

asv_agp <- asv_agp %>%
  dplyr::select(colnames(.)[colSums(.)>0]) # 8197 SVs

agp_biom <- make_biom(t(asv_agp), 
                      sample_metadata = NULL,
                      observation_metadata = NULL,
                      id = "AGP",
                      matrix_element_type = "int") ; write_biom(agp_biom, "picrust2/AGP.biom")


# 2) GGMP
# asv_chi <- read.csv("input/OTU_GGMP_4548x16070_asv.csv", 
#                     row.names = 1, check.names = F)
# meta_chi <- read.csv("input/meta_GGMP_4547.csv", row.names = 1)
# 
# asv_chi <- asv_chi %>%
#   filter(rownames(.) %in% meta_chi$sample_name) %>%
#   dplyr::select(colnames(.)[colSums(.)>0])
#   rownames(asv_chi) <- meta_chi$sample_name
#   write.csv(asv_chi, "input/OTU_GGMP_4547x14454_asv.csv")

asv_chi <- read.csv("input/OTU_GGMP_4547x14454_asv.csv", 
                    row.names = 1, check.names = F)

set.seed(123) ; asv_chi <- 
  as.data.frame(rrarefy(asv_chi, sample = rarefy.thres)) 

asv_chi <- asv_chi %>%
  dplyr::select(colnames(.)[colSums(.)>0]) # 8349 SVs

chi_biom <- make_biom(t(asv_chi), 
                      sample_metadata = NULL,
                      observation_metadata = NULL,
                      id = "GGMP",
                      matrix_element_type = "int") ; write_biom(chi_biom, "picrust2/GGMP.biom")

# Representative sequence
# I will use greengene 13_5 97% OTU sequence as a representative sequence
library(seqinr)

fa <- read.fasta("input/97_otus.fasta")

agp_fa <- colnames(asv_agp) ; chi_fa <- colnames(asv_chi)
agp_fa <- fa[agp_fa] ; chi_fa <- fa[chi_fa]

write.fasta(names = attr(agp_fa, which = "names"),
            sequences = agp_fa,
            file.out = "picrust2/AGP.fasta")

write.fasta(names = attr(chi_fa, which = "names"),
            sequences = chi_fa,
            file.out = "picrust2/GGMP.fasta")

# Combined file

asv_all <- bind_rows(asv_agp, asv_chi) # 6967 samples x 13000 ASV
asv_all[is.na(asv_all)] <- 0
rownames(asv_all) <- c(as.character(meta_agp$SampleID),
                       as.character(meta_chi$sample_name))
# write.csv(asv_all, "input/OTU_both_6967x13000_asv.csv")

asv_biom <- make_biom(t(asv_all), 
                      sample_metadata = NULL,
                      observation_metadata = NULL,
                      id = "Both",
                      matrix_element_type = "int") ; write_biom(asv_biom, "picrust2/aging_both.biom")

asv_names <- colnames(asv_all) ; fa <- fa[asv_names]
write.fasta(names=attr(fa, which="names"),
            sequences=fa,
            file.out="picrust2/aging_both.fasta")



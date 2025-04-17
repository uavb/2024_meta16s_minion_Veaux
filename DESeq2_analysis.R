#!/usr/bin/Rscript



###########################################
######  DESeq2 Analysis on 16S data  ######  
###########################################

# Clean environnement
ls()
rm(list=ls())
ls()


####################
# Install.packages #
####################

#######################
# microbiomeutilities #
#######################

#install.packages("devtools")
#devtools::install_github("microsud/microbiomeutilities")


###############
##  Library  ##
###############

library("microbiomeutilities")
library(microbiome)
library(knitr)
library(tibble)
library(dplyr)
library(RColorBrewer)
library(patchwork)
library(ggpubr)
library(vroom)
library(plyr); library(dplyr) # plyr need to be loaded before dplyr
library(stringr)
library(phyloseq)
library(tibble)
library(vegan)
library(ggplot2)
library(tidyr)
library(devtools)
library(microbiome)
#library(microbiomeSeq)
library(upstartr)
library("easystats")
library(kableExtra)
library("network")
library("hrbrthemes")
library("readxl")
library("microbiomeutilities")
library("microbiomeMarker")
library("DESeq2")
library(openxlsx)

session_info = utils::sessionInfo()
writeLines(capture.output(sessionInfo()), "sessionInfo.txt")



###########################
## Set working directory ##
###########################

setwd("K:/Personnel_non_permanent/Thibault/12_metagenomique16_minion_Morgane/04_Epi2me_241114")



###########################################
##  Import abondance table and taxonomy  ##
###########################################

Abundance_table <- read.table("abundance_table_species_clean.csv", sep = '\t', header = TRUE)
Abundance_table = Abundance_table[,-length(Abundance_table)]



#####################
## Import metadata ##
#####################

metadata      = read_excel("Metadata_16S.xlsx")
metadata_Race = read_excel("NumVeau_Race.xlsx")

metadata = join(metadata, metadata_Race)

# Ajoute une colonne contenant l'information Temps 
metadata$Time <- gsub("_.*", "", metadata$Sample)
metadata <- metadata %>%
  mutate(Sample_name = paste(Time, Calves, sep = "_"))
head(metadata)
tail(metadata)

metadata <- metadata %>%
  mutate(Time_Condition = paste(Time, Condition, sep = "_"))
head(metadata)
tail(metadata)

# Renomme pour avoir la colonne Condition en Anglais
metadata <- metadata %>%
  mutate(Condition = recode(Condition,
                            "Contrôle" = "Untreated", 
                            "Traité" = "Treated"))

metadata <- metadata %>%
  mutate(Time_Condition = recode(Time_Condition,
                                 "T0_Contrôle" = "T0_Untreated",
                                 "T1_Contrôle" = "T1_Untreated",
                                 "T2_Contrôle" = "T2_Untreated",
                                 "T3_Contrôle" = "T3_Untreated",
                                 "T0_Traité" = "T0_Treated",
                                 "T1_Traité" = "T1_Treated",
                                 "T2_Traité" = "T2_Treated",
                                 "T3_Traité" = "T3_Treated"))


###########################
## Format taxonomy table ##
###########################

Abundance_table$ASV_ID = rownames(Abundance_table)

# import column of interest
taxo_table = Abundance_table[,c("ASV_ID","tax")]
colnames(taxo_table) = c("ASV_ID", "Taxon")

# Remoove useless info : d__Bacteria => Bacteria;
#taxo_table$Taxon = gsub(" *.__", "", taxo_table$Taxon)

taxo_table_split = cbind(taxo_table$ASV_ID, data.frame(str_split_fixed(taxo_table$Taxon, pattern = ";", n = 8)))
taxo_table_split = taxo_table_split[,-3] # enleve les bacteria_none
colnames(taxo_table_split) = c("ASV_ID", "Domain", "Phylum", "Class", "Order", "Family", "Genus", "Specie")

taxo_table = taxo_table_split


# Remplacer les valeurs "unclassified" par "Unknown_domain" dans la colonne "Domain"
taxo_table <- mutate(taxo_table, Domain = ifelse(Domain == "unclassified", "Unknown_domain", Domain))

# Remplacer les valeurs vides par "Unknown_phylum", "Unknown_class", etc. dans les colonnes correspondantes
taxo_table <- mutate(taxo_table,
                     Phylum = ifelse(Phylum == "Unknown", "Unknown_phylum", Phylum),
                     Class = ifelse(Class == "Unknown", "Unknown_class", Class),
                     Order = ifelse(Order == "Unknown", "Unknown_order", Order),
                     Family = ifelse(Family == "Unknown", "Unknown_family", Family),
                     Genus = ifelse(Genus == "Unknown", "Unknown_genus", Genus),
                     Specie = ifelse(Specie == "Unknown", "Unknown_specie", Specie))
head(taxo_table)

# Filtrer les lignes où le Domaine est "Unknown"
taxo_table_clean <- taxo_table %>%
  filter(Domain != "Unknown_domain")
taxo_table = taxo_table_clean

# Afficher le tableau filtré
head(taxo_table_clean)


######################
## Format ASV table ##
######################

ASV_table = Abundance_table
#ASV_table$tax = NULL
#ASV_table$ASV_ID = rownames(ASV_table) # correspondance avec la taxonomie se fait par le numero de ligne
colnames(ASV_table)
rownames(ASV_table) = ASV_table$tax
ASV_table$ASV_ID = NULL
ASV_table$tax = NULL
colnames(ASV_table)

# remoove Controle_positifs
ASV_table_c = select(ASV_table, -contains("Controle"))
# remoove ASV_ID column
ASV_table_c = select(ASV_table_c, -contains("ASV_ID"))

## Enlever les échantillons qui ont été repassés ##
samples_to_remove <- c("T2_5bis", "T2_22", "T2_7bis", "T2_23bis")
ASV_table_c_prune <- ASV_table_c %>%
  select(-all_of(samples_to_remove))

# Renommer la colonne "T2_22bis" en "T2_22"
colnames(ASV_table_c_prune)[colnames(ASV_table_c_prune) == "T2_22bis"] <- "T2_22"

ASV_table_c = ASV_table_c_prune

head(ASV_table)
head(ASV_table_c)



#########################
## CHECK before DESeq2 ##
#########################

cts = ASV_table_c
coldata = metadata
rownames(coldata) = coldata$Sample

# Réarrange les colonnes de cts pour quels soient dans le meme ordre que les rownames de coldata
cts_reordered <- cts[, rownames(coldata)]
cts = cts_reordered

head(cts)
head(coldata)
dim(coldata)
dim(cts)

# check : verifie si le nombre d'echantillon est le meme entre coldata et cts
all(rownames(coldata) %in% colnames(cts)) # must be "TRUE"
all(rownames(coldata) == colnames(cts))   # must be "TRUE"
coldata$Time <- as.factor(coldata$Time) 
coldata$Condition <- as.factor(coldata$Condition) 



############
## DESeq2 ##
############

# DESeqDataSet is a subclass of RangedSummarizedExperiment, used to store the input values, intermediate calculations and results of an analysis of differential expression
# Cree un objet DESeqDataSet a partir d une matrice de comptage d expression genique et d un fichier de metadonnees decrivant les echantillons.
dds <- DESeqDataSetFromMatrix(countData = cts,
                              colData = coldata,
                              design = ~Time_Condition) 
#~Prelevement + Treatment + Prelevement:Treatment # but gives only few comparison. Contrast on condition only allow to play with contrast parameter => Here Condition is Time_Condition
# Name_filtered <-rownames(dds)
colSums(cts)

## Normalization
# this function estimates the size factors using the "median ratio method" described by Equation 5 in Anders and Huber (2010)
dds <- estimateSizeFactors(dds) 
dds$sizeFactor
# Permet de compenser les difference de profondeur de sequencage et rendre les comptages comparables entre les echantillons




#################################
###  Filtration des donnees  ###
#################################

# Here we perform a minimal pre-filtering to keep only rows that have at least 10 reads total.
# A popular filter is to ensure at least X (here: X=3) samples with a count of 10 or more.
keep  <- rowSums(counts(dds,normalized=TRUE) >= 10 ) >= 3    # Pour chaque colonne, verifie si le nombre de reads est >=10 + il faut que cela soit vrai au moins 3fois (par ligne).
dds   <- dds[keep,]                                          # Filtre sur les valeurs precedentes
dim(dds)




###########
##  PCA  ##
###########

# Data transformation : quickly estimate dispersion trend and apply a variance stabilizing transformation : remove the dependence of the variance on the mean
#vsd.fast <- vst(dds, fitType='local',blind=FALSE)
vsd.fast <- varianceStabilizingTransformation(dds, fitType='local',blind=FALSE)

### Sample PCA plot for transformed data
plotData<-plotPCA(vsd.fast,intgroup=c("Time_Condition"), returnData=TRUE)
percentVar <- round(100 * attr(plotData, "percentVar"))


PCA<-ggplot(plotData, aes(x = PC1, y = PC2)) + 
  geom_point(aes(fill=factor(Time_Condition)),colour="black", shape=21, alpha=0.9, stroke=1, size=4) +
  scale_fill_manual(values=c('#addd8e','#31a354', '#9ecae1', '#3182bd', '#ffeda0','#feb24c','#c994c7','#dd1c77'), breaks=c("T0_Untreated", "T0_Treated", "T1_Untreated", "T1_Treated", "T2_Untreated", "T2_Treated", "T3_Untreated", "T3_Treated")) +
  xlab(paste0("PC1 (",percentVar[1],"%)")) +
  ylab(paste0("PC2 (",percentVar[2],"%)")) +
  labs(title = "PCA") +
  theme(axis.text.x = element_text(size=14,color="black"), 
        axis.text.y =element_text(size=14,color="black"), 
        axis.title.y = element_text(face="bold", size = 15), 
        axis.title.x = element_text(face="bold", size = 15), 
        panel.background = element_rect(fill = "white",color="black",size=1), 
        legend.position=c(0.1,0.8), 
        legend.background = element_rect(size = 0.5, linetype = "solid", color = "black", fill = "white")) +
  labs(fill = "Time_Condition") + 
  geom_hline(yintercept=0, linetype="dashed") +
  geom_vline(xintercept=0, linetype="dashed")
PCA


######################
## Diff. expression ##
######################

dds <- DESeq(dds, parallel=TRUE)

# Test with design=~Time_Condition
res_T0_UntreatedvsTreated <- results(dds, contrast = c("Time_Condition", "T0_Untreated", "T0_Treated"), alpha = 0.01)
res_T1_UntreatedvsTreated <- results(dds, contrast = c("Time_Condition", "T1_Untreated", "T1_Treated"), alpha = 0.01)
res_T2_UntreatedvsTreated <- results(dds, contrast = c("Time_Condition", "T2_Untreated", "T2_Treated"), alpha = 0.01)
res_T3_UntreatedvsTreated <- results(dds, contrast = c("Time_Condition", "T3_Untreated", "T3_Treated"), alpha = 0.01)

res_Treated_T0vsT1 <- results(dds, contrast = c("Time_Condition", "T0_Treated", "T1_Treated"), alpha = 0.01)
res_Treated_T0vsT2 <- results(dds, contrast = c("Time_Condition", "T0_Treated", "T2_Treated"), alpha = 0.01)
res_Treated_T0vsT3 <- results(dds, contrast = c("Time_Condition", "T0_Treated", "T3_Treated"), alpha = 0.01)
res_Treated_T1vsT2 <- results(dds, contrast = c("Time_Condition", "T1_Treated", "T2_Treated"), alpha = 0.01)
res_Treated_T1vsT3 <- results(dds, contrast = c("Time_Condition", "T1_Treated", "T3_Treated"), alpha = 0.01)
res_Treated_T2vsT3 <- results(dds, contrast = c("Time_Condition", "T2_Treated", "T3_Treated"), alpha = 0.01)

res_Untreated_T0vsT1 <- results(dds, contrast = c("Time_Condition", "T0_Untreated", "T1_Untreated"), alpha = 0.01)
res_Untreated_T0vsT2 <- results(dds, contrast = c("Time_Condition", "T0_Untreated", "T2_Untreated"), alpha = 0.01)
res_Untreated_T0vsT3 <- results(dds, contrast = c("Time_Condition", "T0_Untreated", "T3_Untreated"), alpha = 0.01)
res_Untreated_T1vsT2 <- results(dds, contrast = c("Time_Condition", "T1_Untreated", "T2_Untreated"), alpha = 0.01)
res_Untreated_T1vsT3 <- results(dds, contrast = c("Time_Condition", "T1_Untreated", "T3_Untreated"), alpha = 0.01)
res_Untreated_T2vsT3 <- results(dds, contrast = c("Time_Condition", "T2_Untreated", "T3_Untreated"), alpha = 0.01)


# Filters
summary(res_T0_UntreatedvsTreated)
summary(res_T1_UntreatedvsTreated)
summary(res_T2_UntreatedvsTreated)
summary(res_T3_UntreatedvsTreated)

summary(res_Treated_T0vsT1)
summary(res_Treated_T0vsT2)
summary(res_Treated_T0vsT3)
summary(res_Treated_T1vsT2)
summary(res_Treated_T1vsT3)
summary(res_Treated_T2vsT3)

summary(res_Untreated_T0vsT1)
summary(res_Untreated_T0vsT2)
summary(res_Untreated_T0vsT3)
summary(res_Untreated_T1vsT2)
summary(res_Untreated_T1vsT3)
summary(res_Untreated_T2vsT3)

res_T0_UntreatedvsTreated_0.01 = res_T0_UntreatedvsTreated[which(res_T0_UntreatedvsTreated$padj < 0.01),]
res_T1_UntreatedvsTreated_0.01 = res_T1_UntreatedvsTreated[which(res_T1_UntreatedvsTreated$padj < 0.01),]
res_T2_UntreatedvsTreated_0.01 = res_T2_UntreatedvsTreated[which(res_T2_UntreatedvsTreated$padj < 0.01),]
res_T3_UntreatedvsTreated_0.01 = res_T3_UntreatedvsTreated[which(res_T3_UntreatedvsTreated$padj < 0.01),]

res_Treated_T0vsT1_0.01 = res_Treated_T0vsT1[which(res_Treated_T0vsT1$padj < 0.01),]
res_Treated_T0vsT2_0.01 = res_Treated_T0vsT2[which(res_Treated_T0vsT2$padj < 0.01),]
res_Treated_T0vsT3_0.01 = res_Treated_T0vsT3[which(res_Treated_T0vsT3$padj < 0.01),]
res_Treated_T1vsT2_0.01 = res_Treated_T1vsT2[which(res_Treated_T1vsT2$padj < 0.01),]
res_Treated_T1vsT3_0.01 = res_Treated_T1vsT3[which(res_Treated_T1vsT3$padj < 0.01),]
res_Treated_T2vsT3_0.01 = res_Treated_T2vsT3[which(res_Treated_T2vsT3$padj < 0.01),]

res_Untreated_T0vsT1_0.01 = res_Untreated_T0vsT1[which(res_Untreated_T0vsT1$padj < 0.01),]
res_Untreated_T0vsT2_0.01 = res_Untreated_T0vsT2[which(res_Untreated_T0vsT2$padj < 0.01),]
res_Untreated_T0vsT3_0.01 = res_Untreated_T0vsT3[which(res_Untreated_T0vsT3$padj < 0.01),]
res_Untreated_T1vsT2_0.01 = res_Untreated_T1vsT2[which(res_Untreated_T1vsT2$padj < 0.01),]
res_Untreated_T1vsT3_0.01 = res_Untreated_T1vsT3[which(res_Untreated_T1vsT3$padj < 0.01),]
res_Untreated_T2vsT3_0.01 = res_Untreated_T2vsT3[which(res_Untreated_T2vsT3$padj < 0.01),]

res_Treated_T0vsT1_0.01


# Créer un objet de type workbook pour stockage dans un fichier excel => regroupe toutes les infos pour chaque comparaison
wb <- createWorkbook()

tablexlsx_T0_UntreatedvsTreated_0.01 = cbind(data.frame(res_T0_UntreatedvsTreated_0.01), rownames(data.frame(res_T0_UntreatedvsTreated_0.01)))
tablexlsx_T1_UntreatedvsTreated_0.01 = cbind(data.frame(res_T1_UntreatedvsTreated_0.01), rownames(data.frame(res_T1_UntreatedvsTreated_0.01)))
tablexlsx_T2_UntreatedvsTreated_0.01 = cbind(data.frame(res_T2_UntreatedvsTreated_0.01), rownames(data.frame(res_T2_UntreatedvsTreated_0.01)))
tablexlsx_T3_UntreatedvsTreated_0.01 = cbind(data.frame(res_T3_UntreatedvsTreated_0.01), rownames(data.frame(res_T3_UntreatedvsTreated_0.01)))

tablexlsx_Treated_T0vsT1_0.01 = cbind(data.frame(res_Treated_T0vsT1_0.01), rownames(data.frame(res_Treated_T0vsT1_0.01)))
tablexlsx_Treated_T0vsT2_0.01 = cbind(data.frame(res_Treated_T0vsT2_0.01), rownames(data.frame(res_Treated_T0vsT2_0.01)))
tablexlsx_Treated_T0vsT3_0.01 = cbind(data.frame(res_Treated_T0vsT3_0.01), rownames(data.frame(res_Treated_T0vsT3_0.01)))
tablexlsx_Treated_T1vsT2_0.01 = cbind(data.frame(res_Treated_T1vsT2_0.01), rownames(data.frame(res_Treated_T1vsT2_0.01)))
tablexlsx_Treated_T1vsT3_0.01 = cbind(data.frame(res_Treated_T1vsT3_0.01), rownames(data.frame(res_Treated_T1vsT3_0.01)))
tablexlsx_Treated_T2vsT3_0.01 = cbind(data.frame(res_Treated_T2vsT3_0.01), rownames(data.frame(res_Treated_T2vsT3_0.01)))

tablexlsx_Untreated_T0vsT1_0.01 = cbind(data.frame(res_Untreated_T0vsT1_0.01), rownames(data.frame(res_Untreated_T0vsT1_0.01)))
tablexlsx_Untreated_T0vsT2_0.01 = cbind(data.frame(res_Untreated_T0vsT2_0.01), rownames(data.frame(res_Untreated_T0vsT2_0.01)))
tablexlsx_Untreated_T0vsT3_0.01 = cbind(data.frame(res_Untreated_T0vsT3_0.01), rownames(data.frame(res_Untreated_T0vsT3_0.01)))
tablexlsx_Untreated_T1vsT2_0.01 = cbind(data.frame(res_Untreated_T1vsT2_0.01), rownames(data.frame(res_Untreated_T1vsT2_0.01)))
tablexlsx_Untreated_T1vsT3_0.01 = cbind(data.frame(res_Untreated_T1vsT3_0.01), rownames(data.frame(res_Untreated_T1vsT3_0.01)))
tablexlsx_Untreated_T2vsT3_0.01 = cbind(data.frame(res_Untreated_T2vsT3_0.01), rownames(data.frame(res_Untreated_T2vsT3_0.01)))



addWorksheet(wb, "T0_UntreatedvsTreated_0.01")
addWorksheet(wb, "T1_UntreatedvsTreated_0.01")
addWorksheet(wb, "T2_UntreatedvsTreated_0.01")
addWorksheet(wb, "T3_UntreatedvsTreated_0.01")

addWorksheet(wb, "Treated_T0vsT1_0.01")
addWorksheet(wb, "Treated_T0vsT2_0.01")
addWorksheet(wb, "Treated_T0vsT3_0.01")
addWorksheet(wb, "Treated_T1vsT2_0.01")
addWorksheet(wb, "Treated_T1vsT3_0.01")
addWorksheet(wb, "Treated_T2vsT3_0.01")

addWorksheet(wb, "Untreated_T0vsT1_0.01")
addWorksheet(wb, "Untreated_T0vsT2_0.01")
addWorksheet(wb, "Untreated_T0vsT3_0.01")
addWorksheet(wb, "Untreated_T1vsT2_0.01")
addWorksheet(wb, "Untreated_T1vsT3_0.01")
addWorksheet(wb, "Untreated_T2vsT3_0.01")



writeDataTable(wb, sheet = "T0_UntreatedvsTreated_0.01", x = tablexlsx_T0_UntreatedvsTreated_0.01, startRow = 1, startCol = 1)
writeDataTable(wb, sheet = "T1_UntreatedvsTreated_0.01", x = tablexlsx_T1_UntreatedvsTreated_0.01, startRow = 1, startCol = 1)
writeDataTable(wb, sheet = "T2_UntreatedvsTreated_0.01", x = tablexlsx_T2_UntreatedvsTreated_0.01, startRow = 1, startCol = 1)
writeDataTable(wb, sheet = "T3_UntreatedvsTreated_0.01", x = tablexlsx_T3_UntreatedvsTreated_0.01, startRow = 1, startCol = 1)

writeDataTable(wb, sheet = "Treated_T0vsT1_0.01", x = tablexlsx_Treated_T0vsT1_0.01, startRow = 1, startCol = 1)
writeDataTable(wb, sheet = "Treated_T0vsT2_0.01", x = tablexlsx_Treated_T0vsT2_0.01, startRow = 1, startCol = 1)
writeDataTable(wb, sheet = "Treated_T0vsT3_0.01", x = tablexlsx_Treated_T0vsT3_0.01, startRow = 1, startCol = 1)
writeDataTable(wb, sheet = "Treated_T1vsT2_0.01", x = tablexlsx_Treated_T1vsT2_0.01, startRow = 1, startCol = 1)
writeDataTable(wb, sheet = "Treated_T1vsT3_0.01", x = tablexlsx_Treated_T1vsT3_0.01, startRow = 1, startCol = 1)
writeDataTable(wb, sheet = "Treated_T2vsT3_0.01", x = tablexlsx_Treated_T2vsT3_0.01, startRow = 1, startCol = 1)

writeDataTable(wb, sheet = "Untreated_T0vsT1_0.01", x = tablexlsx_Untreated_T0vsT1_0.01, startRow = 1, startCol = 1)
writeDataTable(wb, sheet = "Untreated_T0vsT2_0.01", x = tablexlsx_Untreated_T0vsT2_0.01, startRow = 1, startCol = 1)
writeDataTable(wb, sheet = "Untreated_T0vsT3_0.01", x = tablexlsx_Untreated_T0vsT3_0.01, startRow = 1, startCol = 1)
writeDataTable(wb, sheet = "Untreated_T1vsT2_0.01", x = tablexlsx_Untreated_T1vsT2_0.01, startRow = 1, startCol = 1)
writeDataTable(wb, sheet = "Untreated_T1vsT3_0.01", x = tablexlsx_Untreated_T1vsT3_0.01, startRow = 1, startCol = 1)
writeDataTable(wb, sheet = "Untreated_T2vsT3_0.01", x = tablexlsx_Untreated_T2vsT3_0.01, startRow = 1, startCol = 1)

# Enregistrer le fichier Excel
saveWorkbook(wb, "DESeq2_analysis_16S_alpha0.01.xlsx", overwrite = TRUE)





################################################################################





###########################
## Build phyloseq object ##
###########################

otu_mat <- ASV_table %>%
  tibble::column_to_rownames("ASV_ID") %>%
  as.matrix()
head(otu_mat)

tax_mat <- taxo_table %>% 
  tibble::column_to_rownames("ASV_ID") %>%
  as.matrix()
head(tax_mat)

samples_df <- metadata %>% 
  tibble::column_to_rownames("Sample") 
head(samples_df)

OTU = otu_table(otu_mat, taxa_are_rows = TRUE)
TAX = tax_table(tax_mat)
samples = sample_data(samples_df)

# Transform to phyloseq objects
pseq <- phyloseq(OTU, TAX, samples)
pseq
pseq@otu_table

sample_names(pseq)
rank_names(pseq)
sample_variables(pseq)

pseq_raw = pseq

<<<<<<< HEAD:DESeq2_analysis.R
## TEST ##
# Normalize number of reads in each sample using median sequencing depth.
total = median(sample_sums(pseq))
standf = function(x, t=total) round(t * (x / sum(x)))
pseq_norm1 = transform_sample_counts(pseq, standf)
## END TEST ##

# Quelques infos utiles avant de se lancer dans les figures
ntaxa(pseq)
nsamples(pseq)
sample_names(pseq)[1:5] 
rank_names(pseq)  
sample_variables(pseq)  
otu_table(pseq)[1:5, 1:5]  
tax_table(pseq)[1:5, 1:4]


# Rarefy the phyloseq object to even depth prior various analysis
set.seed(1)
physeq_rarefy <- rarefy_even_depth(pseq, rngseed=1, sample.size=0.9*min(sample_sums(pseq)), replace=F)
plot_taxa_prevalence(pseq, "Phylum")
plot_taxa_prevalence(physeq_rarefy, "Phylum")

# Transforme en abondance relative
physeq_AR = normalize(pseq, "TSS")

physeq_rarefy = physeq_AR



#################
# Dominant taxa #
#################

Filtre_Time = "T1"

# Filtrer les échantillons en fonction de la colonne "Wash" des métadonnées
physeq_filtered <- subset_samples(physeq_rarefy, Time %in% Filtre_Time)

# Calculer les taxons agrégés au niveau du genre
physeq_rarefy.gen <- aggregate_taxa(physeq_filtered, "Genus")

# Calculer les taxons dominants
x.d <- dominant_taxa(physeq_filtered, level = "Genus", group = "Condition")
x.d



#############################
# Get Taxa summary by group #
#############################

grp_abund <- get_group_abundances(physeq_filtered, level = "Phylum", group="Condition")

mycols <- c("steelblue", "brown3")

mean.plot <- grp_abund %>% # input data
  ggplot(aes(x= reorder(OTUID, mean_abundance), # reroder based on mean abundance
             y= mean_abundance,
             fill=Condition)) + # x and y axis 
  geom_bar(     stat = "identity", 
                position=position_dodge()) + 
  scale_fill_manual("Condition", values=mycols) + # manually specify colors
  theme_bw() + # add a widely used ggplot2 theme
  ylab("Mean Relative Abundance") + # label y axis
  xlab("Phylum") + # label x axis
  coord_flip() + # rotate plot
  ggtitle(paste0("Taxa summary by group / ",Filtre_Time))

mean.plot



###################
# alpha diversity #
###################

mycols <- c("brown3", "steelblue")

p.m <- plot_diversity_stats(physeq_rarefy, group = "Condition", 
                            index = "diversity_shannon", 
                            group.order = c("Traité", "Contrôle"), 
                            group.colors = mycols,
                            label.format="p.format",
                            stats = TRUE)
#> Warning: `guides(<scale> = FALSE)` is deprecated. Please use `guides(<scale> =
#> "none")` instead.
p.m + ylab("Shannon Diversity") + 
  xlab("") + 
  ggtitle(paste0("Alpha diversity Shannon / Wilcox test"))


mycols <- c("steelblue", "brown3", "green")

p.m <- plot_diversity_stats(physeq_rarefy, group = "Time", 
                            index = "diversity_shannon", 
                            group.order = c("T0", "T1", "T2"), 
                            group.colors = mycols,
                            label.format="p.format",
                            stats = TRUE)
#> Warning: `guides(<scale> = FALSE)` is deprecated. Please use `guides(<scale> =
#> "none")` instead.
p.m + ylab("Shannon Diversity") + 
  xlab("") + 
  ggtitle(paste0("Alpha diversity Shannon / Wilcox test"))



mycols <- c("brown1", "cyan", "brown3", "blue", "brown", "darkblue")

p.m <- plot_diversity_stats(pseq, 
                            group = "Time_Condition", 
                            index = "diversity_shannon", 
                            group.order = c("T0_Traité", "T0_Contrôle", "T1_Traité", "T1_Contrôle", "T2_Traité", "T2_Contrôle"), 
                            group.colors = mycols,
                            label.format="p.format",
                            stats = TRUE)
#> Warning: `guides(<scale> = FALSE)` is deprecated. Please use `guides(<scale> =
#> "none")` instead.
p.m + ylab("Shannon Diversity") + 
  xlab("") + 
  ggtitle(paste0("Alpha diversity Shannon / Wilcox test"))



# Pour ajouter les noms des échantillons sur la figure d'alpha diversité
# il faut travailler sur les données non trimmées et donc pseq (le resultat est le meme sur la figure finale)

# Plot initial avec plot_diversity_stats
p.m <- plot_diversity_stats(pseq, 
                            group = "Time_Condition", 
                            index = "diversity_shannon", 
                            group.order = c("T0_Traité", "T0_Contrôle", "T1_Traité", "T1_Contrôle", "T2_Traité", "T2_Contrôle"), 
                            group.colors = mycols,
                            label.format="p.format",
                            stats = TRUE)


# Extraire les données des échantillons et des valeurs de diversité
sample_data <- as.data.frame(sample_data(pseq))
sample_data$Shannon_Diversity <- estimate_richness(pseq, split = TRUE)$Shannon

# Ajouter les noms des échantillons aux points du plot
p.m <- p.m + 
  geom_text(data = sample_data, 
            aes(x = Time_Condition, 
                y = Shannon_Diversity, 
                label = rownames(sample_data)),
            vjust = -0.5, 
            hjust = 0.5, 
            size = 3, 
            check_overlap = TRUE) +
  ylab("Shannon Diversity") + 
  xlab("") + 
  ggtitle(paste0("Alpha diversity Shannon / Wilcox test"))

# Afficher le plot
print(p.m)


################### ESSAI affichage valeur significatives

# Chargez ggpubr pour utiliser stat_compare_means
library(ggpubr)

# Créez le plot initial de la diversité
p.m <- plot_diversity_stats(pseq, 
                            group = "Time_Condition", 
                            index = "diversity_shannon", 
                            group.order = c("T0_Traité", "T0_Contrôle", "T1_Traité", "T1_Contrôle", "T2_Traité", "T2_Contrôle"), 
                            group.colors = mycols,
                            stats = FALSE)  # Désactiver les stats automatiques

# Ajouter les noms des échantillons aux points du plot
p.m <- p.m + 
  geom_text(data = sample_data, 
            aes(x = Time_Condition, 
                y = Shannon_Diversity, 
                label = rownames(sample_data)),
            vjust = -0.5, 
            hjust = 0.5, 
            size = 3, 
            check_overlap = TRUE) +
  ylab("Shannon Diversity") + 
  xlab("") + 
  ggtitle(paste0("Alpha diversity Shannon / Wilcox test"))

# Afficher le plot
print(p.m)


# Ajouter les stats avec filtrage des p-values significatives (< 0.05)
p.m + ylab("Shannon Diversity") + 
  xlab("") + 
  ggtitle("Alpha diversity Shannon / Wilcox test") +
  stat_compare_means(method = "wilcox.test", 
                     label = "p.signif",  # Utiliser les labels significatifs "p.signif" sinon "p.format" pour les valeurs
                     comparisons = list(c("T0_Traité", "T1_Traité"),
                                        c("T0_Traité", "T1_Contrôle"),
                                        c("T0_Contrôle", "T1_Traité"),
                                        c("T0_Contrôle", "T1_Contrôle"),
                                        c("T1_Traité", "T2_Traité"),
                                        c("T1_Traité", "T2_Contrôle"),
                                        c("T1_Contrôle", "T2_Traité"), 
                                        c("T1_Contrôle", "T2_Contrôle")),
                     hide.ns = TRUE)  # Cacher les tests non significatifs




#####################
# plot taxa boxplot #
#####################

mycols <- c("brown1", "cyan", "brown3", "blue", "brown", "darkblue")

# Créer la boxplot avec les familles sélectionnées
pn <- plot_taxa_boxplot(physeq_rarefy,
                        taxonomic.level = "Genus",
                        top.otu = 16, 
                        group = "Time_Condition",
                        add.violin= FALSE,
                        title = paste0("Selected Genus"), 
                        keep.other = FALSE,
                        group.order = c("T0_Contrôle", "T0_Traité", "T1_Contrôle", "T1_Traité", "T2_Contrôle", "T2_Traité"),
                        group.colors = mycols,
                        dot.size = 1)+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Afficher la boxplot
print(pn)

## Avec test stat ##

# Créer la boxplot avec les genus sélectionnées
pn <- plot_taxa_boxplot(physeq_rarefy,
                        taxonomic.level = "Genus",
                        top.otu = 16, 
                        group = "Time_Condition",
                        add.violin= FALSE,
                        title = paste0("Selected Genus"), 
                        keep.other = FALSE,
                        group.order = c("T0_Contrôle", "T0_Traité", "T1_Contrôle", "T1_Traité", "T2_Contrôle", "T2_Traité"),
                        group.colors = mycols,
                        dot.size = 1) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  stat_compare_means(method = "wilcox.test", label = "p.signif", 
                     comparisons = list(
                       c("T0_Contrôle", "T0_Traité"), 
                       c("T1_Contrôle", "T1_Traité"),
                       c("T0_Contrôle", "T1_Traité"), 
                       c("T0_Contrôle", "T1_Contrôle"),
                       c("T2_Contrôle", "T2_Traité"),
                       c("T2_Contrôle", "T1_Traité"),
                       c("T2_Contrôle", "T0_Traité"),
                       c("T0_Contrôle", "T2_Traité"),
                       c("T1_Contrôle", "T2_Traité"),
                       c("T2_Traité", "T0_Traité"),
                       c("T2_Traité", "T1_Traité")
                     ))

# Afficher la boxplot
print(pn)



###########
# HEATMAP #
###########

# create a gradient color palette for abundance
grad_ab <- colorRampPalette(c("#faf3dd","#f7d486" ,"#5e6472"))
grad_ab_pal <- grad_ab(10)

# create a color palette for varaibles of interest
meta_colors <- list(c("T0" = "yellow", 
                      "T1" = "orange", 
                      "T2" = "red",
                      c("Contrôle" = "brown3", 
                        "Traité" = "steelblue")))

# add labels for pheatmap to detect
#names(meta_colors) <- c("Time", "Condition")

# Top X species to keep
subset_top_taxa = 100

# Agréger les données au niveau du genre
physeq_genus <- tax_glom(physeq_rarefy, taxrank = "Genus")


p <- plot_taxa_heatmap(physeq_genus,
                       subset.top = subset_top_taxa,
                       VariableA = c("Time","Condition"),
                       heatcolors = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100), 
                       transformation = "log10",
                       cluster_rows = T,
                       cluster_cols = T,
                       show_colnames = T,
                       main = paste0("Heatmap with subset of top ", subset_top_taxa, " taxa"))


# Extraire le composant ggplot de la liste retournée par plot_taxa_heatmap()
heatmap_plot <- p$plot$tree_row

# Ajouter un titre dynamique au composant ggplot
heatmap_plot_title <- heatmap_plot + ggtitle(paste0("Heatmap with subset of top ", subset_top_taxa, " taxa"))

#> Top X OTUs selected
#> log10, if zeros in data then log10(1+x) will be used
#> First top taxa were selected and 
#> then abundances tranformed to log10(1+X)
#> Warning in transform(phyobj1, "log10"): OTU table contains zeroes. Using log10(1
#> + x) transform.



########
# NMDS #
########


##############
# Ordination #
##############

pseq.ord <- ordinate(physeq_rarefy, "NMDS", "bray")
plot_ordination(physeq_rarefy, pseq.ord, color="Condition", shape = "Time", title="OTUs")

GPr  = transform_sample_counts(physeq_rarefy, function(x) x / sum(x) )
GPfr = filter_taxa(GPr, function(x) mean(x) > 1e-5, TRUE)

pseq.ord <- ordinate(GPfr, "NMDS", "bray")
plot_ordination(GPfr, pseq.ord, color="Condition", shape = "Time", title="OTUs")

plot_ordination(GPfr, pseq.ord, color="Condition", shape="Time", title="OTUs") +
  geom_point(size=3) + # Ajouter des points pour les OTUs
  #stat_ellipse(aes(group=Time), type="norm", level=0.95) + # Ajouter les ellipses
  theme_bw() + # Utiliser un thème blanc pour le fond
  labs(color="Condition", shape="Time") # Légendes pour les couleurs et les formes

plot_ordination(GPfr, pseq.ord, color="Race", shape="Time", title="OTUs") +
  geom_point(size=3) + # Ajouter des points pour les OTUs
  #stat_ellipse(aes(group=Time), type="norm", level=0.95) + # Ajouter les ellipses
  theme_bw() + # Utiliser un thème blanc pour le fond
  labs(color="Race", shape="Time") # Légendes pour les couleurs et les formes

plot_ordination(GPfr, pseq.ord, color="Age_jours", shape="Time", title="OTUs") +
  geom_point(size=3) + # Ajouter des points pour les OTUs
  #stat_ellipse(aes(group=Time), type="norm", level=0.95) + # Ajouter les ellipses
  theme_bw() + # Utiliser un thème blanc pour le fond
  labs(color="Age_jours", shape="Time") # Légendes pour les couleurs et les formes



=======
>>>>>>> 53724252046b33843e4fc90a69658a0b9a49bbbf:Phyloseq.R
##########
# DESeq2 #
##########


####################
## Sur tous les Temps ##
####################

# Phyloseq to DESeq2
Condsdds = phyloseq_to_deseq2(pseq, ~ Condition)
Condsdds = DESeq(Condsdds, test="Wald", fitType="parametric") 
# choices = "parametric", "local", "mean"


res = results(Condsdds, cooksCutoff = FALSE)
alpha = 0.01
sigtab = res[which(res$padj < alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(pseq)[rownames(sigtab), ], "matrix"))
dim(sigtab)

theme_set(theme_bw())
scale_fill_discrete <- function(palname = "Set1", ...) {
  scale_fill_brewer(palette = palname, ...)
}
# Phylum order
x = tapply(sigtab$log2FoldChange, sigtab$Phylum, function(x) max(x))
x = sort(x, TRUE)
sigtab$Phylum = factor(as.character(sigtab$Phylum), levels=names(x))
# Genus order
x = tapply(sigtab$log2FoldChange, sigtab$Genus, function(x) max(x))
x = sort(x, TRUE)
sigtab$Genus = factor(as.character(sigtab$Genus), levels=names(x))
ggplot(sigtab, aes(x=Genus, y=log2FoldChange, color=Phylum)) + geom_point(size=6) + 
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5))

# verifier les niveau des conditions 
levels(Condsdds$Condition) # Supposons que cela renvoie ["A", "B"], où "A" est la condition de référence.
head(sigtab)

# Si log2FoldChange est négatif pour un taxon, cela signifie que ce taxon est moins abondant dans la condition d'intérêt (B) par rapport à la condition de référence (A).
# Si log2FoldChange est positif, le taxon est plus abondant dans la condition d'intérêt (B) par rapport à la condition de référence (A).



#########################
## Sur la condition traité ##
#########################

# Vérifier les colonnes disponibles dans les métadonnées
sample_data(pseq)

Cond_choisie = "treated"

# Sous-échantillonner pour ne conserver que les échantillons avec la condition "Traité"
pseq_subset <- subset_samples(pseq, Condition == Cond_choisie)
pseq_subset <- subset_samples(pseq_subset, Time != "T2")


# Vérifier les échantillons sélectionnés
sample_data(pseq_subset)

print(sample_data(pseq_subset))

# Phyloseq to DESeq2
Condsdds = phyloseq_to_deseq2(pseq_subset, ~ Time)
Condsdds = DESeq(Condsdds, test="Wald", fitType="parametric")


res = results(Condsdds, cooksCutoff = FALSE)
alpha = 0.01
sigtab = res[which(res$padj < alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(pseq_subset)[rownames(sigtab), ], "matrix"))
dim(sigtab)

theme_set(theme_bw())
scale_fill_discrete <- function(palname = "Set1", ...) {
  scale_fill_brewer(palette = palname, ...)
}
# Phylum order
x = tapply(sigtab$log2FoldChange, sigtab$Phylum, function(x) max(x))
x = sort(x, TRUE)
sigtab$Phylum = factor(as.character(sigtab$Phylum), levels=names(x))
# Genus order
x = tapply(sigtab$log2FoldChange, sigtab$Genus, function(x) max(x))
x = sort(x, TRUE)
sigtab$Genus = factor(as.character(sigtab$Genus), levels=names(x))
ggplot(sigtab, aes(x=Genus, y=log2FoldChange, color=Phylum)) + geom_point(size=6) + 
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5))  +
  ggtitle(paste0(Cond_choisie, " fixe / T0 vs T1 / alpha = ", alpha))

# verifier les niveau des conditions 
levels(Condsdds$Time) # Supposons que cela renvoie ["A", "B"], où "A" est la condition de référence.
head(sigtab)


# Si log2FoldChange est négatif pour un taxon, cela signifie que ce taxon est moins abondant dans la condition d'intérêt (B) par rapport à la condition de référence (A).
# Si log2FoldChange est positif, le taxon est plus abondant dans la condition d'intérêt (B) par rapport à la condition de référence (A).



################
## Sur le Time T1 ##
################

# Vérifier les colonnes disponibles dans les métadonnées
sample_data(pseq)

# Sous-échantillonner pour ne conserver que les échantillons avec le temps "T1"
pseq_subset <- subset_samples(pseq, Time == "T1")


# Vérifier les échantillons sélectionnés
sample_data(pseq_subset)

print(sample_data(pseq_subset))

# Phyloseq to DESeq2
Condsdds = phyloseq_to_deseq2(pseq_subset, ~ Condition)
Condsdds = DESeq(Condsdds, test="Wald", fitType="parametric")


res = results(Condsdds, cooksCutoff = FALSE)
alpha = 0.05
sigtab = res[which(res$padj < alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(pseq_subset)[rownames(sigtab), ], "matrix"))
dim(sigtab)

theme_set(theme_bw())
scale_fill_discrete <- function(palname = "Set1", ...) {
  scale_fill_brewer(palette = palname, ...)
}
# Phylum order
x = tapply(sigtab$log2FoldChange, sigtab$Phylum, function(x) max(x))
x = sort(x, TRUE)
sigtab$Phylum = factor(as.character(sigtab$Phylum), levels=names(x))
# Genus order
x = tapply(sigtab$log2FoldChange, sigtab$Genus, function(x) max(x))
x = sort(x, TRUE)
sigtab$Genus = factor(as.character(sigtab$Genus), levels=names(x))
ggplot(sigtab, aes(x=Genus, y=log2FoldChange, color=Phylum)) + geom_point(size=6) + 
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5)) +
  ggtitle(paste0("T1 fixe / Contrôle vs traité / alpha = ", alpha))

# verifier les niveau des conditions 
levels(Condsdds$Condition) # Supposons que cela renvoie ["A", "B"], où "A" est la condition de référence.
head(sigtab)


# Si log2FoldChange est négatif pour un taxon, cela signifie que ce taxon est moins abondant dans la condition d'intérêt (B) par rapport à la condition de référence (A).
# Si log2FoldChange est positif, le taxon est plus abondant dans la condition d'intérêt (B) par rapport à la condition de référence (A).






##############
# Rare curve #
##############

tab <- otu_table(pseq)
class(tab) <- "matrix" # as.matrix() will do nothing
## you get a warning here, but this is what we need to have
tab <- t(tab) # transpose observations to rows
raremax = min(rowSums(tab))
rare <- rarecurve(tab, step=100, lwd=2, ylab="Richness", sample = raremax, col = "blue", label=F)
# Si on utilise la méthode de rarefaction pour normaliser les données : la barre vertical indique ou sera coupé nos échantillons, on ne devrait pas perdre de diversité pour aucun des échantillons.










################################################################################



###########
# barplot #
###########

Mygrid = "Wash"
Mygrid = "Condition"

# Composition plot phylum
plot_bar(physeq_rarefy, fill = "Phylum") + 
  geom_bar(aes(color=Phylum, fill=Phylum), stat="identity", position="stack")

# Composition plot genus
plot_bar(physeq_rarefy, fill = "Genus") + 
  geom_bar(aes(color=Genus, fill=Genus), stat="identity", position="stack")

# Composition plot - Firmicutes subset and focus on Genus
physeq_rarefy_firmi <- subset_taxa(physeq_rarefy, Phylum %in% c("Firmicutes"))
plot_bar(physeq_rarefy_firmi, x="Genus", fill = "Genus", facet_grid = Mygrid) +
  geom_bar(aes(color=Genus, fill=Genus), stat="identity", position="stack") +
  ggtitle(paste0("Firmicutes subset - focus Genus - Condition:", Mygrid)) +
  theme(legend.position = "none")  # Cette ligne supprime la légende

# Composition plot - Bacteroidota subset and focus on Genus
physeq_rarefy_bacte <- subset_taxa(physeq_rarefy, Phylum %in% c("Bacteroidota"))
plot_bar(physeq_rarefy_bacte, x="Genus", fill = "Genus", facet_grid = Mygrid) +
  geom_bar(aes(color=Genus, fill=Genus), stat="identity", position="stack") +
  ggtitle(paste0("Bacteroidota subset - focus Genus - Condition:", Mygrid)) +
  theme(legend.position = "none")  # Cette ligne supprime la légende

# Composition plot - Synergistota subset and focus on Genus
physeq_rarefy_syne <- subset_taxa(physeq_rarefy, Phylum %in% c("Synergistota"))
plot_bar(physeq_rarefy_syne, x="Genus", fill = "Genus", facet_grid = Mygrid) +
  geom_bar(aes(color=Genus, fill=Genus), stat="identity", position="stack") +
  ggtitle(paste0("Synergistota subset - focus Genus - Condition:", Mygrid)) +
  theme(legend.position = "none")  # Cette ligne supprime la légende

# Composition plot - Actinobacteriota subset and focus on Genus
physeq_rarefy_actino <- subset_taxa(physeq_rarefy, Phylum %in% c("Actinobacteriota"))
plot_bar(physeq_rarefy_actino, x="Genus", fill = "Genus", facet_grid = Mygrid) +
  geom_bar(aes(color=Genus, fill=Genus), stat="identity", position="stack") +
  ggtitle(paste0("Actinobacteriota subset - focus Genus - Condition:", Mygrid)) +
  theme(legend.position = "none")  # Cette ligne supprime la légende

# Composition plot - Proteobacteria subset and focus on Genus
physeq_rarefy_p <- subset_taxa(physeq_rarefy, Phylum %in% c("Proteobacteria"))
plot_bar(physeq_rarefy_p, x="Genus", fill = "Genus", facet_grid = Mygrid) +
  geom_bar(aes(color=Genus, fill=Genus), stat="identity", position="stack") +
  ggtitle(paste0("Proteobacteria subset - focus Genus - Condition:", Mygrid)) +
  theme(legend.position = "none")  # Cette ligne supprime la légende



# Filter on metadata
filtre = "L3"
metadata_filtered <- metadata %>%
  filter(grepl(filtre, Wash))

sample_names <- sample_names(physeq_rarefy)

# Filtrer physeq_rarefy en fonction de metadata_filtered
metadata_filtered_samples <- metadata_filtered$Sample
physeq_filtered <- subset_samples(physeq_rarefy, sample_names(physeq_rarefy) %in% metadata_filtered_samples)

# Créer le graphe en utilisant les données filtrées
plot_bar(physeq_filtered, fill = "Phylum") + 
  geom_bar(aes(color=Phylum, fill=Phylum), stat="identity", position="stack") + 
  ggtitle(paste0("Abondance relative par échantillons : wash ", filtre))


## Species
# Voir pour supprimer les légendes qui n'apportent rien...
# Composition plot - Proteobacteria subset and focus on Specie
physeq_rarefy_p <- subset_taxa(physeq_rarefy, Phylum %in% c("Proteobacteria"))
plot_bar(physeq_rarefy_p, x="Genus", fill = "Genus", facet_grid = Mygrid) +
  geom_bar(aes(color=Genus, fill=Genus), stat="identity", position="stack") +
  ggtitle(paste0("Proteobacteria subset - focus Genus - Condition:", Mygrid)) +
  theme(legend.position = "none")  # Cette ligne supprime la légende

plot_bar(physeq_rarefy_p, x="Specie", fill = "Specie", facet_grid = Mygrid) +
  geom_bar(aes(color=Specie, fill=Specie), stat="identity", position="stack") +
  ggtitle(paste0("Proteobacteria subset - focus Specie - Condition:", Mygrid)) +
  theme(legend.position = "none")  # Cette ligne supprime la légende

# Composition plot - Firmicutes subset and focus on Specie
physeq_rarefy_firmi <- subset_taxa(physeq_rarefy, Phylum %in% c("Firmicutes"))
plot_bar(physeq_rarefy_firmi, x = "Specie", fill = "Specie", facet_grid = Mygrid) +
  geom_bar(aes(color = Specie, fill = Specie), stat = "identity", position = "stack") +
  ggtitle(paste0("Firmicutes subset - focus Specie - Condition:", Mygrid)) +
  theme(legend.position = "none")  # Cette ligne supprime la légende

## Genus
plot_bar(physeq_rarefy_firmi, x = "Genus", fill = "Genus", facet_grid = Mygrid) +
  geom_bar(aes(color = Genus, fill = Genus), stat = "identity", position = "stack") +
  ggtitle(paste0("Firmicutes subset - focus Genus - Condition:", Mygrid)) +
  theme(legend.position = "none")  # Cette ligne supprime la légende

plot_bar(physeq_rarefy_firmi, x = "Genus", fill = "Genus", facet_grid = Mygrid) +
  geom_bar(aes(color = Genus, fill = Genus), stat = "identity", position = "stack") +
  ggtitle(paste0("Firmicutes subset - focus Genus - Condition:", Mygrid)) +
  theme(legend.position = "none")  # Cette ligne supprime la légende




###################
### Add filters ###
###################

# Filtres préalables
filtre_Cond = "controle"
Filtre_Time = "L1"
metadata_filtered_w  <- metadata             %>% filter(grepl(Filtre_Time, Wash))
metadata_filtered_wc <- metadata_filtered_w  %>% filter(grepl(filtre_Cond, Condition))


# Ordonnées sur la figure:
Mygrid = "Condition"

# Subset
physeq_rarefy_firmi <- subset_taxa(physeq_rarefy, Phylum %in% c("Firmicutes"))





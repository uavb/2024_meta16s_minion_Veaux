#!/usr/bin/Rscript

# Lien DESeq2 explication Log2FC
# https://nbisweden.github.io/workshop-RNAseq/2111/lab_dge.html

###########################################
#######  Metagenomic with phyloseq  #######  
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

install.packages("devtools")
devtools::install_github("microsud/microbiomeutilities")

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("Deseq2")


install_github("umerijaz/microbiomeSeq") 
library(microbiomeSeq)



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
library(microbial)
library(microbiome)
library(upstartr)
library("easystats")
library(kableExtra)
library("network")
library("hrbrthemes")
library("readxl")
library("microbiomeutilities")
library("microbiomeMarker") 
library("DESeq2") 
library(tidyverse)
library(ggpubr)
library(openxlsx)


session_info = utils::sessionInfo()
writeLines(capture.output(sessionInfo()), "sessionInfo.txt")



###########################
## Set working directory ##
###########################

setwd("K:/Personnel_non_permanent/Thibault/12_metagenomique16_minion_Morgane/05_Epi2me_250211")



###########################################
##  Import abondance table and taxonomy  ##
###########################################

Abundance_table <- read.csv("abundance_table_species_clean.csv", sep = '\t')


## Enlever les échantillons qui ont été repassés ##
samples_to_remove <- c("T2_5bis", "T2_22", "T2_7bis", "T2_23bis")
abundance_table_pruned <- Abundance_table %>%
  select(-all_of(samples_to_remove))
Abundance_table = abundance_table_pruned

### Regrouper les contrôles positifs entre eux ###
control_positif_table <- Abundance_table %>%
  select(contains("Controle_positif"))

### Remettre le reste dans Abundance_table ### 
Abundance_table <- Abundance_table %>%
  select(-contains("Control_positif"))

### Enlever le "bis" de T2_22bis ###
Abundance_table <- Abundance_table %>%
  rename_with(~ str_replace(., "T2_22bis", "T2_22"))

Abundance_table = Abundance_table[,-length(Abundance_table)]


# Lire la matrice en ignorant la première ligne (tax)
matrice <- Abundance_table[,-1]
matrice_numeric <- apply(matrice, 2, function(x) as.numeric(as.character(x)))


# Calculer le total des reads par échantillon (somme des lignes pour chaque colonne)
reads_par_echantillon <- colSums(matrice_numeric)

# Afficher le résultat
print(reads_par_echantillon)

# Sauvegarder le résultat dans un fichier
write.table(reads_par_echantillon, file = "reads_par_echantillon.txt", sep = "\t", quote = FALSE, col.names = NA)



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

metadata <- metadata %>%
  mutate(Time_Condition = paste(Time, Condition, sep = "_"))
head(metadata)
tail(metadata)

## Changer "Contrôle" par Untreated et "Traité" par Treated
# Replace specific values using dplyr's recode function
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

metadata <- metadata %>%
  mutate(Condition = recode(Condition,
                            "Contrôle" = "Untreated", 
                            "Traité" = "Treated"))



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
print(taxo_table_clean)



######################
## Format ASV table ##
######################

ASV_table         = Abundance_table
ASV_table$tax     = NULL
ASV_table$ASV_ID  = rownames(ASV_table)
colnames(ASV_table)
head(ASV_table)



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


# Quelques infos utiles avant de se lancer dans les figures
ntaxa(pseq)
nsamples(pseq)
sample_names(pseq)[1:5] 
rank_names(pseq)  
sample_variables(pseq)  
otu_table(pseq)[1:5, 1:5]  
tax_table(pseq)[1:5, 1:4]

# Retire les samples "Controle_positif" de la suite de l'analyse
pseq <- subset_samples(pseq, !grepl("Controle_positif", sample_names(pseq)))



###################
## Normalization ##
###################

## Rarefaction ##
# Rarefy the phyloseq object to even depth prior various analysis
set.seed(1)
physeq_rarefy <- rarefy_even_depth(pseq, rngseed=1, sample.size=0.9*min(sample_sums(pseq)), replace=F)

plot_taxa_prevalence(pseq, "Phylum")
plot_taxa_prevalence(physeq_rarefy, "Phylum")

## transform ##
# Transforme en abondance relative
#physeq_AR = microbial::normalize(pseq, "TSS") # does not work => use microbiom package instead
physeq_transform = microbiome::transform(pseq, "compositional")



####################################################

  #######################
## CHOOSE one method :  ##
  #######################

physeq_normalize = physeq_transform

####################################################



########
# PCoA #
########

# La PCoA effectue une transformation linéaire de la matrice de distance pour produire des axes principaux (de façon semblable à l'ACP), tandis que le NMDS est un algorithme qui déplace les points dans un nombre de dimensions choisis au hasard afin de s'approcher le plus possible des distances réelles.

# Sur des données normalisées
dist = phyloseq::distance(physeq_normalize, method="bray")
ordination = ordinate(physeq_normalize, method="PCoA", distance=dist)
plot_ordination(physeq_normalize, ordination, color="Condition", shape = "Time") + 
  theme_classic() +
  theme(strip.background = element_blank()) +
  stat_ellipse(aes(group=Time, color=Time),      # Ajouter des ellipses
               type="norm", level=0.95) +
  ggtitle(paste0("PCoA - Bray-Curtis method on normalize data"))
  


########
# NMDS #
########

#################################
# NMDS focus on a specific Time #
#################################

Diff_Time = c("T0", "T1", "T2", "T3")

for (Filtre_Time in Diff_Time) {
  print(Filtre_Time)
  # Filtrer les échantillons en fonction de la colonne "Wash" des métadonnées
  pseq_filtre_time <- subset_samples(pseq, Time %in% Filtre_Time)
  
  pseq.ord <- ordinate(pseq_filtre_time, "NMDS", "bray")
  plot_ordination(pseq_filtre_time, pseq.ord, color="Condition", shape = "Time", title="OTUs")
  
  GPr  = transform_sample_counts(pseq_filtre_time, function(x) x / sum(x) )
  GPfr = filter_taxa(GPr, function(x) mean(x) > 1e-5, TRUE)
  
  pseq.ord <- ordinate(GPfr, "NMDS", "bray")
  plot_ordination(GPfr, pseq.ord, color="Condition", shape = "Time", title="NMDS")
  
  # Création du fichier PNG
  output_file <- paste0(Filtre_Time, "_NMDS.png")
  png(filename = output_file, width = 800, height = 600)
  
  
  # Tracé du graphique
  print(
    plot_ordination(GPfr, pseq.ord, color="Condition", shape="Time", title=paste0("NMDS_",Filtre_Time,".png")) +
      geom_point(size=2) +                     # Ajouter des points pour les OTUs
      stat_ellipse(aes(group=Condition),      # Ajouter des ellipses
                   type="norm", level=0.95) +
      theme_bw() +                            # Utiliser un thème blanc
      labs(color="Condition", shape="Time")  # Légendes
  )
  
  dev.off() # Fermeture du fichier PNG
}
# Results will be written on your working directory


##############
# Ordination #
##############

# Sur des données non normalisées : 
pseq.ord <- ordinate(pseq, "NMDS", "bray")
plot_ordination(pseq, pseq.ord, color="Condition", shape = "Time", title="ASVs")

GPr  = transform_sample_counts(pseq, function(x) x / sum(x) )
GPfr = filter_taxa(GPr, function(x) mean(x) > 1e-5, TRUE)

pseq.ord <- ordinate(GPfr, "NMDS", "bray")
plot_ordination(GPfr, pseq.ord, color="Condition", shape = "Time", title="ASVs")

# Extraire la table d'OTU (matrix ou data.frame)
otu_matrix <- as(otu_table(GPfr), "matrix")

# Calculer le stress pour différents nombres de dimensions
dims <- 1:6  # Tester de 1 à 6 dimensions
stress_values <- sapply(dims, function(k) {
  mds <- metaMDS(otu_matrix, distance = "bray", k = k, trymax = 20, autotransform = FALSE, trace = FALSE)
  mds$stress
})

# Tracer le stress plot
plot(dims, stress_values, type = "b", pch = 19,
     xlab = "Nombre de dimensions", ylab = "Stress", 
     main = "Stress Plot for NMDS")


# Vérifier si les taxons sont en lignes, si oui, les transposer pour avoir les échantillons en lignes
if(taxa_are_rows(pseq)) {
  otu_matrix <- t(otu_matrix)
}

# Lancer metaMDS en utilisant la distance de Bray-Curtis
nmds_result <- metaMDS(otu_matrix, distance = "bray", k = 2, trymax = 20)

print(nmds_result)
stressplot(nmds_result, main = "Shepard plot")


# Time & Condition
plot_ordination(GPfr, pseq.ord, color="Condition", shape="Time") +
  geom_point(size=3) + # Ajouter des points pour les OTUs
  stat_ellipse(aes(group=Time), type="norm", level=0.95) + # Ajouter les ellipses
  theme_bw() + # Utiliser un thème blanc pour le fond
  labs(color="Condition", shape="Time") + # Légendes pour les couleurs et les formes 
  ggtitle(paste0("NMDS - Bray-Curtis - Time & Condition"))

# Time and Race
plot_ordination(GPfr, pseq.ord, color="Breed", shape="Time", title="OTUs") +
  geom_point(size=3) + # Ajouter des points pour les OTUs
  stat_ellipse(aes(group=Time), type="norm", level=0.95) + # Ajouter les ellipses
  theme_bw() + # Utiliser un thème blanc pour le fond
  labs(color="Race", shape="Time") + # Légendes pour les couleurs et les formes 
  ggtitle(paste0("NMDS - Bray-Curtis - Time & Race"))

plot_ordination(GPfr, pseq.ord, color="Age_jours", shape="Time", title="OTUs") +
  geom_point(size=3) + # Ajouter des points pour les OTUs
  stat_ellipse(aes(group=Time), type="norm", level=0.95) + # Ajouter les ellipses
  theme_bw() + # Utiliser un thème blanc pour le fond
  labs(color="Age_jours", shape="Time") + # Légendes pour les couleurs et les formes 
  ggtitle(paste0("NMDS - Bray-Curtis - Time & Age (days)"))

plot_ordination(GPfr, pseq.ord, color="Diarrhea", shape="Time", title="OTUs") +
  geom_point(size=3) + # Ajouter des points pour les OTUs
  stat_ellipse(aes(group=Time), type="norm", level=0.95) + # Ajouter les ellipses
  theme_bw() + # Utiliser un thème blanc pour le fond
  labs(color="Diarrhea", shape="Time") + # Légendes pour les couleurs et les formes 
  ggtitle(paste0("NMDS - Bray-Curtis - Time & Diarrhea"))

mycols <- c("deeppink4", "maroon3", "cornflowerblue", "lightskyblue", "lightseagreen", "darkseagreen","khaki2","khaki3")

plot_ordination(GPfr, pseq.ord, color="Time_Condition", shape="Time") +
  geom_point(size=3) + # Ajouter des points pour les OTUs
  stat_ellipse(aes(group=Time), type="norm", level=0.95) + # Ajouter les ellipses
  theme_bw() + # Utiliser un thème blanc pour le fond
  labs(color="Time_Condition", shape="Time") + # Légendes pour les couleurs et les formes 
  ggtitle(paste0("NMDS - Bray-Curtis - Time & Time_Condition"))
                       
plot_ordination(GPfr, pseq.ord, color="Time_Condition", shape="Time") +
  geom_point(size=3) +                            # Ajouter des points pour les OTUs
  #stat_ellipse(aes(group = Time_Condition, color = Time_Condition), type = "norm", level = 0.95) +
  theme_bw() +                                    # Utiliser un thème blanc pour le fond
  labs(color="Time_Condition", shape="Time") +    # Légendes pour les couleurs et les formes 
  ggtitle("NMDS - Bray-Curtis - Time & Time_Condition") +
  scale_color_manual(values = mycols)             # Appliquer vos couleurs personnalisées

# Ajout des noms à la figure
ordination_plot <- plot_ordination(GPfr, pseq.ord, color="Time_Condition", shape="Time")
ordination_df <- ordination_plot$data  # Récupérer les données du plot

ordination_df$Sample <- rownames(ordination_df)  # Ajouter noms des échantillons

ggplot(ordination_df, aes(x = NMDS1, y = NMDS2, color = Time_Condition, shape = Time)) +
  geom_point(size=3) +
  geom_text(aes(label = Sample), vjust = -1, size = 3) +
  theme_bw() +
  labs(color="Time_Condition", shape="Time") +    
  ggtitle("NMDS - Bray-Curtis - Time & Time_Condition") +
  scale_color_manual(values = mycols)





###################
# Figure générale #
###################

# Sur les données normalisées
ps.rel = physeq_normalize

# agglomerate taxa by phylum
glom <- tax_glom(ps.rel, taxrank = 'Phylum', NArm = FALSE)
ps.melt <- psmelt(glom)
ps.melt$Abundance = ps.melt$Abundance*100

# change to character for easy-adjusted level
ps.melt$Phylum <- as.character(ps.melt$Phylum)

# Selectionne les phylum les plus abondant (Cette information est donnée par epi2me => fichier wf-metagenomic-report.htlm)
Phylum_s = c("Bacillota", "Pseudomonadota", "Bacteroidota", "Fusobacteriota", "Actinomycetota", "Cyanobacteriota", "Mycoplasmatota", "Thermodesulfobacteriota")

# Regrouper les phylums non sélectionnés sous "Other"
ps.melt_F <- ps.melt %>%
  mutate(Phylum = ifelse(Phylum %in% Phylum_s, Phylum, "Other"))

# Créer le barplot empilé
ggplot(ps.melt_F, aes(x = Sample, y = Abundance, fill = Phylum)) + 
  geom_bar(stat = "identity", position = "stack") + 
  labs(x = "Samples", y = "Abondance relative (%)", fill = "Phylum") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size=5)) +
  scale_fill_manual(values = c(
    "Bacillota" = "#1f78b4", "Pseudomonadota" = "#33a02c", "Bacteroidota" = "#e31a1c",
    "Fusobacteriota" = "#ff7f00", "Actinomycetota" = "#6a3d9a", "Cyanobacteriota" = "#b15928",
    "Mycoplasmatota" = "#a6cee3", "Thermodesulfobacteriota" = "#b2df8a", "Other" = "grey50"
  )) +
  ggtitle("Barplot of Relative Abundance for Each Sample") +
  facet_grid(~ Time+Condition, scale = "free")  # Séparation en fonction de "Condition"


# agglomerate taxa by Family
glom <- tax_glom(ps.rel, taxrank = 'Family', NArm = FALSE)
ps.melt <- psmelt(glom)
ps.melt$Abundance = ps.melt$Abundance*100

# change to character for easy-adjusted level
ps.melt$Family <- as.character(ps.melt$Family)

# Selectionne les phylum les plus abondant (Cette information est donnée par epi2me => fichier wf-metagenomic-report.htlm)
Family_s = c("Lachnospiraceae", "Oscillospiraceae", "Lactobacillaceae", "Clostridiaceae", "Enterobacteriaceae", "Enterococcaceae", "Selenomonadaceae", "Peptostreptococcaceae", "Streptococcaceae" )

# Regrouper les phylums non sélectionnés sous "Other"
ps.melt_F <- ps.melt %>%
  mutate(Family = ifelse(Family %in% Family_s, Family, "Other"))

# Créer le barplot empilé
ggplot(ps.melt_F, aes(x = Sample, y = Abundance, fill = Family)) + 
  geom_bar(stat = "identity", position = "stack") + 
  labs(x = "Samples", y = "Abondance relative (%)", fill = "Family") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size=5)) +
  scale_fill_manual(values = c(
    "Lachnospiraceae" = "#1f78b4", "Oscillospiraceae" = "#33a02c", "Lactobacillaceae" = "#e31a1c",
    "Clostridiaceae" = "#ff7f00", "Enterobacteriaceae" = "#6a3d9a", "Enterococcaceae" = "#b15928",
    "Selenomonadaceae" = "#a6cee3", "Peptostreptococcaceae" = "#b2df8a", "Streptococcaceae" = "yellow", "Other" = "grey50"
  )) +
  ggtitle("Barplot of Relative Abundance for Each Sample") +
  facet_grid(~ Time+Condition, scale = "free")  # Séparation en fonction de "Condition"


# agglomerate taxa by Genus
glom <- tax_glom(ps.rel, taxrank = 'Genus', NArm = FALSE)
ps.melt <- psmelt(glom)
ps.melt$Abundance = ps.melt$Abundance*100

# change to character for easy-adjusted level
ps.melt$Genus <- as.character(ps.melt$Genus)

# Selectionne les phylum les plus abondant (Cette information est donnée par epi2me => fichier wf-metagenomic-report.htlm)
Genus_s = c("Blautia", "Faecalibacterium", "Mediterraneeibacter", "Ruminococcus", "Clostridium", "Lachnoclostridium", "Coprococcus", "Enterococcus", "Escherichia" )

# Regrouper les phylums non sélectionnés sous "Other"
ps.melt_F <- ps.melt %>%
  mutate(Genus = ifelse(Genus %in% Genus_s, Genus, "Other"))

# Créer le barplot empilé
ggplot(ps.melt_F, aes(x = Sample, y = Abundance, fill = Genus)) + 
  geom_bar(stat = "identity", position = "stack") + 
  labs(x = "Samples", y = "Abondance relative (%)", fill = "Genus") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size=5)) +
  scale_fill_manual(values = c(
    "Blautia" = "#1f78b4", "Faecalibacterium" = "#33a02c", "Mediterraneeibacter" = "#e31a1c",
    "Ruminococcus" = "#ff7f00", "Clostridium" = "#6a3d9a", "Lachnoclostridium" = "#b15928",
    "Coprococcus" = "#a6cee3", "Enterococcus" = "#b2df8a", "Escherichia" = "yellow", "Other" = "grey50"
  )) +
  ggtitle("Barplot of Relative Abundance for Each Sample") +
  facet_grid(~ Time+Condition, scale = "free")  # Séparation en fonction de "Condition"



#################
# Dominant taxa #
#################

# VARIABLE A MODIFIER :
Filtre_Time = "T1"

# Filtrer les échantillons en fonction de la colonne "Wash" des métadonnées
physeq_filtered <- subset_samples(physeq_normalize, Time %in% Filtre_Time)

# Calculer les taxons agrégés au niveau du genre
physeq_normalize.gen <- aggregate_taxa(physeq_filtered, "Genus")

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

mycols <- c("deeppink4", "maroon3", "cornflowerblue", "lightskyblue", "lightseagreen", "darkseagreen","khaki2","khaki3")


# Affichage valeur significatives
# Chargez ggpubr pour utiliser stat_compare_means
symnum.args <- list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, Inf), symbols = c("****", "***", "**", "*", "ns"))
symnum.args <- list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, Inf), symbols = c("***", "**", "*", "ns", "ns"))

# Créez le plot initial de la diversité
p.m <- plot_diversity_stats(physeq_normalize, 
                            group = "Time_Condition", 
                            index = "diversity_shannon", 
                            group.order = c("T0_Untreated", "T0_Treated", "T1_Untreated", "T1_Treated", "T2_Untreated","T2_Treated", "T3_Untreated", "T3_Treated"),
                            group.colors = mycols,
                            stats = FALSE)  # Désactiver les stats automatiques

# Extraire les données des échantillons et des valeurs de diversité
sample_data <- as.data.frame(sample_data(pseq))
sample_data$Shannon_Diversity <- estimate_richness(pseq, split = TRUE)$Shannon

# Ajouter les noms des échantillons aux points du plot
#p.m <- p.m + 
#  geom_text(data = sample_data, 
#            aes(x = Time_Condition, 
#                y = Shannon_Diversity, 
#                label = rownames(sample_data)),
#            vjust = -0.5, 
#            hjust = 0.5, 
#            size = 3, 
#            check_overlap = TRUE) +
#  ylab("Shannon Diversity") + 
#  xlab("") + 
#  ggtitle(paste0("Alpha diversity Shannon / Wilcox test"))

# Ajouter les stats avec filtrage des p-values significatives (< 0.05)
p.m + ylab("Shannon Diversity") + 
  xlab("") + 
  ggtitle("Alpha diversity Shannon / Wilcoxon test") +
  stat_compare_means(method = "wilcox.test", 
                     symnum = symnum.args,
                     label = "p.signif",  # Utiliser les labels significatifs "p.signif" sinon "p.format" pour les valeurs
                     comparisons = list(c("T0_Untreated", "T1_Untreated"),
                                        c("T0_Untreated", "T1_Treated"),
                                        c("T0_Untreated", "T2_Treated"),
                                        c("T0_Untreated", "T3_Treated"),
                                        c("T0_Treated", "T1_Treated"),
                                        c("T0_Treated", "T2_Treated"),
                                        c("T0_Treated", "T3_Treated"),
                                        c("T1_Untreated", "T1_Treated"),
                                        c("T1_Untreated", "T2_Untreated"),
                                        c("T1_Untreated", "T2_Treated"),
                                        c("T1_Untreated", "T3_Untreated"), 
                                        c("T1_Untreated", "T3_Treated"),
                                        c("T1_Treated", "T2_Untreated"),
                                        c("T1_Treated", "T2_Treated"),
                                        c("T1_Treated", "T3_Untreated"), 
                                        c("T1_Treated", "T3_Treated"),
                                        c("T2_Untreated", "T3_Treated"), 
                                        c("T3_Untreated", "T3_Treated")),
                     hide.ns = FALSE)  # Cacher les tests non significatifs



##########################
# Top 10 taxa => boxplot #
##########################

mycols <- c("deeppink4", "maroon3", "cornflowerblue", "lightskyblue", "lightseagreen", "darkseagreen","khaki2","khaki3")

# Créer la boxplot avec les familles sélectionnées
pn <- plot_taxa_boxplot(physeq_normalize,
                        taxonomic.level = "Phylum",
                        top.otu = 10, 
                        group = "Time_Condition",
                        add.violin= FALSE,
                        title = paste0("Selected Phylum"), 
                        keep.other = FALSE,
                        group.order = c("T0_Untreated", "T0_Treated", "T1_Untreated", "T1_Treated", "T2_Untreated","T2_Treated", "T3_Untreated", "T3_Treated"),
                        group.colors = mycols,
                        dot.size = 1)+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Afficher la boxplot
print(pn)


####################################
### PLOT SPECIFIC TAXA RANK ####
####################################

  #######################################
## Filtre sur les 5 phylums majoritaires ##
  #######################################

# Ensure tax_table is properly formatted
tax_table(physeq_normalize) <- tax_table(physeq_normalize)

# Get the top 5 phyla
top_5_phyla <- names(sort(tapply(taxa_sums(physeq_normalize), tax_table(physeq_normalize)[, "Phylum"], sum), decreasing = TRUE)[1:5])

# Get OTUs corresponding to top 5 phyla
top_otus <- taxa_names(physeq_normalize)[tax_table(physeq_normalize)[, "Phylum"] %in% top_5_phyla]

# Subset using prune_taxa
physeq_top5 <- prune_taxa(top_otus, physeq_normalize)

# Print the phylum names
print(unique(tax_table(physeq_top5)[, "Phylum"]))

# Melt the phyloseq object
melted_physeq <- psmelt(physeq_top5)


##################
  ### Phylum ###
##################

# Get unique phyla
unique_phyla <- unique(melted_physeq$Phylum)

# Filter melted_physeq for the selected phyla
filtered_data <- melted_physeq %>%
  filter(Phylum %in% unique_phyla)%>% 
  group_by(Time_Condition, Phylum, Calves) %>%
  dplyr::summarise(Rel_Abundance = sum(Abundance*100))%>%
  ungroup()


# Check if filtered_data has the expected rows
print(unique(filtered_data$Phylum))
print(unique(filtered_data$Time_Condition))

filtered_data<- melted_physeq %>%
  group_by(Time_Condition,Phylum,Calves) %>%
  dplyr::summarise(Rel_Abundance = sum(Abundance*100)) %>%
  ungroup()

filtered_data$Time_Condition <- factor(filtered_data$Time_Condition, levels = c("T0_Untreated", "T0_Treated", "T1_Untreated", "T1_Treated", "T2_Untreated","T2_Treated", "T3_Untreated", "T3_Treated")) 
levels=c("T0_Untreated", "T0_Treated", "T1_Untreated", "T1_Treated", "T2_Untreated","T2_Treated", "T3_Untreated", "T3_Treated")


   ############################
### Nom du phylum à modifier : ### 
   ############################

specific_Phylum = "Bacteroidota"

# Autres possibilités : Bacillota, Actinomycetota, Bacteroidota, Pseudomonadota


# On filtre sur le phylum 
data <- filtered_data %>% filter(Phylum == specific_Phylum)

# Wilcox test sur les données filtrée
test_results <- compare_means(
  Rel_Abundance ~ Time_Condition,
  data = data,
  method = "wilcox.test"
)

# On garde seulement les comparaison significativement différentes
significant_comparisons <- test_results %>%
  filter(p.adj < 0.05) %>%
  select(group1, group2) %>%
  split(., seq(nrow(.))) %>%
  lapply(function(x) as.character(unlist(x)))

# Plot specific taxa rank
plot_Specfic_taxa_rank <- ggplot(data, aes(x = Time_Condition, y = Rel_Abundance)) +
  geom_boxplot(aes(fill = Time_Condition), outlier.shape = NA, alpha = 0.5) +
  scale_fill_manual(values = mycols) +
  geom_jitter(aes(color = Time_Condition), alpha = 0.5) +
  scale_color_manual(values = mycols) +
  theme_minimal(base_size = 15) +
  theme(
    panel.grid.major = element_line(color = "grey90"),
    panel.grid.minor = element_line(color = "white"),
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) +
  labs(y = "Relative Abundance (%)", x = "Time_Condition") +
  scale_y_continuous(
    labels = scales::percent_format(scale = 1), 
    limits = c(0, NA),  # Pas de limite supérieure explicite pour conserver les données
    breaks = seq(0, 100, by = 20),  # Afficher uniquement les ticks jusqu'à 100%
  ) +
  stat_compare_means(
    comparisons = significant_comparisons,
    method = "wilcox.test",
    label = "p.signif"
  ) +
  ggtitle(paste0("Distribution of Relative Abundance by Time_Condition - Phylum : ", specific_Phylum, sep=" "))

# Afficher le graphique
plot_Specfic_taxa_rank



###################
  ### Familly ###
###################

# Get unique Familly
unique_Family <- unique(melted_physeq$Family)

# Filter melted_physeq for the selected Familly
filtered_data <- melted_physeq %>%
  filter(Family %in% unique_Family)%>% 
  group_by(Time_Condition, Family, Calves) %>%
  dplyr::summarise(Rel_Abundance = sum(Abundance*100))%>%
  ungroup()

# Check if filtered_data has the expected rows
print(unique(filtered_data$Family))
print(unique(filtered_data$Time_Condition))

filtered_data<- melted_physeq %>%
  group_by(Time_Condition,Family,Calves) %>%
  dplyr::summarise(Rel_Abundance = sum(Abundance*100)) %>%
  ungroup()

filtered_data$Time_Condition <- factor(filtered_data$Time_Condition, levels = c("T0_Untreated", "T0_Treated", "T1_Untreated", "T1_Treated", "T2_Untreated","T2_Treated", "T3_Untreated", "T3_Treated")) 
levels=c("T0_Untreated", "T0_Treated", "T1_Untreated", "T1_Treated", "T2_Untreated","T2_Treated", "T3_Untreated", "T3_Treated")

## Plot with top 15 Family ##
# Ensure tax_table is properly formatted
tax_table(physeq_normalize) <- tax_table(physeq_normalize)

# Get the top 15 family
top_15_family <- names(sort(tapply(taxa_sums(physeq_normalize), tax_table(physeq_normalize)[, "Family"], sum), decreasing = TRUE)[1:15])

# Get OTUs corresponding to top 15 families
top_otus_family <- taxa_names(physeq_normalize)[tax_table(physeq_normalize)[, "Family"] %in% top_15_family]

# Subset using prune_taxa
physeq_top15_family <- prune_taxa(top_otus_family, physeq_normalize)

# Print the family names
print(unique(tax_table(physeq_top15_family)[, "Family"]))


    ###############################
### Nom de la family à modifier :  ###
    ###############################

# Apperçu top15 Family : Oscillospiraceae, Lachnospiraceae, Enterobacteriaceae, Enterococcaceae, Selenomonadaceae, Lactobacillaceae, Clostridiaceae, Peptostreptococcaceae, Streptococcaceae, Prevotellaceae, Acidaminococcaceae, Pasteurellaceae, Anaerotignaceae, Eubacteriales_Incertae_sedis, Bacteroidaceae


for (i in top_15_family) {
  print(i)
  
  # On filtre sur la family 
  data <- filtered_data %>% filter(Family == i)
  
  # Wilcox test sur les données filtrée
  test_results <- compare_means(
    Rel_Abundance ~ Time_Condition,
    data = data,
    method = "wilcox.test"
  )
  
  # On garde seulement les comparaison significativement différentes
  significant_comparisons <- test_results %>%
    filter(p.adj < 0.05) %>%
    select(group1, group2) %>%
    split(., seq(nrow(.))) %>%
    lapply(function(x) as.character(unlist(x)))
  
  # Création du fichier PNG
  output_file <- paste0("Rplot_boxplot_", i,".png", sep="")
  png(filename = output_file, width = 1200, height = 703)
  
  print(
    # Plot specific taxa rank
    plot_Specfic_taxa_rank <- ggplot(data, aes(x = Time_Condition, y = Rel_Abundance)) +
      geom_boxplot(aes(fill = Time_Condition), outlier.shape = NA, alpha = 0.5) +
      scale_fill_manual(values = mycols) +
      geom_jitter(aes(color = Time_Condition), alpha = 0.5) +
      scale_color_manual(values = mycols) +
      theme_minimal(base_size = 15) +
      theme(
        panel.grid.major = element_line(color = "grey90"),
        panel.grid.minor = element_line(color = "white"),
        axis.text.x = element_text(angle = 45, hjust = 1)
      ) +
      labs(y = "Relative Abundance (%)", x = "Time_Condition") +
      scale_y_continuous(
        labels = scales::percent_format(scale = 1), 
        limits = c(0, NA),  # Pas de limite supérieure explicite pour conserver les données
        breaks = seq(0, 100, by = 10),  # Afficher uniquement les ticks jusqu'à 100%
      ) +
      stat_compare_means(
        comparisons = significant_comparisons,
        method = "wilcox.test",
        label = "p.signif"
      ) +
      ggtitle(paste0("Distribution of Relative Abundance by Time_Condition - Family : ", i, sep=" "))
  )
  
  dev.off()
}



###################
   ### Genus ###
###################

# Get unique Genus
unique_Genus <- unique(melted_physeq$Genus)

# Filter melted_physeq for the selected Genus
filtered_data <- melted_physeq %>%
  filter(Genus %in% unique_Genus)%>% 
  group_by(Time_Condition, Genus, Calves) %>%
  dplyr::summarise(Rel_Abundance = sum(Abundance*100))%>%
  ungroup()


# Check if filtered_data has the expected rows
print(unique(filtered_data$Genus))
print(unique(filtered_data$Time_Condition))

filtered_data<- melted_physeq %>%
  group_by(Time_Condition,Genus,Calves) %>%
  dplyr::summarise(Rel_Abundance = sum(Abundance*100)) %>%
  ungroup()

filtered_data$Time_Condition <- factor(filtered_data$Time_Condition, levels = c("T0_Untreated", "T0_Treated", "T1_Untreated", "T1_Treated", "T2_Untreated","T2_Treated", "T3_Untreated", "T3_Treated")) 
levels=c("T0_Untreated", "T0_Treated", "T1_Untreated", "T1_Treated", "T2_Untreated","T2_Treated", "T3_Untreated", "T3_Treated")

## Plot with top 10 Genus ##
# Ensure tax_table is properly formatted
tax_table(physeq_normalize) <- tax_table(physeq_normalize)

# Get the top 10 Genus
top_15_Genus <- names(sort(tapply(taxa_sums(physeq_normalize), tax_table(physeq_normalize)[, "Genus"], sum), decreasing = TRUE)[1:15])

# Get OTUs corresponding to top 10 families
top_otus_Genus <- taxa_names(physeq_normalize)[tax_table(physeq_normalize)[, "Genus"] %in% top_15_Genus]

# Subset using prune_taxa
physeq_top15_Genus <- prune_taxa(top_otus_Genus, physeq_normalize)

# Print the Genus names
print(unique(tax_table(physeq_top15_Genus)[, "Genus"]))


##################################
# Nom de la Genus à modifier :  #
##################################


# Apperçu top15 Genus : Faecalibacterium, Blautia, Mediterraneibacter, Coprococcus, Ruminococcus, Lachnoclostridium,Megamonas, Escherichia, Enterococcus, Ligilactobacillus, Lactobacillus, Anaerostipes, Clostridium, Romboutsia, Streptococcus

for (i in top_15_Genus) {
  # On filtre sur le Genus 
  data <- filtered_data %>% filter(Genus == i)
  
  # Wilcox test sur les données filtrée
  test_results <- compare_means(
    Rel_Abundance ~ Time_Condition,
    data = data,
    method = "wilcox.test"
  )
  test_results
  
  # On garde seulement les comparaison significativement différentes
  significant_comparisons <- test_results %>%
    filter(p.adj < 0.05) %>%
    select(group1, group2) %>%
    split(., seq(nrow(.))) %>%
    lapply(function(x) as.character(unlist(x)))
  
  output_file <- paste0("Rplot_boxplot_", i,".png", sep="")
  png(filename = output_file, width = 1200, height = 703)
  
  print(# Plot specific taxa rank
    plot_Specfic_taxa_rank <- ggplot(data, aes(x = Time_Condition, y = Rel_Abundance)) +
      geom_boxplot(aes(fill = Time_Condition), outlier.shape = NA, alpha = 0.5) +
      scale_fill_manual(values = mycols) +
      geom_jitter(aes(color = Time_Condition), alpha = 0.5) +
      scale_color_manual(values = mycols) +
      theme_minimal(base_size = 15) +
      theme(
        panel.grid.major = element_line(color = "grey90"),
        panel.grid.minor = element_line(color = "white"),
        axis.text.x = element_text(angle = 45, hjust = 1)
      ) +
      labs(y = "Relative Abundance (%)", x = "Time_Condition") +
      scale_y_continuous(
        labels = scales::percent_format(scale = 1), 
        limits = c(0, NA),  # Pas de limite supérieure explicite pour conserver les données
        breaks = seq(0, 100, by = 10),  # Afficher uniquement les ticks jusqu'à 100%
      ) +
      stat_compare_means(
        comparisons = significant_comparisons,
        method = "wilcox.test",
        label = "p.signif"
      ) +
      ggtitle(paste0("Distribution of Relative Abundance by Time_Condition - Genus : ", i, sep=" "))
    )
  dev.off()

}

# Afficher le graphique
plot_Specfic_taxa_rank



############
## DESeq2 ##
############

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
                      "T3" = "green",
                    c("Contrôle" = "brown3", 
                      "Traité" = "steelblue")))

# add labels for pheatmap to detect
#names(meta_colors) <- c("Time", "Condition")

# Top X species to keep
subset_top_taxa = 20

# Agréger les données au niveau du genre
physeq_genus <- tax_glom(physeq_normalize, taxrank = "Genus")


p <- plot_taxa_heatmap(physeq_genus,
                       subset.top = subset_top_taxa,
                       VariableA = c("Time","Condition"),
                       heatcolors = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100), 
                       transformation = "log10",
                       cluster_rows = T,
                       cluster_cols = T,
                       show_colnames = T,
                       main = paste0("Heatmap with subset of top ", subset_top_taxa, " taxa"))



##########
# DESeq2 #
##########

# https://www.yanh.org/2021/01/01/microbiome-r/#differential-abundance-analysis

sample_data(pseq)$Time_Condition <- as.factor(sample_data(pseq)$Time_Condition) # factorize for DESeq2
pseq.taxa <- tax_glom(pseq, taxrank = 'Genus', NArm = FALSE)

# pairwise comparison between Treated and Untreated
pseq.taxa.sub <- subset_samples(pseq.taxa, Time_Condition %in% c("T1_Untreated", "T1_Treated"))

# filter sparse features, with > 90% zeros
pseq.taxa.pse.sub <- prune_taxa(rowSums(otu_table(pseq.taxa.sub) == 0) < ncol(otu_table(pseq.taxa.sub)) * 0.9, pseq.taxa.sub)
pseq_ds = phyloseq_to_deseq2(pseq.taxa.pse.sub, ~ Time_Condition)

# use alternative estimator on a condition of "every gene contains a sample with a zero"
ds <- estimateSizeFactors(pseq_ds, type="poscounts")
ds = DESeq(ds, test="Wald", fitType="parametric")
alpha = 0.05 
res = results(ds, alpha=alpha)
res = res[order(res$padj, na.last=NA), ]
taxa_sig = rownames(res[1:15, ]) # select bottom 20 with lowest p.adj values

pseq.taxa.rel <- transform_sample_counts(pseq, function(x) x/sum(x)*100)
pseq.taxa.rel.sig <- prune_taxa(taxa_sig, pseq.taxa.rel)

# Only keep gut and tongue samples
pseq.taxa.rel.sig <- prune_samples(colnames(otu_table(pseq.taxa.pse.sub)), pseq.taxa.rel.sig)

## Visualisation des 20##

matrix <- as.matrix(data.frame(otu_table(pseq.taxa.rel.sig)))
rownames(matrix) <- as.character(tax_table(pseq.taxa.rel.sig)[, "Genus"])
metadata_sub <- data.frame(sample_data(pseq.taxa.rel.sig))
# Define the annotation color for columns and rows
annotation_col = data.frame(
  Condition = as.factor(metadata_sub$Condition), 
  Breed = as.factor(metadata_sub$Breed), 
  check.names = FALSE
)
rownames(annotation_col) = rownames(metadata_sub)

annotation_row = data.frame(
  Phylum = as.factor(tax_table(pseq.taxa.rel.sig)[, "Phylum"])
)
rownames(annotation_row) = rownames(matrix)

ComplexHeatmap::pheatmap(matrix, scale= "row", 
                         annotation_col = annotation_col, 
                         annotation_row = annotation_row)

## Pour modifier les couleurs
# ann_color should be named vectors
#phylum_col = RColorBrewer::brewer.pal(length(levels(annotation_row$Phylum)), "Paired")
#names(phylum_col) = levels(annotation_row$Phylum)
#ann_colors = list(
#  Condition = c(Untreated = "green", Treated = "red"),
#  Phylum = phylum_col
#)

#ComplexHeatmap::pheatmap(matrix, scale= "row", 
#                         annotation_col = annotation_col, 
#                         annotation_row = annotation_row, 
#                         annotation_colors = ann_colors)




  #########################
## Sur la condition traité ##
  #########################

# Vérifier les colonnes disponibles dans les métadonnées
sample_data(pseq)


# taxaglom
pseq.taxa_Genus <- tax_glom(pseq, taxrank = 'Genus', NArm = FALSE)
pseq.taxa_Family <- tax_glom(pseq, taxrank = 'Family', NArm = FALSE)

# Variables
pseq_taxa = pseq.taxa_Genus
taxa_rank = "Genus"         ### WARNING : il faut changer le code
Cond_choisie = "Treated"
Time_selected <- c("T2", "T3")
Time_selected_name_fig = "T2 vs T3"

pseq_subset <- subset_samples(pseq_taxa, Condition == Cond_choisie)
pseq_subset <- subset_samples(pseq_subset, Time %in% Time_selected)

# Vérif
unique(pseq_subset@sam_data$Time_Condition)


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
# Family order
x = tapply(sigtab$log2FoldChange, sigtab$Family, function(x) max(x))
x = sort(x, TRUE)
sigtab$Family = factor(as.character(sigtab$Family), levels=names(x))
write.xlsx(sigtab, paste0(Cond_choisie, "Fixed", Time_selected_name_fig, "_alpha=", alpha, "_at_taxa_rank_", taxa_rank,".xlsx"))

ggplot(sigtab, aes(x=Family, y=log2FoldChange, color=Phylum)) + geom_point(size=2) + 
  theme(axis.text.x = element_text(angle = 0, hjust = 0, vjust=0.5, size = 8))  +
  theme(axis.text.y = element_text(size = 8))  +
  coord_flip() +
  ggtitle(paste0(Cond_choisie, " fixed / ", Time_selected_name_fig, " / alpha = ", alpha, " at taxa rank : ", taxa_rank))

# Genus order
x = tapply(sigtab$log2FoldChange, sigtab$Genus, function(x) max(x))
x = sort(x, TRUE)
sigtab$Genus = factor(as.character(sigtab$Genus), levels=names(x))

write.xlsx(sigtab, paste0(Cond_choisie, "Fixed", Time_selected_name_fig, "_alpha=", alpha, "_at_taxa_rank_", taxa_rank,".xlsx"))

ggplot(sigtab, aes(x=Genus, y=log2FoldChange, color=Phylum)) + geom_point(size=2) + 
  theme(axis.text.x = element_text(angle = 0, hjust = 0, vjust=0.5, size = 8))  +
  theme(axis.text.y = element_text(size = 8))  +
  coord_flip() +
  ggtitle(paste0(Cond_choisie, " fixed / ", Time_selected_name_fig, " / alpha = ", alpha, " at taxa rank : ", taxa_rank))

# verifier les niveau des conditions 
levels(Condsdds$Time) # Supposons que cela renvoie ["A", "B"], où "A" est la condition de référence.
head(sigtab)


# Si log2FoldChange est négatif pour un taxon, cela signifie que ce taxon est moins abondant dans la condition d'intérêt (B) par rapport à la condition de référence (A).
# Si log2FoldChange est positif, le taxon est plus abondant dans la condition d'intérêt (B) par rapport à la condition de référence (A).



  #############
## Sur le Time ##
  #############

# taxaglom
pseq.taxa_Genus <- tax_glom(pseq, taxrank = 'Genus', NArm = FALSE)
pseq.taxa_Family <- tax_glom(pseq, taxrank = 'Family', NArm = FALSE)

# Variables
pseq_taxa = pseq.taxa_Genus
taxa_rank = "Genus"
Time_selected = "T1"
setwd(dir = "K:/Personnel_non_permanent/Thibault/12_metagenomique16_minion_Morgane/05_Epi2me_250211/DESeq2")

# Vérifier les colonnes disponibles dans les métadonnées
sample_data(pseq)

# Sous-échantillonner pour ne conserver que les échantillons avec le temps "X"
pseq_subset <- subset_samples(pseq_taxa, Time %in% Time_selected)

# Vérifier les échantillons sélectionnés
sample_data(pseq_subset)

# Phyloseq to DESeq2
Condsdds = phyloseq_to_deseq2(pseq_subset, ~ Condition)
Condsdds = DESeq(Condsdds, sfType = "poscounts" ,test="Wald", fitType="parametric")


res = results(Condsdds, cooksCutoff = FALSE)

resultsNames(Condsdds)

alpha = 0.01
sigtab = res[which(res$padj < alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(pseq_subset)[rownames(sigtab), ], "matrix"))
dim(sigtab)


# Family order
x = tapply(sigtab$log2FoldChange, sigtab$Family, function(x) max(x))
x = sort(x, TRUE)
sigtab$Family = factor(as.character(sigtab$Family), levels=names(x))

output_file <- paste0(Time_selected, "Fixed_TreatedvsUntreated_Taxa", taxa_rank, "_alpha=", alpha, ".xlsx")
write.xlsx(sigtab, file = output_file)

output_file <- paste0(Time_selected, "Fixed_TreatedvsUntreated_Taxa", taxa_rank, "_alpha=", alpha, ".png")
png(filename = output_file, width = 1200, height = 703)
print(
  ggplot(sigtab[sigtab$Family,], aes(x=Family, y=log2FoldChange, color=Phylum)) + geom_point(size=6) + 
    coord_flip() +
    theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5)) +
    ggtitle(paste0(Time_selected," fixed / Treated vs Untreated / taxa rank : ", taxa_rank, " alpha = ", alpha))
)
dev.off()


# Genus order
x = tapply(sigtab$log2FoldChange, sigtab$Genus, function(x) max(x))
x = sort(x, TRUE)
sigtab$Genus = factor(as.character(sigtab$Genus), levels=names(x))

output_file <- paste0(Time_selected, "Fixed_TreatedvsUntreated_Taxa", taxa_rank, "_alpha=", alpha, ".xlsx")
write.xlsx(sigtab, file = output_file)

output_file <- paste0(Time_selected, "Fixed_TreatedvsUntreated_Taxa", taxa_rank, "_alpha=", alpha, ".png")
png(filename = output_file, width = 1200, height = 703)
print(
  ggplot(sigtab[sigtab$Genus,], aes(x=Genus, y=log2FoldChange, color=Phylum)) + geom_point(size=6) + 
    coord_flip() +
    theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5)) +
    ggtitle(paste0(Time_selected," fixed / Treated vs Untreated / taxa rank : ", taxa_rank, " alpha = ", alpha))
)
dev.off()

### Avec filtre sur les top 15

## Rappel des top 15 family ##
top_15_family <- names(sort(tapply(taxa_sums(physeq_normalize), tax_table(physeq_normalize)[, "Family"], sum), decreasing = TRUE)[1:15])
top_15_family
theme_set(theme_bw())
scale_fill_discrete <- function(palname = "Set1", ...) {
  scale_fill_brewer(palette = palname, ...)
}

# Family order
x = tapply(sigtab$log2FoldChange, sigtab$Family, function(x) max(x))
x = sort(x, TRUE)
sigtab$Family = factor(as.character(sigtab$Family), levels=names(x))
ggplot(sigtab[sigtab$Family %in% top_15_family,], aes(x=Family, y=log2FoldChange, color=Phylum)) + geom_point(size=6) + 
  coord_flip() +
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5)) +
  ggtitle(paste0(Time_selected," fixed / Treated vs Untreated / taxa rank : ", taxa_rank, " alpha = ", alpha))


## Rappel des top 15 Genus ##
top_15_Genus <- names(sort(tapply(taxa_sums(physeq_normalize), tax_table(physeq_normalize)[, "Genus"], sum), decreasing = TRUE)[1:15])
top_15_Genus
theme_set(theme_bw())
scale_fill_discrete <- function(palname = "Set1", ...) {
  scale_fill_brewer(palette = palname, ...)
}

# Genus order
x = tapply(sigtab$log2FoldChange, sigtab$Genus, function(x) max(x))
x = sort(x, TRUE)
sigtab$Genus = factor(as.character(sigtab$Genus), levels=names(x))
ggplot(sigtab[sigtab$Genus %in% top_15_Genus,], aes(x=Genus, y=log2FoldChange, color=Phylum)) + geom_point(size=6) + 
  coord_flip() +
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5)) +
  ggtitle(paste0(Time_selected," fixed / Treated vs Untreated / taxa rank : ", taxa_rank, " alpha = ", alpha))

# verifier les niveau des conditions 
levels(Condsdds$Condition) # Supposons que cela renvoie ["A", "B"], où "A" est la condition de référence.
head(sigtab)

## Download sur excel le fichier sigtab pour chaque comparaison 
write_xlsx(sigtab, paste0())


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
# Si on utilise la méthode de rarefaction pour normaliser les données : la barre vertical indique ou sera coupé nos échantillons.
# On attend pas vraiment le plateau pour nos échantillons

# Plot les alpha div "Observed", "Chao1" et "Shannon" pour voir si on a un effet différents en fontion du nombre de reads.
plot_richness(pseq, measures = c("Observed", "Shannon", "Chao1")) + theme(axis.text.x = element_text(size = 5))



# Adonis #

# Charger les librairies nécessaires
library(phyloseq)
library(vegan)

# Calculer une matrice de distance (ex : Bray-Curtis)
dist_matrix <- phyloseq::distance(pseq, method = "bray")

# Si metadata est un objet phyloseq (sample_data), le convertir en data.frame
metadata <- as.data.frame(as.matrix(sample_data(pseq)))

adonis2(dist_matrix ~ Condition + Time, data = metadata, permutations = 999)



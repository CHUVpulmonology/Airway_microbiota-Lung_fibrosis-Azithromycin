# General information on the study and analyses ----

# This is the script used to perform the bacterial community analyses for the manuscript:
# "Azithromycin alters spatial and temporal dynamics of airway microbiota in
# idiopathic pulmonary fibrosis"

# Pieter-Jan Gijs (a), Cecile Daccord (a), Eric Bernasconi (a), Martin Brutsche (c), Christian Clarenbach (d),
# Katrin Hostettler (e), Sabina A. Guler (b), Louis Mercier (a), Niki Ubags (a),
# Manuela Funke-Chambour (b) and Christophe von Garnier (a) 

# Affiliations:
# (a) Division of Pulmonary Medicine, Department of Medicine, CHUV, Lausanne University Hospital, Lausanne, Switzerland
# (b) Department of Pulmonary Medicine, Inselspital, Bern University Hospital, Bern, Switzerland
# (c) Lung center, Kantonsspital St. Gallen, St. Gallen, Switzerland
# (d) Division of Pulmonary Medicine, University Hospital of Zurich, Zurich, Switzerland
# (e) Clinics of Respiratory Medicine, University Hospital Basel, Basel, Switzerland

# Script authors: Pieter-Jan Gijs and Eric Bernasconi

# Data type: 16S rRNA amplicon sequencing + clinical metadata
# Sample type: Sputum + Oropharyngeal swab
# Sequencing location: Lausanne Genomic Facility, 1015 Lausanne, Switwerland
# Date: December 13, 2022

# Loading packages ----

library(tidyverse) # Collection of packages for importing/exporting (readr),
# manipulating/analyzing (dplyr) and plotting (ggplot2) data
library(rlang) # Toolbox for working with core 'Tidyverse' features
library(plyr) # For facilitating data separation, processing and re-composition
library(magrittr) # For working with pipe operators
library(hablar) # For converting columns to new data types
library(lubridate) # For working with dates
library(ggpubr) # For data visualization
library(viridis) # Color maps designed to improve graph readability
# for readers with common forms of color blindness and/or color vision deficiency
library(RColorBrewer) # For working with ready-to-use color palettes
library(scales) # For controlling the appearance of axis labels and legends on graphics
library(ggrepel) # For avoiding overlapping text labels with ggplot2
library(phyloseq) # For analyzing microbiome profiles
library(vegan) # For studying community ecology (e.g. ordination, diversity analysis)
library(picante) # For measuring phylogenetic diversity
library(microbiome) # For analyzing microbiome profiles
library(microbiomeutilities) # Tools for the processing and visualisation of microbiome data
library(usedist) # For computing between-group distances
library(ranacapa) # For making rarefaction curves using ggplot2
library(ape) # For working with phylogenetic trees
library(ggordiplots) # For adding features to ordination graphs made with ggplot2
library(ggVennDiagram) # For group data visualization as Venn diagram 
library(rstatix) # Pipe-friendly framework for basic statistical tests
library(decontam) # For contaminant identification
library(flextable) # For tabular reporting


# Importing files ----
# Raw sequencing data, processed sequencing data (i.e. ASV table and taxonomy table) and
# secondary data are available at https://doi.org/10.5281/zenodo.7065053
# Raw sequencing data were processed using a dada2-based pipeline available at https://github.com/chuvpne/dada2-pipeline

# Opening the ASV table previously saved in R data format
# (available at https://github.com/CHUVpulmonology/Airway_microbiota-Lung_fibrosis-Azithromycin)
ASV <- readRDS(file = "asv_table.rds")
# Opening the taxonomy table previously saved in R data format
# (available at https://github.com/CHUVpulmonology/Airway_microbiota-Lung_fibrosis-Azithromycin)
TaxSpecies <- readRDS(file = "Taxspecies.rds")
# Opening the secondary data table previously saved in R data format
# (available at https://github.com/CHUVpulmonology/Airway_microbiota-Lung_fibrosis-Azithromycin)
SecDataMerge <- readRDS(file = "SecDataMerge.rds")

# Cleaning imported files in view of preparing a phyloseq object ----
# The number of samples and taxa and their names must match between the 3 file types
# Sample names in ASV table
ASV %>% colnames
# Sample names in Secondary data file
SecDataMerge$SampleID_PNE_Seq
# Shared Sample names to be kept
intersect(ASV %>% colnames, SecDataMerge$SampleID_PNE_Seq)
# Sample names in ASV but not in Secondary data file
SamplesToBeRemovedFromASV <- setdiff(ASV %>% colnames, SecDataMerge$SampleID_PNE_Seq)
# ASV table without samples to be removed
ASV2 <- ASV %>% select(-all_of(SamplesToBeRemovedFromASV))
# Sample names in Secondary data file but not in ASV2 
SamplesToBeRemovedFromSecData <- setdiff(SecDataMerge$SampleID_PNE_Seq, ASV2 %>% colnames)
# Secondary data without sample to be removed
SecDataMerge2 <- SecDataMerge %>% filter(!SampleID_PNE_Seq %in% SamplesToBeRemovedFromSecData)

# Fixing inconsistencies in secondary data.
# Renaming, reordering and converting variables.
# A warning message is expected because mutate() will modify variables containing NAs
SecDataMerge2 %>% str()
SecDataMerge3 <- SecDataMerge2 %>%
  select(-c(TTTgroupSynth, TTTgroupDet)) %>%
  dplyr::rename(SampleID = SampleID_PNE_Seq,
         PatientID = PatientID_Bern,
         Load = "16S_count_per_ul",
         TLC_percent = "TLC_%",
         FVC_percent = "FVC_%",
         FEV1_percent = "FEV1_%",
         FEV1_FVC_ratio_percent = "FEV1_FVC_ratio_%",
         DLCO_unco_percent = "DLCO_unco_%",
         DLCO_co_percent = "DLCO_co_%") %>% 
  mutate(OriginType = paste(Origin, Type, sep = "_"), 
         OriginType3 = case_when(Origin == "Patient" & Type == "Swab" ~ "Swab",
                                 Origin == "Patient" & Type == "Sputum" ~ "Sputum",
                                 Origin == "Control" ~ "Control")) %>%
  relocate(Origin, .before = Visit) %>% 
  relocate(Type, .after = Origin) %>% 
  relocate(OriginType, .after = Type) %>% 
  relocate(OriginType3, .after = OriginType) %>%
  convert(fct(SampleID:TTTPrePost7gr, DLCO_name),
          dbl(Load:DLCO_co_percent, tetW_crude_20200610:msrE_normal))

SecDataMerge3$tetW_normal <- round(SecDataMerge3$tetW_normal, digit = 4)
SecDataMerge3$mel_normal <- round(SecDataMerge3$mel_normal, digit = 4)
SecDataMerge3$tetM_normal <- round(SecDataMerge3$tetM_normal, digit = 4)
SecDataMerge3$ermB_normal <- round(SecDataMerge3$ermB_normal, digit = 4)
SecDataMerge3$ermF_normal <- round(SecDataMerge3$ermF_normal, digit = 4)
SecDataMerge3$mef_normal <- round(SecDataMerge3$mef_normal, digit = 4)
SecDataMerge3$msrE_normal <- round(SecDataMerge3$msrE_normal, digit = 4)

# Examining the order of levels
SecDataMerge3 %>% select(SampleID:TTTPrePost7gr, DLCO_name) %>% purrr::map(., levels)
# Reordering levels
SecDataMerge3 %<>%
  mutate(TTTPrePost4gr = factor(TTTPrePost4gr, levels = c("NegCtrl", "NegCTRL", "PreAZT", 
                                                          "StartAZT", "EndAZT"))) %>%
  mutate(TTTPrePost5gr = factor(TTTPrePost5gr, levels = c("NegCtrl", "NegCTRL", "PreAZT", 
                                                          "StartAZT", "EndAZT", 
                                                          "PostAZT_1mo", "PostAZT>=3mo"))) %>%
  mutate(TTTPrePost7gr = factor(TTTPrePost7gr, levels = c("NegCtrl", "NegCTRL", "PreAZT_4mo", 
                                                          "PreAZT_1mo", "StartAZT", "EndAZT", 
                                                          "PostAZT_1mo", "PostAZT_4mo", 
                                                          "PostAZT_5mo" ))) %>%
  mutate(DLCO_name = factor(DLCO_name, levels = c("normal", "mild", "moderate", "severe")))

# Converting to date format
SecDataMerge3 %<>% mutate(HarvestDDMMYY = lubridate::dmy(HarvestDDMMYY))

# Counting the number of samples in Controls, Sputum and Oropharyngeal swab (OPS) samples
SecDataMerge3 %>% dplyr::count(OriginType3)

# Preparing phyloseq objects ----
# 1) Preparing an ASV abundance table in a matrix format to create a phyloseq object
ASVmat <- as.matrix(ASV2)
class(ASVmat)

ASVps <- otu_table(ASVmat, taxa_are_rows = TRUE)
ASVps[1:5, 1:5]
sample_names(ASVps)[1:5]
taxa_names(ASVps)[1:5]

# 2) Preparing a taxonomy table in a matrix format to create a phyloseq object
TaxSpecies_mat <- as.matrix(TaxSpecies)

TaxSpeciesPS <- tax_table(TaxSpecies_mat)

# 3) Preparing a data frame with sample variables to create a phyloseq object
SecDataMerge3DF <- as.data.frame(SecDataMerge3)
# Checking that sample names are the same in the data frame of secondary data (as rownames) and 
# in the ASV table (as column names)
rownames(SecDataMerge3DF)
intersect(sample_names(ASVps), SecDataMerge3DF$SampleID)
rownames(SecDataMerge3DF) <- SecDataMerge3DF$SampleID
intersect(rownames(SecDataMerge3DF), sample_names(ASVps))
rownames(SecDataMerge3DF)

SecDataMerge3SD <- sample_data(SecDataMerge3DF)

# Merging 1) ASV abundance table, 2) Taxonomy table
# and 3) Sample variables in a phyloseq object
physeqATS <- phyloseq(ASVps, TaxSpeciesPS, SecDataMerge3SD)
# 4) Adding a phylogenetic tree to the phyloseq object
# Using ape::rtree to generate a phylogenetic tree (ape was loaded with GUniFrac)
tree1 <-  ape::rtree(ntaxa(physeqATS), rooted = TRUE, tip.label = taxa_names(physeqATS))
plot(tree1)
# Merging secondary data and tree
physeq1 <- merge_phyloseq(physeqATS, tree1)
summarize_phyloseq(physeq1)

# Controlling colours ----
# Colorblind-friendly palette of 8 colours
colorBlindGrey8   <- c("#999999", "#E69F00", "#56B4E9", "#009E73", 
                       "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
scales::show_col(colorBlindGrey8)

# Phyla
MyBactCol5BlFr <- c(Firmicutes = "#CC79A7", Fusobacteria = "#999999",
                    Actinobacteria = "#56B4E9", Bacteroidetes = "#E69F00",
                    Proteobacteria = "#009E73")

MyBactCol6BlFr <- c(Firmicutes = "#CC79A7", Fusobacteria = "#999999",
                    Actinobacteria = "#56B4E9", Bacteroidetes = "#E69F00",
                    Proteobacteria = "#009E73", Epsilonbacteraeota = "0072B2")

MyBactCol8BlFr <- c(Firmicutes = "#CC79A7", Fusobacteria = "#999999",
                    Actinobacteria = "#56B4E9", Bacteroidetes = "#E69F00",
                    Proteobacteria = "#009E73", Epsilonbacteraeota = "#0072B2",
                    Synergistetes = "#F0E442", Spirochaete = "#D55E00")

# Sample types
PatientCtrlCol2 <- c(Patient = "#D55E00", Control = "#009E73")

SampleTypeCol2 <- c(Sputum = "#f19046", Swab = "#2a9ddf")

SampleTypeCol3 <- c(Sputum = "#f0E442", Swab = "#2a9ddf", Control = "#a3a3a3")

# Treatment phases
StartEndTTTCol2 <- c(StartAZT = "#f19046", EndAZT = "#2a9ddf")

TTTphaseCol5 <- c(PreAZT = "#999999", StartAZT = "#E69F00", EndAZT = "#56B4E9", 
                  PostAZT_1mo = "#009E73", "PostAZT>=3mo" = "#CC79A7")

TTTphaseCol7 <- c(PreAZT_4mo = "#D55E00", PreAZT_1mo = "#999999", StartAZT = "#E69F00",
                  EndAZT = "#56B4E9", PostAZT_1mo = "#009E73", PostAZT_4mo = "#CC79A7",
                  PostAZT_5mo = "#0072B2")

# Resistance status
ResistStatusCol2 <- c(StableResist = "#009E73", IncreasedResist = "#D55E00")

# Lists of statistical comparisons to be performed ----
# Sample types
CompCtrlOROSwab <- list(c("Control", "Swab"))

CompCtrlSputum <- list(c("Control", "Sputum"))

# Treatment phases
CompTTT1 <- list(c("PreAZT", "StartAZT"), c("StartAZT", "EndAZT"), c("StartAZT", "PostAZT_1mo"),
                 c("StartAZT", "PostAZT>=3mo"))

CompTTT2 <- list(c("PreAZT", "StartAZT"), c("StartAZT", "EndAZT"), c("StartAZT", "PostAZT_1mo"),
                 c("StartAZT", "PostAZT>=3mo"), c("PostAZT_1mo", "PostAZT>=3mo"))

CompTTT3 <- list(c("PreAZT_4mo", "EndAZT"), c("PreAZT_1mo", "EndAZT"),
                 c("StartAZT", "EndAZT"), c("StartAZT", "PostAZT_1mo"),
                 c("StartAZT", "PostAZT_4mo"), c("StartAZT", "PostAZT_5mo"))

CompTTT4 <- list(c("PreAZT_4mo", "EndAZT"), c("PreAZT_1mo", "EndAZT"),
                 c("StartAZT", "EndAZT"), c("StartAZT", "PostAZT_1mo"),
                 c("EndAZT", "PostAZT_1mo"))

CompTTT5 <- list(c("StartAZT", "EndAZT"), c("StartAZT", "PostAZT_1mo"),
                 c("StartAZT", "PostAZT_4mo"), c("StartAZT", "PostAZT_5mo"))

CompTTT7 <- list(c("PreAZT_4mo", "EndAZT"), c("PreAZT_1mo", "EndAZT"),
                 c("StartAZT", "EndAZT"), c("StartAZT", "PostAZT_1mo"),
                 c("StartAZT", "PostAZT_4mo"), c("StartAZT", "PostAZT_5mo"),
                 c("EndAZT", "PostAZT_5mo"))

# Resistance status
CompResistStatus <- list(c("StableResist", "IncreasedResist"))

# Preliminary actions on the phyloseq object including filtering ----
# Replacing "sp" labels with "ASV"
taxa_names(physeq1) %<>% gsub("sp", "ASV", .)

# Removing taxa whose total number of reads did not exceed 1
physeqNoLowRead <- prune_taxa(taxa_sums(physeq1) > 1, physeq1)
# Number of taxa removed during abundance filtering of absolute read counts 
physeq1 %>% ntaxa() - physeqNoLowRead %>% ntaxa()

# Keeping only ASVs assigned to Bacteria
physeqBact <- subset_taxa(physeqNoLowRead, Kingdom == "Bacteria")
# Number of taxa removed during "taxonomy assigment filtering"
physeqNoLowRead %>% ntaxa() - physeqBact %>% ntaxa()

physeqBact %>% sample_sums()

sample_data(physeqBact)$ReadDepth <- sample_sums(physeqBact)

# Supplementary Figure E4
LoadSeqDepth_Dotplot <- meta(physeqBact) %>% 
  ggplot(aes(x = ReadDepth,
             y = Load,
             fill = Origin)) +
  geom_smooth(aes(fill = Origin),
              method = lm,
              color = "white",
              alpha = .3) +
  geom_point(shape = 21,
             size = 1.5,
             alpha = .4) +
  scale_x_continuous(trans = "log10",
                     labels = scales::scientific) +
  scale_y_continuous(trans = "log10",
                     labels = scales::scientific) +
  scale_fill_manual(values = PatientCtrlCol2) +
  stat_cor() +
  theme_classic()

ggsave("LoadSeqDepth_Dotplot.pdf",
       width = 12, height = 10, dpi = 200, units = "cm")

# Preparing a boxplot for the top 15 ASVs in Controls (Supplementary Figure E7, panel c) ----
# with separate display for Controls and Patient samples.
# Working with relative read counts
physeqBactRel <- transform_sample_counts(physeqBact, 
                                         function(x) {return(100 * x / sum(x))})

physeqBactRelControls <- subset_samples(physeqBactRel, OriginType3 == "Control")

# Sorting 1:15 top abundant ASVs in Controls in descending order
TopAb15ASVControls <- sort(taxa_sums(physeqBactRelControls), decreasing = TRUE)[1:15]
# Computing the names of the most abundant ASVs in Controls
TopAb15ASVControlName  <-  names(TopAb15ASVControls)
# Keeping samples from Controls AND Patients but only the top 15 most abundant ASVs in Controls
physeqTopAb15ControlTaxa_AllSamples  <-  prune_taxa(TopAb15ASVControlName, physeqBactRel)
# Obtaining a data frame from phyloseq object
TopAb15ControlTaxa_AllSamplesDF <- psmelt(physeqTopAb15ControlTaxa_AllSamples)

TopAb15ControlTaxa_AllSamplesDF2 <- TopAb15ControlTaxa_AllSamplesDF %>%
  mutate(ASV = OTU) %>% 
  unite(ASV_Genus_Species, c(OTU, Genus, Species), sep ="_") %>%
  convert(fct(ASV, ASV_Genus_Species, SampleID, Origin, Phylum)) %>%
  dplyr::rename(RelAbundance = Abundance) %>%
  arrange(desc(RelAbundance)) %>%
  group_by(ASV_Genus_Species, Phylum, RelAbundance)

# Ranking the 15 most abundant ASVs in Controls in descending order of medians
physeqBactRelControlsTopAb15ASVControl <- prune_taxa(TopAb15ASVControlName, physeqBactRelControls)

psMeltControlsTopAb15ASVControlDF <- psmelt(physeqBactRelControlsTopAb15ASVControl)

ControlsTopAb15ASVControlGroupedDF <- psMeltControlsTopAb15ASVControlDF %>%
  mutate(ASV = OTU) %>%
  unite(ASV_Genus_Species, c(ASV, Genus, Species), sep ="_") %>%
  convert(fct(ASV_Genus_Species, Phylum)) %>%
  dplyr::group_by(ASV_Genus_Species) %>%
  dplyr::summarise(MedianAbundance = median(Abundance)) %>%
  dplyr::arrange(desc(MedianAbundance))

ControlsTopAb15ASVControlGroupedDF$ASV_Genus_Species %>% levels()
Top15ControlTaxaNamesDescMed <- ControlsTopAb15ASVControlGroupedDF$ASV_Genus_Species

# Boxplot for the 15 most abundant ASVs in Controls,
# displayed for Controls and Patient samples separately (Supplementary Figure S7, panel C)
Top15ControlTax_boxplot <- TopAb15ControlTaxa_AllSamplesDF2 %>%
  ggplot(aes(x = factor(ASV_Genus_Species, level = Top15ControlTaxaNamesDescMed), 
             y = RelAbundance,
             fill = Phylum)) +
  geom_boxplot(width = .7,
               size = .3,
               alpha = .4,
               position = position_dodge(.8),
               outlier.shape = NA,
               show.legend = TRUE) +
  geom_jitter(position = position_jitterdodge(.8),
              shape = 21,
              colour = "black",
              size = 1,
              alpha = .7) +
  ylim(0, 60) +
  scale_fill_manual(values = MyBactCol5BlFr) +
  facet_wrap(~factor(Origin, levels = c("Control", "Patient")),
             nrow = 2) +
  labs(x = "", y = "Relative abundance") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  theme(strip.background = element_blank())

ggsave("Top15ControlTax_boxplot.pdf",
       width = 15, height = 18, dpi = 200, units = "cm")

# Preparing a boxplot for the top 15 ASVs in Patient samples (Supplementary Figure E7, panel c) ----
# with separate display for Patient samples and Controls.
physeqBactRelPatients <- subset_samples(physeqBactRel, Origin == "Patient")

# Removing samples with zero read count
physeqBactRelPatients <- prune_samples(sample_sums(physeqBactRelPatients) > 0, physeqBactRelPatients)
# Removing the taxa that are no longer represented after filtering
physeqBactRelPatients2 <- prune_taxa(taxa_sums(physeqBactRelPatients) > 0, physeqBactRelPatients)
# Sorting 1:15 top abundant ASVs in Patient samples in descending order
TopAb15ASVPatients <- sort(taxa_sums(physeqBactRelPatients2), decreasing = TRUE)[1:15]
# Computing the names of the most abundant ASVs in Patient samples
TopAb15ASVPatientsName  <-  names(TopAb15ASVPatients)
# Keeping samples from Patient samples AND Controls but only the top 15 most abundant ASVs in Patient samples
physeqTopAb15PatientsTaxa_AllSamples  <-  prune_taxa(TopAb15ASVPatientsName,
                                                     physeqBactRel)
# Obtaining a data frame from phyloseq object
TopAb15PatientsTaxa_AllSamplesDF <- psmelt(physeqTopAb15PatientsTaxa_AllSamples)

TopAb15PatientsTaxa_AllSamplesDF2 <- TopAb15PatientsTaxa_AllSamplesDF %>%
  mutate(ASV = OTU) %>% 
  unite(ASV_Genus_Species, c(OTU, Genus, Species), sep ="_") %>%
  convert(fct(ASV, ASV_Genus_Species, SampleID, Origin, Phylum)) %>%
  dplyr::rename(RelAbundance = Abundance) %>%
  arrange(desc(RelAbundance)) %>%
  group_by(ASV_Genus_Species, Phylum, RelAbundance)

# Ranking the 15 most abundant ASVs in Patient samples in descending order of medians
physeqBactRelPatientsTopAb15ASVPatients <- prune_taxa(TopAb15ASVPatientsName, physeqBactRelPatients)

psMeltPatientsTopAb15ASVPatientsDF <- psmelt(physeqBactRelPatientsTopAb15ASVPatients)

PatientsTopAb15ASVPatientsGroupedDF <- psMeltPatientsTopAb15ASVPatientsDF %>%
  mutate(ASV = OTU) %>%
  unite(ASV_Genus_Species, c(ASV, Genus, Species), sep ="_") %>%
  convert(fct(ASV_Genus_Species, Phylum)) %>%
  dplyr::group_by(ASV_Genus_Species) %>%
  dplyr::summarise(MedianAbundance = median(Abundance)) %>%
  dplyr::arrange(desc(MedianAbundance))

PatientsTopAb15ASVPatientsGroupedDF$ASV_Genus_Species %>% levels()
Top15PatientsTaxaNamesDescMed <- PatientsTopAb15ASVPatientsGroupedDF$ASV_Genus_Species

# Boxplot for the 15 most abundant ASVs in Patient samples
# displayed for Patient samples and Controls separately (Supplementary Figure E7, panel c)
Top15PatientsTax_boxplot <- TopAb15PatientsTaxa_AllSamplesDF2 %>%
  ggplot(aes(x = factor(ASV_Genus_Species, level = Top15PatientsTaxaNamesDescMed), 
             y = RelAbundance,
             fill = Phylum)) +
  geom_boxplot(width = .7,
               size = .3,
               alpha = .4,
               position = position_dodge(.8),
               outlier.shape = NA,
               show.legend = TRUE) +
  geom_jitter(position = position_jitterdodge(.8),
              shape = 21,
              colour = "black",
              size = 1,
              alpha = .7) +
  ylim(0, 60) +
  scale_fill_manual(values = MyBactCol5BlFr) +
  facet_wrap(~factor(Origin, levels = c("Patient", "Control")),
             nrow = 2) +
  labs(x = "", y = "Relative abundance") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  theme(strip.background = element_blank())

ggsave("Top15PatientsTax_boxplot.pdf",
       width = 15, height = 15, dpi = 200, units = "cm")

# Analysis of contaminants (Supplementary Figure E2) ----
# Frequency-based analysis
# Frequency distribution of each ASV as a function of 16S rRNA gene copy number
physeqBactForContam <- physeqBact %>% prune_samples(sample_sums(.) > 0, .)

Contam_Freq_physeqBactDF <- decontam::isContaminant(physeqBactForContam, 
                                                    method = "frequency",
                                                    threshold = 0.1,
                                                    conc = "Load")

# Number of contaminants based on frequency
table(Contam_Freq_physeqBactDF$contaminant)

# Taxonomy table of frequency-based contaminants
Contam_Freq_physeqBactDF2 <- Contam_Freq_physeqBactDF %>% 
  rownames_to_column("ASV")

FreqContam <- Contam_Freq_physeqBactDF2 %>%
  filter(contaminant == "TRUE") %>%
  arrange(p) %>%
  pull(ASV)

FreqContamTaxDF <- tax_table(physeqBactForContam)[FreqContam[1:15], c("Phylum", "Family", "Genus", "Species")] %>% 
  as.data.frame()

FreqContamTaxDF[is.na(FreqContamTaxDF)] <- ""
FreqContamTaxDF2  <- FreqContamTaxDF %>% 
  rownames_to_column("ASV")

# Table E3
FreqContamTable <- flextable(FreqContamTaxDF2)
save_as_docx("Frequency-based contaminants" = FreqContamTable,
             path = "FrequencyContamTaxa.docx")

# Supplementary Figure E2, panel b
FrequencyContamAllPlot <- plot_frequency(physeqBactForContam,
               taxa_names(physeqBactForContam)[FreqContamIndex],
               conc = "Load") + 
  xlab("16S rRNA gene copies") +
  theme(axis.text.x = element_text(angle = 60, hjust = 1))

ggsave("FrequencyContamAllPlot.pdf",
       width = 30, height = 30, dpi = 200, units = "cm")

# Prevalence-based analysis
# Prevalence (presence/absence across samples) of each ASV in Patients samples versus Controls
sample_data(physeqBactForContam)$IsCtrl <- sample_data(physeqBactForContam)$OriginType3 == "Control"
Contam_Prev_physeqBact <- isContaminant(physeqBactForContam,
                                        method = "prevalence", 
                                        neg = "IsCtrl")

table(Contam_Prev_physeqBact$contaminant)

# Phyloseq object of presence-absence in negative controls and patient samples
physeqBactForContam2 <- physeqBactForContam
# Replacing positive read counts with 1
otu_table(physeqBactForContam2)[otu_table(physeqBactForContam) > 0] <- 1

Presence_Absence_CrtlPS <- prune_samples(sample_data(physeqBactForContam2)$Origin == "Control",
                                         physeqBactForContam2)

Presence_Absence_PatientPS <- prune_samples(sample_data(physeqBactForContam2)$Origin == "Patient",
                                            physeqBactForContam2)

# Data frame of prevalence in Patients samples and Controls
Presence_Absence_DF <- data.frame(PresAbsPat = taxa_sums(Presence_Absence_PatientPS), 
                                  PresAbsCtrl = taxa_sums(Presence_Absence_CrtlPS),
                                  Contaminant = Contam_Prev_physeqBact$contaminant) %>% 
  rownames_to_column("ASV")

# Supplementary Figure E2, panel a
PrevalencePlot <- Presence_Absence_DF %>% 
  mutate(HighPrevPat = ifelse(PresAbsPat > 100, ASV, ""),
         HighPrevCtrl = ifelse(PresAbsCtrl > 5, ASV, "")) %>% 
  ggplot(aes(x = PresAbsCtrl,
             y = PresAbsPat,
             fill = Contaminant)) +
  geom_point(shape = 21,
             color = "black",
             size = 2,
             alpha = .6) +
  labs(x = "Prevalence (Controls)",
       y = "Prevalence (Patient samples)") +
  scale_fill_manual(values = c("#D55E00", "#009E73")) +
  geom_text_repel(aes(label = HighPrevPat), size = 2.5, max.overlaps = Inf) +
  geom_text_repel(aes(label = HighPrevCtrl), size = 2.5, max.overlaps = Inf) +
  theme_classic()

ggsave("PrevalenceContam_Plot.pdf",
       width = 12, height = 12, dpi = 200, units = "cm")

# Taxonomy table of prevalence-based contaminants
PrevContam <- Presence_Absence_DF %>%
  filter(Contaminant == "TRUE") %>% 
  arrange(desc(PresAbsCtrl)) %>% 
  pull(ASV)

PrevContamTaxDF <- tax_table(physeqBactForContam2)[PrevContam[1:15], c("Phylum", "Family", "Genus", "Species")] %>% 
  as.data.frame()

PrevContamTaxDF[is.na(PrevContamTaxDF)] <- ""
PrevContamTaxDF2  <- PrevContamTaxDF %>% 
  rownames_to_column("ASV")

# Table E2
PrevContamTable <- flextable(PrevContamTaxDF2)
save_as_docx("Prevalence-based contaminants" = PrevContamTable,
             path = "PrevContamTaxa.docx")

# Taxonomy table of most prevalent ASVs in patient samples
PrevPatients <- Presence_Absence_DF %>%
  filter(Contaminant == FALSE) %>%
  arrange(desc(PresAbsPat)) %>% 
  pull(ASV)

PrevPatientTaxDF <- tax_table(physeqBactForContam2)[PrevPatients[1:15], c("Phylum", "Family", "Genus", "Species")] %>% 
  as.data.frame()

PrevPatientTaxDF[is.na(PrevPatientTaxDF)] <- ""
PrevPatientTaxDF2  <- PrevPatientTaxDF %>% 
  rownames_to_column("ASV")

PrevPatientTable <- flextable(PrevPatientTaxDF2)
save_as_docx("Prevalence-based patient taxa" = PrevPatientTable,
             path = "PrevPatientTaxa.docx")

# ASVs identified as contaminants by Decontam's "frequency" and/or "prevalence" methods
FreqAndPrevContam <- intersect(FreqContam, PrevContam)

setdiff(FreqContam, PrevContam) %>% length()
setdiff(PrevContam, FreqContam) %>% length()

FreqOrPrevContam <- c(FreqContam, PrevContam)
FreqOrPrevContam %>% unique() %>% length()

# Filtering of ASVs identified as contaminants by Decontam's prevalence and/or frequency methods
physeqBactDecontam <- prune_taxa(!(taxa_names(physeqBact) %in% FreqOrPrevContam), physeqBact)

# Comparison of bacterial load between Controls and Patient samples (Supplementary Figure E7, panels a and b) ----
# Controls vs. OPS (Supplementary Figure E7, panel b)
CtrlOROSwabLoad_boxplot <- meta(physeqBact) %>%
  filter(OriginType3 != "Sputum") %>% 
  ggplot(aes(x = OriginType3,
             y = Load,
             fill = OriginType3)) +
  geom_boxplot(width = .7,
               size = .5,
               alpha = .3,
               position = position_dodge(.7),
               outlier.shape = NA,
               show.legend = TRUE) +
  geom_jitter(shape = 21,
              colour = "black",
              size = 1,
              alpha = .6,
              width = .2) +
  scale_y_continuous(breaks = c(10^2, 10^4, 10^6, 10^8), 
                     limits = c(10^1, 10^8), 
                     trans = "log10") +
  scale_fill_manual(values = c(Control = "#F0E442", Swab = "#2a9ddf")) +
  stat_compare_means(label.y = 8) + # Adding Wilcoxon test result
  labs(x = "", y = "16S rRNA gene copies\n per oropharyngeal swab") +
  theme_classic() +
  theme(strip.background = element_blank(),
        axis.text.x = element_blank())

ggsave("CtrlOROSwabLoad_boxplot.pdf",
       width = 8, height = 10, dpi = 200, units = "cm")

# Two sample Wilcoxon test
meta(physeqBact) %>%
  filter(OriginType3 != "Sputum") %>%
  droplevels() %>%
  rstatix::wilcox_test(data =., Load ~ OriginType3) %>%
  adjust_pvalue(method = "holm") %>%
  add_significance("p.adj")

# Controls vs. Sputum (Supplementary Figure E7, panel a)
CtrlSputumLoad_boxplot <- meta(physeqBact) %>%
  filter(OriginType3 != "Swab") %>% 
  ggplot(aes(x = factor(OriginType3, levels = c("Control", "Sputum")),
             y = Load,
             fill = OriginType3)) +
  geom_boxplot(width = .7,
               size = .5,
               alpha = .3,
               position = position_dodge(.7),
               outlier.shape = NA,
               show.legend = TRUE) +
  geom_jitter(shape = 21,
              colour = "black",
              size = 1,
              alpha = .6,
              width = .2) +
  scale_y_continuous(breaks = c(10^2, 10^4, 10^6, 10^8), 
                     limits = c(10^1, 10^8), 
                     trans = "log10") +
  scale_fill_manual(values = c(Control = "#f0E442", Sputum = "#f19046")) +
  stat_compare_means(label.y = 8) + # Adding Wilcoxon test result
  labs(x = "", y = "16S rRNA gene copies\n per microlitre of sputum") +
  theme_classic() +
  theme(strip.background = element_blank(),
        axis.text.x = element_blank())

ggsave("CtrlSputumLoad_boxplot.pdf",
       width = 8, height = 10, dpi = 200, units = "cm")

# Two sample Wilcoxon test
meta(physeqBact) %>%
  filter(OriginType3 != "Swab") %>%
  droplevels() %>%
  rstatix::wilcox_test(data =., Load ~ OriginType3) %>%
  adjust_pvalue(method = "holm") %>%
  add_significance("p.adj")

# Summary statistics
meta(physeqBact) %>%
  dplyr::group_by(OriginType3) %>%
  dplyr::summarize(MedianLoad = median(Load),
                   IQRLoad = IQR(Load))

meta(physeqBact) %>%
  filter(OriginType3 == "Control") %>%
  pull(Load) %>% 
  quantile() %>% 
  format(scientific = TRUE)

meta(physeqBact) %>%
  filter(OriginType3 == "Sputum") %>% 
  pull(Load) %>% 
  quantile() %>% 
  format(scientific = TRUE)

meta(physeqBact) %>% 
  filter(OriginType3 == "Swab") %>% 
  pull(Load) %>% 
  quantile() %>% 
  format(scientific = TRUE)

# Lollipop chart of load in Samples and Controls
Load_Lolli <- meta(physeqBact) %>% 
  ggpubr::ggdotchart(x = "SampleID",
                     y = "Load",
                     color = "OriginType3",
                     palette = c(Sputum = "#D55E00",
                                 Swab = "#56B4E9", 
                                 Control = "#F0E442"),
                     sorting = "descending",
                     dot.size = 1.5,
                     add = "segments",
                     xlab = "Samples",
                     ylab = "16S rRNA gene copies",
                     ggtheme = theme_pubr()) +
  rremove(c("x.text")) +
  rremove(c("x.ticks")) +
  yscale("log10")

ggsave("Load_Lolli_20221108.pdf",
       width = 30, height = 10, dpi = 200, units = "cm")

# Beta diversity analysis in Samples and Controls BEFORE Decontam filtering
physeqBact2 <- physeqBact %>% prune_samples(sample_sums(.) > 0, .)

physeqBrayOrd <- physeqBact2 %>% ordinate(., 
                                          "PCoA", 
                                          "bray")
physeqBact2 %>% plot_ordination(., 
                                physeqBrayOrd,
                                type = "samples", 
                                color = "OriginType3")

# Plotting ellipses and spiders
PCoABrayOriginType3Plot1 <- gg_ordiplot(physeqBrayOrd$vectors, 
                                        groups = sample_data(physeqBact2)$OriginType3,
                                        spiders = TRUE, 
                                        ellipse = TRUE,
                                        pt.size = 1.5, 
                                        plot = TRUE)

# Modifying the plot using ggplot2
PCoABrayOriginType3Plot2 <- PCoABrayOriginType3Plot1$plot +
  geom_segment(data = PCoABrayOriginType3Plot1$df_spiders, 
               aes(x = cntr.x, 
                   xend = x, 
                   y = cntr.y, 
                   yend = y, 
                   color = Group), 
               size = .1,
               show.legend = FALSE) +
  geom_label(data = PCoABrayOriginType3Plot1$df_mean.ord, 
             aes(x = x, 
                 y = y, 
                 label = Group),
             colour = "black",
             size = 2.7,
             label.padding = unit(.2, "lines"),
             label.r = unit(0, "lines"),
             label.size = .1,
             show.legend = FALSE) +
  scale_color_manual(values = c(Sputum = "#D55E00",
                                Swab = "#56B4E9", 
                                Control = "#F0E442")) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  theme(legend.position = "none") +
  labs(x = "PC1", y = "PC2") +
  theme(axis.title.x = element_text(size = 10),
        axis.title.y = element_text(size = 10),
        axis.text.x = element_text(size = 8),
        axis.text.y = element_text(size = 8))

ggsave(filename = "PCoABrayOriginType3Plot2.pdf",
       width = 20, height = 15, dpi = 200, units = "cm", device='pdf')

# Performing PERMANOVA to test dissimilarity between the three groups (i.e. different centroids).
# The null hypothesis is that the sample groups have the same centroid.
# Groups to be compared
sample_data(physeqBact2)$OriginType3
# Seed was set to ensure reproducibility of the analysis
set.seed(9284)
# Bray distance matrix
physeqBactBray <- phyloseq::distance(physeqBact2, method = "bray")
# Adonis (column R2 indicates how much of the variance can be explained by the groups)
SputumSwabCtrl_Bray_adonis2 <- vegan::adonis2(physeqBactBray ~ OriginType3, 
                                              data = meta(physeqBact2))

# Computing the average distance of group members to the group centroid 
dist_Bray_OriginType3_betadisper <- betadisper(physeqBactBray, 
                                               meta(physeqBact2)$OriginType3, 
                                               type = c("median","centroid"))
# Performing an ANOVA to compare group dispersions. A non-significant ANOVA p-value
# would mean that the dispersions of the groups are homogeneous.
anova(dist_Bray_OriginType3_betadisper)
# Alternative: Run a permutation test for homogeneity of multivariate dispersions
permutest(dist_Bray_OriginType3_betadisper)
# Non-significant permutest results mean the null hypothesis that the groups 
# have the same dispersions cannot be rejected.
# This means we can be more confident that the subsequent adonis result is real,
# and not due to differences in group dispersions.

# Where appropriate, performing a Tukey's Honest test to assess which groups differ in their variances
# When comparing multiple groups, this test returns which pairs of comparisons are significant
TukeyHSD(dist_Bray_OriginType3_betadisper)

# Beta diversity analysis in Samples and Controls AFTER Decontam filtering
physeqBactDecontam2 <- physeqBactDecontam %>% prune_samples(sample_sums(.) > 0, .) %>% 
  prune_samples(sample_sums(.) > 0, .)

physeqDecontamBrayOrd <- physeqBactDecontam2 %>% ordinate(., 
                                                  "PCoA", 
                                                  "bray")
physeqBactDecontam2 %>% plot_ordination(., 
                                        physeqDecontamBrayOrd,
                                        type = "samples", 
                                        color = "OriginType3")

# Plotting ellipses and spiders
PCoABrayOriginType3DecontamPlot1 <- gg_ordiplot(physeqDecontamBrayOrd$vectors, 
                                                groups = sample_data(physeqBactDecontam2)$OriginType3,
                                                spiders = TRUE, 
                                                ellipse = TRUE,
                                                pt.size = 1.5, 
                                                plot = TRUE)

# Modifying the plot using ggplot2
PCoABrayOriginType3DecontamPlot2 <- PCoABrayOriginType3DecontamPlot1$plot +
  geom_segment(data = PCoABrayOriginType3DecontamPlot1$df_spiders, 
               aes(x = cntr.x, 
                   xend = x, 
                   y = cntr.y, 
                   yend = y, 
                   color = Group), 
               size = .1,
               show.legend = FALSE) +
  geom_label(data = PCoABrayOriginType3DecontamPlot1$df_mean.ord, 
             aes(x = x, 
                 y = y, 
                 label = Group),
             colour = "black",
             size = 2.7,
             label.padding = unit(.2, "lines"),
             label.r = unit(0, "lines"),
             label.size = .1,
             show.legend = FALSE) +
  scale_color_manual(values = c(Sputum = "#D55E00",
                                Swab = "#56B4E9", 
                                Control = "#F0E442")) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  theme(legend.position = "none") +
  labs(x = "PC1", y = "PC2") +
  theme(axis.title.x = element_text(size = 10),
        axis.title.y = element_text(size = 10),
        axis.text.x = element_text(size = 8),
        axis.text.y = element_text(size = 8))

ggsave(filename = "PCoABrayOriginType3DecontamPlot2.pdf",
       width = 20, height = 15, dpi = 200, units = "cm", device='pdf')

# Performing PERMANOVA to test dissimilarity between the three groups (i.e. different centroids).
# The null hypothesis is that the sample groups have the same centroid.
# Groups to be compared
sample_data(physeqBactDecontam2)$OriginType3
# Seed was set to ensure reproducibility of the analysis
set.seed(9834)
# Bray distance matrix
physeqBactDecontamBray <- phyloseq::distance(physeqBactDecontam2, method = "bray")
# Adonis (column R2 indicates how much of the variance can be explained by the groups)
SputumSwabCtrl_Bray_adonis2 <- vegan::adonis2(physeqBactDecontamBray ~ OriginType3, 
                                              data = meta(physeqBactDecontam2))

# Computing the average distance of group members to the group centroid 
dist_Bray_OriginType3_Decontam_betadisper <- betadisper(physeqBactDecontamBray, 
                                                        meta(physeqBactDecontam2)$OriginType3, 
                                                        type = c("median","centroid"))
# Performing an ANOVA to compare group dispersions. A non-significant ANOVA p-value
# would mean that the dispersions of the groups are homogeneous.
anova(dist_Bray_OriginType3_Decontam_betadisper)
# Alternative: Run a permutation test for homogeneity of multivariate dispersions
permutest(dist_Bray_OriginType3_Decontam_betadisper)
# Non-significant permutest results mean the null hypothesis that the groups 
# have the same dispersions cannot be rejected.
# This means we can be more confident that the subsequent adonis result is real,
# and not due to differences in group dispersions.

# Where appropriate, performing a Tukey's Honest test to assess which groups differ in their variances
# When comparing multiple groups, this test returns which pairs of comparisons are significant
TukeyHSD(dist_Bray_OriginType3_Decontam_betadisper)

# Computing, inspecting and plotting sequencing depth per sample (Supplementary Figure E3) ----
# after initial filtering
SampleSumsAll <- sample_sums(physeqBactDecontam)
# Displaying sequencing depth as a tibble
SampleSumsAllTib <- tibble(SampleID = names(SampleSumsAll),
                           ReadsPerSample = unname(SampleSumsAll))

# Adding Read count per sample in phyloseq object
identical(SampleSumsAllTib$SampleID, sample_names(physeqBactDecontam))

sample_data(physeqBactDecontam)$ReadsPerSample <- SampleSumsAllTib$ReadsPerSample
# Inspecting read counts
SampleData1DF <- sample_data(physeqBactDecontam) %>% 
  as.matrix() %>% 
  as.data.frame() %>% 
  select(SampleID, ReadsPerSample, OriginType3) %>% 
  convert(fct(SampleID, OriginType3), dbl(ReadsPerSample)) %>% 
  arrange(desc(ReadsPerSample))

SampleData1DF %>% print.data.frame()
SampleData1DF %>% summary

# Lollipop chart of sequencing depth per sample (Supplementary Figure E3)
SeqDepth_Lolli <- SampleData1DF %>% ggpubr::ggdotchart(x = "SampleID",
                                                       y = "ReadsPerSample",
                                                       color = "OriginType3",
                                                       palette = c(Sputum = "#D55E00",
                                                                   Swab = "#56B4E9", 
                                                                   Control = "#F0E442"),
                                                       sorting = "descending",
                                                       dot.size = 1.5,
                                                       add = "segments",
                                                       xlab = "Samples",
                                                       ylab = "Sequencing depth",
                                                       ggtheme = theme_pubr()) +
  rremove(c("x.text")) +
  rremove(c("x.ticks")) +
  yscale("log10") +
  geom_hline(yintercept = 10^4, linetype = "dashed", color = "black", size = .4)

ggsave("SeqDepth_Lolli.pdf",
       width = 30, height = 10, dpi = 200, units = "cm")

# Inspecting samples with read counts above a threshold that broadly distinguishes
# Patient samples from Controls 
SampleData1AboveDF <- SampleData1DF %>% 
  arrange(desc(ReadsPerSample)) %>% 
  filter(ReadsPerSample > 10^4)
# Count per sample type
SampleData1AboveDF %>% dplyr::count(OriginType3)

# Inspecting samples with a number of reads below threshold
SampleData1BelowDF <- SampleData1DF %>% 
  dplyr::arrange(desc(ReadsPerSample)) %>% 
  filter(ReadsPerSample < 10^4)
# Count per sample type
SampleData1BelowDF %>% dplyr::count(OriginType3)

# Taxa count before filtering based on threshold
physeqBactDecontam2 %>% ntaxa()
# Removing from PS object samples with a number of total reads below threshold
physeqHighDepth <- prune_samples(sample_sums(physeqBactDecontam2) > 10^4, physeqBactDecontam2)
# Keeping only the taxa represented after filtering and counting them
physeqHighDepth %<>% prune_taxa(taxa_sums(.) > 0, .)
physeqHighDepth %>% ntaxa()
# Number of taxa removed during last filtering
physeqBactDecontam2 %>% ntaxa() - physeqHighDepth %>% ntaxa()
# Number of samples remaining
sample_data(physeqHighDepth) %>% pull(OriginType) %>% table
# Reference phyloseq objects - Rarefaction - Transformation (Supplementary Figure E5) ----
# 1) Reference phyloseq objects - Rarefaction - Transformation for Total samples including Controls
# Absolute read counts (i.e. no Hellinger transformation):
# Reference object with Absolute read counts (above 10^4 per sample) and No rarefaction
# for the Total number of samples including Controls
AbsNoRfPs <- physeqHighDepth
# Rarefaction with Absolute read counts (i.e. no Hellinger transformation).
# rngseed() was used to initialize repeatable random subsampling.
# Reference object with Absolute read counts (above 10^4 per sample) and Rarefaction
# for the Total number of samples including Controls
AbsRfPs <- rarefy_even_depth(AbsNoRfPs,
                             rngseed = 567765,
                             sample.size = min(sample_sums(AbsNoRfPs)))
# Inspecting sequencing depth after rarefaction                                   
sample_sums(AbsRfPs)[1:5]
# Number of ASVs removed during rarefaction
AbsNoRfPs %>% ntaxa() - AbsRfPs %>% ntaxa()

# Hellinger transformation
# Reference object with Hellinger-transformed ASV data and No rarefaction
# for the Total number of samples including Controls
HelNoRfPs <- microbiome::transform(AbsNoRfPs, transform = "hellinger", target = "OTU")

# Reference object with Hellinger-transformed ASV data and Rarefaction
# for the Total number of samples including controls
HelRfPs <- microbiome::transform(AbsRfPs, transform = "hellinger", target = "OTU")

# 2) Reference phyloseq objects - Rarefaction - Transformation for Patient samples (Sputum + Swab)
# Absolute read counts (i.e. no Hellinger transformation):
# Reference object with Absolute read counts and No rarefaction for Patient samples
AbsNoRfPatientPs <- AbsNoRfPs %>% subset_samples(. , Origin == "Patient")
summarize_phyloseq(AbsNoRfPatientPs)
# Rarefaction with Absolute read counts (i.e. no Hellinger transformation).
# rngseed() was used to initialize repeatable random subsampling.
# Reference object with Absolute read counts and Rarefaction
# for Patient samples
AbsRfPatientPs <- rarefy_even_depth(AbsNoRfPatientPs, rngseed = 456654)
# ASVs were removed because they were no longer present in any sample
# after random subsampling.
# Compute the selected sequencing depth
sample_sums(AbsRfPatientPs)[1:5]
# Number of ASVs removed during rarefaction
AbsNoRfPatientPs %>% ntaxa() - AbsRfPatientPs %>% ntaxa()

# Plot of rarefaction curves
RareCurvesPatients <- ranacapa::ggrare(AbsNoRfPatientPs,
                                       step = 100, 
                                       color = "Type",
                                       label = "SampleID",
                                       se = FALSE)
# Used to construct Supplementary Figure E5, panels a and b
RareCurvesPatients2 <- RareCurvesPatients + 
  facet_wrap(~factor(Type, levels = c("Sputum", "Swab")), nrow = 2) +
  geom_vline(xintercept = 11741, linetype="dotted") +
  scale_x_continuous(labels = scales::scientific) +
  labs(x = "Sample size", y = "ASV richness") +
  scale_color_manual(values = c(Sputum = "#D55E00", Swab = "#56B4E9")) +
  theme_classic() +
  theme(strip.background = element_blank(),
        legend.position = "none")

ggsave("RareCurvesPatients10055.pdf",
       width = 15, height = 15, dpi = 200, units = "cm")

# Hellinger transformation
# Reference object with Hellinger-transformed ASV data and No rarefaction for Patient samples
HelNoRfPatientPs <- microbiome::transform(AbsNoRfPatientPs, transform = "hellinger", target = "OTU")
# Reference object with Hellinger-transformed ASV data and Rarefaction for Patient samples
HelRfPatientPs <- microbiome::transform(AbsRfPatientPs, transform = "hellinger", target = "OTU")

# Plotting bacterial load per treatment phase (Figure 1) ----
PatientLoadBoxplot <- microbiome::meta(AbsNoRfPatientPs) %>%
  ggplot(aes(x = TTTPrePost5gr, 
             y = Load,
             fill = TTTPrePost5gr)) +
  geom_boxplot(width = .7,
               size = .3,
               alpha = .4,
               position = position_dodge(.8),
               outlier.shape = NA,
               show.legend = TRUE) +
  geom_jitter(position = position_jitterdodge(.8),
              shape = 21,
              colour = "black",
              size = 1,
              alpha = .7) +
  scale_fill_manual(values = TTTphaseCol5) +
  stat_compare_means(comparisons = CompTTT1) + 
  stat_compare_means(label.y = 10.9) +
  facet_wrap(~Type, ncol = 2) +
  scale_y_continuous(breaks = c(10^2, 10^4, 10^6, 10^8), 
                     limits = c(10^2, 10^11), 
                     trans = "log10") +
  labs(fill = "Treatment phase",
       y = "16S rRNA gene copies\n per microlitre of sputum or total OPS sample") +
  theme_classic() +
  theme(strip.background = element_blank()) +
  rotate_x_text(angle = 60, hjust = 1)

ggsave(filename = "PatientLoadBoxplot.pdf",
       width = 15, height = 15, dpi = 200, units = "cm", device='pdf')

# Kruskal-Wallis test for Sputum samples
meta(AbsNoRfPatientPs) %>%
  filter(OriginType3 == "Sputum") %>% 
  droplevels() %>% 
  rstatix::kruskal_test(data = .,
                        Load ~ TTTPrePost5gr) %>%
  adjust_pvalue(method = "holm") %>%
  add_significance("p.adj")

# Kruskal-Wallis test for OPS samples
meta(AbsNoRfPatientPs) %>%
  filter(OriginType3 == "Swab") %>% 
  droplevels() %>% 
  rstatix::kruskal_test(data = .,
                        Load ~ TTTPrePost5gr) %>%
  adjust_pvalue(method = "holm") %>%
  add_significance("p.adj")

# Alpha diversity: Evenness - Richness (Figure 2, panels a, c and d; Supplementary Figure E8, panel a: Supplementary Figure E10) ----
# Alpha diversity based on a selection of metrics:
AlphaDivDF <- microbiome::alpha(AbsRfPatientPs, index = c("chao1", 
                                                          "evenness_camargo",
                                                          "diversity_shannon",
                                                          "dominance_dmn",
                                                          "dominance_core_abundance"))

# Moving the rownames to a new column in diversity data frame
AlphaDivDF$SampleID <- rownames(AlphaDivDF)

# Extracting the metadata from phyloseq object and merging them with diversity data
AbsRfPatientMetaDF <- microbiome::meta(AbsRfPatientPs)
AbsRfPatientMetaDivDF <- left_join(AbsRfPatientMetaDF, AlphaDivDF, by = "SampleID")

# Plotting Chao1 for Sputum (LRT) and OPS (URT) samples (Figure 2a)
Chao1SputumSwabBoxplot <- AbsRfPatientMetaDivDF %>%
  ggplot(aes(x = TTTPrePost5gr, 
             y = chao1,
             fill = TTTPrePost5gr)) +
  geom_boxplot(width = .7,
               size = .3,
               alpha = .4,
               position = position_dodge(.8),
               outlier.shape = NA,
               show.legend = TRUE) +
  geom_jitter(position = position_jitterdodge(.8),
              shape = 21,
              colour = "black",
              size = 1,
              alpha = .7) +
  scale_fill_manual(values = TTTphaseCol5) +
  stat_compare_means(comparisons = CompTTT1) +
  stat_compare_means(label.y = 700) +
  facet_wrap(~Type, ncol = 2) +
  labs(fill = "Treatment phase",
       y = "Chao1 richness") +
  theme_classic() +
  theme(strip.background = element_blank()) +
  rotate_x_text(angle = 60, hjust = 1)

ggsave(filename = "Chao1SputumSwabBoxplot.pdf",
       width = 15, height = 15, dpi = 200, units = "cm", device='pdf')

# Kruskal-Wallis test for Sputum samples
AbsRfPatientMetaDivDF %>%
  filter(OriginType3 == "Sputum") %>% 
  droplevels() %>% 
  kruskal_test(data = .,
               chao1 ~ TTTPrePost5gr) %>%
  adjust_pvalue(method = "BH") %>%
  add_significance("p.adj")
# Dunn's test for pairwise multiple comparisons
AbsRfPatientMetaDivDF %>%
  filter(OriginType3 == "Sputum") %>% 
  droplevels() %>%
  dunn_test(data = .,
            chao1 ~ TTTPrePost5gr,
            p.adjust.method = "BH")

# Kruskal-Wallis test for OPS samples
AbsRfPatientMetaDivDF %>%
  filter(OriginType3 == "Swab") %>% 
  droplevels() %>% 
  kruskal_test(data = .,
               chao1 ~ TTTPrePost5gr) %>%
  adjust_pvalue(method = "BH") %>%
  add_significance("p.adj")
# Dunn's test for pairwise multiple comparisons
AbsRfPatientMetaDivDF %>%
  filter(OriginType3 == "Swab") %>% 
  droplevels() %>%
  dunn_test(data = .,
            chao1 ~ TTTPrePost5gr,
            p.adjust.method = "BH")

# Intra-individual difference in Chao1 richness between StartAZT and EndAZT in sputum samples (LRT)
EndAZTStartAZTChao1DiffSputumDF <- AbsRfPatientMetaDivDF %>%
  select(SampleID, PatientID, OriginType3, TTTPrePost5gr, chao1) %>% 
  filter(OriginType3 == "Sputum",
         TTTPrePost5gr %in% c("StartAZT", "EndAZT")) %>% 
  dplyr::group_by(PatientID) %>%
  arrange(PatientID,
          desc(TTTPrePost5gr)) %>% 
  dplyr::mutate(nSamplesPerPatient = n_distinct(TTTPrePost5gr)) %>% 
  filter(nSamplesPerPatient == 2) %>%
  mutate(EndAZTStartAZTChao1Diff = diff(chao1),
         EndAZTStartAZTChao1Diff = round(EndAZTStartAZTChao1Diff, 0)) %>% 
  select(PatientID, TTTPrePost5gr, chao1, EndAZTStartAZTChao1Diff) %>%
  filter(TTTPrePost5gr == "StartAZT") %>% 
  mutate(EndAZTStartAZTChao1DiffPercentStartAZT = round(EndAZTStartAZTChao1Diff / chao1 * 100, 1))

EndAZTStartAZTChao1DiffSputumDF$EndAZTStartAZTChao1Diff %>% sort()

EndAZTStartAZTChao1DiffSputumDF$EndAZTStartAZTChao1DiffPercentStartAZT %>% sort()

# Intra-individual difference in Chao1 richness between StartAZT and EndAZT in OPS samples (URT)
EndAZTStartAZTChao1DiffOPSDF <- AbsRfPatientMetaDivDF %>%
  select(SampleID, PatientID, OriginType3, TTTPrePost5gr, chao1) %>% 
  filter(OriginType3 == "Swab",
         TTTPrePost5gr %in% c("StartAZT", "EndAZT")) %>% 
  dplyr::group_by(PatientID) %>%
  arrange(PatientID,
          desc(TTTPrePost5gr)) %>% 
  dplyr::mutate(nSamplesPerPatient = n_distinct(TTTPrePost5gr)) %>% 
  filter(nSamplesPerPatient == 2) %>%
  mutate(EndAZTStartAZTChao1Diff = diff(chao1),
         EndAZTStartAZTChao1Diff = round(EndAZTStartAZTChao1Diff, 0)) %>% 
  select(PatientID, TTTPrePost5gr, chao1, EndAZTStartAZTChao1Diff) %>%
  filter(TTTPrePost5gr == "StartAZT") %>% 
  mutate(EndAZTStartAZTChao1DiffPercentStartAZT = round(EndAZTStartAZTChao1Diff / chao1 * 100, 1))

EndAZTStartAZTChao1DiffOPSDF$EndAZTStartAZTChao1Diff %>% sort()

EndAZTStartAZTChao1DiffOPSDF$EndAZTStartAZTChao1DiffPercentStartAZT %>% sort()

# Paired boxplot of Chao1 richness between StartAZT and EndAZT in sputum samples (LRT) (Supplementary Figure E8, panel a)
Chao1StartAZTEndAZTSputum_BoxPaired <- AbsRfPatientMetaDivDF %>% 
  filter(OriginType3 == "Sputum",
         TTTPrePost5gr %in% c("StartAZT", "EndAZT"),
         PatientID %in% as.character(EndAZTStartAZTChao1DiffSputumDF$PatientID)) %>%
ggplot(aes(x = TTTPrePost5gr,
           y = chao1,
           fill = TTTPrePost5gr)) +
  geom_boxplot(width = .5,
               size = .5,
               alpha = .4,
               position = position_dodge(.7),
               outlier.shape = NA,
               show.legend = FALSE) +
  geom_line(aes(group = PatientID),
            colour = "#939696",
            size = .5) +
  geom_point(shape = 21,
             colour = "black",
             size = 1.5,
             alpha = .6) +
  scale_fill_manual(values = StartEndTTTCol2) +
  stat_compare_means(label.y = 475) +
  labs(x = "", y = "Chao1 richness") +
  theme_classic() +
  theme(strip.background = element_blank(),
        axis.text.x = element_blank())

ggsave(filename = "Chao1StartAZTEndAZTSputum_BoxPaired.pdf",
       width = 8, height = 10, dpi = 200, units = "cm", device='pdf')

# Paired boxplot of Chao1 richness between StartAZT and EndAZT in OPS samples (URT) (Supplementary Figure E8, panel a)
Chao1StartAZTEndAZTOPS_BoxPaired <- AbsRfPatientMetaDivDF %>% 
  filter(OriginType3 == "Swab",
         TTTPrePost5gr %in% c("StartAZT", "EndAZT"),
         PatientID %in% as.character(EndAZTStartAZTChao1DiffOPSDF$PatientID)) %>%
  ggplot(aes(x = TTTPrePost5gr,
             y = chao1,
             fill = TTTPrePost5gr)) +
  geom_boxplot(width = .5,
               size = .5,
               alpha = .4,
               position = position_dodge(.7),
               outlier.shape = NA,
               show.legend = FALSE) +
  geom_line(aes(group = PatientID),
            colour = "#939696",
            size = .5) +
  geom_point(shape = 21,
             colour = "black",
             size = 1.5,
             alpha = .6) +
  scale_fill_manual(values = StartEndTTTCol2) +
  stat_compare_means(label.y = 475) +
  labs(x = "", y = "Chao1 richness") +
  theme_classic() +
  theme(strip.background = element_blank(),
        axis.text.x = element_blank())

ggsave(filename = "Chao1StartAZTEndAZTOPS_BoxPaired.pdf",
       width = 8, height = 10, dpi = 200, units = "cm", device='pdf')

# Plotting Camargo evenness for Sputum (LRT) and Swab (URT) (Figure 2c)
CamargoEvenSputumSwabBoxplot <- AbsRfPatientMetaDivDF %>% 
  ggplot(aes(x = TTTPrePost5gr, 
             y = evenness_camargo,
             fill = TTTPrePost5gr)) +
  geom_boxplot(width = .7,
               size = .3,
               alpha = .4,
               position = position_dodge(.8),
               outlier.shape = NA,
               show.legend = TRUE) +
  geom_jitter(position = position_jitterdodge(.8),
              shape = 21,
              colour = "black",
              size = 1,
              alpha = .7) +
  scale_fill_manual(values = TTTphaseCol5) +
  scale_y_continuous(limits=c(.7, 1.3),
                     breaks = seq(.7, 1.3, .1)) +
  stat_compare_means(comparisons = CompTTT1) +
  stat_compare_means(label.y = 1.3) +
  facet_wrap(~Type, ncol = 2) +
  labs(fill = "Treatment phase",
       y = "Camargo evenness") +
  theme_classic() +
  theme(strip.background = element_blank()) +
  rotate_x_text(angle = 60, hjust = 1)

ggsave(filename = "CamargoEvenSputumSwabBoxplot.pdf",
       width = 15, height = 15, dpi = 200, units = "cm", device='pdf')

# Kruskal-Wallis test for Sputum samples
AbsRfPatientMetaDivDF %>%
  filter(OriginType3 == "Sputum") %>% 
  droplevels() %>% 
  kruskal_test(data = .,
               evenness_camargo ~ TTTPrePost5gr) %>%
  adjust_pvalue(method = "BH") %>%
  add_significance("p.adj")

# Kruskal-Wallis test for OPS samples
AbsRfPatientMetaDivDF %>%
  filter(OriginType3 == "Swab") %>% 
  droplevels() %>% 
  kruskal_test(data = .,
               evenness_camargo ~ TTTPrePost5gr) %>%
  adjust_pvalue(method = "BH") %>%
  add_significance("p.adj")
# Dunn's test for pairwise multiple comparisons
AbsRfPatientMetaDivDF %>%
  filter(OriginType3 == "Swab") %>% 
  droplevels() %>%
  dunn_test(data = .,
            evenness_camargo ~ TTTPrePost5gr,
            p.adjust.method = "BH")

# Plotting Shannon diversity for Sputum (LRT) and Swab (URT) (Figure 2d)
ShannonSputumSwabBoxplot <- AbsRfPatientMetaDivDF %>%
  ggplot(aes(x = TTTPrePost5gr, 
             y = diversity_shannon,
             fill = TTTPrePost5gr)) +
  geom_boxplot(width = .7,
               size = .3,
               alpha = .4,
               position = position_dodge(.8),
               outlier.shape = NA,
               show.legend = TRUE) +
  geom_jitter(position = position_jitterdodge(.8),
              shape = 21,
              colour = "black",
              size = 1,
              alpha = .7) +
  scale_fill_manual(values = TTTphaseCol5) +
  stat_compare_means(comparisons = CompTTT1) +
  stat_compare_means(label.y = 6.7) +
  facet_wrap(~Type, ncol = 2) +
  labs(fill = "Treatment phase",
       y = "Shannon diversity") +
  theme_classic() +
  theme(strip.background = element_blank()) +
  rotate_x_text(angle = 60, hjust = 1)

ggsave(filename = "ShannonSputumSwabBoxplot.pdf",
       width = 15, height = 15, dpi = 200, units = "cm", device='pdf')

# Kruskal-Wallis test for Sputum samples
AbsRfPatientMetaDivDF %>%
  filter(OriginType3 == "Sputum") %>% 
  droplevels() %>% 
  kruskal_test(data = .,
               diversity_shannon ~ TTTPrePost5gr) %>%
  adjust_pvalue(method = "BH") %>%
  add_significance("p.adj")
# Dunn's test for pairwise multiple comparisons
AbsRfPatientMetaDivDF %>%
  filter(OriginType3 == "Sputum") %>% 
  droplevels() %>%
  dunn_test(data = .,
            diversity_shannon ~ TTTPrePost5gr,
            p.adjust.method = "BH")

# Kruskal-Wallis test for OPS samples
AbsRfPatientMetaDivDF %>%
  filter(OriginType3 == "Swab") %>% 
  droplevels() %>% 
  kruskal_test(data = .,
               diversity_shannon ~ TTTPrePost5gr) %>%
  adjust_pvalue(method = "BH") %>%
  add_significance("p.adj")

# Plotting Core dominance for Sputum (LRT) and Swab (URT) (Supplementary Figure E10)
CoreDomSputumSwabBoxplot <- AbsRfPatientMetaDivDF %>%
  ggplot(aes(x = TTTPrePost5gr, 
             y = dominance_core_abundance,
             fill = TTTPrePost5gr)) +
  geom_boxplot(width = .7,
               size = .3,
               alpha = .4,
               position = position_dodge(.8),
               outlier.shape = NA,
               show.legend = TRUE) +
  geom_jitter(position = position_jitterdodge(.8),
              shape = 21,
              colour = "black",
              size = 1,
              alpha = .7) +
  scale_fill_manual(values = TTTphaseCol5) +
  stat_compare_means(comparisons = CompTTT1) +
  stat_compare_means(label.y = .95) +
  facet_wrap(~Type, ncol = 2) +
  labs(fill = "Treatment phase",
       y = "Core dominance") +
  theme_classic() +
  theme(strip.background = element_blank()) +
  rotate_x_text(angle = 60, hjust = 1)

ggsave(filename = "CoreDomSputumSwabBoxplot.pdf",
       width = 15, height = 15, dpi = 200, units = "cm", device='pdf')

# Kruskal-Wallis test for Sputum samples
AbsRfPatientMetaDivDF %>%
  filter(OriginType3 == "Sputum") %>% 
  droplevels() %>% 
  kruskal_test(data = .,
               dominance_core_abundance ~ TTTPrePost5gr) %>%
  adjust_pvalue(method = "BH") %>%
  add_significance("p.adj")

# Kruskal-Wallis test for OPS samples
AbsRfPatientMetaDivDF %>%
  filter(OriginType3 == "Swab") %>% 
  droplevels() %>% 
  kruskal_test(data = .,
               dominance_core_abundance ~ TTTPrePost5gr) %>%
  adjust_pvalue(method = "BH") %>%
  add_significance("p.adj")

# Plotting McNaughton dominance (dmn) for Sputum (LRT) and Swab (URT)
DMNSputumSwabBoxplot <- AbsRfPatientMetaDivDF %>%
  ggplot(aes(x = TTTPrePost5gr, 
             y = dominance_dmn,
             fill = TTTPrePost5gr)) +
  geom_boxplot(width = .7,
               size = .3,
               alpha = .4,
               position = position_dodge(.8),
               outlier.shape = NA,
               show.legend = TRUE) +
  geom_jitter(position = position_jitterdodge(.8),
              shape = 21,
              colour = "black",
              size = 1,
              alpha = .7) +
  scale_fill_manual(values = TTTphaseCol5) +
  stat_compare_means(comparisons = CompTTT1) +
  stat_compare_means(label.y = 1.45) +
  facet_wrap(~Type, ncol = 2) +
  labs(fill = "Treatment phase",
       y = "McNaughton dominance") +
  theme_classic() +
  theme(strip.background = element_blank()) +
  rotate_x_text(angle = 60, hjust = 1)

ggsave(filename = "DMNSputumSwabBoxplot.pdf",
       width = 15, height = 15, dpi = 200, units = "cm", device='pdf')

# Kruskal-Wallis test for Sputum samples
AbsRfPatientMetaDivDF %>%
  filter(OriginType3 == "Sputum") %>% 
  droplevels() %>% 
  kruskal_test(data = .,
               dominance_dmn ~ TTTPrePost5gr) %>%
  adjust_pvalue(method = "BH") %>%
  add_significance("p.adj")

# Kruskal-Wallis test for OPS samples
AbsRfPatientMetaDivDF %>%
  filter(OriginType3 == "Swab") %>% 
  droplevels() %>% 
  kruskal_test(data = .,
               dominance_dmn ~ TTTPrePost5gr) %>%
  adjust_pvalue(method = "BH") %>%
  add_significance("p.adj")

# Inspecting Faith's phylogenetic diversity (Figure 2b; Supplementary Figure E8, panel b) ----
# Corresponds to the sum of the branch lengths of phylogenetic tree.
# Obtaining a DF with transposed data from the ASV table
AbsRfPatientASVtabTranspDF <- as.data.frame(otu_table(AbsRfPatientPs)) %>% t()
# Phylogenetic diversity will be measured from a rooted tree
AbsRfPatientTree <- phy_tree(AbsRfPatientPs)
# Using picante::pd to measure phylogenetic diversity
PhyloDivSputumSwabDF <- pd(AbsRfPatientASVtabTranspDF, AbsRfPatientTree, include.root = T)

PhyloDivSputumSwabDF %<>% dplyr::rename(FaithPhylogenDiv = PD, SpeciesRichness = SR)
PhyloDivSputumSwabDF$SampleID <- rownames(PhyloDivSputumSwabDF)

AbsRfPatientMetaDivPhyloDivDF <- left_join(AbsRfPatientMetaDivDF, 
                                           PhyloDivSputumSwabDF, by = "SampleID")

# Plotting Faith's phylogenetic diversity for Sputum (LRT) and Swab (URT) (Figure 2b)
PhyloDivSputumSwabBoxplot <- AbsRfPatientMetaDivPhyloDivDF %>%
  ggplot(aes(x = TTTPrePost5gr, 
             y = FaithPhylogenDiv,
             fill = TTTPrePost5gr)) +
  geom_boxplot(width = .7,
               size = .3,
               alpha = .4,
               position = position_dodge(.8),
               outlier.shape = NA,
               show.legend = TRUE) +
  geom_jitter(position = position_jitterdodge(.8),
              shape = 21,
              colour = "black",
              size = 1,
              alpha = .7) +
  scale_fill_manual(values = TTTphaseCol5) +
  stat_compare_means(comparisons = CompTTT1) +
  stat_compare_means(label.y = 1500) +
  facet_wrap(~Type, ncol = 2) +
  labs(fill = "Treatment phase",
       y = "Faith's phylogenetic diversity") +
  theme_classic() +
  theme(strip.background = element_blank()) +
  rotate_x_text(angle = 60, hjust = 1)

ggsave(filename = "PhyloDivSputumSwabBoxplot.pdf",
       width = 15, height = 15, dpi = 200, units = "cm", device='pdf')

# Kruskal-Wallis test for Sputum samples
AbsRfPatientMetaDivPhyloDivDF %>%
  filter(OriginType3 == "Sputum") %>% 
  droplevels() %>% 
  kruskal_test(data = .,
               FaithPhylogenDiv ~ TTTPrePost5gr) %>%
  adjust_pvalue(method = "BH") %>%
  add_significance("p.adj")
# Dunn's test for pairwise multiple comparisons
AbsRfPatientMetaDivPhyloDivDF %>%
  filter(OriginType3 == "Sputum") %>% 
  droplevels() %>%
  dunn_test(data = .,
            FaithPhylogenDiv ~ TTTPrePost5gr,
            p.adjust.method = "BH")

# Kruskal-Wallis test for OPS samples
AbsRfPatientMetaDivPhyloDivDF %>%
  filter(OriginType3 == "Swab") %>% 
  droplevels() %>% 
  kruskal_test(data = .,
               FaithPhylogenDiv ~ TTTPrePost5gr) %>%
  adjust_pvalue(method = "BH") %>%
  add_significance("p.adj")
# Dunn's test for pairwise multiple comparisons
AbsRfPatientMetaDivPhyloDivDF %>%
  filter(OriginType3 == "Swab") %>% 
  droplevels() %>%
  dunn_test(data = .,
            FaithPhylogenDiv ~ TTTPrePost5gr,
            p.adjust.method = "BH")

# Intra-individual difference in Faith's phylogenetic diversity between StartAZT and EndAZT in sputum samples (LRT)
EndAZTStartAZTFaithDiffSputumDF <- AbsRfPatientMetaDivPhyloDivDF %>%
  select(SampleID, PatientID, OriginType3, TTTPrePost5gr, FaithPhylogenDiv) %>% 
  filter(OriginType3 == "Sputum",
         TTTPrePost5gr %in% c("StartAZT", "EndAZT")) %>% 
  dplyr::group_by(PatientID) %>%
  arrange(PatientID,
          desc(TTTPrePost5gr)) %>% 
  dplyr::mutate(nSamplesPerPatient = n_distinct(TTTPrePost5gr)) %>% 
  filter(nSamplesPerPatient == 2) %>%
  mutate(EndAZTStartAZTFaithDiff = diff(FaithPhylogenDiv),
         EndAZTStartAZTFaithDiff = round(EndAZTStartAZTFaithDiff, 0)) %>% 
  select(PatientID, TTTPrePost5gr, FaithPhylogenDiv, EndAZTStartAZTFaithDiff) %>%
  filter(TTTPrePost5gr == "StartAZT") %>% 
  mutate(EndAZTStartAZTFaithDiffPercentStartAZT = round(EndAZTStartAZTFaithDiff / FaithPhylogenDiv * 100, 1))

EndAZTStartAZTFaithDiffSputumDF$EndAZTStartAZTFaithDiff %>% sort()

EndAZTStartAZTFaithDiffSputumDF$EndAZTStartAZTFaithDiffPercentStartAZT %>% sort()

# Intra-individual difference in Faith's phylogenetic diversity between StartAZT and EndAZT in OPS samples (URT)
EndAZTStartAZTFaithDiffOPSDF <- AbsRfPatientMetaDivPhyloDivDF %>%
  select(SampleID, PatientID, OriginType3, TTTPrePost5gr, FaithPhylogenDiv) %>% 
  filter(OriginType3 == "Swab",
         TTTPrePost5gr %in% c("StartAZT", "EndAZT")) %>% 
  dplyr::group_by(PatientID) %>%
  arrange(PatientID,
          desc(TTTPrePost5gr)) %>% 
  dplyr::mutate(nSamplesPerPatient = n_distinct(TTTPrePost5gr)) %>% 
  filter(nSamplesPerPatient == 2) %>%
  mutate(EndAZTStartAZTFaithDiff = diff(FaithPhylogenDiv),
         EndAZTStartAZTFaithDiff = round(EndAZTStartAZTFaithDiff, 0)) %>% 
  select(PatientID, TTTPrePost5gr, FaithPhylogenDiv, EndAZTStartAZTFaithDiff) %>%
  filter(TTTPrePost5gr == "StartAZT") %>% 
  mutate(EndAZTStartAZTFaithDiffPercentStartAZT = round(EndAZTStartAZTFaithDiff / FaithPhylogenDiv * 100, 1))

EndAZTStartAZTFaithDiffOPSDF$EndAZTStartAZTFaithDiff %>% sort()

EndAZTStartAZTFaithDiffOPSDF$EndAZTStartAZTFaithDiffPercentStartAZT %>% sort()

# Paired boxplot of Faith's phylogenetic diversity between StartAZT and EndAZT in sputum samples (LRT) (Supplementary Figure E8, panel b)
FaithStartAZTEndAZTSputum_BoxPaired <- AbsRfPatientMetaDivPhyloDivDF %>% 
  filter(OriginType3 == "Sputum",
         TTTPrePost5gr %in% c("StartAZT", "EndAZT"),
         PatientID %in% as.character(EndAZTStartAZTFaithDiffSputumDF$PatientID)) %>%
  ggplot(aes(x = TTTPrePost5gr,
             y = FaithPhylogenDiv,
             fill = TTTPrePost5gr)) +
  geom_boxplot(width = .5,
               size = .5,
               alpha = .4,
               position = position_dodge(.7),
               outlier.shape = NA,
               show.legend = FALSE) +
  geom_line(aes(group = PatientID),
            colour = "#939696",
            size = .5) +
  geom_point(shape = 21,
             colour = "black",
             size = 1.5,
             alpha = .6) +
  scale_fill_manual(values = StartEndTTTCol2) +
  stat_compare_means(label.y = 1050) +
  labs(x = "", y = "Faith's phylogenetic diversity") +
  theme_classic() +
  theme(strip.background = element_blank(),
        axis.text.x = element_blank())

ggsave(filename = "FaithStartAZTEndAZTSputum_BoxPaired.pdf",
       width = 8, height = 10, dpi = 200, units = "cm", device='pdf')

# Paired boxplot of Faith's phylogenetic diversity between StartAZT and EndAZT in OPS samples (URT) (Supplementary Figure E8, panel b)
FaithStartAZTEndAZTOPS_BoxPaired <- AbsRfPatientMetaDivPhyloDivDF %>% 
  filter(OriginType3 == "Swab",
         TTTPrePost5gr %in% c("StartAZT", "EndAZT"),
         PatientID %in% as.character(EndAZTStartAZTFaithDiffOPSDF$PatientID)) %>%
  ggplot(aes(x = TTTPrePost5gr,
             y = FaithPhylogenDiv,
             fill = TTTPrePost5gr)) +
  geom_boxplot(width = .5,
               size = .5,
               alpha = .4,
               position = position_dodge(.7),
               outlier.shape = NA,
               show.legend = FALSE) +
  geom_line(aes(group = PatientID),
            colour = "#939696",
            size = .5) +
  geom_point(shape = 21,
             colour = "black",
             size = 1.5,
             alpha = .6) +
  scale_fill_manual(values = StartEndTTTCol2) +
  stat_compare_means(label.y = 1050) +
  labs(x = "", y = "Faith's phylogenetic diversity") +
  theme_classic() +
  theme(strip.background = element_blank(),
        axis.text.x = element_blank())

ggsave(filename = "FaithStartAZTEndAZTOPS_BoxPaired.pdf",
       width = 8, height = 10, dpi = 200, units = "cm", device='pdf')

# PCoA analysis on Hellinger-transformed ASV data, comparing Sputum with Swab (Figure 3a) ----
# Obtaining a PCoA plot based on UniFrac distance
# Computing the distance matrix and ordinating
distUniFracPatient <- phyloseq::distance(HelRfPatientPs, method = "unifrac")
ordPCoAPatient <- ordinate(HelRfPatientPs, method = "PCoA", distance = distUniFracPatient)
# Plotting ellipses and spiders
PCoA_unifracHellingerPatientPlot1 <- ggordiplots::gg_ordiplot(ordPCoAPatient$vectors, 
                                                              groups = sample_data(HelRfPatientPs)$OriginType3,
                                                              spiders = TRUE, 
                                                              ellipse = TRUE,
                                                              pt.size = 1.5, 
                                                              plot = TRUE)

# Modifying the plot using ggplot2 (Figure 3a)
PCoA_UniFrac_Hellinger_SputumSwab_plot <- PCoA_unifracHellingerPatientPlot1$plot +
  geom_segment(data = PCoA_unifracHellingerPatientPlot1$df_spiders, 
               aes(x = cntr.x, 
                   xend = x, 
                   y = cntr.y, 
                   yend = y, 
                   color = Group), 
               size = .1,
               show.legend = FALSE) +
  geom_label(data = PCoA_unifracHellingerPatientPlot1$df_mean.ord, 
             aes(x = x, 
                 y = y, 
                 label = Group),
             colour = "black",
             size = 2.7,
             label.padding = unit(.2, "lines"),
             label.r = unit(0, "lines"),
             label.size = .1,
             show.legend = FALSE) +
  scale_color_manual(values = SampleTypeCol2) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  theme(legend.position = "none") +
  labs(x = "PC1", y = "PC2") +
  theme(axis.title.x = element_text(size = 10),
        axis.title.y = element_text(size = 10),
        axis.text.x = element_text(size = 8),
        axis.text.y = element_text(size = 8))

ggsave(filename = "PCoA_UniFrac_Hellinger_SputumSwab_plot.pdf",
       width = 20, height = 15, dpi = 200, units = "cm", device='pdf')

# Performing PERMANOVA to test dissimilarity between the two groups (i.e. different centroids).
# The null hypothesis is that the sample groups have the same centroid.
# Groups to be compared
sample_data(HelRfPatientPs)$Type
# Seed was set to ensure reproducibility of the analysis
set.seed(5243)
# Calling the pre-calculated UniFrac distance matrix
distUniFracPatient
# Making a data frame from the sample_data
SputumSwabSecDataDF <- sample_data(HelRfPatientPs) %>% data.frame()
# Adonis (column R2 indicates how much of the variance can be explained by the groups)
SputumSwab_unifrac_adonis2 <- vegan::adonis2(distUniFracPatient ~ Type, 
                                             data = SputumSwabSecDataDF)
# Analysis of similarity (ANOSIM)
SputumSwab_unifrac_anosim <- vegan::anosim(distUniFracPatient, 
                                           SputumSwabSecDataDF$Type)

# Distance between Swab samples and the centroid of Sputum samples for each TTT phase (Figure 3b) ----
# Computing the distance matrix 
TTTPatients_UniFrac <- phyloseq::distance(HelRfPatientPs, method = "unifrac")
# Data frame with Secondary data
PatientsSecDataDF <- sample_data(HelRfPatientPs) %>% data.frame()
# Adding a "TTTphase_Type" column
PatientsSecData_TTT_TypeDF <- PatientsSecDataDF %>% 
  mutate(TTTPrePost5grCopy = TTTPrePost5gr, TypeCopy = Type) %>% 
  unite(TTTphase_Type, c(TTTPrePost5grCopy, TypeCopy), sep = "_") %>% 
  convert(fct(TTTphase_Type))
# Computing distance between each sample and each centroid
Patients_centr_distDF <- usedist::dist_to_centroids(TTTPatients_UniFrac, 
                                                    PatientsSecData_TTT_TypeDF$TTTphase_Type)
# Adding TTT group variable in data frame with distances to centroids
Patients_centr_distDF %<>% dplyr::rename(SampleID = Item)
Patients_centr_distDF2 <- left_join(Patients_centr_distDF, 
                                    PatientsSecData_TTT_TypeDF, 
                                    by = "SampleID") %>% 
  select(SampleID, 
         CentroidGroup, 
         CentroidDistance, 
         TTTphase_Type, 
         TTTPrePost5gr, 
         Type)

# Data frame with Swab sample distances to Sputum centroid for each TTT phase
SwabSamples_to_SputumCentroidsDF <- Patients_centr_distDF2 %>%
  filter(Type == "Swab") %>%
  filter(CentroidGroup %in% c("PreAZT_Sputum", "StartAZT_Sputum", "EndAZT_Sputum", 
                              "PostAZT_1mo_Sputum", "PostAZT>=3mo_Sputum")) %>% 
  mutate(CentroidGroup = factor(CentroidGroup, 
                                levels = c("PreAZT_Sputum", "StartAZT_Sputum", "EndAZT_Sputum", 
                                           "PostAZT_1mo_Sputum", "PostAZT>=3mo_Sputum"),
                                labels = c("PreAZT", "StartAZT", "EndAZT", 
                                           "PostAZT_1mo", "PostAZT>=3mo"))) %>% 
  arrange(CentroidGroup)

# Plotting Swab sample distances to Sputum centroid for each TTT phase (Figure 3b)
SwabSamples_to_SputumCentroidsBoxplot <- SwabSamples_to_SputumCentroidsDF %>%
  ggplot(aes(x = CentroidGroup, 
             y = CentroidDistance,
             fill = CentroidGroup)) +
  geom_boxplot(width = .7,
               size = .3,
               alpha = .4,
               position = position_dodge(.8),
               outlier.shape = NA,
               show.legend = TRUE) +
  geom_jitter(position = position_jitterdodge(.8),
              shape = 21,
              colour = "black",
              size = 1,
              alpha = .7) +
  scale_fill_manual(values = TTTphaseCol5) +
  stat_compare_means(comparisons = CompTTT2) +
  stat_compare_means(label.y = .87) +
  labs(fill = "Treatment phase",
       y = "OPS to sputum distance") +
  theme_classic() +
  theme(strip.background = element_blank()) +
  rotate_x_text(angle = 60, hjust = 1)

ggsave(filename = "SwabSamples_to_SputumCentroidsBoxplot.pdf",
       width = 12, height = 12, dpi = 200, units = "cm", device='pdf')

# Kruskal-Wallis test
SwabSamples_to_SputumCentroidsDF %>%
  kruskal_test(data = .,
               CentroidDistance ~ CentroidGroup) %>%
  adjust_pvalue(method = "BH") %>%
  add_significance("p.adj")
# Dunn's test for pairwise multiple comparisons
SwabSamples_to_SputumCentroidsDF %>%
  dunn_test(data = .,
            CentroidDistance ~ CentroidGroup,
            p.adjust.method = "BH")

# wUniFrac - Distance between Swab samples and the centroid of Sputum samples for each TTT phase (Supplementary Figure E12) ----
# Computing the distance matrix 
TTTPatients_wUniFrac <- phyloseq::distance(HelRfPatientPs, method = "wunifrac")
# Data frame with Secondary data
PatientsSecDataDF <- sample_data(HelRfPatientPs) %>% data.frame()
# Adding a "TTTphase_Type" column
PatientsSecData_TTT_TypeDF <- PatientsSecDataDF %>% 
  mutate(TTTPrePost5grCopy = TTTPrePost5gr, TypeCopy = Type) %>% 
  unite(TTTphase_Type, c(TTTPrePost5grCopy, TypeCopy), sep = "_") %>% 
  convert(fct(TTTphase_Type))
# Computing distance between each sample and each centroid
Patients_centr_wUniFracDistDF <- usedist::dist_to_centroids(TTTPatients_wUniFrac, 
                                                            PatientsSecData_TTT_TypeDF$TTTphase_Type)
# Adding TTT group variable in data frame with distances to centroids
Patients_centr_wUniFracDistDF %<>% dplyr::rename(SampleID = Item)
Patients_centr_wUniFracDistDF2 <- left_join(Patients_centr_wUniFracDistDF, 
                                    PatientsSecData_TTT_TypeDF, 
                                    by = "SampleID") %>% 
  select(SampleID, 
         CentroidGroup, 
         CentroidDistance, 
         TTTphase_Type, 
         TTTPrePost5gr, 
         Type)

# Data frame with Swab sample distances to Sputum centroid for each TTT phase
SwabSamples_to_SputumCentroids_wUniFracDF <- Patients_centr_wUniFracDistDF2 %>%
  filter(Type == "Swab") %>%
  filter(CentroidGroup %in% c("PreAZT_Sputum", "StartAZT_Sputum", "EndAZT_Sputum", 
                              "PostAZT_1mo_Sputum", "PostAZT>=3mo_Sputum")) %>% 
  mutate(CentroidGroup = factor(CentroidGroup, 
                                levels = c("PreAZT_Sputum", "StartAZT_Sputum", "EndAZT_Sputum", 
                                           "PostAZT_1mo_Sputum", "PostAZT>=3mo_Sputum"),
                                labels = c("PreAZT", "StartAZT", "EndAZT", 
                                           "PostAZT_1mo", "PostAZT>=3mo"))) %>% 
  arrange(CentroidGroup)

# Plotting Swab sample distances to Sputum centroid for each TTT phase (Supplementary Figure E12)
SwabSamples_to_SputumCentroids_wUniFracBoxplot <- SwabSamples_to_SputumCentroids_wUniFracDF %>%
  ggplot(aes(x = CentroidGroup, 
             y = CentroidDistance,
             fill = CentroidGroup)) +
  geom_boxplot(width = .7,
               size = .3,
               alpha = .4,
               position = position_dodge(.8),
               outlier.shape = NA,
               show.legend = TRUE) +
  geom_jitter(position = position_jitterdodge(.8),
              shape = 21,
              colour = "black",
              size = 1,
              alpha = .7) +
  scale_fill_manual(values = TTTphaseCol5) +
  stat_compare_means(comparisons = CompTTT2) +
  stat_compare_means(label.y = .87) +
  labs(fill = "Treatment phase",
       y = "OPS to sputum distance") +
  theme_classic() +
  theme(strip.background = element_blank()) +
  rotate_x_text(angle = 60, hjust = 1)

ggsave(filename = "SwabSamples_to_SputumCentroids_wUniFracBoxplot.pdf",
       width = 12, height = 12, dpi = 200, units = "cm", device='pdf')

# Kruskal-Wallis test
SwabSamples_to_SputumCentroids_wUniFracDF %>%
  kruskal_test(data = .,
               CentroidDistance ~ CentroidGroup) %>%
  adjust_pvalue(method = "BH") %>%
  add_significance("p.adj")
# Dunn's test for pairwise multiple comparisons
SwabSamples_to_SputumCentroids_wUniFracDF %>%
  dunn_test(data = .,
            CentroidDistance ~ CentroidGroup,
            p.adjust.method = "BH")

# Bray - Distance between Swab samples and the centroid of Sputum samples for each TTT phase ----
# Computing the distance matrix 
TTTPatientsBray <- phyloseq::distance(HelRfPatientPs, method = "bray")
# Data frame with Secondary data
PatientsSecDataDF <- sample_data(HelRfPatientPs) %>% data.frame()
# Adding a "TTTphase_Type" column
PatientsSecData_TTT_TypeDF <- PatientsSecDataDF %>% 
  mutate(TTTPrePost5grCopy = TTTPrePost5gr, TypeCopy = Type) %>% 
  unite(TTTphase_Type, c(TTTPrePost5grCopy, TypeCopy), sep = "_") %>% 
  convert(fct(TTTphase_Type))
# Computing distance between each sample and each centroid
Patients_centrBrayDistDF <- usedist::dist_to_centroids(TTTPatientsBray, 
                                                            PatientsSecData_TTT_TypeDF$TTTphase_Type)
# Adding TTT group variable in data frame with distances to centroids
Patients_centrBrayDistDF %<>% dplyr::rename(SampleID = Item)
Patients_centrBrayDistDF2 <- left_join(Patients_centrBrayDistDF, 
                                            PatientsSecData_TTT_TypeDF, 
                                            by = "SampleID") %>% 
  select(SampleID, 
         CentroidGroup, 
         CentroidDistance, 
         TTTphase_Type, 
         TTTPrePost5gr, 
         Type)

# Data frame with Swab sample distances to Sputum centroid for each TTT phase
SwabSamples_to_SputumCentroidsBrayDF <- Patients_centrBrayDistDF2 %>%
  filter(Type == "Swab") %>%
  filter(CentroidGroup %in% c("PreAZT_Sputum", "StartAZT_Sputum", "EndAZT_Sputum", 
                              "PostAZT_1mo_Sputum", "PostAZT>=3mo_Sputum")) %>% 
  mutate(CentroidGroup = factor(CentroidGroup, 
                                levels = c("PreAZT_Sputum", "StartAZT_Sputum", "EndAZT_Sputum", 
                                           "PostAZT_1mo_Sputum", "PostAZT>=3mo_Sputum"),
                                labels = c("PreAZT", "StartAZT", "EndAZT", 
                                           "PostAZT_1mo", "PostAZT>=3mo"))) %>% 
  arrange(CentroidGroup)

# Plotting Swab sample distances to Sputum centroid for each TTT phase
SwabSamples_to_SputumCentroidsBrayBoxplot <- SwabSamples_to_SputumCentroidsBrayDF %>%
  ggplot(aes(x = CentroidGroup, 
             y = CentroidDistance,
             fill = CentroidGroup)) +
  geom_boxplot(width = .7,
               size = .3,
               alpha = .4,
               position = position_dodge(.8),
               outlier.shape = NA,
               show.legend = TRUE) +
  geom_jitter(position = position_jitterdodge(.8),
              shape = 21,
              colour = "black",
              size = 1,
              alpha = .7) +
  scale_fill_manual(values = TTTphaseCol5) +
  stat_compare_means(comparisons = CompTTT2) +
  stat_compare_means(label.y = .87) +
  labs(fill = "Treatment phase",
       y = "OPS to sputum distance") +
  theme_classic() +
  theme(strip.background = element_blank()) +
  rotate_x_text(angle = 60, hjust = 1)

ggsave(filename = "SwabSamples_to_SputumCentroidsBrayBoxplot.pdf",
       width = 12, height = 12, dpi = 200, units = "cm", device='pdf')

# Kruskal-Wallis test
SwabSamples_to_SputumCentroidsBrayDF %>%
  kruskal_test(data = .,
               CentroidDistance ~ CentroidGroup) %>%
  adjust_pvalue(method = "BH") %>%
  add_significance("p.adj")
# Dunn's test for pairwise multiple comparisons
SwabSamples_to_SputumCentroidsBrayDF %>%
  dunn_test(data = .,
            CentroidDistance ~ CentroidGroup,
            p.adjust.method = "BH")

# PCoA analysis on Hellinger-transformed ASV data: Sputum as per TTT phase (Supplementary Figure E11, panel a) ----
# Obtaining a PCoA plot based on unweigthed UniFrac distance
HelRfSputumPs <- subset_samples(HelRfPatientPs, OriginType3 == "Sputum")
# Computing the distance matrix 
distUniFracSputum <- phyloseq::distance(HelRfSputumPs, method = "unifrac")
# Ordination
ordSputum <- phyloseq::ordinate(HelRfSputumPs, method = "PCoA", distance = distUniFracSputum)

# Plot with ellipses and spiders
PCoA_unifracSputumTTT5grPlot <- gg_ordiplot(ordSputum$vectors,
                                            groups = sample_data(HelRfSputumPs)$TTTPrePost5gr,
                                            ellipse = TRUE,
                                            pt.size = 1.5,
                                            kind = "sd", # Ellipses show one standard deviation around centroids
                                            plot = TRUE)

# Modifying the plot using ggplot2 (Supplementary Figure E11, panel a)
PCoA_unifrac_SputumTTT5gr_plot <- PCoA_unifracSputumTTT5grPlot$plot +
  geom_point(data = PCoA_unifracSputumTTT5grPlot$df_ord, 
             aes(x = x, 
                 y = y, 
                 fill = Group),
             shape = 21, 
             colour = "Black", 
             size = 3) +
  geom_label(data = PCoA_unifracSputumTTT5grPlot$df_mean.ord, 
             aes(x = x, 
                 y = y, 
                 label = Group),
             size = 2.7,
             label.padding = unit(.2, "lines"),
             label.r = unit(0, "lines"),
             label.size = .1,
             show.legend = FALSE) +
  scale_fill_manual(values = TTTphaseCol5) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  labs(x = "PC1", y = "PC2") +
  theme(axis.title.x = element_text(size = 11),
        axis.title.y = element_text(size = 11),
        axis.text.x = element_text(size = 9),
        axis.text.y = element_text(size = 9))

ggsave(filename = "PCoA_unifrac_SputumTTT5gr_plot.pdf",
       width = 20, height = 15, dpi = 200, units = "cm", device='pdf')

# Betadisper and adonis/ANOSIM in Sputum as per TTT phase ----
# Testing for homogeneity of group dispersions and compositional dissimilarity.
HelRfSputumPs <- subset_samples(HelRfPatientPs, OriginType3 == "Sputum")
# Computing the distance matrix 
distUniFracSputum <- phyloseq::distance(HelRfSputumPs, method = "unifrac")
# Computing the average distance of group members to the group centroid 
UniFrac_TTTSputum_betadisper <- betadisper(distUniFracSputum, 
                                           sample_data(HelRfSputumPs)$TTTPrePost5gr, 
                                           type = c("median","centroid"))
# Permutation test for homogeneity of multivariate dispersions
permutest(UniFrac_TTTSputum_betadisper)
# Non-significant permutest results mean the null hypothesis that the groups 
# have the same dispersions cannot be rejected.
# This means we can be more confident that the subsequent adonis result is real,
# and not due to differences in group dispersions.

# PERMANOVA (i.e multivariate analysis of variance starting from distance matrix)
# to test whether two or more groups have similar compositions (i.e. centroids).
# The null hypothesis is that the sample groups have the same centroid.
# Groups to be compared
sample_data(HelRfSputumPs)$TTTPrePost5gr
# Seed was set to ensure reproducibility of the analysis
set.seed(124421)
# Calculating UniFrac distance matrix
TTTSputum_unifrac <- phyloseq::distance(HelRfSputumPs, method = "unifrac")
# Making a data frame from the sample_data
SputumSecDataDF <- sample_data(HelRfSputumPs) %>% data.frame()
# Adonis (column R2 indicates how much of the variance can be explained by the groups)
TTTSputum_unifrac_TTTPrePost5gr_adonis2 <- vegan::adonis2(TTTSputum_unifrac ~ TTTPrePost5gr, 
                                                          data = SputumSecDataDF)
# Analysis of Similarity (ANOSIM)
TTTSputum_unifrac_TTTPrePost5gr_anosim <- vegan::anosim(TTTSputum_unifrac, 
                                                        SputumSecDataDF$TTTPrePost5gr)

# wUniFrac - PCoA analysis on Hellinger-transformed ASV data: Sputum as per TTT phase ----
# Obtaining a PCoA plot based on wUniFrac distance
# Computing the distance matrix 
dist_wUniFracSputum <- phyloseq::distance(HelRfSputumPs, method = "wunifrac")
# Ordination
ord_wUniFracSputum <- phyloseq::ordinate(HelRfSputumPs, method = "PCoA", distance = dist_wUniFracSputum)

# Plot with ellipses and spiders
PCoA_wUniFracSputumTTT5grPlot <- gg_ordiplot(ord_wUniFracSputum$vectors,
                                            groups = sample_data(HelRfSputumPs)$TTTPrePost5gr,
                                            ellipse = TRUE,
                                            pt.size = 1.5,
                                            kind = "sd", # Ellipses show one standard deviation around centroids
                                            plot = TRUE)

# Modifying the plot using ggplot2
PCoA_wUniFrac_SputumTTT5gr_plot <- PCoA_wUniFracSputumTTT5grPlot$plot +
  geom_point(data = PCoA_wUniFracSputumTTT5grPlot$df_ord, 
             aes(x = x, 
                 y = y, 
                 fill = Group),
             shape = 21, 
             colour = "Black", 
             size = 3) +
  geom_label(data = PCoA_wUniFracSputumTTT5grPlot$df_mean.ord, 
             aes(x = x, 
                 y = y, 
                 label = Group),
             size = 2.7,
             label.padding = unit(.2, "lines"),
             label.r = unit(0, "lines"),
             label.size = .1,
             show.legend = FALSE) +
  scale_fill_manual(values = TTTphaseCol5) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  labs(x = "PC1", y = "PC2") +
  theme(axis.title.x = element_text(size = 11),
        axis.title.y = element_text(size = 11),
        axis.text.x = element_text(size = 9),
        axis.text.y = element_text(size = 9))

ggsave(filename = "PCoA_wUniFrac_SputumTTT5gr_plot.pdf",
       width = 20, height = 15, dpi = 200, units = "cm", device='pdf')

# wUniFrac - Betadisper and adonis/ANOSIM in Sputum as per TTT phase ----
# Computing the average distance of group members to the group centroid 
wUniFrac_TTTSputum_betadisper <- betadisper(dist_wUniFracSputum, 
                                           sample_data(HelRfSputumPs)$TTTPrePost5gr, 
                                           type = c("median","centroid"))
# Permutation test for homogeneity of multivariate dispersions
permutest(wUniFrac_TTTSputum_betadisper)
# Non-significant permutest results mean the null hypothesis that the groups 
# have the same dispersions cannot be rejected.
# This means we can be more confident that the subsequent adonis result is real,
# and not due to differences in group dispersions.

# PERMANOVA (i.e multivariate analysis of variance starting from distance matrix)
# to test whether two or more groups have similar compositions (i.e. centroids).
# The null hypothesis is that the sample groups have the same centroid.
# Groups to be compared
sample_data(HelRfSputumPs)$TTTPrePost5gr
# Seed was set to ensure reproducibility of the analysis
set.seed(367976)
# Calculating UniFrac distance matrix
wUniFrac_TTTSputum <- phyloseq::distance(HelRfSputumPs, method = "wunifrac")
# Making a data frame from the sample_data
SputumSecDataDF <- meta(HelRfSputumPs)
# Adonis (column R2 indicates how much of the variance can be explained by the groups)
TTTSputum_wUniFrac_TTTPrePost5gr_adonis2 <- vegan::adonis2(wUniFrac_TTTSputum ~ TTTPrePost5gr, 
                                                          data = SputumSecDataDF)
# Analysis of Similarity (ANOSIM)
TTTSputum_unifrac_TTTPrePost5gr_anosim <- vegan::anosim(wUniFrac_TTTSputum, 
                                                        SputumSecDataDF$TTTPrePost5gr)

# PCoA analysis on Hellinger-transformed ASV data > Swab as per TTT phase (Supplementary Figure E11, panel b) ----
# Obtaining a PCoA plot based on UniFrac distance
HelRfSwabPs <- HelRfPatientPs %>% subset_samples(., OriginType3 == "Swab")
# Computing the distance matrix 
distUniFracSwab <- phyloseq::distance(HelRfSwabPs, method = "unifrac")
# Ordination
ordSwab <- phyloseq::ordinate(HelRfSwabPs, method = "PCoA", distance = distUniFracSwab)

# Plot with ellipses and spiders
PCoA_unifracSwabTTT5grPlot <- gg_ordiplot(ordSwab$vectors,
                                          groups = sample_data(HelRfSwabPs)$TTTPrePost5gr,
                                          ellipse = TRUE, 
                                          pt.size = 1.5, 
                                          kind = "sd", # Ellipses show one standard deviation around centroids
                                          plot = TRUE)

# Modifying the plot using ggplot2 (Supplementary Figure E11, panel b)
PCoA_unifrac_SwabTTT5gr_plot <- PCoA_unifracSwabTTT5grPlot$plot +
  geom_point(data = PCoA_unifracSwabTTT5grPlot$df_ord, 
             aes(x = x, 
                 y = y, 
                 fill = Group),
             shape = 21, 
             colour = "Black", 
             size = 3) +
  geom_label(data = PCoA_unifracSwabTTT5grPlot$df_mean.ord, 
             aes(x = x, 
                 y = y, 
                 label = Group),
             size = 2.7,
             label.padding = unit(.2, "lines"),
             label.r = unit(0, "lines"),
             label.size = .1,
             show.legend = FALSE) +
  scale_fill_manual(values = TTTphaseCol5) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  labs(x = "PC1", y = "PC2") +
  theme(axis.title.x = element_text(size = 11),
        axis.title.y = element_text(size = 11),
        axis.text.x = element_text(size = 9),
        axis.text.y = element_text(size = 9))

ggsave(filename = "PCoA_unifrac_SwabTTT5gr_plot.pdf",
       width = 20, height = 15, dpi = 200, units = "cm", device='pdf')

# Betadisper and adonis/ANOSIM in Swab as per TTT phase ----
# Testing for homogeneity of group dispersions and compositional dissimilarity.
HelRfSwabPs <- HelRfPatientPs %>% subset_samples(., OriginType3 == "Swab")
# Computing the distance matrix 
distUniFracSwab <- phyloseq::distance(HelRfSwabPs, method = "unifrac")
# Computing the average distance of group members to the group centroid 
UniFrac_TTTSwab_betadisper <- betadisper(distUniFracSwab, 
                                         sample_data(HelRfSwabPs)$TTTPrePost5gr, 
                                         type = c("median","centroid"))
# Permutation test for homogeneity of multivariate dispersions
permutest(UniFrac_TTTSwab_betadisper)
# Non-significant permutest results mean the null hypothesis that the groups 
# have the same dispersions cannot be rejected.
# This means we can be more confident that the subsequent adonis result is real,
# and not due to differences in group dispersions.

# PERMANOVA (i.e multivariate analysis of variance starting from distance matrix)
# to test whether two or more groups have similar compositions (i.e. centroids).
# The null hypothesis is that the sample groups have the same centroid.
# Groups to be compared
sample_data(HelRfSwabPs)$TTTPrePost5gr
# Seed was set to ensure reproducibility of the analysis
set.seed(135531)
# Calculating UniFrac distance matrix
TTTSwab_unifrac <- phyloseq::distance(HelRfSwabPs, method = "unifrac")
# Making a data frame from the sample_data
SwabSecDataDF <- sample_data(HelRfSwabPs) %>% data.frame()
# Adonis (column R2 indicates how much of the variance can be explained by the groups)
TTTSwab_unifrac_TTTPrePost5gr_adonis2 <- vegan::adonis2(TTTSwab_unifrac ~ TTTPrePost5gr, 
                                                        data = SwabSecDataDF)
# Analysis of similarity (ANOSIM)
TTTSwab_unifrac_TTTPrePost5gr_anosim <- vegan::anosim(TTTSwab_unifrac, 
                                                      SwabSecDataDF$TTTPrePost5gr)

# CCA in sputum focusing on antibiotic resistance genes (Figure 5, panels a and b) ----
# 1) Data frame with ASV information
# First removing taxa with zero count
HelRfSputumPsNoZeroBact <- HelRfSputumPs %>% 
  prune_taxa(taxa_sums(.) > 0, .)
# Obtaining a data frame with samples in rows
ASVSputumInfoDF <- HelRfSputumPsNoZeroBact %>% 
  otu_table() %>% 
  as.matrix() %>% 
  as.data.frame()
# Transposing the data frame
ASVSputumInfoDF2 <- ASVSputumInfoDF %>% 
  t() %>% 
  as.data.frame()
# 2) Data frame with Secondary data
SecDataSputumDF <- microbiome::meta(HelRfSputumPs)
# Renaming variables to improve readability
SecDataResistGenesSputumDF <- SecDataSputumDF %>% 
  dplyr::rename(tetW = tetW_normal,
                mel = mel_normal,
                tetM = tetM_normal,
                ermB = ermB_normal,
                ermF = ermF_normal,
                mef = mef_normal,
                msrE = msrE_normal) %>% 
  select(SampleID,
         TTTPrePost5gr, 
         mel, 
         mef, 
         ermB, 
         ermF, 
         msrE,
         tetW, 
         tetM)

# Computing Canonical Correspondence Analysis
CCA1 <- SecDataResistGenesSputumDF
ASVSecDataResistGenesSputum_cca <- cca(ASVSputumInfoDF2 ~
                                         CCA1$mel +
                                         CCA1$mef +
                                         CCA1$ermB +
                                         CCA1$ermF +
                                         CCA1$msrE +
                                         CCA1$tetW +
                                         CCA1$tetM)
summary(ASVSecDataResistGenesSputum_cca)

set.seed(866442)
ASVSecDataResistGenesSputum_cca %>% anova.cca(.)

# Plotting the results with "sp" for species scores, "wa" for sample scores and
# "bp" for biplot arrows (For Figure 5a)
plot(ASVSecDataResistGenesSputum_cca, display = c("sp", "wa", "bp"))

# CCA based on carriage of resistance genes: Plotting samples as per TTT phase
CCASputumResistGenesTTTphase_ggordiplot <- gg_ordiplot(ASVSecDataResistGenesSputum_cca,
                                                       groups = SecDataResistGenesSputumDF$TTTPrePost5gr,
                                                       ellipse = TRUE, 
                                                       pt.size = 1.5, 
                                                       plot = TRUE)
# Modifying the plot using ggplot2 (Figure 5b)
CCA_SputumResistGenesTTT5gr_plot <- CCASputumResistGenesTTTphase_ggordiplot$plot +
  geom_point(data = CCASputumResistGenesTTTphase_ggordiplot$df_ord, 
             aes(x = x, 
                 y = y, 
                 fill = Group),
             shape = 21, 
             colour = "Black", 
             size = 3) +
  geom_label(data = CCASputumResistGenesTTTphase_ggordiplot$df_mean.ord, 
             aes(x = x, 
                 y = y, 
                 label = Group),
             size = 2.7,
             label.padding = unit(.2, "lines"),
             label.r = unit(0, "lines"),
             label.size = .1,
             show.legend = FALSE) +
  scale_fill_manual(values = TTTphaseCol5) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  labs(x = "CCA1 (21.0%)", y = "CCA2 (18.3)") +
  theme(legend.position = "top") +
  theme(axis.title.x = element_text(size = 11),
        axis.title.y = element_text(size = 11),
        axis.text.x = element_text(size = 9),
        axis.text.y = element_text(size = 9))

ggsave(filename = "CCA_SputumResistGenesTTT5gr_plot.pdf",
       width = 16, height = 15, dpi = 200, units = "cm", device = "pdf")

# Performing a CCA with inclusion of treatment group
ASVResistGenesTTT5grSputum_cca2 <- cca(ASVSputumInfoDF2 ~
                                         CCA1$TTTPrePost5gr +
                                         CCA1$mel +
                                         CCA1$mef +
                                         CCA1$ermB +
                                         CCA1$ermF +
                                         CCA1$msrE +
                                         CCA1$tetW +
                                         CCA1$tetM)

plot(ASVResistGenesTTT5grSputum_cca2, display = c("sp", "wa", "bp"))
summary(ASVResistGenesTTT5grSputum_cca2)
set.seed(977553)
ASVResistGenesTTT5grSputum_cca2 %>% anova.cca(.)

# Total ARG carriage in Sputum: Longitudinal analyses (Figure 5C; Supplementary Figure E14, panels a and b) ----
# Completing the ARG carriage data frame with patient and treatment phase information
SecDataSputumResistGenesPatientDF <- SecDataResistGenesSputumDF %>%
  mutate(SampleID = rownames(.)) %>% 
  left_join(., 
            SecDataSputumDF %>% select(SampleID, 
                                       PatientID, 
                                       TTTPrePost5gr, 
                                       TTTPrePost7gr), 
            by = "SampleID") %>% 
  relocate(SampleID, .before = mel) %>% 
  relocate(PatientID, .after = SampleID) %>% 
  relocate(TTTPrePost7gr, .after = PatientID)

# Summing carriage of the 7 ARG
Sputum7ResistGenesDF <- SecDataSputumResistGenesPatientDF %>%
  dplyr::rowwise() %>%
  mutate(Total7ResistGenes = sum(across(mel:tetM)))

# Summing carriage of the 7 ARG with rescaling (0 to 1 range)
Sputum_7ResistGenesScale01DF <- SecDataSputumResistGenesPatientDF %>% 
  mutate(across(mel:tetM, ~ scales::rescale(.))) 

Sputum_7ResistGenesScale01DF2 <- Sputum_7ResistGenesScale01DF %>% 
  dplyr::rowwise() %>%
  mutate(Total7ResistGenes = sum(across(mel:tetM)))

# Line plot showing the scaled cumulative carriage of 7 ARG in total patients (Figure 5C)
Sputum_7ResistGenesScale01_LinePlot <- Sputum_7ResistGenesScale01DF2 %>% 
  ggplot(aes(x = TTTPrePost7gr, 
             y = Total7ResistGenes, 
             group = PatientID)) +
  geom_line(aes(color = PatientID), 
            size = .4) +
  geom_point(aes(color = PatientID)) +
  theme_classic() +
  labs(y = "Cumulative carriage of antibiotic\nresistance genes (arbitrary units)") +
  theme(panel.grid = element_blank()) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  theme(axis.title.x = element_blank())

ggsave(filename = "Sputum_7ResistGenesScale01_LinePlot.pdf",
       width = 15, height = 9, dpi = 200, units = "cm", device = "pdf")

# Sputum: Paired analyses between StartAZT and EndAZT
# Selecting patients with StartAZT or EndAZT available
SputumPatientsEndStartDF <- Sputum7ResistGenesDF %>% 
  filter(TTTPrePost7gr == "StartAZT" | TTTPrePost7gr == "EndAZT")
# Selecting patients with the two samples available
SputumPatientsEndStartDF2 <- SputumPatientsEndStartDF %>% 
  group_by(PatientID) %>% 
  filter(n_distinct(SampleID) == 2) %>% 
  ungroup(PatientID) %>% 
  arrange(PatientID)
# Writing a function to identify the difference in total ARG carriage 
# between StartAZT and EndAZT. Input will be a vector of Patient IDs
SelectResistPatients <- function(x) {
  Sputum7ResistGenesEndAZTDF <-  SputumPatientsEndStartDF2 %>% 
    filter(PatientID == x, TTTPrePost7gr == "EndAZT") %>% 
    select(Total7ResistGenes)
  
  Sputum7ResistGenesStartAZTDF <- SputumPatientsEndStartDF2 %>% 
    filter(PatientID == x, TTTPrePost7gr == "StartAZT") %>% 
    select(Total7ResistGenes)
  
  Sputum7ResistGenesEndAZTDF / Sputum7ResistGenesStartAZTDF
}
# Vector of names of patients with StartAZT and EndAZT available
SputumPatientsEndStartVect <- SputumPatientsEndStartDF2 %>%
  distinct(PatientID) %>% 
  pull(PatientID)
# Applying the function to all patients with StartAZT and EndAZT available
EndStartDiff <- SputumPatientsEndStartVect %>% 
  map_df(SelectResistPatients) %>% 
  round(2)

# Including the difference in total ARG carriage into the data frame of secondary data
SputumPatientsEndStartDF3 <- bind_cols(SputumPatientsEndStartVect,
                                       EndStartDiff)

SputumPatientsEndStartScale01DF <- SputumPatientsEndStartDF3 %>% 
  dplyr::rename(PatientID = ...1, 
                ResistGenesEndStartFoldCh = Total7ResistGenes) %>% 
  right_join(.,
             Sputum_7ResistGenesScale01DF2, 
             by = "PatientID") %>% 
  relocate(ResistGenesEndStartFoldCh,
           .after = TTTPrePost7gr) %>% 
  dplyr::rename(TTTPrePost5gr = TTTPrePost5gr.x)

# Adding a new variable related to total ARG carriage based on
# median fold change between "StartAZT" and "EndAZT"
MedFoldChange <- SputumPatientsEndStartScale01DF %>%
  filter(!is.na(ResistGenesEndStartFoldCh)) %>% 
  pull(ResistGenesEndStartFoldCh) %>% median()

SputumResistStatusScale01DF <- SputumPatientsEndStartScale01DF %>% 
  mutate(PatientResistStatus = ifelse(ResistGenesEndStartFoldCh > MedFoldChange,
                                      "IncreasedResist",
                                      "StableResist")) %>% 
  relocate(PatientResistStatus, .after = TTTPrePost7gr) %>% 
  convert(fct(PatientResistStatus))

# Paired box plot of total resistance gene carriage between StartAZT and EndAZT
# with facet by PatientResistStatus (Supplementary Figure E14, panel b)
SputumResistStatusScale01DF2 <- SputumResistStatusScale01DF %>% 
  filter(PatientResistStatus == "IncreasedResist" | PatientResistStatus == "StableResist") %>% 
  filter(TTTPrePost5gr == "StartAZT" | TTTPrePost5gr == "EndAZT")

BoxPaired_StartAZT_EndAZT_TotalResist <- ggpaired(SputumResistStatusScale01DF2, 
                                                  x = "TTTPrePost5gr", 
                                                  y = "Total7ResistGenes",
                                                  id = "PatientID",
                                                  color = "Black",
                                                  alpha = .6,
                                                  fill = "TTTPrePost5gr",
                                                  line.color = "black",
                                                  line.size = .4,
                                                  point.size = 1.2,
                                                  palette = StartEndTTTCol2,
                                                  facet.by = "PatientResistStatus",
                                                  ylab = "Total carriage of antibiotic\nresistance genes (arbitrary units)") +
  stat_compare_means(paired = TRUE) +
  theme(legend.position = "right") +
  theme(panel.grid = element_blank()) +
  theme(axis.text.x = element_blank()) +
  theme(axis.title.x = element_blank())

ggsave("BoxPaired_StartAZT_EndAZT_TotalResist.pdf",
       width = 12, height = 10, dpi = 200, units = "cm", limitsize = FALSE)

# Bar graph showing the difference in total ARG carriage between StartAZT and EndAZT
# in sputum samples, with the median fold increase indicated (Supplementary Figure E14, panel a)
SputumResistDiffMedianBarPlot <- SputumResistStatusScale01DF2 %>%
  filter(PatientID %in% SputumPatientsEndStartVect & TTTPrePost5gr == "EndAZT") %>%
  mutate(PatientID = fct_reorder(PatientID,
                                 desc(ResistGenesEndStartFoldCh))) %>% 
  ggplot(aes(x = PatientID,
             y = ResistGenesEndStartFoldCh)) +
  geom_bar(stat="identity", fill =  "#59aad9") +
  geom_hline(yintercept = MedFoldChange,
             linetype = "dashed",
             color = "red",
             size = .5) +
  theme_classic() +
  labs(y = "Fold-change in carriage of total resistance\ngenes between StartATZ and EndAZT") +
  theme(panel.grid = element_blank()) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  theme(axis.title.x = element_blank())

ggsave("SputumResistDiffMedianBarPlot.pdf",
       width = 7, height = 10, dpi = 200, units = "cm", limitsize = FALSE)

# Comparison of OPS and sputum samples (Figure 3c) ---- 
# 1A) Inspecting taxa present in OPS and sputum samples at 4 or 1 mo before treatment (patients starting with placebo)
# Patients with OPS and sputum samples available at 4 mo before treatment
Placebo4moPreAZTSputumSwabPat <- meta(AbsRfPatientPs) %>%
  filter(TTTPrePost7gr == "PreAZT_4mo") %>% 
  select(PatientID, 
         TTTPrePost7gr, 
         OriginType3) %>% 
  group_by(PatientID) %>%
  filter(n() == 2) %>% 
  distinct(PatientID) %>% 
  pull(PatientID) %>% 
  droplevels()

AbsRfPreAZT4moSputumSwabPs <- subset_samples(AbsRfPatientPs,
                                             PatientID %in% Placebo4moPreAZTSputumSwabPat &
                                               TTTPrePost7gr == "PreAZT_4mo")

# Keeping only the taxa represented after filtering
AbsRfPreAZT4moSputumSwabPs %<>% prune_taxa(taxa_sums(.) > 0, .)

# Patients with OPS and sputum samples available at 1 mo before treatment (end placebo)
Placebo1moPreAZTSputumSwabPat <- meta(AbsRfPatientPs) %>%
  filter(TTTPrePost7gr == "PreAZT_1mo") %>% 
  select(PatientID, TTTPrePost7gr, OriginType3) %>% 
  group_by(PatientID) %>%
  filter(n() == 2) %>% 
  distinct(PatientID) %>% 
  pull(PatientID) %>% 
  droplevels()

AbsRfPreAZT1moSputumSwabPs <- subset_samples(AbsRfPatientPs,
                                             PatientID %in% Placebo1moPreAZTSputumSwabPat &
                                               TTTPrePost7gr == "PreAZT_1mo")

# Keeping only the taxa represented after filtering
AbsRfPreAZT1moSputumSwabPs %<>% prune_taxa(taxa_sums(.) > 0, .)

# Patients with OPS and sputum samples available at 4 or 1 mo before treatment
AbsRfPreAZT4mo1moSputumSwabPs <- merge_phyloseq(otu_table(AbsRfPreAZT4moSputumSwabPs),
                                                otu_table(AbsRfPreAZT1moSputumSwabPs),
                                                tax_table(AbsRfPreAZT4moSputumSwabPs),
                                                tax_table(AbsRfPreAZT1moSputumSwabPs),
                                                sample_data(AbsRfPreAZT4moSputumSwabPs),
                                                sample_data(AbsRfPreAZT1moSputumSwabPs))

# 1B) Inspecting taxa present in OPS and sputum samples at 4 mo or 1 mo before treatment
SputumSwabPreAZT4mo1moNumber <- AbsRfPreAZT4mo1moSputumSwabPs %>% ntaxa()
SputumSwabPreAZT4mo1moNames <- AbsRfPreAZT4mo1moSputumSwabPs %>% taxa_names()
# 1C) Inspecting taxa present in OPS samples at 4 mo or 1 mo before treatment
AbsRfPreAZT4mo1moSwabPs <- subset_samples(AbsRfPreAZT4mo1moSputumSwabPs,
                                          OriginType3 == "Swab")
# Keeping only the taxa represented after filtering, then counting them and listing their names
AbsRfPreAZT4mo1moSwabPs %<>% prune_taxa(taxa_sums(.) > 0, .)
SwabPreAZT4mo1moNumber <- AbsRfPreAZT4mo1moSwabPs %>% ntaxa()
SwabPreAZT4mo1moNames <- AbsRfPreAZT4mo1moSwabPs %>% taxa_names()
# 1D) Inspecting taxa present in sputum samples at 4 mo before or 1 mo before treatment
AbsRfPreAZT4mo1moSputumPs <- subset_samples(AbsRfPreAZT4mo1moSputumSwabPs,
                                            OriginType3 == "Sputum")
# Keeping only the taxa represented after filtering, then counting them and listing their names
AbsRfPreAZT4mo1moSputumPs %<>% prune_taxa(taxa_sums(.) > 0, .)
SputumPreAZT4mo1moNumber <- AbsRfPreAZT4mo1moSputumPs %>% ntaxa()
SputumPreAZT4mo1moNames <- AbsRfPreAZT4mo1moSputumPs %>% taxa_names()

# Preparing a list of vectors with taxa names for visualisation by Venn diagram
SputumSwabPreAZT4mo1moList <- list(OPS = SwabPreAZT4mo1moNames,
                                   Sputum = SputumPreAZT4mo1moNames)
# Calculating the coordinates required for plotting
SputumSwabPreAZT4mo1moVenn <- Venn(SputumSwabPreAZT4mo1moList)
SputumSwabPreAZT4mo1moVennData <- process_data(SputumSwabPreAZT4mo1moVenn)
# Obtaining the Venn diagram (For Figure 3c)
SputumSwabPreAZT4mo1moVennDiagram <- ggVennDiagram(SputumSwabPreAZT4mo1moList,
                                                   set_color = "white", # Hiding original label
                                                   label =  "both", # Showing "count" and "percent"
                                                   label_size = 10,
                                                   label_alpha = .6,
                                                   edge_lty = "blank") +
  # to customize the label layer
  geom_sf_text(aes(label = name),
               fontface = "plain",
               size = 9,
               data = venn_setlabel(SputumSwabPreAZT4mo1moVennData)) +
  # to customize the edge layer
  geom_sf(size = .1,
          color = "#f2f2f2",
          alpha = .6,
          data = venn_setedge(SputumSwabPreAZT4mo1moVennData)) +
  scale_fill_gradient(low = "#f0f0f0", high = "#999999") +
  theme_void() +
  theme(legend.position = "none")

ggsave("SputumSwabPreAZT4mo1moVennDiagram.pdf",
       width = 12, height = 10, dpi = 200, units = "cm", limitsize = FALSE)

# 2A) Inspecting taxa present in OPS and sputum samples at start of treatment
# Patients with OPS and sputum samples available
StartAZTSputumSwabIncrStableARPat <- meta(AbsRfPatientPs) %>%
  filter(TTTPrePost7gr == "StartAZT") %>% 
  select(PatientID, 
         TTTPrePost5gr, 
         OriginType3) %>% 
  group_by(PatientID) %>%
  arrange(PatientID) %>% 
  filter(n() == 2) %>% 
  distinct(PatientID) %>% 
  pull(PatientID) %>% 
  droplevels()

StartAZTSputumSwabIncrStableARPs <- subset_samples(AbsRfPatientPs,
                                                   PatientID %in% StartAZTSputumSwabIncrStableARPat &
                                                     TTTPrePost5gr == "StartAZT")

# Keeping only the taxa represented after filtering
StartAZTSputumSwabIncrStableARPs %<>% prune_taxa(taxa_sums(.) > 0, .)

# Proportion of unique or shared taxa in OPS and sputum samples at the start of treatment
# 2B) Inspecting taxa present in OPS and sputum samples
StartAZTSputumSwabIncrStableNumber <- StartAZTSputumSwabIncrStableARPs %>% ntaxa()
StartAZTSputumSwabIncrStableNames <- StartAZTSputumSwabIncrStableARPs %>% taxa_names()
# 2C) Inspecting taxa present in OPS at the start of treatment
AbsRfStartAZTSwabIncrStableARPs <- subset_samples(StartAZTSputumSwabIncrStableARPs,
                                                  OriginType3 == "Swab")
# Keeping only the taxa represented after filtering, then counting them and listing their names
AbsRfStartAZTSwabIncrStableARPs %<>% prune_taxa(taxa_sums(.) > 0, .)
StartAZTSwabIncrStableARNumber <- AbsRfStartAZTSwabIncrStableARPs %>% ntaxa()
StartAZTSwabIncrStableARNames <- AbsRfStartAZTSwabIncrStableARPs %>% taxa_names()
# 2D) Inspecting taxa present in sputum samples at the start of treatment
AbsRfStartAZTSputumIncrStableARPs <- subset_samples(StartAZTSputumSwabIncrStableARPs,
                                                    OriginType3 == "Sputum")
# Keeping only the taxa represented after filtering, then counting them and listing their names
AbsRfStartAZTSputumIncrStableARPs %<>% prune_taxa(taxa_sums(.) > 0, .)
StartAZTSputumIncrStableARPNumber <- AbsRfStartAZTSputumIncrStableARPs %>% ntaxa()
StartAZTSputumIncrStableARPNames <- AbsRfStartAZTSputumIncrStableARPs %>% taxa_names()

# Preparing a list of vectors with taxa names for visualisation by Venn diagram
SputumSwabStartAZTIncrStableARList <- list(OPS = StartAZTSwabIncrStableARNames,
                                           Sputum = StartAZTSputumIncrStableARPNames)
# Calculating the coordinates required for plotting
SputumSwabStartAZTIncrStableARVenn <- Venn(SputumSwabStartAZTIncrStableARList)
SputumSwabStartAZTIncrStableARVennData <- process_data(SputumSwabStartAZTIncrStableARVenn)
# Obtaining the Venn diagram (For Figure 3c)
SputumSwabStartAZTIncrStableARVennDiagram <- ggVennDiagram(SputumSwabStartAZTIncrStableARList,
                                                           set_color = "white", # Hiding original label
                                                           label =  "both", # Showing "count" and "percent"
                                                           label_size = 10,
                                                           label_alpha = .6,
                                                           edge_lty = "blank") +
  # to customize the label layer
  geom_sf_text(aes(label = name),
               fontface = "plain",
               size = 9,
               data = venn_setlabel(SputumSwabStartAZTIncrStableARVennData)) +
  # to customize the edge layer
  geom_sf(size = .1,
          color = "#EFCA81",
          alpha = .6,
          data = venn_setedge(SputumSwabStartAZTIncrStableARVennData)) +
  scale_fill_gradient(low = "#F7DFAD", high = "#E69F00") +
  theme_void() +
  theme(legend.position = "none")

ggsave("SputumSwabStartAZTIncrStableARVennDiagram.pdf",
       width = 12, height = 10, dpi = 200, units = "cm", limitsize = FALSE)

# 3A) Inspecting taxa present in OPS and sputum samples at the end of treatment
# Patients with OPS and sputum samples available
EndAZTSputumSwabIncrStableARPat <- meta(AbsRfPatientPs) %>%
  filter(TTTPrePost5gr == "EndAZT") %>% 
  select(PatientID, TTTPrePost5gr, OriginType3) %>% 
  group_by(PatientID) %>%
  arrange(PatientID) %>%
  filter(n() == 2) %>% 
  distinct(PatientID) %>% 
  pull(PatientID) %>% 
  droplevels()

EndAZTSputumSwabIncrStableARPs <- subset_samples(AbsRfPatientPs,
                                                 PatientID %in% EndAZTSputumSwabIncrStableARPat &
                                                   TTTPrePost5gr == "EndAZT")

# Keeping only the taxa represented after filtering
EndAZTSputumSwabIncrStableARPs %<>% prune_taxa(taxa_sums(.) > 0, .)

# Proportion of unique or shared taxa in OPS and sputum samples at the end of treatment
# 3B) Inspecting taxa present in OPS and sputum samples
EndAZTSputumSwabIncrStableNumber <- EndAZTSputumSwabIncrStableARPs %>% ntaxa()
EndAZTSputumSwabIncrStableNames <- EndAZTSputumSwabIncrStableARPs %>% taxa_names()
# 3C) Inspecting taxa present in OPS at the end of treatment
AbsRfEndAZTSwabIncrStableARPs <- subset_samples(EndAZTSputumSwabIncrStableARPs,
                                                OriginType3 == "Swab")
# Keeping only the taxa represented after filtering, then counting them and listing their names
AbsRfEndAZTSwabIncrStableARPs %<>% prune_taxa(taxa_sums(.) > 0, .)
EndAZTSwabIncrStableARNumber <- AbsRfEndAZTSwabIncrStableARPs %>% ntaxa()
EndAZTSwabIncrStableARNames <- AbsRfEndAZTSwabIncrStableARPs %>% taxa_names()
# 3D) Inspecting taxa present in sputum samples at the end of treatment
AbsRfEndAZTSputumIncrStableARPs <- subset_samples(EndAZTSputumSwabIncrStableARPs,
                                                  OriginType3 == "Sputum")
# Keeping only the taxa represented after filtering, then counting them and listing their names
AbsRfEndAZTSputumIncrStableARPs %<>% prune_taxa(taxa_sums(.) > 0, .)
EndAZTSputumIncrStableARPNumber <- AbsRfEndAZTSputumIncrStableARPs %>% ntaxa()
EndAZTSputumIncrStableARPNames <- AbsRfEndAZTSputumIncrStableARPs %>% taxa_names()

# Preparing a list of vectors with taxa names for visualisation by Venn diagram
SputumSwabEndAZTIncrStableARList <- list(OPS = EndAZTSwabIncrStableARNames,
                                         Sputum = EndAZTSputumIncrStableARPNames)
# Calculating the coordinates required for plotting
SputumSwabEndAZTIncrStableARVenn <- Venn(SputumSwabEndAZTIncrStableARList)
SputumSwabEndAZTIncrStableARVennData <- process_data(SputumSwabEndAZTIncrStableARVenn)
# Obtaining the Venn diagram (For Figure 3c)
SputumSwabEndAZTIncrStableARVennDiagram <- ggVennDiagram(SputumSwabEndAZTIncrStableARList,
                                                         set_color = "white", # For hiding original label
                                                         label =  "both", # For showing "count" and "percent"
                                                         label_size = 10, # For region label size
                                                         label_alpha = .6, # For adding transparency to region labels
                                                         edge_lty = "blank") +
  # to customize the label layer
  geom_sf_text(aes(label = name),
               fontface = "plain",
               size = 9,
               data = venn_setlabel(SputumSwabEndAZTIncrStableARVennData)) +
  # to customize the edge layer
  geom_sf(size = .1,
          color = "#a5d6f3",
          alpha = .6,
          data = venn_setedge(SputumSwabEndAZTIncrStableARVennData)) +
  scale_fill_gradient(low = "#d2ebf9", high = "#34a4e5") +
  theme_void() +
  theme(legend.position = "none")

ggsave("SputumSwabEndAZTIncrStableARVennDiagram.pdf",
       width = 12, height = 10, dpi = 200, units = "cm", limitsize = FALSE)

# 4A) Inspecting taxa present in OPS and sputum samples at 1 mo post-treatment
# Patients with OPS and sputum samples available
PostAZT1moSputumSwabIncrStableARPat <- meta(AbsRfPatientPs) %>%
  filter(TTTPrePost7gr == "PostAZT_1mo") %>% 
  select(PatientID, TTTPrePost5gr, OriginType3) %>% 
  group_by(PatientID) %>%
  arrange(PatientID) %>% 
  filter(n() == 2) %>% 
  distinct(PatientID) %>% 
  pull(PatientID) %>% 
  droplevels()

PostAZT1moSputumSwabIncrStableARPs <- subset_samples(AbsRfPatientPs,
                                                     PatientID %in% PostAZT1moSputumSwabIncrStableARPat &
                                                       TTTPrePost5gr == "PostAZT_1mo")

# Keeping only the taxa represented after filtering
PostAZT1moSputumSwabIncrStableARPs %<>% prune_taxa(taxa_sums(.) > 0, .)

# Proportion of unique or shared taxa in OPS and sputum samples at 1 mo post-treatment
# 4B) Inspecting taxa present in OPS and sputum samples
PostAZT1moSputumSwabIncrStableNumber <- PostAZT1moSputumSwabIncrStableARPs %>% ntaxa()
PostAZT1moSputumSwabIncrStableNames <- PostAZT1moSputumSwabIncrStableARPs %>% taxa_names()
# 4C) Inspecting taxa present in OPS at 1 mo post-treatment
AbsRfPostAZT1moSwabIncrStableARPs <- subset_samples(PostAZT1moSputumSwabIncrStableARPs,
                                                    OriginType3 == "Swab")
# Keeping only the taxa represented after filtering, then counting them and listing their names
AbsRfPostAZT1moSwabIncrStableARPs %<>% prune_taxa(taxa_sums(.) > 0, .)
PostAZT1moSwabIncrStableARNumber <- AbsRfPostAZT1moSwabIncrStableARPs %>% ntaxa()
PostAZT1moSwabIncrStableARNames <- AbsRfPostAZT1moSwabIncrStableARPs %>% taxa_names()
# 4D) Inspecting taxa present in sputum samples at 1 mo post-treatment
AbsRfPostAZT1moSputumIncrStableARPs <- subset_samples(PostAZT1moSputumSwabIncrStableARPs,
                                                      OriginType3 == "Sputum")
# Keeping only the taxa represented after filtering, then counting them and listing their names
AbsRfPostAZT1moSputumIncrStableARPs %<>% prune_taxa(taxa_sums(.) > 0, .)
PostAZT1moSputumIncrStableARPNumber <- AbsRfPostAZT1moSputumIncrStableARPs %>% ntaxa()
PostAZT1moSputumIncrStableARPNames <- AbsRfPostAZT1moSputumIncrStableARPs %>% taxa_names()

# Preparing a list of vectors with taxa names for visualisation by Venn diagram
SputumSwabPostAZT1moIncrStableARList <- list(OPS = PostAZT1moSwabIncrStableARNames,
                                             Sputum = PostAZT1moSputumIncrStableARPNames)

# Calculating the coordinates required for plotting
SputumSwabPostAZT1moIncrStableARVenn <- Venn(SputumSwabPostAZT1moIncrStableARList)
SputumSwabPostAZT1moIncrStableARVennData <- process_data(SputumSwabPostAZT1moIncrStableARVenn)
# Obtaining the Venn diagram (For Figure 3c)
SputumSwabPostAZT1moIncrStableARVennDiagram <- ggVennDiagram(SputumSwabPostAZT1moIncrStableARList,
                                                             set_color = "white", # Hiding original label
                                                             label =  "both", # Showing "count" and "percent"
                                                             label_size = 10,
                                                             label_alpha = .6,
                                                             edge_lty = "blank") +
  # to customize the label layer
  geom_sf_text(aes(label = name),
               fontface = "plain",
               size = 9,
               data = venn_setlabel(SputumSwabPostAZT1moIncrStableARVennData)) +
  # to customize the edge layer
  geom_sf(size = .1,
          color = "#99ffe4",
          alpha = .6,
          data = venn_setedge(SputumSwabPostAZT1moIncrStableARVennData)) +
  scale_fill_gradient(low = "#ccfff1", high = "#00cc96") +
  theme_void() +
  theme(legend.position = "none")

ggsave("SputumSwabPostAZT1moIncrStableARVennDiagram.pdf",
       width = 12, height = 10, dpi = 200, units = "cm", limitsize = FALSE)

# 5A) Inspecting taxa present in OPS and sputum samples at 4 mo post-treatment
# Patients with OPS and sputum samples available
PostAZT4moSputumSwabIncrStableARPat <- meta(AbsRfPatientPs) %>%
  filter(TTTPrePost7gr == "PostAZT_4mo") %>%
  group_by(PatientID) %>%
  arrange(PatientID) %>% 
  filter(n() == 2) %>% 
  distinct(PatientID) %>% 
  pull(PatientID) %>% 
  droplevels()

PostAZT4moSputumSwabIncrStableARPs <- subset_samples(AbsRfPatientPs,
                                                     PatientID %in% PostAZT4moSputumSwabIncrStableARPat &
                                                       TTTPrePost7gr == "PostAZT_4mo")

# Keeping only the taxa represented after filtering
PostAZT4moSputumSwabIncrStableARPs %<>% prune_taxa(taxa_sums(.) > 0, .)

# Proportion of unique or shared taxa in OPS and sputum samples at 4 mo post-treatment
# 5B) Inspecting taxa present in OPS and sputum samples
PostAZT4moSputumSwabIncrStableNumber <- PostAZT4moSputumSwabIncrStableARPs %>% ntaxa()
PostAZT4moSputumSwabIncrStableNames <- PostAZT4moSputumSwabIncrStableARPs %>% taxa_names()
# 5C) Inspecting taxa present in OPS at 4 mo post-treatment
AbsRfPostAZT4moSwabIncrStableARPs <- subset_samples(PostAZT4moSputumSwabIncrStableARPs,
                                                    OriginType3 == "Swab")
# Keeping only the taxa represented after filtering, then counting them and listing their names
AbsRfPostAZT4moSwabIncrStableARPs %<>% prune_taxa(taxa_sums(.) > 0, .)
PostAZT4moSwabIncrStableARNumber <- AbsRfPostAZT4moSwabIncrStableARPs %>% ntaxa()
PostAZT4moSwabIncrStableARNames <- AbsRfPostAZT4moSwabIncrStableARPs %>% taxa_names()
# 5D) Inspecting taxa present in sputum samples at 4 mo post-treatment
AbsRfPostAZT4moSputumIncrStableARPs <- subset_samples(PostAZT4moSputumSwabIncrStableARPs,
                                                      OriginType3 == "Sputum")
# Keeping only the taxa represented after filtering, then counting them and listing their names
AbsRfPostAZT4moSputumIncrStableARPs %<>% prune_taxa(taxa_sums(.) > 0, .)
PostAZT4moSputumIncrStableARPNumber <- AbsRfPostAZT4moSputumIncrStableARPs %>% ntaxa()
PostAZT4moSputumIncrStableARPNames <- AbsRfPostAZT4moSputumIncrStableARPs %>% taxa_names()

# Preparing a list of vectors with taxa names for visualisation by Venn diagram
SputumSwabPostAZT4moIncrStableARList <- list(OPS = PostAZT4moSwabIncrStableARNames,
                                             Sputum = PostAZT4moSputumIncrStableARPNames)

# Calculating the coordinates required for plotting
SputumSwabPostAZT4moIncrStableARVenn <- Venn(SputumSwabPostAZT4moIncrStableARList)
SputumSwabPostAZT4moIncrStableARVennData <- process_data(SputumSwabPostAZT4moIncrStableARVenn)
# Obtaining the Venn diagram (For Figure 3c)
SputumSwabPostAZT4moIncrStableARVennDiagram <- ggVennDiagram(SputumSwabPostAZT4moIncrStableARList,
                                                             set_color = "white", # Hiding original label
                                                             label =  "both", # Showing "count" and "percent"
                                                             label_size = 10,
                                                             label_alpha = .6,
                                                             edge_lty = "blank") +
  # to customize the label layer
  geom_sf_text(aes(label = name),
               fontface = "plain",
               size = 9,
               data = venn_setlabel(SputumSwabPostAZT4moIncrStableARVennData)) +
  # to customize the edge layer
  geom_sf(size = .1,
          color = "#eac8da",
          alpha = .6,
          data = venn_setedge(SputumSwabPostAZT4moIncrStableARVennData)) +
  scale_fill_gradient(low = "#f1dae7", high = "#cc79a7") +
  theme_void() +
  theme(legend.position = "none")

ggsave("SputumSwabPostAZT3moIncrStableARVennDiagram.pdf",
       width = 12, height = 10, dpi = 200, units = "cm", limitsize = FALSE)

# Statistics for comparisons between OPS and Sputum per patient (For Figure 3c) ----
# 1) Comparisons between OPS and Sputum per patient before AZM treatment
# Patients with a sample pair (OPS swab + Sputum) available at 4 and/or 1 mo before AZM treatment
AbsRfPreAZT4mo1moSputumSwabPs2 <- AbsRfPreAZT4mo1moSputumSwabPs

# Adding a "Comparisons" column in sample_data object for later use in a function
sample_data(AbsRfPreAZT4mo1moSputumSwabPs2)$Comparisons <- meta(AbsRfPreAZT4mo1moSputumSwabPs2) %>%
  mutate(Comparisons = paste(PatientID, TTTPrePost7gr, sep = "_")) %>% 
  pull(Comparisons)

# Adding a "Conditions" column in sample_data object for later use in a function
sample_data(AbsRfPreAZT4mo1moSputumSwabPs2)$Conditions <- meta(AbsRfPreAZT4mo1moSputumSwabPs2) %>%
  mutate(Conditions = ifelse(OriginType3 == "Swab", "Condition1", "Condition2")) %>% 
  pull(Conditions)

# Obtaining a data frame with all the data
SputumSwabCond1Cond2DF <- psmelt(AbsRfPreAZT4mo1moSputumSwabPs2) %>%
  arrange(Comparisons, Conditions, OTU, Abundance) %>% 
  select(Comparisons, Conditions, OTU, Abundance)

# Vector of comparisons to be made (i.e. distinct Comparisons)
ComparisonsPreAZT4mo1moSputumSwabVect <- SputumSwabCond1Cond2DF %>% 
  distinct(Comparisons) %>% 
  pull(Comparisons)

# Function to determine the distribution of ASVs between sputum and OPS (i.e. conditions)
# for each patient (i.e. comparisons). The function requires that the data be contained in a data frame
# called "SputumSwabCond1Cond2DF" and containing a "Comparisons" column and a "Conditions" column
ASVSputumSwabComparison <- function(Comparison){
  # Listing the names and counting ASVs represented in condition 1
  Cond1ASVNames <- SputumSwabCond1Cond2DF %>%
    filter(Comparisons == Comparison & Conditions == "Condition1" & Abundance > 0) %>% 
    pull(OTU)
  
  Cond1ASVNumber <- length(Cond1ASVNames)
  
  # Listing the names and counting ASVs represented in condition 2
  Cond2ASVNames <- SputumSwabCond1Cond2DF %>%
    filter(Comparisons == Comparison & Conditions == "Condition2" & Abundance > 0) %>% 
    pull(OTU)
  
  Cond2ASVNumber <- length(Cond2ASVNames)
  
  # Listing the names and counting ASVs represented in condition 1 ONLY
  Cond1OnlyASVNames <- setdiff(Cond1ASVNames, Cond2ASVNames)
  Cond1OnlyASVNumber <- length(Cond1OnlyASVNames)
  
  # Listing the names and counting ASVs represented in condition 2 ONLY
  Cond2OnlyASVNames <- setdiff(Cond2ASVNames, Cond1ASVNames)
  Cond2OnlyASVNumber <- length(Cond2OnlyASVNames)
  
  # Listing the names and counting ASVs represented in both condition 1 and condition 2
  SharedASVNames <- intersect(Cond1ASVNames, Cond2ASVNames)
  SharedASVNumber <- length(SharedASVNames)
  
  # Counting total ASVs per comparison
  Cond1Cond2TotalASVNumber <- Cond1OnlyASVNumber + Cond2OnlyASVNumber + SharedASVNumber
  
  # Calculating the percentage of ASVs in condition 1 ONLY
  Cond1OnlyPerCentASVNumber <- round(Cond1OnlyASVNumber / Cond1Cond2TotalASVNumber * 100, digits = 1)
  
  # Calculating the percentage of ASVs in condition 2 ONLY
  Cond2OnlyPerCentASVNumber <- round(Cond2OnlyASVNumber / Cond1Cond2TotalASVNumber * 100, digits = 1)
  
  # Calculating the percentage of shared ASVs
  PerCentSharedASVNumber <- round(SharedASVNumber / Cond1Cond2TotalASVNumber * 100, digits = 1)
  
  ComparisonOutput <- c(Cond1OnlyASVNumber,
                        SharedASVNumber,
                        Cond2OnlyASVNumber,
                        Cond1OnlyPerCentASVNumber,
                        PerCentSharedASVNumber,
                        Cond2OnlyPerCentASVNumber)
  
  return(ComparisonOutput)
}

# Obtaining a data frame containing the output of the function call 
ComparisonVect <- ComparisonsPreAZT4mo1moSputumSwabVect
ComparisonList <- as.list(ComparisonVect)

PreAZT4mo1moSputumSwabComparisonList <- map(ComparisonList, ASVSputumSwabComparison)

PreAZT4mo1moSputumSwabComparisonDF <- PreAZT4mo1moSputumSwabComparisonList %>%
  as.data.frame(col.names = ComparisonVect,
                row.names = c("Cond1OnlyASVNumber",
                              "SharedASVNumber",
                              "Cond2OnlyASVNumber",
                              "Cond1OnlyPerCentASVNumber",
                              "PerCentSharedASVNumber",
                              "Cond2OnlyPerCentASVNumber")) %>% 
  # Transposing and coercing to a data frame
  t() %>% 
  as.data.frame()

# Testing normality and plotting the distribution
shapiro.test(PreAZT4mo1moSputumSwabComparisonDF$Cond1OnlyASVNumber)

ggdensity(PreAZT4mo1moSputumSwabComparisonDF$Cond1OnlyASVNumber, 
          main = "Density plot of Cond1OnlyASVNumber",
          xlab = "Number of ASVs")

shapiro.test(PreAZT4mo1moSputumSwabComparisonDF$SharedASVNumber)

ggdensity(PreAZT4mo1moSputumSwabComparisonDF$SharedASVNumber, 
          main = "Density plot of SharedASVNumber",
          xlab = "Number of ASVs")

shapiro.test(PreAZT4mo1moSputumSwabComparisonDF$Cond2OnlyASVNumber)

ggdensity(PreAZT4mo1moSputumSwabComparisonDF$Cond2OnlyASVNumber, 
          main = "Density plot of Cond2OnlyASVNumber",
          xlab = "Number of ASVs")

# Obtaining a long data frame to compare medians
PreAZT4mo1moSputumSwabComparisonLongDF <- PreAZT4mo1moSputumSwabComparisonDF %>%
  select("Cond1OnlyASVNumber", "SharedASVNumber", "Cond2OnlyASVNumber") %>% 
  gather(., key = ASVDistribution, value = ASVNumber) %>% 
  convert(fct(ASVDistribution))

# Computing Kruskal-Wallis test
kruskal.test(ASVNumber ~ ASVDistribution,
             data = PreAZT4mo1moSputumSwabComparisonLongDF)

# Multiple pairwise-comparison between groups, with Benjamini & Hochberg correction for multiple testing
pairwise.wilcox.test(PreAZT4mo1moSputumSwabComparisonLongDF$ASVNumber,
                     PreAZT4mo1moSputumSwabComparisonLongDF$ASVDistribution,
                     p.adjust.method = "BH")

# Computing summary statistics per condition
PreAZT4mo1moSputumSwabComparisonLongPerCentDF <- PreAZT4mo1moSputumSwabComparisonDF %>%
  select("Cond1OnlyPerCentASVNumber", "PerCentSharedASVNumber", "Cond2OnlyPerCentASVNumber") %>% 
  gather(., key = ASVDistribution, value = ASVPerCent) %>% 
  convert(fct(ASVDistribution)) %>%
  mutate(ASVDistribution = factor(ASVDistribution,
                                  levels = c("Cond1OnlyPerCentASVNumber",
                                             "PerCentSharedASVNumber",
                                             "Cond2OnlyPerCentASVNumber")))

PreAZT4mo1moSputumSwabSummaryStatDF <- PreAZT4mo1moSputumSwabComparisonLongPerCentDF %>% 
  group_by(ASVDistribution) %>% 
  dplyr::summarize(medianASVPerCent = median(ASVPerCent),
            IQR = IQR(ASVPerCent)) %>% 
  ungroup()

# 2) Comparisons between OPS and Sputum per patient at start of AZM treatment
# Patients with a sample pair (OPS swab + Sputum) available at Start of AZM treatment
AbsRfStartSputumSwabPs2 <- StartAZTSputumSwabIncrStableARPs

# Adding a "Comparisons" column in sample_data object for later use in a function
sample_data(AbsRfStartSputumSwabPs2)$Comparisons <- meta(AbsRfStartSputumSwabPs2) %>%
  mutate(Comparisons = paste(PatientID, TTTPrePost7gr, sep = "_")) %>% 
  pull(Comparisons)

# Adding a "Conditions" column in sample_data object for later use in a function
sample_data(AbsRfStartSputumSwabPs2)$Conditions <- meta(AbsRfStartSputumSwabPs2) %>%
  mutate(Conditions = ifelse(OriginType3 == "Swab", "Condition1", "Condition2")) %>% 
  pull(Conditions)

# Obtaining a data frame with all the data
SputumSwabCond1Cond2StartDF <- psmelt(AbsRfStartSputumSwabPs2) %>%
  arrange(Comparisons, Conditions, OTU, Abundance) %>% 
  select(Comparisons, Conditions, OTU, Abundance)

# Vector of comparisons to be made (i.e. distinct Comparisons)
ComparisonsStartSputumSwabVect <- SputumSwabCond1Cond2StartDF %>% 
  distinct(Comparisons) %>% 
  pull(Comparisons)

# Function to determine the distribution of ASVs between sputum and OPS (i.e. conditions)
# for each patient (i.e. comparisons). The function requires that the data be contained in a data frame
# called "SputumSwabCond1Cond2StartDF" and containing a "Comparisons" column and a "Conditions" column
ASVSputumSwabComparisonStart <- function(Comparison){
  # Listing the names and counting ASVs represented in condition 1
  Cond1ASVNames <- SputumSwabCond1Cond2StartDF %>%
    filter(Comparisons == Comparison & Conditions == "Condition1" & Abundance > 0) %>% 
    pull(OTU)
  
  Cond1ASVNumber <- length(Cond1ASVNames)
  
  # Listing the names and counting ASVs represented in condition 2
  Cond2ASVNames <- SputumSwabCond1Cond2StartDF %>%
    filter(Comparisons == Comparison & Conditions == "Condition2" & Abundance > 0) %>% 
    pull(OTU)
  
  Cond2ASVNumber <- length(Cond2ASVNames)
  
  # Listing the names and counting ASVs represented in condition 1 ONLY
  Cond1OnlyASVNames <- setdiff(Cond1ASVNames, Cond2ASVNames)
  Cond1OnlyASVNumber <- length(Cond1OnlyASVNames)
  
  # Listing the names and counting ASVs represented in condition 2 ONLY
  Cond2OnlyASVNames <- setdiff(Cond2ASVNames, Cond1ASVNames)
  Cond2OnlyASVNumber <- length(Cond2OnlyASVNames)
  
  # Listing the names and counting ASVs represented in both condition 1 and condition 2
  SharedASVNames <- intersect(Cond1ASVNames, Cond2ASVNames)
  SharedASVNumber <- length(SharedASVNames)
  
  # Counting total ASVs per comparison
  Cond1Cond2TotalASVNumber <- Cond1OnlyASVNumber + Cond2OnlyASVNumber + SharedASVNumber
  
  # Calculating the percentage of ASVs in condition 1 ONLY
  Cond1OnlyPerCentASVNumber <- round(Cond1OnlyASVNumber / Cond1Cond2TotalASVNumber * 100, digits = 1)
  
  # Calculating the percentage of ASVs in condition 2 ONLY
  Cond2OnlyPerCentASVNumber <- round(Cond2OnlyASVNumber / Cond1Cond2TotalASVNumber * 100, digits = 1)
  
  # Calculating the percentage of shared ASVs
  PerCentSharedASVNumber <- round(SharedASVNumber / Cond1Cond2TotalASVNumber * 100, digits = 1)
  
  ComparisonOutput <- c(Cond1OnlyASVNumber,
                        SharedASVNumber,
                        Cond2OnlyASVNumber,
                        Cond1OnlyPerCentASVNumber,
                        PerCentSharedASVNumber,
                        Cond2OnlyPerCentASVNumber)
  
  return(ComparisonOutput)
}

# Obtaining a data frame containing the output of the function call 
ComparisonVect <- ComparisonsStartSputumSwabVect
ComparisonList <- as.list(ComparisonVect)

StartSputumSwabComparisonList <- map(ComparisonList, ASVSputumSwabComparisonStart)

StartSputumSwabComparisonDF <- StartSputumSwabComparisonList %>%
  as.data.frame(col.names = ComparisonVect,
                row.names = c("Cond1OnlyASVNumber",
                              "SharedASVNumber",
                              "Cond2OnlyASVNumber",
                              "Cond1OnlyPerCentASVNumber",
                              "PerCentSharedASVNumber",
                              "Cond2OnlyPerCentASVNumber")) %>% 
  # Transposing and coercing to a data frame
  t() %>% 
  as.data.frame()

# Testing normality and plotting the distribution
shapiro.test(StartSputumSwabComparisonDF$Cond1OnlyASVNumber)

ggdensity(StartSputumSwabComparisonDF$Cond1OnlyASVNumber, 
          main = "Density plot of Cond1OnlyASVNumber",
          xlab = "Number of ASVs")

shapiro.test(StartSputumSwabComparisonDF$SharedASVNumber)

ggdensity(StartSputumSwabComparisonDF$SharedASVNumber, 
          main = "Density plot of SharedASVNumber",
          xlab = "Number of ASVs")

shapiro.test(StartSputumSwabComparisonDF$Cond2OnlyASVNumber)

ggdensity(StartSputumSwabComparisonDF$Cond2OnlyASVNumber, 
          main = "Density plot of Cond2OnlyASVNumber",
          xlab = "Number of ASVs")

# Obtaining a long data frame to compare medians
StartSputumSwabComparisonLongDF <- StartSputumSwabComparisonDF %>%
  select("Cond1OnlyASVNumber", "SharedASVNumber", "Cond2OnlyASVNumber") %>% 
  gather(., key = ASVDistribution, value = ASVNumber) %>% 
  convert(fct(ASVDistribution))

# Computing Kruskal-Wallis test
kruskal.test(ASVNumber ~ ASVDistribution,
             data = StartSputumSwabComparisonLongDF)

# Multiple pairwise-comparison between groups, with Benjamini & Hochberg correction for multiple testing
pairwise.wilcox.test(StartSputumSwabComparisonLongDF$ASVNumber,
                     StartSputumSwabComparisonLongDF$ASVDistribution,
                     p.adjust.method = "BH")

# Computing summary statistics per condition
StartSputumSwabComparisonLongPerCentDF <- StartSputumSwabComparisonDF %>%
  select("Cond1OnlyPerCentASVNumber", "PerCentSharedASVNumber", "Cond2OnlyPerCentASVNumber") %>% 
  gather(., key = ASVDistribution, value = ASVPerCent) %>% 
  convert(fct(ASVDistribution)) %>%
  mutate(ASVDistribution = factor(ASVDistribution,
                                  levels = c("Cond1OnlyPerCentASVNumber",
                                             "PerCentSharedASVNumber",
                                             "Cond2OnlyPerCentASVNumber")))

StartSputumSwabSummaryStatDF <- StartSputumSwabComparisonLongPerCentDF %>% 
  group_by(ASVDistribution) %>% 
  dplyr::summarize(medianASVPerCent = median(ASVPerCent),
            IQR = IQR(ASVPerCent)) %>% 
  ungroup()

# 3) Comparisons between OPS and Sputum per patient at the end of AZM treatment
# Patients with a sample pair (OPS swab + Sputum) available at the end of AZM treatment
AbsRfEndSputumSwabPs2 <- EndAZTSputumSwabIncrStableARPs

# Adding a "Comparisons" column in sample_data object for later use in a function
sample_data(AbsRfEndSputumSwabPs2)$Comparisons <- meta(AbsRfEndSputumSwabPs2) %>%
  mutate(Comparisons = paste(PatientID, TTTPrePost7gr, sep = "_")) %>% 
  pull(Comparisons)

# Adding a "Conditions" column in sample_data object for later use in a function
sample_data(AbsRfEndSputumSwabPs2)$Conditions <- meta(AbsRfEndSputumSwabPs2) %>%
  mutate(Conditions = ifelse(OriginType3 == "Swab", "Condition1", "Condition2")) %>% 
  pull(Conditions)

# Obtaining a data frame with all the data
SputumSwabCond1Cond2EndDF <- psmelt(AbsRfEndSputumSwabPs2) %>%
  arrange(Comparisons, Conditions, OTU, Abundance) %>% 
  select(Comparisons, Conditions, OTU, Abundance)

# Vector of comparisons to be made (i.e. distinct Comparisons)
ComparisonsEndSputumSwabVect <- SputumSwabCond1Cond2EndDF %>% 
  distinct(Comparisons) %>% 
  pull(Comparisons)

# Function to determine the distribution of ASVs between sputum and OPS (i.e. conditions)
# for each patient (i.e. comparisons). The function requires that the data be contained in a data frame
# called "SputumSwabCond1Cond2EndDF" and containing a "Comparisons" column and a "Conditions" column
ASVSputumSwabComparisonEnd <- function(Comparison){
  # Listing the names and counting ASVs represented in condition 1
  Cond1ASVNames <- SputumSwabCond1Cond2EndDF %>%
    filter(Comparisons == Comparison & Conditions == "Condition1" & Abundance > 0) %>% 
    pull(OTU)
  
  Cond1ASVNumber <- length(Cond1ASVNames)
  
  # Listing the names and counting ASVs represented in condition 2
  Cond2ASVNames <- SputumSwabCond1Cond2EndDF %>%
    filter(Comparisons == Comparison & Conditions == "Condition2" & Abundance > 0) %>% 
    pull(OTU)
  
  Cond2ASVNumber <- length(Cond2ASVNames)
  
  # Listing the names and counting ASVs represented in condition 1 ONLY
  Cond1OnlyASVNames <- setdiff(Cond1ASVNames, Cond2ASVNames)
  Cond1OnlyASVNumber <- length(Cond1OnlyASVNames)
  
  # Listing the names and counting ASVs represented in condition 2 ONLY
  Cond2OnlyASVNames <- setdiff(Cond2ASVNames, Cond1ASVNames)
  Cond2OnlyASVNumber <- length(Cond2OnlyASVNames)
  
  # Listing the names and counting ASVs represented in both condition 1 and condition 2
  SharedASVNames <- intersect(Cond1ASVNames, Cond2ASVNames)
  SharedASVNumber <- length(SharedASVNames)
  
  # Counting total ASVs per comparison
  Cond1Cond2TotalASVNumber <- Cond1OnlyASVNumber + Cond2OnlyASVNumber + SharedASVNumber
  
  # Calculating the percentage of ASVs in condition 1 ONLY
  Cond1OnlyPerCentASVNumber <- round(Cond1OnlyASVNumber / Cond1Cond2TotalASVNumber * 100, digits = 1)
  
  # Calculating the percentage of ASVs in condition 2 ONLY
  Cond2OnlyPerCentASVNumber <- round(Cond2OnlyASVNumber / Cond1Cond2TotalASVNumber * 100, digits = 1)
  
  # Calculating the percentage of shared ASVs
  PerCentSharedASVNumber <- round(SharedASVNumber / Cond1Cond2TotalASVNumber * 100, digits = 1)
  
  ComparisonOutput <- c(Cond1OnlyASVNumber,
                        SharedASVNumber,
                        Cond2OnlyASVNumber,
                        Cond1OnlyPerCentASVNumber,
                        PerCentSharedASVNumber,
                        Cond2OnlyPerCentASVNumber)
  
  return(ComparisonOutput)
}

# Obtaining a data frame containing the output of the function call 
ComparisonVect <- ComparisonsEndSputumSwabVect
ComparisonList <- as.list(ComparisonVect)

EndSputumSwabComparisonList <- map(ComparisonList, ASVSputumSwabComparisonEnd)

EndSputumSwabComparisonDF <- EndSputumSwabComparisonList %>%
  as.data.frame(col.names = ComparisonVect,
                row.names = c("Cond1OnlyASVNumber",
                              "SharedASVNumber",
                              "Cond2OnlyASVNumber",
                              "Cond1OnlyPerCentASVNumber",
                              "PerCentSharedASVNumber",
                              "Cond2OnlyPerCentASVNumber")) %>% 
  # Transposing and coercing to a data frame
  t() %>% 
  as.data.frame()

# Testing normality and plotting the distribution
shapiro.test(EndSputumSwabComparisonDF$Cond1OnlyASVNumber)

ggdensity(EndSputumSwabComparisonDF$Cond1OnlyASVNumber, 
          main = "Density plot of Cond1OnlyASVNumber",
          xlab = "Number of ASVs")

shapiro.test(EndSputumSwabComparisonDF$SharedASVNumber)

ggdensity(EndSputumSwabComparisonDF$SharedASVNumber, 
          main = "Density plot of SharedASVNumber",
          xlab = "Number of ASVs")

shapiro.test(EndSputumSwabComparisonDF$Cond2OnlyASVNumber)

ggdensity(EndSputumSwabComparisonDF$Cond2OnlyASVNumber, 
          main = "Density plot of Cond2OnlyASVNumber",
          xlab = "Number of ASVs")

# Obtaining a long data frame to compare medians
EndSputumSwabComparisonLongDF <- EndSputumSwabComparisonDF %>%
  select("Cond1OnlyASVNumber", "SharedASVNumber", "Cond2OnlyASVNumber") %>% 
  gather(., key = ASVDistribution, value = ASVNumber) %>% 
  convert(fct(ASVDistribution))

# Computing Kruskal-Wallis test
kruskal.test(ASVNumber ~ ASVDistribution,
             data = EndSputumSwabComparisonLongDF)

# Multiple pairwise-comparison between groups, with Benjamini & Hochberg correction for multiple testing
pairwise.wilcox.test(EndSputumSwabComparisonLongDF$ASVNumber,
                     EndSputumSwabComparisonLongDF$ASVDistribution,
                     p.adjust.method = "BH")

# Computing summary statistics per condition
EndSputumSwabComparisonLongPerCentDF <- EndSputumSwabComparisonDF %>%
  select("Cond1OnlyPerCentASVNumber", "PerCentSharedASVNumber", "Cond2OnlyPerCentASVNumber") %>% 
  gather(., key = ASVDistribution, value = ASVPerCent) %>% 
  convert(fct(ASVDistribution)) %>%
  mutate(ASVDistribution = factor(ASVDistribution,
                                  levels = c("Cond1OnlyPerCentASVNumber",
                                             "PerCentSharedASVNumber",
                                             "Cond2OnlyPerCentASVNumber")))

EndSputumSwabSummaryStatDF <- EndSputumSwabComparisonLongPerCentDF %>% 
  group_by(ASVDistribution) %>% 
  dplyr::summarize(medianASVPerCent = median(ASVPerCent),
            IQR = IQR(ASVPerCent)) %>% 
  ungroup()

# 4) Comparisons between OPS and Sputum per patient at 1 mo after AZM treatment
# Patients with a sample pair (OPS swab + Sputum) available at 1 mo after AZM treatment
AbsRfPostAZT1moSputumSwabPs2 <- PostAZT1moSputumSwabIncrStableARPs

# Adding a "Comparisons" column in sample_data object for later use in a function
sample_data(AbsRfPostAZT1moSputumSwabPs2)$Comparisons <- meta(AbsRfPostAZT1moSputumSwabPs2) %>%
  mutate(Comparisons = paste(PatientID, TTTPrePost7gr, sep = "_")) %>% 
  pull(Comparisons)

# Adding a "Conditions" column in sample_data object for later use in a function
sample_data(AbsRfPostAZT1moSputumSwabPs2)$Conditions <- meta(AbsRfPostAZT1moSputumSwabPs2) %>%
  mutate(Conditions = ifelse(OriginType3 == "Swab", "Condition1", "Condition2")) %>% 
  pull(Conditions)

# Obtaining a data frame with all the data
SputumSwabCond1Cond2PostAZT1moDF <- psmelt(AbsRfPostAZT1moSputumSwabPs2) %>%
  arrange(Comparisons, Conditions, OTU, Abundance) %>% 
  select(Comparisons, Conditions, OTU, Abundance)

# Vector of comparisons to be made (i.e. distinct Comparisons)
ComparisonsPostAZT1moSputumSwabVect <- SputumSwabCond1Cond2PostAZT1moDF %>% 
  distinct(Comparisons) %>% 
  pull(Comparisons)

# Function to determine the distribution of ASVs between sputum and OPS (i.e. conditions)
# for each patient (i.e. comparisons). The function requires that the data be contained in a data frame
# called "SputumSwabCond1Cond2PostAZT1moDF" and containing a "Comparisons" column and a "Conditions" column
ASVSputumSwabComparisonPostAZT1mo <- function(Comparison){
  # Listing the names and counting ASVs represented in condition 1
  Cond1ASVNames <- SputumSwabCond1Cond2PostAZT1moDF %>%
    filter(Comparisons == Comparison & Conditions == "Condition1" & Abundance > 0) %>% 
    pull(OTU)
  
  Cond1ASVNumber <- length(Cond1ASVNames)
  
  # Listing the names and counting ASVs represented in condition 2
  Cond2ASVNames <- SputumSwabCond1Cond2PostAZT1moDF %>%
    filter(Comparisons == Comparison & Conditions == "Condition2" & Abundance > 0) %>% 
    pull(OTU)
  
  Cond2ASVNumber <- length(Cond2ASVNames)
  
  # Listing the names and counting ASVs represented in condition 1 ONLY
  Cond1OnlyASVNames <- setdiff(Cond1ASVNames, Cond2ASVNames)
  Cond1OnlyASVNumber <- length(Cond1OnlyASVNames)
  
  # Listing the names and counting ASVs represented in condition 2 ONLY
  Cond2OnlyASVNames <- setdiff(Cond2ASVNames, Cond1ASVNames)
  Cond2OnlyASVNumber <- length(Cond2OnlyASVNames)
  
  # Listing the names and counting ASVs represented in both condition 1 and condition 2
  SharedASVNames <- intersect(Cond1ASVNames, Cond2ASVNames)
  SharedASVNumber <- length(SharedASVNames)
  
  # Counting total ASVs per comparison
  Cond1Cond2TotalASVNumber <- Cond1OnlyASVNumber + Cond2OnlyASVNumber + SharedASVNumber
  
  # Calculating the percentage of ASVs in condition 1 ONLY
  Cond1OnlyPerCentASVNumber <- round(Cond1OnlyASVNumber / Cond1Cond2TotalASVNumber * 100, digits = 1)
  
  # Calculating the percentage of ASVs in condition 2 ONLY
  Cond2OnlyPerCentASVNumber <- round(Cond2OnlyASVNumber / Cond1Cond2TotalASVNumber * 100, digits = 1)
  
  # Calculating the percentage of shared ASVs
  PerCentSharedASVNumber <- round(SharedASVNumber / Cond1Cond2TotalASVNumber * 100, digits = 1)
  
  ComparisonOutput <- c(Cond1OnlyASVNumber,
                        SharedASVNumber,
                        Cond2OnlyASVNumber,
                        Cond1OnlyPerCentASVNumber,
                        PerCentSharedASVNumber,
                        Cond2OnlyPerCentASVNumber)
  
  return(ComparisonOutput)
}

# Obtaining a data frame containing the output of the function call 
ComparisonVect <- ComparisonsPostAZT1moSputumSwabVect
ComparisonList <- as.list(ComparisonVect)

PostAZT1moSputumSwabComparisonList <- map(ComparisonList, ASVSputumSwabComparisonPostAZT1mo)

PostAZT1moSputumSwabComparisonDF <- PostAZT1moSputumSwabComparisonList %>%
  as.data.frame(col.names = ComparisonVect,
                row.names = c("Cond1OnlyASVNumber",
                              "SharedASVNumber",
                              "Cond2OnlyASVNumber",
                              "Cond1OnlyPerCentASVNumber",
                              "PerCentSharedASVNumber",
                              "Cond2OnlyPerCentASVNumber")) %>% 
  # Transposing and coercing to a data frame
  t() %>% 
  as.data.frame()

# Testing normality and plotting the distribution
shapiro.test(PostAZT1moSputumSwabComparisonDF$Cond1OnlyASVNumber)

ggdensity(PostAZT1moSputumSwabComparisonDF$Cond1OnlyASVNumber, 
          main = "Density plot of Cond1OnlyASVNumber",
          xlab = "Number of ASVs")

shapiro.test(PostAZT1moSputumSwabComparisonDF$SharedASVNumber)

ggdensity(PostAZT1moSputumSwabComparisonDF$SharedASVNumber, 
          main = "Density plot of SharedASVNumber",
          xlab = "Number of ASVs")

shapiro.test(PostAZT1moSputumSwabComparisonDF$Cond2OnlyASVNumber)

ggdensity(PostAZT1moSputumSwabComparisonDF$Cond2OnlyASVNumber, 
          main = "Density plot of Cond2OnlyASVNumber",
          xlab = "Number of ASVs")

# Obtaining a long data frame to compare medians
PostAZT1moSputumSwabComparisonLongDF <- PostAZT1moSputumSwabComparisonDF %>%
  select("Cond1OnlyASVNumber", "SharedASVNumber", "Cond2OnlyASVNumber") %>% 
  gather(., key = ASVDistribution, value = ASVNumber) %>% 
  convert(fct(ASVDistribution))
PostAZT1moSputumSwabComparisonLongDF %>% print.data.frame()
# Computing Kruskal-Wallis test
kruskal.test(ASVNumber ~ ASVDistribution,
             data = PostAZT1moSputumSwabComparisonLongDF)

# Multiple pairwise-comparison between groups, with Benjamini & Hochberg correction for multiple testing
pairwise.wilcox.test(PostAZT1moSputumSwabComparisonLongDF$ASVNumber,
                     PostAZT1moSputumSwabComparisonLongDF$ASVDistribution,
                     p.adjust.method = "BH")

# Computing summary statistics per condition
PostAZT1moSputumSwabComparisonLongPerCentDF <- PostAZT1moSputumSwabComparisonDF %>%
  select("Cond1OnlyPerCentASVNumber", "PerCentSharedASVNumber", "Cond2OnlyPerCentASVNumber") %>% 
  gather(., key = ASVDistribution, value = ASVPerCent) %>% 
  convert(fct(ASVDistribution)) %>%
  mutate(ASVDistribution = factor(ASVDistribution,
                                  levels = c("Cond1OnlyPerCentASVNumber",
                                             "PerCentSharedASVNumber",
                                             "Cond2OnlyPerCentASVNumber")))

PostAZT1moSputumSwabSummaryStatDF <- PostAZT1moSputumSwabComparisonLongPerCentDF %>% 
  group_by(ASVDistribution) %>% 
  dplyr::summarize(medianASVPerCent = median(ASVPerCent),
            IQR = IQR(ASVPerCent)) %>% 
  ungroup()

# 5) Comparisons between OPS and Sputum per patient at >=4 mo after AZM treatment
# Patients with a sample pair (OPS swab + Sputum) available at >=4 mo after AZM treatment
AbsRfPostAZT4moSputumSwabPs2 <- PostAZT4moSputumSwabIncrStableARPs

# Adding a "Comparisons" column in sample_data object for later use in a function
sample_data(AbsRfPostAZT4moSputumSwabPs2)$Comparisons <- meta(AbsRfPostAZT4moSputumSwabPs2) %>%
  mutate(Comparisons = paste(PatientID, TTTPrePost7gr, sep = "_")) %>% 
  pull(Comparisons)

# Adding a "Conditions" column in sample_data object for later use in a function
sample_data(AbsRfPostAZT4moSputumSwabPs2)$Conditions <- meta(AbsRfPostAZT4moSputumSwabPs2) %>%
  mutate(Conditions = ifelse(OriginType3 == "Swab", "Condition1", "Condition2")) %>% 
  pull(Conditions)

# Obtaining a data frame with all the data
SputumSwabCond1Cond2PostAZT4moDF <- psmelt(AbsRfPostAZT4moSputumSwabPs2) %>%
  arrange(Comparisons, Conditions, OTU, Abundance) %>% 
  select(Comparisons, Conditions, OTU, Abundance)

# Vector of comparisons to be made (i.e. distinct Comparisons)
ComparisonsPostAZT4moSputumSwabVect <- SputumSwabCond1Cond2PostAZT4moDF %>% 
  distinct(Comparisons) %>% 
  pull(Comparisons)

# Function to determine the distribution of ASVs between sputum and OPS (i.e. conditions)
# for each patient (i.e. comparisons). The function requires that the data be contained in a data frame
# called "SputumSwabCond1Cond2PostAZT4moDF" and containing a "Comparisons" column and a "Conditions" column
ASVSputumSwabComparisonPostAZT4mo <- function(Comparison){
  # Listing the names and counting ASVs represented in condition 1
  Cond1ASVNames <- SputumSwabCond1Cond2PostAZT4moDF %>%
    filter(Comparisons == Comparison & Conditions == "Condition1" & Abundance > 0) %>% 
    pull(OTU)
  
  Cond1ASVNumber <- length(Cond1ASVNames)
  
  # Listing the names and counting ASVs represented in condition 2
  Cond2ASVNames <- SputumSwabCond1Cond2PostAZT4moDF %>%
    filter(Comparisons == Comparison & Conditions == "Condition2" & Abundance > 0) %>% 
    pull(OTU)
  
  Cond2ASVNumber <- length(Cond2ASVNames)
  
  # Listing the names and counting ASVs represented in condition 1 ONLY
  Cond1OnlyASVNames <- setdiff(Cond1ASVNames, Cond2ASVNames)
  Cond1OnlyASVNumber <- length(Cond1OnlyASVNames)
  
  # Listing the names and counting ASVs represented in condition 2 ONLY
  Cond2OnlyASVNames <- setdiff(Cond2ASVNames, Cond1ASVNames)
  Cond2OnlyASVNumber <- length(Cond2OnlyASVNames)
  
  # Listing the names and counting ASVs represented in both condition 1 and condition 2
  SharedASVNames <- intersect(Cond1ASVNames, Cond2ASVNames)
  SharedASVNumber <- length(SharedASVNames)
  
  # Counting total ASVs per comparison
  Cond1Cond2TotalASVNumber <- Cond1OnlyASVNumber + Cond2OnlyASVNumber + SharedASVNumber
  
  # Calculating the percentage of ASVs in condition 1 ONLY
  Cond1OnlyPerCentASVNumber <- round(Cond1OnlyASVNumber / Cond1Cond2TotalASVNumber * 100, digits = 1)
  
  # Calculating the percentage of ASVs in condition 2 ONLY
  Cond2OnlyPerCentASVNumber <- round(Cond2OnlyASVNumber / Cond1Cond2TotalASVNumber * 100, digits = 1)
  
  # Calculating the percentage of shared ASVs
  PerCentSharedASVNumber <- round(SharedASVNumber / Cond1Cond2TotalASVNumber * 100, digits = 1)
  
  ComparisonOutput <- c(Cond1OnlyASVNumber,
                        SharedASVNumber,
                        Cond2OnlyASVNumber,
                        Cond1OnlyPerCentASVNumber,
                        PerCentSharedASVNumber,
                        Cond2OnlyPerCentASVNumber)
  
  return(ComparisonOutput)
}

# Obtaining a data frame containing the output of the function call 
ComparisonVect <- ComparisonsPostAZT4moSputumSwabVect
ComparisonList <- as.list(ComparisonVect)

PostAZT4moSputumSwabComparisonList <- map(ComparisonList, ASVSputumSwabComparisonPostAZT4mo)

PostAZT4moSputumSwabComparisonDF <- PostAZT4moSputumSwabComparisonList %>%
  as.data.frame(col.names = ComparisonVect,
                row.names = c("Cond1OnlyASVNumber",
                              "SharedASVNumber",
                              "Cond2OnlyASVNumber",
                              "Cond1OnlyPerCentASVNumber",
                              "PerCentSharedASVNumber",
                              "Cond2OnlyPerCentASVNumber")) %>% 
  # Transposing and coercing to a data frame
  t() %>% 
  as.data.frame()

# Testing normality and plotting the distribution
shapiro.test(PostAZT4moSputumSwabComparisonDF$Cond1OnlyASVNumber)

ggdensity(PostAZT4moSputumSwabComparisonDF$Cond1OnlyASVNumber, 
          main = "Density plot of Cond1OnlyASVNumber",
          xlab = "Number of ASVs")

shapiro.test(PostAZT4moSputumSwabComparisonDF$SharedASVNumber)

ggdensity(PostAZT4moSputumSwabComparisonDF$SharedASVNumber, 
          main = "Density plot of SharedASVNumber",
          xlab = "Number of ASVs")

shapiro.test(PostAZT4moSputumSwabComparisonDF$Cond2OnlyASVNumber)

ggdensity(PostAZT4moSputumSwabComparisonDF$Cond2OnlyASVNumber, 
          main = "Density plot of Cond2OnlyASVNumber",
          xlab = "Number of ASVs")

# Obtaining a long data frame to compare medians
PostAZT4moSputumSwabComparisonLongDF <- PostAZT4moSputumSwabComparisonDF %>%
  select("Cond1OnlyASVNumber", "SharedASVNumber", "Cond2OnlyASVNumber") %>% 
  gather(., key = ASVDistribution, value = ASVNumber) %>% 
  convert(fct(ASVDistribution))

# Computing Kruskal-Wallis test
kruskal.test(ASVNumber ~ ASVDistribution,
             data = PostAZT4moSputumSwabComparisonLongDF)

# Multiple pairwise-comparison between groups, with Benjamini & Hochberg correction for multiple testing
pairwise.wilcox.test(PostAZT4moSputumSwabComparisonLongDF$ASVNumber,
                     PostAZT4moSputumSwabComparisonLongDF$ASVDistribution,
                     p.adjust.method = "BH")

# Computing summary statistics per condition
PostAZT4moSputumSwabComparisonLongPerCentDF <- PostAZT4moSputumSwabComparisonDF %>%
  select("Cond1OnlyPerCentASVNumber", "PerCentSharedASVNumber", "Cond2OnlyPerCentASVNumber") %>% 
  gather(., key = ASVDistribution, value = ASVPerCent) %>% 
  convert(fct(ASVDistribution)) %>%
  mutate(ASVDistribution = factor(ASVDistribution,
                                  levels = c("Cond1OnlyPerCentASVNumber",
                                             "PerCentSharedASVNumber",
                                             "Cond2OnlyPerCentASVNumber")))

PostAZT4moSputumSwabSummaryStatDF <- PostAZT4moSputumSwabComparisonLongPerCentDF %>% 
  group_by(ASVDistribution) %>% 
  dplyr::summarize(medianASVPerCent = median(ASVPerCent),
            IQR = IQR(ASVPerCent)) %>% 
  ungroup()

# Microbiota dynamics in OPS (Figure 3d) ----
# Analyses on all patients with available pairs of OPS samples independent of AR gene carriage
# 1A) Determining the proportion of taxa gained/conserved/lost between StartATZ and EndAZT
# Patients with a pair of OPS samples available
PatientsStartEndSwab <- meta(AbsRfPatientPs) %>% 
  filter(Type == "Swab") %>%
  filter(TTTPrePost5gr == "StartAZT" | TTTPrePost7gr == "EndAZT") %>% 
  group_by(PatientID) %>% 
  dplyr::summarize(NumberSamplesAvailable = n()) %>% 
  filter(NumberSamplesAvailable == 2) %>% 
  select(PatientID)

PatientsStartEndSwab$PatientID %<>% droplevels()

AbsRfSwabStartEndPs <- subset_samples(AbsRfPatientPs,
                                      PatientID %in% PatientsStartEndSwab$PatientID &
                                        Type == "Swab" &
                                        TTTPrePost5gr %in% c("StartAZT",
                                                             "EndAZT"))

# Keeping only the taxa represented after filtering, then counting them and listing their names
AbsRfSwabStartEndPs %<>% prune_taxa(taxa_sums(.) > 0, .)
SwabStartEndNumber <- AbsRfSwabStartEndPs %>% ntaxa()
SwabStartEndNames <- AbsRfSwabStartEndPs %>% taxa_names()
# 1B) Inspecting taxa present in ORO swabs at StartAZT
AbsRfSwabStartAZTPs <- subset_samples(AbsRfSwabStartEndPs,
                                      TTTPrePost5gr == "StartAZT")
# Keeping only the taxa represented after filtering, then counting them and listing their names
AbsRfSwabStartAZTPs %<>% prune_taxa(taxa_sums(.) > 0, .)
SwabStartNumber <- AbsRfSwabStartAZTPs %>% ntaxa()
SwabStartNames <- AbsRfSwabStartAZTPs %>% taxa_names()
# 1C) Inspecting taxa present in ORO swabs at EndAZT
AbsRfSwabEndAZTPs <- subset_samples(AbsRfSwabStartEndPs,
                                    TTTPrePost5gr == "EndAZT")
# Keeping only the taxa represented after filtering, then counting them and listing their names
AbsRfSwabEndAZTPs %<>% prune_taxa(taxa_sums(.) > 0, .)
SwabEndNumber <- AbsRfSwabEndAZTPs %>% ntaxa()
SwabEndNames <- AbsRfSwabEndAZTPs %>% taxa_names()

# Preparing a list of vectors with taxa names for visualisation by Venn diagram
SwabStartEndAZTList <- list(StartAZT = SwabStartNames,
                            EndAZT = SwabEndNames)
# Calculating the coordinates required for plotting a Venn diagram
SwabStartEndAZTVenn <- Venn(SwabStartEndAZTList)
SwabStartEndAZTVennData <- process_data(SwabStartEndAZTVenn)
# Obtaining the Venn diagram (For Figure 3d)
SwabStartEndAZTVennDiagram <- ggVennDiagram(SwabStartEndAZTList,
                                            set_color = "white", # Hiding original label
                                            label =  "both", # Showing "count" and "percent"
                                            label_size = 10,
                                            label_alpha = .6,
                                            edge_lty = "blank") +
  # to customize the label layer
  geom_sf_text(aes(label = name),
               fontface = "plain",
               size = 9,
               data = venn_setlabel(SwabStartEndAZTVennData)) +
  # to customize the edge layer
  geom_sf(size = .1,
          color = "#ff7d1a",
          alpha = .6,
          data = venn_setedge(SwabStartEndAZTVennData)) +
  scale_fill_gradient(low = "#d2ebf9", high = "#E69F00") +
  theme_void() +
  theme(legend.position = "none")

ggsave("SwabStartEndAZTVennDiagram.pdf",
       width = 12, height = 10, dpi = 200, units = "cm", limitsize = FALSE)

# 2A) Determining the proportion of taxa gained/conserved/lost between PostAZT_1mo and PostAZT>=3mo
# Patients with a pair of OPS samples available
PatientsPost1moPost4moSwab <- meta(AbsRfPatientPs) %>% 
  filter(Type == "Swab") %>%
  filter(TTTPrePost7gr == "PostAZT_1mo" | TTTPrePost7gr == "PostAZT_4mo") %>% 
  group_by(PatientID) %>% 
  dplyr::summarize(NumberSamplesAvailable = n()) %>% 
  filter(NumberSamplesAvailable >= 2) %>% 
  select(PatientID)

PatientsPost1moPost4moSwab$PatientID %<>% droplevels()

AbsRfSwabPost1moPost4moPs <- subset_samples(AbsRfPatientPs,
                                            PatientID %in% PatientsPost1moPost4moSwab$PatientID &
                                              Type == "Swab" &
                                              TTTPrePost7gr %in% c("PostAZT_1mo",
                                                                   "PostAZT_4mo"))

# Keeping only the taxa represented after filtering, then counting them and listing their names
AbsRfSwabPost1moPost4moPs %<>% prune_taxa(taxa_sums(.) > 0, .)
SwabPost1moPost4moNumber <- AbsRfSwabPost1moPost4moPs %>% ntaxa()
SwabPost1moPost4moNames <- AbsRfSwabPost1moPost4moPs %>% taxa_names()
# 2B) Inspecting taxa present in ORO swabs at PostAZT_1mo
AbsRfSwabPost1moPs <- subset_samples(AbsRfSwabPost1moPost4moPs,
                                     TTTPrePost7gr == "PostAZT_1mo")
# Keeping only the taxa represented after filtering, then counting them and listing their names
AbsRfSwabPost1moPs %<>% prune_taxa(taxa_sums(.) > 0, .)
SwabPost1moNumber <- AbsRfSwabPost1moPs %>% ntaxa()
SwabPost1moNames <- AbsRfSwabPost1moPs %>% taxa_names()
# 2C) Inspecting taxa present in ORO swabs at PostAZT_4mo
AbsRfSwabPost4moPs <- subset_samples(AbsRfSwabPost1moPost4moPs,
                                     TTTPrePost7gr == "PostAZT_4mo")
# Keeping only the taxa represented after filtering, then counting them and listing their names
AbsRfSwabPost4moPs %<>% prune_taxa(taxa_sums(.) > 0, .)
SwabPost4moNumber <- AbsRfSwabPost4moPs %>% ntaxa()
SwabPost4moNames <- AbsRfSwabPost4moPs %>% taxa_names()

# Preparing a list of vectors with taxa names for visualisation by Venn diagram
SwabPost1moPost4moList <- list(PostAZT_1mo = SwabPost1moNames,
                               PostAZT_4mo = SwabPost4moNames)
# Calculating the coordinates required for plotting a Venn diagram
SwabPost1moPost4moVenn <- Venn(SwabPost1moPost4moList)
SwabPost1moPost4moVennData <- process_data(SwabPost1moPost4moVenn)
# Obtaining the Venn diagram (For Figure 3d)
SwabPost1moPost4moVennDiagram <- ggVennDiagram(SwabPost1moPost4moList,
                                               set_color = "white", # Hiding original label
                                               label =  "both", # Showing "count" and "percent"
                                               label_size = 10,
                                               label_alpha = .6,
                                               edge_lty = "blank") +
  # to customize the label layer
  geom_sf_text(aes(label = name),
               fontface = "plain",
               size = 9,
               data = venn_setlabel(SwabPost1moPost4moVennData)) +
  # to customize the edge layer
  geom_sf(size = .1,
          color = "#ff7d1a",
          alpha = .6,
          data = venn_setedge(SwabPost1moPost4moVennData)) +
  scale_fill_gradient(low = "#00CC96", high = "#CC79A7") +
  theme_void() +
  theme(legend.position = "none")

ggsave("SwabPost1moPost4moVennDiagram.pdf",
       width = 12, height = 10, dpi = 200, units = "cm", limitsize = FALSE)

# 3A) Determining the proportion of taxa gained/conserved/lost between PreAZT_4mo and PreAZT_1mo
# in patients receiving placebo first.
# Keeping only patients with a pair of OPS samples available
PatientsPlaceboFirstSwab <- meta(AbsRfPatientPs) %>% 
  filter(Type == "Swab") %>%
  filter(TTTPrePost7gr == "PreAZT_4mo" | TTTPrePost7gr == "PreAZT_1mo") %>% 
  group_by(PatientID) %>% 
  dplyr::summarize(NumberSamplesAvailable = n()) %>% 
  filter(NumberSamplesAvailable == 2) %>% 
  select(PatientID)

PatientsPlaceboFirstSwab$PatientID %<>% droplevels()

AbsRfSwabPlaceboFirst4Mo1MoPs <- subset_samples(AbsRfPatientPs,
                                                PatientID %in% PatientsPlaceboFirstSwab$PatientID &
                                                  Type == "Swab" &
                                                  TTTPrePost7gr %in% c("PreAZT_1mo",
                                                                       "PreAZT_4mo"))
# Keeping only the taxa represented after filtering, then counting them and listing their names
AbsRfSwabPlaceboFirst4Mo1MoPs %<>% prune_taxa(taxa_sums(.) > 0, .)
SwabPlaceboFirst4Mo1MoNumber <- AbsRfSwabPlaceboFirst4Mo1MoPs %>% ntaxa()
SwabPlaceboFirst4Mo1ModNames <- AbsRfSwabPlaceboFirst4Mo1MoPs %>% taxa_names()
# 3B) Inspecting taxa present at PreAZT_4mo in the sputum of patients receiving placebo first
AbsRfSwabPlaceboFirstPre4MoPs <- subset_samples(AbsRfSwabPlaceboFirst4Mo1MoPs,
                                                TTTPrePost7gr == "PreAZT_4mo")
# Keeping only the taxa represented after filtering, then counting them and listing their names
AbsRfSwabPlaceboFirstPre4MoPs %<>% prune_taxa(taxa_sums(.) > 0, .)
SwabPlaceboFirst4MoNumber <- AbsRfSwabPlaceboFirstPre4MoPs %>% ntaxa()
SwabPlaceboFirst4MoNames <- AbsRfSwabPlaceboFirstPre4MoPs %>% taxa_names()
# 3C) Inspecting taxa present at PreAZT_1mo in the sputum of patients receiving placebo first
AbsRfSwabPlaceboFirstPre1MoPs <- subset_samples(AbsRfSwabPlaceboFirst4Mo1MoPs,
                                                TTTPrePost7gr == "PreAZT_1mo")
# Keeping only the taxa represented after filtering, then counting them and listing their names
AbsRfSwabPlaceboFirstPre1MoPs %<>% prune_taxa(taxa_sums(.) > 0, .)
SwabPlaceboFirst1MoNumber <- AbsRfSwabPlaceboFirstPre1MoPs %>% ntaxa()
SwabPlaceboFirst1MoNames <- AbsRfSwabPlaceboFirstPre1MoPs %>% taxa_names()

# Preparing a list of vectors with taxa names for visualisation by Venn diagram
SwabPlaceboFirst4Mo1MoList <- list("PreTTT_4mo" = SwabPlaceboFirst4MoNames,
                                   "PreTTT_1moEndAZT" = SwabPlaceboFirst1MoNames)
# Calculating the coordinates required for plotting
SwabPlaceboFirst4Mo1MoVenn <- Venn(SwabPlaceboFirst4Mo1MoList)
SwabPlaceboFirst4Mo1MoVennData <- process_data(SwabPlaceboFirst4Mo1MoVenn)
# Obtaining the Venn diagram (For Figure 3d)
SwabPlaceboFirst4Mo1MoVennDiagram <- ggVennDiagram(SwabPlaceboFirst4Mo1MoList,
                                                   set_color = "white", # For hiding original label
                                                   label =  "both", # For showing "count" and "percent"
                                                   label_size = 10, # For region label size
                                                   label_alpha = .6, # For adding transparency to region labels
                                                   edge_lty = "blank") +
  # to customize the label layer
  geom_sf_text(aes(label = name),
               fontface = "plain",
               size = 9,
               data = venn_setlabel(SwabPlaceboFirst4Mo1MoVennData)) +
  # to customize the edge layer
  geom_sf(size = .1,
          color = "#f2f2f2",
          alpha = .6,
          data = venn_setedge(SwabPlaceboFirst4Mo1MoVennData)) +
  scale_fill_gradient(low = "#f0f0f0", high = "#999999") +
  theme_void() +
  theme(legend.position = "none")

ggsave("SwabPlaceboFirst4Mo1MoVennDiagram.pdf",
       width = 12, height = 10, dpi = 200, units = "cm", limitsize = FALSE)

# Microbiota dynamics in Sputum (Figure 3d) ----
# Analyses on all patients with available pairs of sputum samples independent of AR gene carriage
# 1A) Determining the proportion of taxa gained/conserved/lost between StartATZ and EndAZT
# Patients with a pair of sputum samples available
PatientsStartEndSputum <- meta(AbsRfPatientPs) %>% 
  filter(Type == "Sputum") %>%
  filter(TTTPrePost5gr == "StartAZT" | TTTPrePost7gr == "EndAZT") %>% 
  group_by(PatientID) %>% 
  dplyr::summarize(NumberSamplesAvailable = n()) %>% 
  filter(NumberSamplesAvailable == 2) %>% 
  select(PatientID)

PatientsStartEndSputum$PatientID %<>% droplevels()

AbsRfSputumStartEndPs <- subset_samples(AbsRfPatientPs,
                                        PatientID %in% PatientsStartEndSputum$PatientID &
                                          Type == "Sputum" &
                                          TTTPrePost5gr %in% c("StartAZT",
                                                               "EndAZT"))
# Keeping only the taxa represented after filtering, then counting them and listing their names
AbsRfSputumStartEndPs %<>% prune_taxa(taxa_sums(.) > 0, .)
SputumStartEndNumber <- AbsRfSputumStartEndPs %>% ntaxa()
SputumStartEndNames <- AbsRfSputumStartEndPs %>% taxa_names()
# 1B) Inspecting taxa present in sputum samples at StartAZT
AbsRfSputumStartAZTPs <- subset_samples(AbsRfSputumStartEndPs,
                                        TTTPrePost5gr == "StartAZT")
# Keeping only the taxa represented after filtering, then counting them and listing their names
AbsRfSputumStartAZTPs %<>% prune_taxa(taxa_sums(.) > 0, .)
SputumStartNumber <- AbsRfSputumStartAZTPs %>% ntaxa()
SputumStartNames <- AbsRfSputumStartAZTPs %>% taxa_names()
# 1C) Inspecting taxa present in sputum samples at EndAZT
AbsRfSputumEndAZTPs <- subset_samples(AbsRfSputumStartEndPs,
                                      TTTPrePost5gr == "EndAZT")
# Keeping only the taxa represented after filtering, then counting them and listing their names
AbsRfSputumEndAZTPs %<>% prune_taxa(taxa_sums(.) > 0, .)
SputumEndNumber <- AbsRfSputumEndAZTPs %>% ntaxa()
SputumEndNames <- AbsRfSputumEndAZTPs %>% taxa_names()

# Preparing a list of vectors with taxa names for visualisation by Venn diagram
SputumStartEndAZTList <- list(StartAZT = SputumStartNames,
                              EndAZT = SputumEndNames)
# Calculating the coordinates required for plotting a Venn diagram
SputumStartEndAZTVenn <- ggVennDiagram::Venn(SputumStartEndAZTList)
SputumStartEndAZTVennData <- ggVennDiagram::process_data(SputumStartEndAZTVenn)
# Obtaining the Venn diagram (For Figure 3d)
SputumStartEndAZTVennDiagram <- ggVennDiagram(SputumStartEndAZTList,
                                              set_color = "white", # Hiding original label
                                              label =  "both", # Showing "count" and "percent"
                                              label_size = 10,
                                              label_alpha = .6,
                                              edge_lty = "blank") +
  # to customize the label layer
  geom_sf_text(aes(label = name),
               fontface = "plain",
               size = 9,
               data = venn_setlabel(SputumStartEndAZTVennData)) +
  # to customize the edge layer
  geom_sf(size = .1,
          color = "#ff7d1a",
          alpha = .6,
          data = venn_setedge(SputumStartEndAZTVennData)) +
  scale_fill_gradient(low = "#d2ebf9", high = "#E69F00") +
  theme_void() +
  theme(legend.position = "none")

ggsave("SputumStartEndAZTVennDiagram.pdf",
       width = 12, height = 10, dpi = 200, units = "cm", limitsize = FALSE)

# 2A) Determining the proportion of taxa gained/conserved/lost between PostAZT_1mo and PostAZT>=3mo
# Patients with a pair of sputum samples available
PatientsPost1moPost4moSputum <- meta(AbsRfPatientPs) %>% 
  filter(Type == "Sputum") %>%
  filter(TTTPrePost7gr == "PostAZT_1mo" | TTTPrePost7gr == "PostAZT_4mo") %>% 
  group_by(PatientID) %>% 
  dplyr::summarize(NumberSamplesAvailable = n()) %>% 
  filter(NumberSamplesAvailable >= 2) %>% 
  select(PatientID)

PatientsPost1moPost4moSputum$PatientID %<>% droplevels()

AbsRfSputumPost1moPost4moPs <- subset_samples(AbsRfPatientPs,
                                              PatientID %in% PatientsPost1moPost4moSputum$PatientID &
                                                Type == "Sputum" &
                                                TTTPrePost7gr %in% c("PostAZT_1mo",
                                                                     "PostAZT_4mo"))

# Keeping only the taxa represented after filtering, then counting them and listing their names
AbsRfSputumPost1moPost4moPs %<>% prune_taxa(taxa_sums(.) > 0, .)
SputumPost1moPost4moNumber <- AbsRfSputumPost1moPost4moPs %>% ntaxa()
SputumPost1moPost4moNames <- AbsRfSputumPost1moPost4moPs %>% taxa_names()
# 2B) Inspecting taxa present in sputum samples at PostAZT_1mo
AbsRfSputumPost1moPs <- subset_samples(AbsRfSputumPost1moPost4moPs,
                                       TTTPrePost7gr == "PostAZT_1mo")
# Keeping only the taxa represented after filtering, then counting them and listing their names
AbsRfSputumPost1moPs %<>% prune_taxa(taxa_sums(.) > 0, .)
SputumPost1moNumber <- AbsRfSputumPost1moPs %>% ntaxa()
SputumPost1moNames <- AbsRfSputumPost1moPs %>% taxa_names()
# 2C) Inspecting taxa present in sputum samples at PostAZT_4mo
AbsRfSputumPost4moPs <- subset_samples(AbsRfSputumPost1moPost4moPs,
                                       TTTPrePost7gr == "PostAZT_4mo")
# Keeping only the taxa represented after filtering, then counting them and listing their names
AbsRfSputumPost4moPs %<>% prune_taxa(taxa_sums(.) > 0, .)
SputumPost4moNumber <- AbsRfSputumPost4moPs %>% ntaxa()
SputumPost4moNames <- AbsRfSputumPost4moPs %>% taxa_names()

# Preparing a list of vectors with taxa names for visualisation by Venn diagram
SputumPost1moPost4moList <- list(PostAZT_1mo = SputumPost1moNames,
                                 PostAZT_4mo = SputumPost4moNames)
# Calculating the coordinates required for plotting a Venn diagram
SputumPost1moPost4moVenn <- Venn(SputumPost1moPost4moList)
SputumPost1moPost4moVennData <- process_data(SputumPost1moPost4moVenn)
# Obtaining the Venn diagram (For Figure 3d)
SputumPost1moPost4moVennDiagram <- ggVennDiagram(SputumPost1moPost4moList,
                                                 set_color = "white", # Hiding original label
                                                 label =  "both", # Showing "count" and "percent"
                                                 label_size = 10,
                                                 label_alpha = .6,
                                                 edge_lty = "blank") +
  # to customize the label layer
  geom_sf_text(aes(label = name),
               fontface = "plain",
               size = 9,
               data = venn_setlabel(SputumPost1moPost4moVennData)) +
  # to customize the edge layer
  geom_sf(size = .1,
          color = "#ff7d1a",
          alpha = .6,
          data = venn_setedge(SputumPost1moPost4moVennData)) +
  scale_fill_gradient(low = "#00CC96", high = "#CC79A7") +
  theme_void() +
  theme(legend.position = "none")

ggsave("SputumPost1moPost4moVennDiagram.pdf",
       width = 12, height = 10, dpi = 200, units = "cm", limitsize = FALSE)

# 3A) Determining the proportion of taxa gained/conserved/lost between PreAZT_4mo and PreAZT_1mo
# in patients receiving placebo first.
# Keeping only patients with a pair of sputum samples available
PatientsPlaceboFirst <- meta(AbsRfPatientPs) %>% 
  filter(Type == "Sputum") %>%
  filter(TTTPrePost7gr == "PreAZT_4mo" | TTTPrePost7gr == "PreAZT_1mo") %>% 
  group_by(PatientID) %>% 
  dplyr::summarize(NumberSamplesAvailable = n()) %>% 
  filter(NumberSamplesAvailable == 2) %>% 
  select(PatientID)

PatientsPlaceboFirst$PatientID %<>% droplevels()

AbsRfSputumPlaceboFirst4Mo1MoPs <- subset_samples(AbsRfPatientPs,
                                                  PatientID %in% PatientsPlaceboFirst$PatientID &
                                                    Type == "Sputum" &
                                                    TTTPrePost7gr %in% c("PreAZT_1mo",
                                                                         "PreAZT_4mo"))
# Keeping only the taxa represented after filtering, then counting them and listing their names
AbsRfSputumPlaceboFirst4Mo1MoPs %<>% prune_taxa(taxa_sums(.) > 0, .)
SputumPlaceboFirst4Mo1MoNumber <- AbsRfSputumPlaceboFirst4Mo1MoPs %>% ntaxa()
SputumPlaceboFirst4Mo1ModNames <- AbsRfSputumPlaceboFirst4Mo1MoPs %>% taxa_names()
# 3B) Inspecting taxa present at PreAZT_4mo in the sputum of patients receiving placebo first
AbsRfSputumPlaceboFirstPre4MoPs <- subset_samples(AbsRfSputumPlaceboFirst4Mo1MoPs,
                                                  TTTPrePost7gr == "PreAZT_4mo")
# Keeping only the taxa represented after filtering, then counting them and listing their names
AbsRfSputumPlaceboFirstPre4MoPs %<>% prune_taxa(taxa_sums(.) > 0, .)
SputumPlaceboFirst4MoNumber <- AbsRfSputumPlaceboFirstPre4MoPs %>% ntaxa()
SputumPlaceboFirst4MoNames <- AbsRfSputumPlaceboFirstPre4MoPs %>% taxa_names()
# 3C) Inspecting taxa present at PreAZT_1mo in the sputum of patients receiving placebo first
AbsRfSputumPlaceboFirstPre1MoPs <- subset_samples(AbsRfSputumPlaceboFirst4Mo1MoPs,
                                                  TTTPrePost7gr == "PreAZT_1mo")
# Keeping only the taxa represented after filtering, then counting them and listing their names
AbsRfSputumPlaceboFirstPre1MoPs %<>% prune_taxa(taxa_sums(.) > 0, .)
SputumPlaceboFirst1MoNumber <- AbsRfSputumPlaceboFirstPre1MoPs %>% ntaxa()
SputumPlaceboFirst1MoNames <- AbsRfSputumPlaceboFirstPre1MoPs %>% taxa_names()

# Preparing a list of vectors with taxa names for visualisation by Venn diagram
SputumPlaceboFirst4Mo1MoList <- list("PreTTT_4mo" = SputumPlaceboFirst4MoNames,
                                     "PreTTT_1moEndAZT" = SputumPlaceboFirst1MoNames)
# Calculating the coordinates required for plotting a Venn diagram
SputumPlaceboFirst4Mo1MoVenn <- Venn(SputumPlaceboFirst4Mo1MoList)
SputumPlaceboFirst4Mo1MoVennData <- process_data(SputumPlaceboFirst4Mo1MoVenn)
# Obtaining the Venn diagram (For Figure 3d)
SputumPlaceboFirst4Mo1MoVennDiagram <- ggVennDiagram(SputumPlaceboFirst4Mo1MoList,
                                                     set_color = "white", # Hiding original label
                                                     label =  "both", # Showing "count" and "percent"
                                                     label_size = 10,
                                                     label_alpha = .6,
                                                     edge_lty = "blank") +
  # to customize the label layer
  geom_sf_text(aes(label = name),
               fontface = "plain",
               size = 9,
               data = venn_setlabel(SputumPlaceboFirst4Mo1MoVennData)) +
  # to customize the edge layer
  geom_sf(size = .1,
          color = "#f2f2f2",
          alpha = .6,
          data = venn_setedge(SputumPlaceboFirst4Mo1MoVennData)) +
  scale_fill_gradient(low = "#f0f0f0", high = "#999999") +
  theme_void() +
  theme(legend.position = "none")

ggsave("SputumPlaceboFirst4Mo1MoVennDiagram.pdf",
       width = 12, height = 10, dpi = 200, units = "cm", limitsize = FALSE)

# Abundance, prevalence and identity of bacteria present at StartAZT and/or EndAZT ----
# Working with Relative abundance
RelRfSputumAllStartEndPs <- microbiome::transform(AbsRfSputumAllStartEndPs, "compositional")
# Removing the phy_tree
RelRfSputumAllStartEndNoTreePs <- RelRfSputumAllStartEndPs
RelRfSputumAllStartEndNoTreePs@phy_tree <- NULL
# Merging at Phylum level
RelRfSputumAllStartEndPhylumPs <- microbiomeutilities::aggregate_top_taxa2(RelRfSputumAllStartEndNoTreePs,
                                                                           "Phylum", top = 6)
# Using microbiome::plot_composition for a bar plot
SputumRelAbPhylumPlotComp <- plot_composition(RelRfSputumAllStartEndPhylumPs,
                                              sample.sort = "TTTPrePost5gr",
                                              x.label = "PatientID") + 
  theme(legend.position = "bottom") +
  scale_fill_brewer("Phylum", palette = "Paired") + 
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90)) +
  ggtitle("Relative abundance") +
  theme(legend.title = element_text(size = 18))
# Extracting data from plot_composition object
SputumRelAbPhylumPlotCompData <- SputumRelAbPhylumPlotComp$data
# Adding variables from Secondary data
SputumRelAbPhylumPlotCompData2 <- SputumRelAbPhylumPlotCompData %>% 
  mutate(SampleID = Sample, PatientID = xlabel) %>% 
  left_join(., meta(AbsRfSputumAllStartEndPs) %>% select(SampleID, 
                                                         TTTPrePost5gr, 
                                                         PatientResistStatus, 
                                                         ResistGenesEndStartFoldCh), 
            by = "SampleID")
# Customizing the plot
SputumRelAbPhylumPlotComp2 <- SputumRelAbPhylumPlotCompData2 %>% ggplot(aes(x = SampleID,
                                                                            y = Abundance, 
                                                                            fill = Tax)) + 
  geom_bar(position = "stack", stat = "identity") +
  scale_x_discrete(labels = SputumRelAbPhylumPlotCompData2$PatientID, 
                   breaks = SputumRelAbPhylumPlotCompData2$SampleID) +
  scale_fill_brewer("Family", palette = "Paired") +
  facet_wrap(PatientResistStatus~PatientID, scales = "free") +
  theme_bw()
# Using data from the plot_composition object to obtain a heat map
SputumRelAbPhylumPlotHeat <- SputumRelAbPhylumPlotCompData2 %>% ggplot(aes(x = Sample,
                                                                           y = Tax)) + 
  geom_tile(aes(fill = Abundance)) +
  scale_fill_distiller("Abundance", palette = "RdYlBu") +
  facet_wrap(PatientResistStatus~PatientID, scales = "free") +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.y = element_blank()) 

# Merging at Family level
RelRfSputumAllStartEndFamilyPs <- microbiomeutilities::aggregate_top_taxa2(RelRfSputumAllStartEndNoTreePs,
                                                                           "Family", top = 11)
# Using microbiome::plot_composition for a bar plot
SputumRelAbFamilyPlotComp <- plot_composition(RelRfSputumAllStartEndFamilyPs,
                                              sample.sort = "TTTPrePost5gr",
                                              x.label = "PatientID")
# Extracting data from plot_composition object
SputumRelAbFamilyPlotCompData <- SputumRelAbFamilyPlotComp$data
# Adding variables from Secondary data
SputumRelAbFamilyPlotCompData2 <- SputumRelAbFamilyPlotCompData %>% 
  mutate(SampleID = Sample, PatientID = xlabel) %>% 
  left_join(., meta(AbsRfSputumAllStartEndPs) %>% select(SampleID, 
                                                         TTTPrePost5gr, 
                                                         PatientResistStatus, 
                                                         ResistGenesEndStartFoldCh), 
            by = "SampleID")
# Customizing the plot
SputumRelAbFamilyPlotComp2 <- SputumRelAbFamilyPlotCompData2 %>% ggplot(aes(x = SampleID,
                                                                            y = Abundance, 
                                                                            fill = Tax)) + 
  geom_bar(position = "stack", stat = "identity") +
  scale_x_discrete(labels = SputumRelAbFamilyPlotCompData2$PatientID, 
                   breaks = SputumRelAbFamilyPlotCompData2$SampleID) +
  scale_fill_brewer("Family", palette = "Paired") +
  facet_wrap(PatientResistStatus~PatientID, scales = "free") + 
  theme_bw()
# Using data from the plot_composition object to obtain a heat map
SputumRelAbFamilyPlotHeat <- SputumRelAbFamilyPlotCompData2 %>% ggplot(aes(x = Sample,
                                                                           y = Tax)) + 
  geom_tile(aes(fill = Abundance)) +
  scale_fill_distiller("Abundance", palette = "RdYlBu") +
  facet_wrap(PatientResistStatus~PatientID, scales = "free") +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.y = element_blank())

# Sputum > Top 25 ASVs present ONLY at StartAZT in patients with Stable resistance ----
# Inspecting taxa present in the sputum of patients with Stable resistance at StartAZT+EndAZT
AbsRfSputumStableResistStartEndPs <- subset_samples(AbsRfResist2PatientPs, 
                                                Type == "Sputum" & 
                                                  TTTPrePost5gr %in% c("StartAZT","EndAZT") &
                                                  ResistGenesEndStartFoldCh < 5)
# Keeping only the taxa represented after filtering, then counting them and listing their names
AbsRfSputumStableResistStartEndPs %<>% prune_taxa(taxa_sums(.) > 0, .)
SputumStableResistStartEndNumber <- AbsRfSputumStableResistStartEndPs %>% ntaxa()
SputumStableResistStartEndNames <- AbsRfSputumStableResistStartEndPs %>% taxa_names()
# Inspecting taxa present in the sputum of patients with No resistance at StartAZT
AbsRfSputumStableResistStartAZTPs <- subset_samples(AbsRfSputumStableResistStartEndPs,
                                                TTTPrePost5gr == "StartAZT")
# Keeping only the taxa represented after filtering, then counting them and listing their names
AbsRfSputumStableResistStartAZTPs %<>% prune_taxa(taxa_sums(.) > 0, .)
SputumStableResistStartNumber <- AbsRfSputumStableResistStartAZTPs %>% ntaxa()
SputumStableResistStartNames <- AbsRfSputumStableResistStartAZTPs %>% taxa_names()
# Inspecting taxa present in the sputum of patients with Stable resistance at EndAZT
AbsRfSputumStableResistEndAZTPs <- subset_samples(AbsRfSputumStableResistStartEndPs,
                                              TTTPrePost5gr == "EndAZT")
# Keeping only the taxa represented after filtering, then counting them and listing their names
AbsRfSputumStableResistEndAZTPs %<>% prune_taxa(taxa_sums(.) > 0, .)
SputumStableResistEndNumber <- AbsRfSputumStableResistEndAZTPs %>% ntaxa()
SputumStableResistEndNames <- AbsRfSputumStableResistEndAZTPs %>% taxa_names()

# Obtaining taxa names
SputumStableResistStartOnlyNames <- setdiff(SputumStableResistStartNames, SputumStableResistEndNames)
# Obtaining a PS object with ONLY ASVs present at StartAZT in patients with Stable resistance
AbsRfSputumStableResistStartAZTASVPs <- prune_taxa(SputumStableResistStartOnlyNames, 
                                               AbsRfSputumStableResistStartAZTPs)
# Sorting the 1:25 top abundant ASVs
TopAb25ASVSputumStableResistStartAZT <- sort(taxa_sums(AbsRfSputumStableResistStartAZTASVPs), 
                                         decreasing = TRUE)[1:25]
# Computing the names of the most abundant ASVs
TopAb25ASVSputumStableResistStartAZTName  <-  names(TopAb25ASVSputumStableResistStartAZT)
# Obtaining a PS object containing all ASVs and only sputum samples from patients 
# at the start of treatment, among those with a pair of samples available at the start
# and end of treatment
AbsRfSputumAllASVStableResistStartAZTSamplesPs <- subset_samples(AbsRfPatientPs, 
                                                             OriginType == "Patient_Sputum" &
                                                               SampleID %in% sample_names(AbsRfSputumStableResistStartAZTPs))
# Working with relative abundance
RelRfSputumAllASVStableResistStartAZTSamplesPs <- microbiome::transform(AbsRfSputumAllASVStableResistStartAZTSamplesPs, "compositional")
# Keeping only the 25 most abundant ASVs lost during treatment
RelRfSputumTop25ASVLostStableResistSamplesPs <- RelRfSputumAllASVStableResistStartAZTSamplesPs %>% 
  prune_taxa(TopAb25ASVSputumStableResistStartAZTName, .)
# Obtaining a DF
RelRfSputumTop25ASVLostStableResistSamplesDF <- psmelt(RelRfSputumTop25ASVLostStableResistSamplesPs)

RelRfSputumTop25ASVLostStableResistSamplesDF2 <- RelRfSputumTop25ASVLostStableResistSamplesDF %>%
  mutate(ASV2 = OTU, Abundance = round(Abundance *100, 1)) %>% 
  dplyr::rename(RelAbundPercent = Abundance) %>% 
  unite(ASV_Genus_Species, c(OTU, Genus, Species), sep ="_") %>%
  convert(fct(ASV_Genus_Species, Phylum, ASV2)) %>%
  select(ASV_Genus_Species, SampleID, PatientID, ASV2, Phylum, RelAbundPercent)

#  Dot plot for the 25 most abundant ASVs present at StartAZT ONLY in patients with Stable resistance
LevelsForPlotVect <- RelRfSputumTop25ASVLostStableResistSamplesDF2$ASV_Genus_Species %>% unique
# Reversing the vector due to the horizontal graph
RevLevelsForPlotVect <- rev(LevelsForPlotVect)

SputumTop25ASVLostStableResist_jitterplot <- RelRfSputumTop25ASVLostStableResistSamplesDF2 %>%
  ggplot(aes(x = factor(ASV_Genus_Species, levels = RevLevelsForPlotVect),
             y = RelAbundPercent,
             fill = Phylum)) +
  geom_jitter(shape = 21,
              colour = "black",
              size = 2.5,
              alpha = .7,
              width = .25) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "black") +
  scale_y_continuous(trans = "log10", 
                     limits = c(.01, 100), 
                     labels = format_format(big.mark = " ", 
                                            decimal.mark = ".", 
                                            scientific = FALSE)) +
  scale_fill_manual(values = MyBactCol6BlFr) +
  labs(x = "", y = "Relative abundance") +
  theme_classic() +
  coord_flip() +
  theme(strip.background = element_blank())

ggsave("Top25SputumStableResistStartAZTASV_jitterplot.pdf",
       width = 15, height = 15, dpi = 200, units = "cm")

# Number of ASVs removed during TTT stratified by phylum
Top20GeneraSputumStableResistStartAZT <- tax_table(AbsRfSputumStableResistStartAZTASVPs)[ , "Genus"] %>%
  table %>%
  as_tibble %>% 
  dplyr::rename(Genus = ".", NumberOfASVs = n) %>% 
  convert(fct(Genus)) %>% 
  arrange(desc(NumberOfASVs)) %>% 
  slice_max(NumberOfASVs, n = 20, with_ties = FALSE)
# Adding phylum information
PhylumInfoStableResDF <- tax_table(AbsRfSputumStableResistStartAZTASVPs)[ , c("Phylum", "Genus")] %>% as.data.frame()

Top20GeneraSputumStableResistStartAZT2 <- Top20GeneraSputumStableResistStartAZT %>%
  left_join(., PhylumInfoNoResDF, by = "Genus") %>% 
  distinct(Genus, .keep_all = TRUE) %>% 
  convert(fct(Genus, Phylum))

# Lollipop chart of the number of ASVs lost per phylum (Supplementary Figure E17, panel a)
SputumNumberASVsLostDuringTTTStableResist_Lolli <- Top20GeneraSputumStableResistStartAZT2 %>% 
  ggpubr::ggdotchart(x = "Genus", 
                     y = "NumberOfASVs",
                     color = "Phylum",
                     palette = MyBactCol8BlFr,
                     sorting = "descending",
                     rotate = TRUE,
                     dot.size = 3,
                     add = "segments",
                     xlab = "Top 20 genera with the greatest number\nof ASVs lost during treatment",
                     ylab = "Number of ASVs",
                     ggtheme = theme_pubr())

ggsave("SputumNumberASVsLostDuringTTTStableResist_Lolli.pdf",
       width = 12, height = 9, dpi = 200, units = "cm")

# Sputum > Top 25 ASVs present ONLY at EndAZT in patients with Stable resistance ----
SputumStableResistEndOnlyNames <- setdiff(SputumStableResistEndNames, SputumStableResistStartNames)
# Obtaining a PS object with ONLY ASVs present at EndAZT in patients with Stable resistance
AbsRfSputumStableResistEndAZTASVPs <- prune_taxa(SputumStableResistEndOnlyNames, 
                                             AbsRfSputumStableResistEndAZTPs)
# Sorting the 1:25 top abundant ASVs
TopAb25ASVSputumStableResistEndAZT <- sort(taxa_sums(AbsRfSputumStableResistEndAZTASVPs), 
                                       decreasing = TRUE)[1:25]
# Computing the names of the most abundant ASVs
TopAb25ASVSputumStableResistEndAZTName  <-  names(TopAb25ASVSputumStableResistEndAZT)
# Obtaining a PS object containing all ASVs and only sputum samples from patients at the end
# of TTT, among those with a pair of samples available at the start and end of TTT
AbsRfSputumAllASVStableResistEndAZTSamplesPs <- subset_samples(AbsRfPatientPs, 
                                                           OriginType == "Patient_Sputum" &
                                                             SampleID %in% sample_names(AbsRfSputumStableResistEndAZTPs))
# Working with relative abundance
RelRfSputumAllASVStableResistEndAZTSamplesPs <- microbiome::transform(AbsRfSputumAllASVStableResistEndAZTSamplesPs, "compositional")
# Keeping only the 25 most abundant ASVs acquired during TTT
RelRfSputumTop25ASVAcquiredStableResistSamplesPs <- prune_taxa(TopAb25ASVSputumStableResistEndAZTName, 
                                                           RelRfSputumAllASVStableResistEndAZTSamplesPs)
# Obtaining a DF
RelRfSputumTop25ASVAcquiredStableResistSamplesDF <- psmelt(RelRfSputumTop25ASVAcquiredStableResistSamplesPs)

RelRfSputumTop25ASVAcquiredStableResistSamplesDF2 <- RelRfSputumTop25ASVAcquiredStableResistSamplesDF %>%
  mutate(ASV2 = OTU, Abundance = round(Abundance *100, 1)) %>% 
  dplyr::rename(RelAbundPercent = Abundance) %>% 
  unite(ASV_Genus_Species, c(OTU, Genus, Species), sep ="_") %>%
  convert(fct(ASV_Genus_Species, Phylum, ASV2)) %>%
  select(ASV_Genus_Species, SampleID, PatientID, ASV2, Phylum, RelAbundPercent)

#  Dot plot for the 25 most abundant ASVs present at EndAZT ONLY in patients with Stable resistance
LevelsForPlotVect2 <- RelRfSputumTop25ASVAcquiredStableResistSamplesDF2$ASV_Genus_Species %>% unique
# Reversing the vector due to the horizontal graph
RevLevelsForPlotVect2 <- rev(LevelsForPlotVect2)

SputumTop25ASVAcquiredStableResist_jitterplot <- RelRfSputumTop25ASVAcquiredStableResistSamplesDF2 %>%
  ggplot(aes(x = factor(ASV_Genus_Species, levels = RevLevelsForPlotVect2),
             y = RelAbundPercent,
             fill = Phylum)) +
  geom_jitter(shape = 21,
              colour = "black",
              size = 2.5,
              alpha = .7,
              width = .25) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "black") +
  scale_y_continuous(trans = "log10", 
                     limits = c(.01, 100), 
                     labels = format_format(big.mark = " ", 
                                            decimal.mark = ".", 
                                            scientific = FALSE)) +
  scale_fill_manual(values = MyBactCol6BlFr) +
  labs(x = "", y = "Relative abundance") +
  theme_classic() +
  coord_flip() +
  theme(strip.background = element_blank())

ggsave("Top25SputumNoResistEndAZTASV_jitterplot.pdf",
       width = 15, height = 15, dpi = 200, units = "cm")

# Number of ASVs acquired during TTT stratified by phylum
Top20GeneraSputumStableResistEndAZT <- tax_table(AbsRfSputumStableResistEndAZTASVPs)[ , "Genus"] %>%
  table %>%
  as_tibble %>% 
  dplyr::rename(Genus = ".", NumberOfASVs = n) %>% 
  convert(fct(Genus)) %>% 
  arrange(desc(NumberOfASVs)) %>% 
  slice_max(NumberOfASVs, n = 20, with_ties = FALSE)
# Adding phylum information
PhylumInfoEndStableResDF <- tax_table(AbsRfSputumStableResistEndAZTASVPs)[ , c("Phylum", "Genus")] %>% as.data.frame()

Top20GeneraSputumStableResistEndAZT2 <- Top20GeneraSputumStableResistEndAZT %>%
  left_join(., PhylumInfoEndNoResDF, by = "Genus") %>% 
  distinct(Genus, .keep_all = TRUE) %>% 
  convert(fct(Genus, Phylum))
# Lollipop chart of the number of ASVs lost per phylum (Supplementary Figure E17, panel a)
SputumNumberASVsAcquiredDuringTTTStableResist_Lolli <- Top20GeneraSputumStableResistEndAZT2 %>% 
  ggdotchart(x = "Genus", 
             y = "NumberOfASVs",
             color = "Phylum",
             palette = MyBactCol8BlFr,
             sorting = "descending",
             rotate = TRUE,
             dot.size = 3,
             add = "segments",
             xlab = "Top 20 genera with the greatest number\nof ASVs acquired during treatment",
             ylab = "Number of ASVs",
             ggtheme = theme_pubr())

ggsave("SputumNumberASVsAcquiredDuringTTTNoResist_Lolli.pdf",
       width = 12, height = 9, dpi = 200, units = "cm")

# Sputum > Top 25 ASVs present ONLY at StartAZT in patients with Increasing resistance ----
# Inspecting taxa present in the sputum of patients with Resistance at StartAZT+EndAZT
AbsRfSputumResistStartEndPs <- subset_samples(AbsRfResist2PatientPs, 
                                              Type == "Sputum" & 
                                                TTTPrePost5gr %in% c("StartAZT","EndAZT") &
                                                ResistGenesEndStartFoldCh > 3.83)
# Keeping only the taxa represented after filtering, then counting them and listing their names
AbsRfSputumResistStartEndPs %<>% prune_taxa(taxa_sums(.) > 0, .)
SputumResistStartEndNumber <- AbsRfSputumResistStartEndPs %>% ntaxa()
SputumResistStartEndNames <- AbsRfSputumResistStartEndPs %>% taxa_names()
# Inspecting taxa present in the sputum of patients with Resistance at StartAZT
AbsRfSputumResistStartAZTPs <- subset_samples(AbsRfSputumResistStartEndPs,
                                              TTTPrePost5gr == "StartAZT")
# Keeping only the taxa represented after filtering, then counting them and listing their names
AbsRfSputumResistStartAZTPs %<>% prune_taxa(taxa_sums(.) > 0, .)
SputumResistStartNumber <- AbsRfSputumResistStartAZTPs %>% ntaxa()
SputumResistStartNames <- AbsRfSputumResistStartAZTPs %>% taxa_names()
# Inspecting taxa present in the sputum of patients with Resistance at EndAZT
AbsRfSputumResistEndAZTPs <- subset_samples(AbsRfSputumResistStartEndPs,
                                            TTTPrePost5gr == "EndAZT")
# Keeping only the taxa represented after filtering, then counting them and listing their names
AbsRfSputumResistEndAZTPs %<>% prune_taxa(taxa_sums(.) > 0, .)
SputumResistEndNumber <- AbsRfSputumResistEndAZTPs %>% ntaxa()
SputumResistEndNames <- AbsRfSputumResistEndAZTPs %>% taxa_names()

SputumResistStartOnlyNames <- setdiff(SputumResistStartNames, SputumResistEndNames)
# Obtaining a PS object with ONLY ASVs present at StartAZT in patients with Increasing resistance
AbsRfSputumResistStartAZTASVPs <- prune_taxa(SputumResistStartOnlyNames, 
                                             AbsRfSputumResistStartAZTPs)
# Sorting the 1:25 top abundant ASVs
TopAb25ASVSputumResistStartAZT <- sort(taxa_sums(AbsRfSputumResistStartAZTASVPs), 
                                       decreasing = TRUE)[1:25]
# Computing the names of the most abundant ASVs
TopAb25ASVSputumResistStartAZTName  <-  names(TopAb25ASVSputumResistStartAZT)
# Obtaining a PS object containing all ASVs and only sputum samples from patients 
# at the start of TTT, among those with a pair of samples available at the start and end of TTT
AbsRfSputumAllASVResistStartAZTSamplesPs <- subset_samples(AbsRfPatientPs, 
                                                           OriginType == "Patient_Sputum" &
                                                             SampleID %in% sample_names(AbsRfSputumResistStartAZTPs))
# Working with relative abundance
RelRfSputumAllASVResistStartAZTSamplesPs <- microbiome::transform(AbsRfSputumAllASVResistStartAZTSamplesPs, "compositional")
# Keeping only the 25 most abundant ASVs lost during treatment
RelRfSputumTop25ASVLostResistSamplesPs <- prune_taxa(TopAb25ASVSputumResistStartAZTName, 
                                                     RelRfSputumAllASVResistStartAZTSamplesPs)
# Obtaining a DF
RelRfSputumTop25ASVLostResistSamplesDF <- psmelt(RelRfSputumTop25ASVLostResistSamplesPs)

RelRfSputumTop25ASVLostResistSamplesDF2 <- RelRfSputumTop25ASVLostResistSamplesDF %>%
  mutate(ASV2 = OTU, Abundance = round(Abundance * 100, 1)) %>% 
  dplyr::rename(RelAbundPercent = Abundance) %>% 
  unite(ASV_Genus_Species, c(OTU, Genus, Species), sep ="_") %>%
  convert(fct(ASV_Genus_Species, Phylum, ASV2)) %>%
  select(ASV_Genus_Species, SampleID, PatientID, ASV2, Phylum, RelAbundPercent)

# Dot plot for the 25 most abundant ASVs present at StartAZT ONLY in patients with Increasing resistance
LevelsForPlotVect3 <- RelRfSputumTop25ASVLostResistSamplesDF2$ASV_Genus_Species %>% unique
# Reversing the vector due to the horizontal graph
RevLevelsForPlotVect3 <- rev(LevelsForPlotVect3)

SputumTop25ASVLostResist_jitterplot <- RelRfSputumTop25ASVLostResistSamplesDF2 %>%
  ggplot(aes(x = factor(ASV_Genus_Species, levels = RevLevelsForPlotVect3),
             y = RelAbundPercent,
             fill = Phylum)) +
  geom_jitter(shape = 21,
              colour = "black",
              size = 2.5,
              alpha = .7,
              width = .25) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "black") +
  scale_y_continuous(trans = "log10", 
                     limits = c(.01, 100), 
                     labels = format_format(big.mark = " ", 
                                            decimal.mark = ".", 
                                            scientific = FALSE)) +
  scale_fill_manual(values = MyBactCol6BlFr) +
  labs(x = "", y = "Relative abundance") +
  theme_classic() +
  coord_flip() +
  theme(strip.background = element_blank())

ggsave("Top25SputumResistStartAZTASV_jitterplot.pdf",
       width = 15, height = 15, dpi = 200, units = "cm")

# Number of ASVs removed during TTT stratified by phylum
Top20GeneraSputumResistStartAZT <- tax_table(AbsRfSputumResistStartAZTASVPs)[ , "Genus"] %>%
  table %>%
  as_tibble %>% 
  dplyr::rename(Genus = ".", NumberOfASVs = n) %>% 
  convert(fct(Genus)) %>% 
  arrange(desc(NumberOfASVs)) %>% 
  slice_max(NumberOfASVs, n = 20, with_ties = FALSE)
# Adding phylum information
PhylumInfoResDF <- tax_table(AbsRfSputumResistStartAZTASVPs)[ , c("Phylum", "Genus")] %>% as.data.frame()

Top20GeneraSputumResistStartAZT2 <- Top20GeneraSputumResistStartAZT %>%
  left_join(., PhylumInfoResDF, by = "Genus") %>% 
  distinct(Genus, .keep_all = TRUE) %>% 
  convert(fct(Genus, Phylum))
# Lollipop chart of the number of ASVs lost per phylum (Supplementary Figure E17, panel b)
SputumNumberASVsLostDuringTTTResist_Lolli <- Top20GeneraSputumResistStartAZT2 %>% 
  ggdotchart(x = "Genus", 
             y = "NumberOfASVs",
             color = "Phylum",
             palette = MyBactCol8BlFr,
             sorting = "descending",
             rotate = TRUE,
             dot.size = 3,
             add = "segments",
             xlab = "Top 20 genera with the greatest number\nof ASVs lost during treatment",
             ylab = "Number of ASVs",
             ggtheme = theme_pubr())
ggsave("SputumNumberASVsLostDuringTTTIncrResist_Lolli.pdf",
       width = 12, height = 9, dpi = 200, units = "cm")

# Sputum > Top 25 ASVs present ONLY at EndAZT in patients with Increasing resistance ----
SputumResistEndOnlyNames <- setdiff(SputumResistEndNames, SputumResistStartNames)
# Obtaining a PS object with ONLY ASVs present at EndAZT in patients with Increasing resistance
AbsRfSputumResistEndAZTASVPs <- prune_taxa(SputumResistEndOnlyNames, 
                                           AbsRfSputumResistEndAZTPs)
# Sorting the 1:25 top abundant ASVs
TopAb25ASVSputumResistEndAZT <- sort(taxa_sums(AbsRfSputumResistEndAZTASVPs), 
                                     decreasing = TRUE)[1:25]
# Computing the names of the most abundant ASVs
TopAb25ASVSputumResistEndAZTName  <-  names(TopAb25ASVSputumResistEndAZT)
# Obtaining a PS object containing all ASVs and only sputum samples from patients at the end
# of TTT, among those with a pair of samples available at the start and end of TTT
AbsRfSputumAllASVResistEndAZTSamplesPs <- subset_samples(AbsRfPatientPs, 
                                                         OriginType == "Patient_Sputum" &
                                                           SampleID %in% sample_names(AbsRfSputumResistEndAZTPs))
# Working with relative abundance
RelRfSputumAllASVResistEndAZTSamplesPs <- microbiome::transform(AbsRfSputumAllASVResistEndAZTSamplesPs, "compositional")
# Keeping only the 25 most abundant ASVs acquired during TTT
RelRfSputumTop25ASVAcquiredResistSamplesPs <- prune_taxa(TopAb25ASVSputumResistEndAZTName, 
                                                         RelRfSputumAllASVResistEndAZTSamplesPs)
# Obtaining a DF
RelRfSputumTop25ASVAcquiredResistSamplesDF <- psmelt(RelRfSputumTop25ASVAcquiredResistSamplesPs)

RelRfSputumTop25ASVAcquiredResistSamplesDF2 <- RelRfSputumTop25ASVAcquiredResistSamplesDF %>%
  mutate(ASV2 = OTU, Abundance = round(Abundance * 100, 1)) %>% 
  dplyr::rename(RelAbundPercent = Abundance) %>% 
  unite(ASV_Genus_Species, c(OTU, Genus, Species), sep ="_") %>%
  convert(fct(ASV_Genus_Species, Phylum, ASV2)) %>%
  select(ASV_Genus_Species, SampleID, PatientID, ASV2, Phylum, RelAbundPercent)

#  Dot plot for the 25 most abundant ASVs present at EndAZT ONLY in patients with Increasing resistance
LevelsForPlotVect4 <- RelRfSputumTop25ASVAcquiredResistSamplesDF2$ASV_Genus_Species %>% unique
# Reversing the vector due to the horizontal graph
RevLevelsForPlotVect4 <- rev(LevelsForPlotVect4)

SputumTop25ASVAcquiredResist_jitterplot <- RelRfSputumTop25ASVAcquiredResistSamplesDF2 %>%
  ggplot(aes(x = factor(ASV_Genus_Species, levels = RevLevelsForPlotVect4),
             y = RelAbundPercent,
             fill = Phylum)) +
  geom_jitter(shape = 21,
              colour = "black",
              size = 2.5,
              alpha = .7,
              width = .25) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "black") +
  scale_y_continuous(trans = "log10", 
                     limits = c(.01, 100), 
                     labels = format_format(big.mark = " ", 
                                            decimal.mark = ".", 
                                            scientific = FALSE)) +
  scale_fill_manual(values = MyBactCol6BlFr) +
  labs(x = "", y = "Relative abundance") +
  theme_classic() +
  coord_flip() +
  theme(strip.background = element_blank())

ggsave("Top25SputumResistEndAZTASV_jitterplot.pdf",
       width = 15, height = 15, dpi = 200, units = "cm")

# Number of ASVs acquired during TTT stratified by phylum
Top20GeneraSputumResistEndAZT <- tax_table(AbsRfSputumResistEndAZTASVPs)[ , "Genus"] %>%
  table %>%
  as_tibble %>% 
  dplyr::rename(Genus = ".", NumberOfASVs = n) %>% 
  convert(fct(Genus)) %>% 
  arrange(desc(NumberOfASVs)) %>% 
  slice_max(NumberOfASVs, n = 20, with_ties = FALSE)
# Adding phylum information
PhylumInfoEndDF <- tax_table(AbsRfSputumResistEndAZTASVPs)[ , c("Phylum", "Genus")] %>% as.data.frame()

Top20GeneraSputumResistEndAZT2 <- Top20GeneraSputumResistEndAZT %>%
  left_join(., PhylumInfoEndDF, by = "Genus") %>% 
  distinct(Genus, .keep_all = TRUE) %>% 
  convert(fct(Genus, Phylum))
# Lollipop chart of the number of ASVs lost per phylum (Supplementary Figure E17, panel b)
SputumNumberASVsAcquiredDuringTTTResist_Lolli <- Top20GeneraSputumResistEndAZT2 %>% 
  ggdotchart(x = "Genus", 
             y = "NumberOfASVs",
             color = "Phylum",
             palette = MyBactCol8BlFr,
             sorting = "descending",
             rotate = TRUE,
             dot.size = 3,
             add = "segments",
             xlab = "Top 20 genera with the greatest number\nof ASVs acquired during treatment",
             ylab = "Number of ASVs",
             ggtheme = theme_pubr())
ggsave("SputumNumberASVsAcquiredDuringTTTIncrResist_Lolli.pdf",
       width = 12.5, height = 9, dpi = 200, units = "cm")

# Inspecting the dynamics of different genera using phyloseq::plot_tree (Figure 7) ----
# Analyzing Streptococcus, Leptotrichia, Fusobacteria, Actinomyces,
# Treponema_2, Fretibacterium and Pseudopropionibacterium together
StreptToPseudoprop_ps <- subset_taxa(AbsRfPatientPs,
                                     Genus %in% c("Streptococcus", "Leptotrichia",
                                                  "Fusobacteria", "Actinomyces",
                                                  "Treponema_2", "Fretibacterium",
                                                  "Pseudopropionibacterium"))
# Obtaining a tree agglomerated at the genus level
StreptToPseudopropGenus_ps <- tax_glom(StreptToPseudoprop_ps, "Genus", NArm = TRUE)

# For Figure 7
StreptToPseudopropGenusPtree <- plot_tree(StreptToPseudopropGenus_ps, 
                                          label.tips = "Genus", 
                                          color = "TTTPrePost7gr", 
                                          size = "Abundance", 
                                          sizebase = 5,
                                          base.spacing = 0.1,
                                          plot.margin = 0, 
                                          ladderize = "right") +
  scale_color_manual(values = colorBlindGrey8) +
  theme(legend.position = "top")

ggsave("StreptToPseudopropGenusPtree.pdf",
       width = 40, height = 6, dpi = 200, units = "cm")

# Genus Treponema_2
Trepon_ps <- subset_taxa(AbsRfPatientPs, Genus == "Treponema_2")
# Obtaining a detailed tree with each ASV belonging to the genus Treponema_2
PtreeTrepon <- plot_tree(Trepon_ps, 
                         label.tips = NULL, 
                         color = "TTTPrePost7gr", 
                         size = "Abundance", 
                         sizebase = 5, 
                         plot.margin = 0, 
                         ladderize = "right") +
  scale_color_manual(values = colorBlindGrey8)

ggsave("PtreeTrepon.pdf",
       width = 23, height = 20, dpi = 200, units = "cm")

# Genus Fretibacterium
Synergi_ps <- subset_taxa(AbsRfPatientPs, Phylum == "Synergistetes")
# Obtaining a detailed tree with each ASV belonging to the genus Treponema_2
PtreeSynergi <- plot_tree(Synergi_ps, 
                          label.tips = "Genus", 
                          color = "TTTPrePost7gr", 
                          size = "Abundance", 
                          sizebase = 10, 
                          plot.margin = 0.5, 
                          ladderize = "right")
# Obtaining a tree agglomerated at the genus level
SynergiGenus_ps <- tax_glom(Synergi_ps, "Genus", NArm = TRUE)

PtreeSynergiGenus <- plot_tree(SynergiGenus_ps, 
                               label.tips = "Genus", 
                               color = "TTTPrePost7gr", 
                               size = "Abundance", 
                               sizebase = 10, 
                               plot.margin = 0.5, 
                               ladderize = "right")

# Genus Campylobacter
Epsilonbact_ps <- subset_taxa(AbsRfPatientPs, Phylum == "Epsilonbacteraeota")
# Obtaining a detailed tree with each ASV belonging to the genus Treponema_2
PtreeEpsilonbact <- plot_tree(Epsilonbact_ps, 
                              label.tips = "Genus", 
                              color = "TTTPrePost7gr", 
                              size = "Abundance", 
                              sizebase = 10, 
                              plot.margin = 0.5, 
                              ladderize = "right")

# Genus Mycoplasma
Teneri_ps <- subset_taxa(AbsRfPatientPs, Phylum == "Tenericutes")
# Obtaining a detailed tree with each ASV belonging to the genus Treponema_2
PtreeTeneri <- plot_tree(Teneri_ps, 
                         label.tips = "Genus", 
                         color = "TTTPrePost7gr", 
                         size = "Abundance", 
                         sizebase = 10, 
                         plot.margin = 0.5, 
                         ladderize = "right")

# Genus Streptococcus
Strept_ps <- subset_taxa(AbsRfPatientPs, Genus == "Streptococcus")
# Obtaining a detailed tree with each ASV belonging to the genus Treponema_2
PtreeStrept <- plot_tree(Strept_ps, 
                         label.tips = "Species", 
                         color = "TTTPrePost7gr",
                         size = "Abundance",
                         sizebase = 10,
                         ladderize = "right") + 
  coord_polar(theta="y")

ggsave("PtreeStrept.pdf",
       width = 50, height = 50, dpi = 200, units = "cm")

# Genus Pseudopropionibacterium
Pseudoprop_ps <- subset_taxa(AbsRfPatientPs, Genus == "Pseudopropionibacterium")
# Obtaining a detailed tree with each ASV belonging to the genus Treponema_2
PseudopropPtree <- plot_tree(Pseudoprop_ps, 
                             label.tips = "Genus", 
                             color = "TTTPrePost7gr", 
                             size = "Abundance", 
                             sizebase = 10, 
                             plot.margin = 0.5, 
                             ladderize = "right")

# Calculating the factor of variation between the start and end of TTT for each ASV in each sample ----
# Starting with a PS object limited to patients with sputum samples available at StartAZT and EndAZT
AbsRfSputumAllStartEndDF <- psmelt(AbsRfSputumAllStartEndPs)
# Calculating the abundance ratio (integer division %/%) between StartAZT and EndAZT for each
# ASV in each patient
AbsRfSputumStartVsEndRatioDF <- AbsRfSputumAllStartEndDF %>%  
  rename(ASV = OTU) %>% 
  unite(ASV_Genus_Species, c(ASV, Genus, Species), sep ="_") %>% 
  convert(fct(ASV_Genus_Species, Phylum)) %>% 
  select(ASV_Genus_Species,Phylum, PatientID, TTTPrePost5gr, Abundance, PatientResistStatus) %>%
  group_by(PatientID, ASV_Genus_Species) %>% 
  mutate(StartVsEndRatio = Abundance[TTTPrePost5gr == "StartAZT"] %/%
           Abundance[TTTPrePost5gr == "EndAZT"]) %>%
  relocate(StartVsEndRatio, .after = Abundance) %>% 
  filter(StartVsEndRatio > 1 & StartVsEndRatio != Inf) %>% 
  arrange(desc(StartVsEndRatio), ASV_Genus_Species, PatientID, TTTPrePost5gr) %>% 
  filter(TTTPrePost5gr == "StartAZT")

# Dot plot for the 25 ASVs reaching the highest ratio between StartAZT and EndAZT 
# in sputum samples from patients with Stable resistance
AbsRfSputumStartVsEndRatioStableResistDF <- AbsRfSputumStartVsEndRatioDF %>% 
  filter(PatientResistStatus == "StableResist")
# ASV names
AbsRfSputumStartVsEndRatioStableResistNamesVect <- AbsRfSputumStartVsEndRatioStableResistDF$ASV_Genus_Species %>% unique()
# Levels for plot
LevelsForPlotVect5 <- AbsRfSputumStartVsEndRatioStableResistNamesVect[1:25]
# Reversing the vector due to the horizontal graph
RevLevelsForPlotVect5 <- rev(LevelsForPlotVect5)

SputumTop25DecreasingASVsStableResist_jitterplot <- AbsRfSputumStartVsEndRatioStableResistDF %>%
  ggplot(aes(x = factor(ASV_Genus_Species, levels = RevLevelsForPlotVect5),
             y = StartVsEndRatio,
             fill = Phylum,
             shape = PatientResistStatus)) +
  geom_jitter(shape = 21,
              colour = "black",
              size = 2.5,
              alpha = .7,
              width = .15) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "black") +
  scale_y_continuous(trans = "log10", limits = c(1, 1000)) +
  scale_fill_manual(values = MyBactCol6BlFr) +
  labs(x = "", y = "StartAZM to EndAZM abundance ratio") +
  theme_classic() +
  coord_flip() +
  theme(strip.background = element_blank())

# Saving as .pdf in the folder containing the R project
ggsave("SputumTop25DecreasingASVsStableResist_jitterplot.pdf",
       width = 15, height = 15, dpi = 200, units = "cm")

# Dot plot for the 25 ASVs reaching the highest ratio between EndAZT and StartAZT 
# in sputum samples from patients with Stable resistance
AbsRfSputumEndVsStartRatioDF <- AbsRfSputumAllStartEndDF %>%  
  rename(ASV = OTU) %>% 
  unite(ASV_Genus_Species, c(ASV, Genus, Species), sep ="_") %>% 
  convert(fct(ASV_Genus_Species, Phylum)) %>% 
  select(ASV_Genus_Species, Phylum, PatientID, TTTPrePost5gr, Abundance, PatientResistStatus) %>%
  group_by(PatientID, ASV_Genus_Species) %>% 
  mutate(EndVsStartRatio = Abundance[TTTPrePost5gr == "EndAZT"] %/%
           Abundance[TTTPrePost5gr == "StartAZT"]) %>%
  relocate(EndVsStartRatio, .after = Abundance) %>% 
  filter(EndVsStartRatio > 1 & EndVsStartRatio != Inf) %>% 
  arrange(desc(EndVsStartRatio), ASV_Genus_Species, PatientID, TTTPrePost5gr) %>% 
  filter(TTTPrePost5gr == "EndAZT")

AbsRfSputumEndVsStartRatioStableResistDF <- AbsRfSputumEndVsStartRatioDF %>% 
  filter(PatientResistStatus == "StableResist")
# ASV names
AbsRfSputumEndVsStartRatioStableResistNamesVect <- AbsRfSputumEndVsStartRatioStableResistDF$ASV_Genus_Species %>% 
  unique()
# Levels for plot
LevelsForPlotVect6 <- AbsRfSputumEndVsStartRatioStableResistNamesVect[1:25]
# Reversing the vector due to the horizontal graph
RevLevelsForPlotVect6 <- rev(LevelsForPlotVect6)

SputumTop25IncreasingASVsStableResist_jitterplot <- AbsRfSputumEndVsStartRatioStableResistDF %>%
  ggplot(aes(x = factor(ASV_Genus_Species, levels = RevLevelsForPlotVect6),
             y = EndVsStartRatio,
             fill = Phylum,
             shape = PatientResistStatus)) +
  geom_jitter(shape = 21,
              colour = "black",
              size = 2.5,
              alpha = .7,
              width = .05) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "black") +
  scale_y_continuous(trans = "log10", limits = c(1, 1000)) +
  scale_fill_manual(values = MyBactCol6BlFr) +
  labs(x = "", y = "EndAZM to StartAZM abundance ratio") +
  theme_classic() +
  coord_flip() +
  theme(strip.background = element_blank())

# Saving as .pdf in the folder containing the R project
ggsave("SputumTop25IncreasingASVsStableResist_jitterplot.pdf",
       width = 15, height = 15, dpi = 200, units = "cm")

# Dot plot for the 25 ASVs reaching the highest ratio between StartAZT and EndAZT 
# in sputum samples from patients with Increasing resistance
AbsRfSputumStartVsEndRatioIncrResistDF <- AbsRfSputumStartVsEndRatioDF %>% 
  filter(PatientResistStatus == "IncreasedResist")
# ASV names
AbsRfSputumStartVsEndRatioIncrResistNamesVect <- AbsRfSputumStartVsEndRatioIncrResistDF$ASV_Genus_Species %>% 
  unique()
# Levels for plot
LevelsForPlotVect7 <- AbsRfSputumStartVsEndRatioIncrResistNamesVect[1:25]
# Reversing the vector due to the horizontal graph
RevLevelsForPlotVect7 <- rev(LevelsForPlotVect7)

SputumTop25DecreasingASVsIncrResist_jitterplot <- AbsRfSputumStartVsEndRatioIncrResistDF %>%
  ggplot(aes(x = factor(ASV_Genus_Species, levels = RevLevelsForPlotVect7),
             y = StartVsEndRatio,
             fill = Phylum,
             shape = PatientResistStatus)) +
  geom_jitter(shape = 21,
              colour = "black",
              size = 2.5,
              alpha = .7,
              width = .1) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "black") +
  scale_y_continuous(trans = "log10", limits = c(1, 1000)) +
  scale_fill_manual(values = MyBactCol6BlFr) +
  labs(x = "", y = "StartAZM to EndAZM abundance ratio") +
  theme_classic() +
  coord_flip() +
  theme(strip.background = element_blank())

# Saving as .pdf in the folder containing the R project
ggsave("SputumTop25DecreasingASVsIncrResist_jitterplot.pdf",
       width = 15, height = 15, dpi = 200, units = "cm")

# Dot plot for the 25 ASVs reaching the highest ratio between EndAZT and StartAZT 
# in sputum samples from patients with Increasing resistance
AbsRfSputumEndVsStartRatioDF <- AbsRfSputumAllStartEndDF %>%  
  rename(ASV = OTU) %>% 
  unite(ASV_Genus_Species, c(ASV, Genus, Species), sep ="_") %>% 
  convert(fct(ASV_Genus_Species, Phylum)) %>% 
  select(ASV_Genus_Species, Phylum, PatientID, TTTPrePost5gr, Abundance, PatientResistStatus) %>%
  group_by(PatientID, ASV_Genus_Species) %>% 
  mutate(EndVsStartRatio = Abundance[TTTPrePost5gr == "EndAZT"] %/%
           Abundance[TTTPrePost5gr == "StartAZT"]) %>%
  relocate(EndVsStartRatio, .after = Abundance) %>% 
  filter(EndVsStartRatio > 1 & EndVsStartRatio != Inf) %>% 
  arrange(desc(EndVsStartRatio), ASV_Genus_Species, PatientID, TTTPrePost5gr) %>% 
  filter(TTTPrePost5gr == "EndAZT")

AbsRfSputumEndVsStartRatioIncrResistDF <- AbsRfSputumEndVsStartRatioDF %>% 
  filter(PatientResistStatus == "IncreasedResist")
# ASV names
AbsRfSputumEndVsStartRatioIncrResistNamesVect <- AbsRfSputumEndVsStartRatioIncrResistDF$ASV_Genus_Species %>% 
  unique()
# Levels for plot
LevelsForPlotVect8 <- AbsRfSputumEndVsStartRatioIncrResistNamesVect[1:25]
# Reversing the vector due to the horizontal graph
RevLevelsForPlotVect8 <- rev(LevelsForPlotVect8)

SputumTop25IncreasingASVsIncrResist_jitterplot <- AbsRfSputumEndVsStartRatioIncrResistDF %>%
  ggplot(aes(x = factor(ASV_Genus_Species, levels = RevLevelsForPlotVect8),
             y = EndVsStartRatio,
             fill = Phylum,
             shape = PatientResistStatus)) +
  geom_jitter(shape = 21,
              colour = "black",
              size = 2.5,
              alpha = .7,
              width = .1) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "black") +
  scale_fill_manual(values = MyBactCol6BlFr) +
  labs(x = "", y = "EndAZM to StartAZM abundance ratio") +
  scale_y_continuous(trans = "log10", limits = c(1, 1000)) +
  theme_classic() +
  coord_flip() +
  theme(strip.background = element_blank())

# Saving as .pdf in the folder containing the R project
ggsave("SputumTop25IncreasingASVsIncrResist_jitterplot.pdf",
       width = 15, height = 15, dpi = 200, units = "cm")





# Assessing in LRT the dynamics of bacteria cleared from LRT during treatment (Figure 4a, left panel) ----
# Sputum: Names of all ASVs present at StartAZT but not EndAZT
SputumStartOnlyNames <- setdiff(SputumStartNames, SputumEndNames)

# Names of patients with a pair of sputum samples available at StartAZT and EndAZT
# and selection of all samples obtained from these patients
SputumStartEndPatientNames <- meta(AbsRfSputumStartAZTPs)$PatientID

AbsRfSputumPatSputumStartEndAllASVPs <- subset_samples(AbsRfPatientPs,
                                                       OriginType == "Patient_Sputum" &
                                                         PatientID %in% SputumStartEndPatientNames)
# Removing ASVs with zero count
AbsRfSputumPatSputumStartEndAllASVPs2 <- prune_taxa(taxa_sums(AbsRfSputumPatSputumStartEndAllASVPs) > 0, 
                                                    AbsRfSputumPatSputumStartEndAllASVPs)
# Merging taxa cleared during treatment
AbsRfSputumAllASVMergeClearedPs2 <- microbiome::merge_taxa2(AbsRfSputumPatSputumStartEndAllASVPs2, 
                                                            taxa = SputumStartOnlyNames,
                                                            name = "Cleared during treatment")
# Working with relative abundance
RelRfSputumAllASVMergeClearedPs2 <- microbiome::transform(AbsRfSputumAllASVMergeClearedPs2, 
                                                          "compositional")
# Data frame for subsequent plotting
RelRfSputumAllASVMergeClearedDF <- psmelt(RelRfSputumAllASVMergeClearedPs2)
# Keeping only taxa removed during treatment
RelRfSputumAllASVMergeClearedDF2 <- RelRfSputumAllASVMergeClearedDF %>%
  dplyr::rename(TaxaGroup = OTU, RelAbundance = Abundance) %>% 
  filter(TaxaGroup == "Cleared during treatment") %>% 
  select(SampleID, PatientID, TTTPrePost5gr, TTTPrePost7gr, RelAbundance, TaxaGroup) %>% 
  convert(fct(SampleID, PatientID, TTTPrePost5gr, TTTPrePost7gr, TaxaGroup))

# Summary statistics
RelRfSputumAllASVMergeClearedDF2 %>% 
  dplyr::group_by(TTTPrePost7gr) %>%
  dplyr::summarize(MedianLoad = median(RelAbundance),
                   IQRLoad = IQR(RelAbundance))

RelRfSputumAllASVMergeClearedDF2 %>%
  filter(TTTPrePost7gr == "StartAZT") %>%
  pull(RelAbundance) %>% 
  quantile() %>% 
  percent()

RelRfSputumAllASVMergeClearedDF2 %>%
  filter(TTTPrePost7gr %in% c("PreAZT_1mo", "PreAZT_4mo")) %>%
  pull(RelAbundance) %>% 
  quantile() %>% 
  percent()

RelRfSputumAllASVMergeClearedDF2 %>%
  filter(TTTPrePost7gr == "PreAZT_4mo") %>%
  pull(RelAbundance) %>% 
  quantile() %>% 
  percent()

# For Figure 4a, left panel
SputumMergeClearedASV_boxplot <- RelRfSputumAllASVMergeClearedDF2 %>%
  ggplot(aes(x = TTTPrePost7gr,
             y = RelAbundance,
             fill = TTTPrePost7gr)) +
  geom_boxplot(width = .5,
               size = .3,
               alpha = .4,
               position = position_dodge(.7),
               outlier.shape = NA,
               show.legend = FALSE) +
  geom_jitter(shape = 21,
              colour = "black",
              size = 1.5,
              alpha = .6,
              width = .2) +
  stat_compare_means(comparisons = CompTTT7) +
  stat_compare_means(label.y = .8) +
  ylim(0, .8) +
  scale_fill_manual(values = TTTphaseCol7) +
  labs(x = "", y = "Relative abundance") +
  theme_classic() +
  theme(strip.background = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1))

ggsave("SputumMergeClearedASV_boxplot.pdf",
       width = 10, height = 12, dpi = 200, units = "cm")

# Assessing in LRT the dynamics of bacteria acquired in LRT during treatment (Figure 4b, left panel) ----
# Sputum: Names of all ASVs present at EndAZT but not StartAZT
SputumEndOnlyNames <- setdiff(SputumEndNames, SputumStartNames)

# Names of patients with a pair of sputum samples available at StartAZT and EndAZT
# and selection of all samples obtained from these patients
SputumStartEndPatientNames <- meta(AbsRfSputumStartAZTPs)$PatientID

AbsRfSputumPatSputumStartEndAllASVPs <- subset_samples(AbsRfPatientPs,
                                                       OriginType == "Patient_Sputum" &
                                                         PatientID %in% SputumStartEndPatientNames)
# Removing ASVs with zero count
AbsRfSputumPatSputumStartEndAllASVPs2 <- prune_taxa(taxa_sums(AbsRfSputumPatSputumStartEndAllASVPs) > 0, 
                                                    AbsRfSputumPatSputumStartEndAllASVPs)
# Merging taxa cleared during treatment
AbsRfSputumAllASVMergeAcquiredPs2 <- microbiome::merge_taxa2(AbsRfSputumPatSputumStartEndAllASVPs2, 
                                                             taxa = SputumEndOnlyNames,
                                                             name = "Acquired during treatment")
# Working with relative abundance
RelRfSputumAllASVMergeAcquiredPs2 <- microbiome::transform(AbsRfSputumAllASVMergeAcquiredPs2, 
                                                           "compositional")
# Data frame for subsequent plotting
RelRfSputumAllASVMergeAcquiredDF <- psmelt(RelRfSputumAllASVMergeAcquiredPs2)
# Keeping only taxa removed during treatment
RelRfSputumAllASVMergeAcquiredDF2 <- RelRfSputumAllASVMergeAcquiredDF %>%
  dplyr::rename(TaxaGroup = OTU, RelAbundance = Abundance) %>% 
  filter(TaxaGroup == "Acquired during treatment") %>% 
  select(SampleID, PatientID, TTTPrePost5gr, TTTPrePost7gr, RelAbundance, TaxaGroup) %>% 
  convert(fct(SampleID, PatientID, TTTPrePost5gr, TTTPrePost7gr, TaxaGroup))

# For Figure 4b, left panel
SputumMergeAcquiredASV_boxplot <- RelRfSputumAllASVMergeAcquiredDF2 %>%
  ggplot(aes(x = TTTPrePost7gr,
             y = RelAbundance,
             fill = TTTPrePost7gr)) +
  geom_boxplot(width = .5,
               size = .3,
               alpha = .4,
               position = position_dodge(.7),
               outlier.shape = NA,
               show.legend = FALSE) +
  geom_jitter(shape = 21,
              colour = "black",
              size = 1.5,
              alpha = .6,
              width = .2) +
  stat_compare_means(comparisons = CompTTT7) +
  stat_compare_means(label.y = .8) +
  ylim(0, .8) +
  scale_fill_manual(values = TTTphaseCol7) +
  labs(x = "", y = "Relative abundance") +
  theme_classic() +
  theme(strip.background = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1))

ggsave("SputumMergeAcquiredASV_boxplot.pdf",
       width = 10, height = 12, dpi = 200, units = "cm")

# Assessing in URT the dynamics of bacteria cleared from LRT during treatment (Figure 4a, right panel) ----
# Only patients who have lost these bacteria in LTR are considered
# Names of patients with a pair of sputum samples available at StartAZT and EndAZT
# and selection of all swab samples obtained from these patients
SputumStartEndPatientNames <- meta(AbsRfSputumStartAZTPs)$PatientID

AbsRfSwabPatSputumStartEndAllASVPs <- subset_samples(AbsRfPatientPs,
                                                     OriginType == "Patient_Swab" &
                                                       PatientID %in% SputumStartEndPatientNames)

# Removing ASVs with zero count
AbsRfSwabPatSputumStartEndAllASVPs2 <- prune_taxa(taxa_sums(AbsRfSwabPatSputumStartEndAllASVPs) > 0, 
                                                  AbsRfSwabPatSputumStartEndAllASVPs)

# Names of all ASVs present at StartAZT but not EndAZT in sputum
SputumStartOnlyNames <- setdiff(SputumStartNames, SputumEndNames)

ClearedSputumPresentSwabNames<- intersect(SputumStartOnlyNames,
                                          taxa_names(AbsRfSwabPatSputumStartEndAllASVPs2))

# Merging taxa cleared during treatment
AbsRfSwabAllASVMergeClearedSputumPs2 <- microbiome::merge_taxa2(AbsRfSwabPatSputumStartEndAllASVPs2, 
                                                                taxa = ClearedSputumPresentSwabNames,
                                                                name = "Cleared from sputum during treatment")
# Working with relative abundance
RelRfSwabAllASVMergeClearedSputumPs2 <- microbiome::transform(AbsRfSwabAllASVMergeClearedSputumPs2, 
                                                              "compositional")
# Data frame for subsequent plotting
RelRfSwabAllASVMergeClearedSputumDF <- psmelt(RelRfSwabAllASVMergeClearedSputumPs2)
# Keeping only taxa removed during treatment
RelRfSwabAllASVMergeClearedSputumDF2 <- RelRfSwabAllASVMergeClearedSputumDF %>%
  dplyr::rename(TaxaGroup = OTU, RelAbundance = Abundance) %>% 
  filter(TaxaGroup == "Cleared from sputum during treatment") %>% 
  select(SampleID, PatientID, TTTPrePost5gr, TTTPrePost7gr, RelAbundance, TaxaGroup) %>% 
  convert(fct(SampleID, PatientID, TTTPrePost5gr, TTTPrePost7gr, TaxaGroup))

# Summary statistics
RelRfSwabAllASVMergeClearedSputumDF2 %>% 
  dplyr::group_by(TTTPrePost7gr) %>%
  dplyr::summarize(MedianLoad = median(RelAbundance),
                   IQRLoad = IQR(RelAbundance))

RelRfSwabAllASVMergeClearedSputumDF2 %>%
  filter(TTTPrePost7gr == "StartAZT") %>%
  pull(RelAbundance) %>% 
  quantile() %>% 
  percent()

RelRfSwabAllASVMergeClearedSputumDF2 %>%
  filter(TTTPrePost7gr %in% c("PreAZT_1mo", "PreAZT_4mo")) %>%
  pull(RelAbundance) %>% 
  quantile() %>% 
  percent()

RelRfSwabAllASVMergeClearedSputumDF2 %>%
  filter(TTTPrePost7gr == "EndAZT") %>%
  pull(RelAbundance) %>% 
  quantile() %>% 
  percent()

# For Figure 4a, right panel
SwabMergeClearedSputumASV_boxplot <- RelRfSwabAllASVMergeClearedSputumDF2 %>%
  ggplot(aes(x = TTTPrePost7gr,
             y = RelAbundance,
             fill = TTTPrePost7gr)) +
  geom_boxplot(width = .5,
               size = .3,
               alpha = .4,
               position = position_dodge(.7),
               outlier.shape = NA,
               show.legend = FALSE) +
  geom_jitter(shape = 21,
              colour = "black",
              size = 1.5,
              alpha = .6,
              width = .2) +
  stat_compare_means(comparisons = CompTTT7) +
  stat_compare_means(label.y = .8) +
  ylim(0, .8) +
  scale_fill_manual(values = TTTphaseCol7) +
  labs(x = "", y = "Relative abundance") +
  theme_classic() +
  theme(strip.background = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1))

ggsave("SwabMergeClearedSputumASV_boxplot.pdf",
       width = 10, height = 12, dpi = 200, units = "cm")

# Assessing in URT the dynamics of bacteria acquired in LRT during treatment (Figure 4b, right panel) ----
# Only patients who have acquired these bacteria in LTR are considered
# Names of patients with a pair of sputum samples available at StartAZT and EndAZT
# and selection of all swab samples obtained from these patients
SputumStartEndPatientNames <- meta(AbsRfSputumStartAZTPs)$PatientID

AbsRfSwabPatSputumStartEndAllASVPs <- subset_samples(AbsRfPatientPs,
                                                     OriginType == "Patient_Swab" &
                                                       PatientID %in% SputumStartEndPatientNames)

# Removing ASVs with zero count
AbsRfSwabPatSputumStartEndAllASVPs2 <- prune_taxa(taxa_sums(AbsRfSwabPatSputumStartEndAllASVPs) > 0, 
                                                  AbsRfSwabPatSputumStartEndAllASVPs)

# Names of all ASVs present at EndAZT but not StartAZT in sputum
SputumEndOnlyNames <- setdiff(SputumEndNames, SputumStartNames)

AcquiredSputumPresentSwabNames<- intersect(SputumEndOnlyNames,
                                           taxa_names(AbsRfSwabPatSputumStartEndAllASVPs2))

# Merging taxa acquired in sputum during treatment
AbsRfSwabAllASVMergeAcquiredSputumPs2 <- microbiome::merge_taxa2(AbsRfSwabPatSputumStartEndAllASVPs2, 
                                                                 taxa = AcquiredSputumPresentSwabNames,
                                                                 name = "Acquired in sputum during treatment")
# Working with relative abundance
RelRfSwabAllASVMergeAcquiredSputumPs2 <- microbiome::transform(AbsRfSwabAllASVMergeAcquiredSputumPs2, 
                                                               "compositional")
# Obtaining a DF for subsequent plotting
RelRfSwabAllASVMergeAcquiredSputumDF <- psmelt(RelRfSwabAllASVMergeAcquiredSputumPs2)
# Keeping only taxa removed during treatment
RelRfSwabAllASVMergeAcquiredSputumDF2 <- RelRfSwabAllASVMergeAcquiredSputumDF %>%
  dplyr::rename(TaxaGroup = OTU, RelAbundance = Abundance) %>% 
  filter(TaxaGroup == "Acquired in sputum during treatment") %>% 
  select(SampleID, PatientID, TTTPrePost5gr, TTTPrePost7gr, RelAbundance, TaxaGroup) %>% 
  convert(fct(SampleID, PatientID, TTTPrePost5gr, TTTPrePost7gr, TaxaGroup))

# Figure 4b, right panel
SwabMergeAcquiredSputumASV_boxplot <- RelRfSwabAllASVMergeAcquiredSputumDF2 %>%
  ggplot(aes(x = TTTPrePost7gr,
             y = RelAbundance,
             fill = TTTPrePost7gr)) +
  geom_boxplot(width = .5,
               size = .3,
               alpha = .4,
               position = position_dodge(.7),
               outlier.shape = NA,
               show.legend = FALSE) +
  geom_jitter(shape = 21,
              colour = "black",
              size = 1.5,
              alpha = .6,
              width = .2) +
  stat_compare_means(comparisons = CompTTT7) +
  stat_compare_means(label.y = .8) +
  ylim(0, .8) +
  scale_fill_manual(values = TTTphaseCol7) +
  labs(x = "", y = "Relative abundance") +
  theme_classic() +
  theme(strip.background = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1))

ggsave("SwabMergeAcquiredSputumASV_boxplot.pdf",
       width = 10, height = 12, dpi = 200, units = "cm")

# Summary statistics
RelRfSwabAllASVMergeAcquiredSputumDF2 %>% 
  dplyr::group_by(TTTPrePost7gr) %>%
  dplyr::summarize(MedianLoad = median(RelAbundance),
                   IQRLoad = IQR(RelAbundance))

RelRfSwabAllASVMergeAcquiredSputumDF2 %>%
  filter(TTTPrePost7gr %in% c("PreAZT_1mo", "PreAZT_4mo", "StartAZT")) %>%
  pull(RelAbundance) %>% 
  quantile() %>% 
  percent()

RelRfSwabAllASVMergeAcquiredSputumDF2 %>%
  filter(TTTPrePost7gr == "EndAZT") %>%
  pull(RelAbundance) %>% 
  quantile() %>% 
  percent()

# Sputum: Asscessing the link between the status of ARG carriage and bacterial load (Figure 6a; Supplementary Figure E15) ----
# Obtaining a phyloseq object with only those patients who acquired Abx resistance 
# between StartATZ and EndAZT  
# Sputum + Swab
AbsRfResistPatientPs <- subset_samples(AbsRfPatientPs, PatientID %in% SputumPatientsEndStartVect)
# Keeping only the taxa represented after filtering
AbsRfResistPatientPs %<>% prune_taxa(taxa_sums(.) > 0, .)
# Including into secondary data the new variables linked to resistance genes
AbsRfResist2PatientDF <- meta(AbsRfResistPatientPs) %>% left_join(.,
                                                                  SputumResistStatusScale01DF %>%
                                                                    select(SampleID, 
                                                                           PatientResistStatus, 
                                                                           ResistGenesEndStartFoldCh),
                                                                  by = "SampleID")
# Checking sample names to update the phyloseq object
# Sample names must be the same in the Secondary Data file (as rownames) and 
# in the ASV table (as column names)
sample_names(otu_table(AbsRfResistPatientPs))
rownames(AbsRfResist2PatientDF)
rownames(AbsRfResist2PatientDF) <- AbsRfResist2PatientDF$SampleID
intersect(rownames(AbsRfResist2PatientDF), sample_names(otu_table(AbsRfResistPatientPs)))


# Using phyloseq::sample_data to make the data frame compatible for
# the subsequent construction of a synthetic phyloseq object
AbsRfResist2Patient3SD <- sample_data(AbsRfResist2PatientDF)
# Updating the phyloseq object with resistance information
AbsRfResist2PatientPs <- merge_phyloseq(AbsRfResistPatientPs, AbsRfResist2Patient3SD)

# 1) Phyloseq object with patients with increased resistance during treatment
AbsRfSputumResistStartEndPs <- subset_samples(AbsRfResist2PatientPs, 
                                              Type == "Sputum" & 
                                                TTTPrePost5gr %in% c("StartAZT","EndAZT") &
                                                ResistGenesEndStartFoldCh > MedFoldChange)



ResistPatientNames <- meta(AbsRfSputumResistStartEndPs)$PatientID

AbsRfSputumAllASVResistPs <- subset_samples(AbsRfPatientPs,
                                            OriginType == "Patient_Sputum" &
                                              PatientID %in% ResistPatientNames)

# 2) Phyloseq object with patients with stable resistance during treatment
AbsRfSputumNoResistStartEndPs <- subset_samples(AbsRfResist2PatientPs, 
                                                Type == "Sputum" & 
                                                  TTTPrePost5gr %in% c("StartAZT","EndAZT") &
                                                  ResistGenesEndStartFoldCh < MedFoldChange)

NoResistPatientNames <- meta(AbsRfSputumNoResistStartEndPs)$PatientID

AbsRfSputumAllASVNoResistPs <- subset_samples(AbsRfPatientPs,
                                            OriginType == "Patient_Sputum" &
                                              PatientID %in% NoResistPatientNames)

PatientResistStatusLoadDF <- meta(AbsRfPatientPs) %>%
  left_join(.,
            AbsRfResist2PatientDF %>% 
              select(SampleID, PatientResistStatus), 
            by = "SampleID") %>% 
  mutate(ResistStatusComplete = case_when(PatientID %in% NoResistPatientNames ~ "StableResist",
                                          PatientID %in% ResistPatientNames ~ "IncreasedResist")) %>% 
  select(SampleID, 
         PatientID,
         OriginType, 
         TTTPrePost5gr, 
         TTTPrePost7gr, 
         Load, 
         ResistStatusComplete) %>% 
  filter(ResistStatusComplete %in% c("StableResist", "IncreasedResist"),
         OriginType == "Patient_Sputum")

# Boxplot per treatment phase, distinguishing between Stable and Increased resistance groups (Supplementary Figure E15)
SputumPatientResistStatusLoadDF_boxplot <- PatientResistStatusLoadDF %>%
  mutate(ResistStatusComplete = factor(ResistStatusComplete, 
                                       levels = c("StableResist", "IncreasedResist"))) %>% 
  ggplot(aes(x = TTTPrePost7gr,
             y = Load,
             fill = TTTPrePost7gr)) +
  geom_boxplot(width = .5,
               size = .3,
               alpha = .4,
               position = position_dodge(.7),
               outlier.shape = NA,
               show.legend = FALSE) +
  geom_point(shape = 21,
              colour = "black",
              size = 1.5,
              alpha = .6) +
  scale_fill_manual(values = TTTphaseCol7) +
  stat_compare_means(comparisons = CompTTT3) + # Add p-value
  stat_compare_means(label.y = 3) + # Add global p-value
  facet_wrap(~ ResistStatusComplete, ncol = 2) +
  scale_y_continuous(breaks = c(10^3, 10^4, 10^5, 10^6, 10^7, 10^8), 
                     limits = c(10^2, 10^8), 
                     trans = "log10") +
  labs(x = "", y = "16S rRNA gene copies\n per microlitre of sputum") +
  theme_classic() +
  theme(strip.background = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1))

ggsave("SputumPatientResistStatusLoad_boxplot.pdf",
       width = 13, height = 12, dpi = 200, units = "cm")

# Comparing the overall bacterial load in patients with Stable vs. Increased resistance (Figure 6a)
PatientResistStatusAllTimesLoad_boxplot <- PatientResistStatusLoadDF %>%
  mutate(ResistStatusComplete = factor(ResistStatusComplete, 
                                       levels = c("StableResist", "IncreasedResist"))) %>% 
  ggplot(aes(x = ResistStatusComplete,
             y = Load,
             fill = ResistStatusComplete)) +
  geom_boxplot(width = .5,
               size = .5,
               alpha = .3,
               position = position_dodge(.7),
               outlier.shape = NA,
               show.legend = FALSE) +
  geom_jitter(shape = 21,
              colour = "black",
              size = 1.5,
              alpha = .6,
              width = .2) +
  stat_compare_means(paired = FALSE) +
  scale_y_continuous(breaks = c(10^3, 10^4, 10^5, 10^6, 10^7), 
                     limits = c(10^3, 10^8), 
                     trans = "log10") +
  scale_fill_manual(values = ResistStatusCol2) +
  labs(x = "", y = "16S rRNA gene copies\n per microlitre of sputum") +
  theme_classic() +
  theme(strip.background = element_blank(),
        axis.text.x = element_blank())

ggsave("SputumPatientResistStatusAllTimesLoad_boxplot.pdf",
       width = 10, height = 10, dpi = 200, units = "cm")

# Sputum: Prevotella_7 and Streptococcus abundance (Supplementary Figure E16, panels a to d) ----
# Comparing the medians (Wilcoxon rank sum test) of the relative abundance of Prevotella_7 or Streptococcus 
# in sputum samples collected at the end and after AZT treatment, in patients with stable or increased AR gene carriage
AbsRfSputumAllStartEndPs <- subset_samples(AbsRfResist2PatientPs, 
                                           Type == "Sputum" & 
                                             TTTPrePost5gr %in% c("StartAZT","EndAZT"))

# Patients with stable AR gene carriage
StableResistPat <- meta(AbsRfSputumAllStartEndPs) %>%
  filter(PatientResistStatus == "StableResist") %>% 
  pull(PatientID) %>% 
  unique() %>% 
  droplevels()
# Patients with increasing AR gene carriage
IncreasedResistPat <- meta(AbsRfSputumAllStartEndPs) %>%
  filter(PatientResistStatus == "IncreasedResist") %>% 
  pull(PatientID) %>% 
  unique() %>% 
  droplevels()
# Patients with stable or increasing AR gene carriage
KnownResistPat <- meta(AbsRfSputumAllStartEndPs) %>% 
  pull(PatientID) %>% 
  unique()

# Phyloseq object focused on patients with known AR gene carriage status and
# samples collected either at the end of treatment or later
AbsRfSputumKnownResistPatEndPostAZTPs <- AbsRfPatientPs %>%
  subset_samples(OriginType3 == "Sputum" &
                   PatientID %in% KnownResistPat &
                   TTTPrePost5gr %in% c("EndAZT", "PostAZT_1mo", "PostAZT>=3mo"))

# Agglomerating at genus level
AbsRfSputumKnownResistPatEndPostAZTGenusPs <- tax_glom(AbsRfSputumKnownResistPatEndPostAZTPs, "Genus")
# Transforming to relative counts
RelRfSputumKnownResistPatEndPostAZTGenusPs <- transform_sample_counts(AbsRfSputumKnownResistPatEndPostAZTGenusPs,
                                                                      function(x) {return(x / sum(x))})
# Keping only the 2 genera of interest
RelRfSputumKnownResistPatEndPostAZTPrevStreptPs <- subset_taxa(RelRfSputumKnownResistPatEndPostAZTGenusPs,
                                                               Genus %in% c("Prevotella_7", "Streptococcus"))

RelRfSputumKnownResistPatEndPostAZTPrevStreptDF <- psmelt(RelRfSputumKnownResistPatEndPostAZTPrevStreptPs) %>% 
  select(PatientID, TTTPrePost5gr, Genus, Abundance) %>% 
  mutate(PatientResistStatus = ifelse(PatientID %in% StableResistPat, "StableResist", "IncreasedResist")) %>%
  mutate(PatientResistStatus = factor(PatientResistStatus, levels = c("StableResist", "IncreasedResist")))

# Prevotella_7 analysis
RelRfSputumKnownResistPatEndPostAZTPrevotellaDF <- RelRfSputumKnownResistPatEndPostAZTPrevStreptDF %>% 
  filter(Genus == "Prevotella_7")

# Box plot showing the difference in Prevotella_7 relative abundance per sputum sample
# between patients with stable or increased AR gene carriage (Supplementary Figure E16, panel c)
RelRfSputumKnownResistPatEndPostAZTPrevotella_boxplot <- RelRfSputumKnownResistPatEndPostAZTPrevotellaDF %>%
  ggplot(aes(x = PatientResistStatus,
             y = Abundance,
             fill = PatientResistStatus)) +
  geom_boxplot(width = .5,
               size = .5,
               alpha = .3,
               position = position_dodge(.7),
               outlier.shape = NA,
               show.legend = FALSE) +
  geom_jitter(shape = 21,
              colour = "black",
              size = 1.5,
              alpha = .6,
              width = .2) +
  stat_compare_means(paired = FALSE) +
  scale_fill_manual(values = ResistStatusCol2) +
  labs(x = "", y = "Relative abundance") +
  theme_classic() +
  theme(strip.background = element_blank(),
        axis.text.x = element_blank())

ggsave("RelRfSputumKnownResistPatEndPostAZTPrevotella_boxplot.pdf",
       width = 10, height = 10, dpi = 200, units = "cm")

# Streptococcus analysis
RelRfSputumKnownResistPatEndPostAZTStreptDF <- RelRfSputumKnownResistPatEndPostAZTPrevStreptDF %>% 
  filter(Genus == "Streptococcus")

# Box plot showing the difference in Streptococcus relative abundance per sputum sample
# between patients with stable or increased AR gene carriage (Supplementary Figure E16, panel d)
RelRfSputumKnownResistPatEndPostAZTStrept_boxplot <- RelRfSputumKnownResistPatEndPostAZTStreptDF %>%
  ggplot(aes(x = PatientResistStatus,
             y = Abundance,
             fill = PatientResistStatus)) +
  geom_boxplot(width = .5,
               size = .5,
               alpha = .3,
               position = position_dodge(.7),
               outlier.shape = NA,
               show.legend = FALSE) +
  geom_jitter(shape = 21,
              colour = "black",
              size = 1.5,
              alpha = .6,
              width = .2) +
  stat_compare_means(paired = FALSE) +
  scale_fill_manual(values = ResistStatusCol2) +
  labs(x = "", y = "Relative abundance") +
  theme_classic() +
  theme(strip.background = element_blank(),
        axis.text.x = element_blank())

ggsave("RelRfSputumKnownResistPatEndPostAZTStrept_boxplot.pdf",
       width = 10, height = 10, dpi = 200, units = "cm")

# Intra-individual genus level analysis
AbsRfSputumAllStartEndGenusAllPs <- tax_glom(AbsRfSputumAllStartEndPs, "Genus")
# Transforming to relative counts
RelRfSputumAllStartEndGenusAllPs <- transform_sample_counts(AbsRfSputumAllStartEndGenusAllPs,
                                                            function(x) {return(x / sum(x))})
# Keping only the 2 genera
RelRfSputumAllStartEndGenusPrevStreptPs <- subset_taxa(RelRfSputumAllStartEndGenusAllPs,
                                                       Genus %in% c("Prevotella_7", "Streptococcus"))

RelRfSputumAllStartEndPrevStreptDF <- psmelt(RelRfSputumAllStartEndGenusPrevStreptPs) %>% 
  select(SampleID, PatientID, TTTPrePost5gr, PatientResistStatus, Genus, Abundance)

# Prevotella_7 analysis
RelRfSputumAllStartEndPrevotellaDF <- RelRfSputumAllStartEndPrevStreptDF %>% 
  filter(Genus == "Prevotella_7")

# Paired box plot of the difference in Prevotella_7 relative abundance per sputum sample
# between StartAZT and EndAZT with facet by PatientResistStatus (Supplementary Figure E16, panel a)
BoxPaired_StartAZT_EndAZT_Prevotella <- ggpaired(RelRfSputumAllStartEndPrevotellaDF, 
                                                 x = "TTTPrePost5gr", 
                                                 y = "Abundance",
                                                 id = "PatientID",
                                                 color = "Black",
                                                 fill = "TTTPrePost5gr",
                                                 line.color = "black",
                                                 line.size = .4,
                                                 point.size = 1.2,
                                                 palette = StartEndTTTCol2,
                                                 facet.by = "PatientResistStatus",
                                                 ylab = "Number of ASVs per sample") +
  stat_compare_means(paired = TRUE) +
  theme(legend.position = "right") +
  theme(panel.grid = element_blank()) +
  theme(axis.text.x = element_blank()) +
  theme(axis.title.x = element_blank())

ggsave("BoxPaired_StartAZT_EndAZT_Prevotella.pdf",
       width = 12, height = 10, dpi = 200, units = "cm", limitsize = FALSE)

# Streptococcus analysis
RelRfSputumAllStartEndStreptDF <- RelRfSputumAllStartEndPrevStreptDF %>% 
  filter(Genus == "Streptococcus")

# Paired box plot of the difference in Streptococcus relative abundance per sputum sample
# between StartAZT and EndAZT with facet by PatientResistStatus (Supplementary Figure E16, panel b)
BoxPaired_StartAZT_EndAZT_Strept <- ggpaired(RelRfSputumAllStartEndStreptDF, 
                                             x = "TTTPrePost5gr", 
                                             y = "Abundance",
                                             id = "PatientID",
                                             color = "Black",
                                             fill = "TTTPrePost5gr",
                                             line.color = "black",
                                             line.size = .4,
                                             point.size = 1.2,
                                             palette = StartEndTTTCol2,
                                             facet.by = "PatientResistStatus",
                                             ylab = "Number of ASVs per sample") +
  stat_compare_means(paired = TRUE) +
  theme(legend.position = "right") +
  theme(panel.grid = element_blank()) +
  theme(axis.text.x = element_blank()) +
  theme(axis.title.x = element_blank())

ggsave("BoxPaired_StartAZT_EndAZT_Strept.pdf",
       width = 12, height = 10, dpi = 200, units = "cm", limitsize = FALSE)

# Sputum: DPCoA with distinction by StartAZT/EndAZT and Resistance status (Figure 6, panels b and c) ----
# Genus agglomeration
AbsRfSputumAllStartEndPs %<>% subset_taxa(., Genus != "Unknown")

AbsRfSputumAllStartEndGenus15Ps <- aggregate_top_taxa2(AbsRfSputumAllStartEndPs,
                                                                  "Genus", top = 15)

taxa_names(AbsRfSputumAllStartEndGenus15Ps)

# Agglomerating at genus level
AbsRfSputumAllStartEndGenusAllPs <- tax_glom(AbsRfSputumAllStartEndPs, "Genus")
# Keeping only the top 15 abundant genera
AbsRfSputumAllStartEndGenus15Ps <- subset_taxa(AbsRfSputumAllStartEndGenusAllPs, 
                                               Genus %in% taxa_names(AbsRfSputumAllStartEndGenus15Ps))
# DPCoA at genus level
AbsRfSputumAllStartEndGenus15DPCoA <- ordinate(AbsRfSputumAllStartEndGenus15Ps, 
                                               "DPCoA")

# Plotting Taxa only at genus level (For Figure 6b)
SputumStartEndGenus15DPCoA_TaxaPlot <- plot_ordination(AbsRfSputumAllStartEndGenus15Ps,
                                                       AbsRfSputumAllStartEndGenus15DPCoA,
                                                       color = "Phylum",
                                                       label = NULL,
                                                       type = "taxa") +
  geom_point(size = 4.5, alpha = .6) +
  scale_color_manual(values = MyBactCol6BlFr) +
  # Taxa names
  geom_text(aes(label = Genus), color = "black", size = 3.2, vjust = 1.5)

ggsave("SputumStartEndGenus15DPCoA_TaxaPlot.pdf",
       width = 17, height = 12, dpi = 200, units = "cm", limitsize = FALSE)

# Plotting Samples only at genus level (For Figure 6c)
SputumStartEndGenus15DPCoA_SamplePlot <- plot_ordination(AbsRfSputumAllStartEndGenus15Ps,
                                                         AbsRfSputumAllStartEndGenus15DPCoA,
                                                         color = "PatientResistStatus",
                                                         label = NULL,
                                                         shape = "TTTPrePost5gr",
                                                         type = "samples") +
  geom_point(size = 4, alpha = .6) +
  scale_color_manual(values = ResistStatusCol2) +
  # Sample names
  geom_text(aes(label = PatientID), color = "black", size = 3.2, vjust = 1.5)

ggsave("SputumStartEndGenus15DPCoA_SamplePlot.pdf",
       width = 17, height = 12, dpi = 200, units = "cm", limitsize = FALSE)

# Analyses integrating lung functions (Supplementary Figure E18, panels a and b) ----
# Including status of ARG carriage in the metadata of the initial phyloseq object
ForJoinResist <- meta(AbsRfSputumAllStartEndPs) %>% select(PatientID,  
                                                           PatientResistStatus, 
                                                           ResistGenesEndStartFoldCh)

MetaHelRfSputumPs2 <- meta(HelRfSputumPs) %>% 
  left_join(., 
            ForJoinResist, 
            by = "PatientID") %>% 
  distinct()

rownames(MetaHelRfSputumPs2) <- MetaHelRfSputumPs2$SampleID

HelRfSputumPs2 <- merge_phyloseq(HelRfSputumPs, 
                                 sample_data(MetaHelRfSputumPs2))

metaLungFunctDF <- meta(HelRfSputumPs2) %>% 
  select(SampleID, 
         PatientID, 
         TTTPrePost7gr, 
         TLC_percent:DLCO_name, 
         PatientResistStatus, 
         ResistGenesEndStartFoldCh)

# Plotting percent FVC per status of ARG carriage (Supplementary Figure E18, panel a)
FVC_percentBoxplot <- metaLungFunctDF %>%
  filter(PatientResistStatus %in% c("StableResist", "IncreasedResist")) %>% 
  ggplot(aes(x = TTTPrePost7gr,
             y = FVC_percent,
             fill = TTTPrePost7gr)) +
  geom_boxplot(width = .5,
               size = .3,
               alpha = .4,
               position = position_dodge(.7),
               outlier.shape = NA,
               show.legend = FALSE) +
  geom_line(aes(group = PatientID),
            colour = "#939696",
            size = .5) +
  geom_point(shape = 21,
             colour = "black",
             size = 1.5,
             alpha = .6) +
  scale_fill_manual(values = TTTphaseCol7) +
  stat_compare_means(label.y = 10) +
  facet_wrap(~ PatientResistStatus, ncol = 2) +
  labs(x = "", y = "FVC (percent)") +
  theme_classic() +
  theme(strip.background = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(filename = "FVC_percentBoxplot.pdf",
       width = 13, height = 12, dpi = 200, units = "cm", device='pdf')

# Plotting percent corrected DLCO per status of ARG carriage (Supplementary Figure E18, panel b)
DLCOco_percentBoxplot <- metaLungFunctDF %>%
  filter(PatientResistStatus %in% c("StableResist", "IncreasedResist")) %>% 
  ggplot(aes(x = TTTPrePost7gr,
             y = DLCO_co_percent,
             fill = TTTPrePost7gr)) +
  geom_boxplot(width = .5,
               size = .3,
               alpha = .4,
               position = position_dodge(.7),
               outlier.shape = NA,
               show.legend = FALSE) +
  geom_line(aes(group = PatientID),
            colour = "#939696",
            size = .5) +
  geom_point(shape = 21,
             colour = "black",
             size = 1.5,
             alpha = .6) +
  scale_fill_manual(values = TTTphaseCol7) +
  stat_compare_means(label.y = 15) +
  facet_wrap(~ PatientResistStatus, ncol = 2) +
  labs(x = "", y = "Corrected DLCO (percent)") +
  theme_classic() +
  theme(strip.background = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(filename = "DLCOco_percentBoxplot.pdf",
       width = 13, height = 12, dpi = 200, units = "cm", device='pdf')

# Plotting percent TLC per status of ARG carriage
TLC_percentBoxplot <- metaLungFunctDF %>%
  filter(PatientResistStatus %in% c("StableResist", "IncreasedResist")) %>% 
  ggplot(aes(x = TTTPrePost7gr,
             y = TLC_percent,
             fill = TTTPrePost7gr)) +
  geom_boxplot(width = .5,
               size = .3,
               alpha = .4,
               position = position_dodge(.7),
               outlier.shape = NA,
               show.legend = FALSE) +
  geom_line(aes(group = PatientID),
            colour = "#939696",
            size = .5) +
  geom_point(shape = 21,
             colour = "black",
             size = 1.5,
             alpha = .6) +
  scale_fill_manual(values = TTTphaseCol7) +
  stat_compare_means(label.y = 15) +
  facet_wrap(~ PatientResistStatus, ncol = 2) +
  labs(x = "", y = "Total Lung Capacity (percent)") +
  theme_classic() +
  theme(strip.background = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(filename = "TLC_percentBoxplot.pdf",
       width = 13, height = 12, dpi = 200, units = "cm", device='pdf')

# Plotting percent FEV1 per status of ARG carriage
FEV1_percentBoxplot <- metaLungFunctDF %>%
  filter(PatientResistStatus %in% c("StableResist", "IncreasedResist")) %>% 
  ggplot(aes(x = TTTPrePost7gr,
             y = FEV1_percent,
             fill = TTTPrePost7gr)) +
  geom_boxplot(width = .5,
               size = .3,
               alpha = .4,
               position = position_dodge(.7),
               outlier.shape = NA,
               show.legend = FALSE) +
  geom_line(aes(group = PatientID),
            colour = "#939696",
            size = .5) +
  geom_point(shape = 21,
             colour = "black",
             size = 1.5,
             alpha = .6) +
  scale_fill_manual(values = TTTphaseCol7) +
  stat_compare_means(label.y = 15) +
  facet_wrap(~ PatientResistStatus, ncol = 2) +
  labs(x = "", y = "FEV1 (percent)") +
  theme_classic() +
  theme(strip.background = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(filename = "FEV1_percentBoxplot.pdf",
       width = 13, height = 12, dpi = 200, units = "cm", device='pdf')

# Genus-level analysis - Alpha diversity: Evenness - Richness (Supplementary Figure E9, panel a) ----
AbsRfPatientGenusPs <- AbsRfPatientPs %>% tax_glom(., "Genus")

# Alpha diversity based on a selection of metrics:
AlphaDivGenusDF <- microbiome::alpha(AbsRfPatientGenusPs,
                                     index = c("chao1", 
                                               "evenness_camargo",
                                               "diversity_shannon",
                                               "dominance_dmn",
                                               "dominance_core_abundance"))

# Moving the rownames to a new column in diversity data frame
AlphaDivGenusDF$SampleID <- rownames(AlphaDivGenusDF)

# Extracting the metadata from phyloseq object and merging them with diversity data
AbsRfPatientGenusMetaDF <- microbiome::meta(AbsRfPatientGenusPs)
AbsRfPatientGenusMetaDivDF <- left_join(AbsRfPatientGenusMetaDF, AlphaDivGenusDF, by = "SampleID")

# Plotting Chao1 for Sputum (LRT) and OPS (URT) samples (Supplementary Figure E9, panel a)
Chao1SputumSwabGenusBoxplot <- AbsRfPatientGenusMetaDivDF %>%
  ggplot(aes(x = TTTPrePost5gr, 
             y = chao1,
             fill = TTTPrePost5gr)) +
  geom_boxplot(width = .7,
               size = .3,
               alpha = .4,
               position = position_dodge(.8),
               outlier.shape = NA,
               show.legend = TRUE) +
  geom_jitter(position = position_jitterdodge(.8),
              shape = 21,
              colour = "black",
              size = 1,
              alpha = .7) +
  scale_fill_manual(values = TTTphaseCol5) +
  stat_compare_means(comparisons = CompTTT1) +
  stat_compare_means(label.y = 70) +
  facet_wrap(~Type, ncol = 2) +
  labs(fill = "Treatment phase",
       y = "Chao1 richness") +
  theme_classic() +
  theme(strip.background = element_blank()) +
  rotate_x_text(angle = 60, hjust = 1)

ggsave(filename = "Chao1SputumSwabGenusBoxplot.pdf",
       width = 15, height = 15, dpi = 200, units = "cm", device='pdf')

# Kruskal-Wallis test for Sputum samples
AbsRfPatientGenusMetaDivDF %>%
  filter(OriginType3 == "Sputum") %>% 
  droplevels() %>% 
  kruskal_test(data = .,
               chao1 ~ TTTPrePost5gr) %>%
  adjust_pvalue(method = "BH") %>%
  add_significance("p.adj")
# Dunn's test for pairwise multiple comparisons
AbsRfPatientGenusMetaDivDF %>%
  filter(OriginType3 == "Sputum") %>% 
  droplevels() %>%
  dunn_test(data = .,
            chao1 ~ TTTPrePost5gr,
            p.adjust.method = "BH")

# Kruskal-Wallis test for OPS samples
AbsRfPatientGenusMetaDivDF %>%
  filter(OriginType3 == "Swab") %>% 
  droplevels() %>% 
  kruskal_test(data = .,
               chao1 ~ TTTPrePost5gr) %>%
  adjust_pvalue(method = "BH") %>%
  add_significance("p.adj")
# Dunn's test for pairwise multiple comparisons
AbsRfPatientGenusMetaDivDF %>%
  filter(OriginType3 == "Swab") %>% 
  droplevels() %>%
  dunn_test(data = .,
            chao1 ~ TTTPrePost5gr,
            p.adjust.method = "BH")

# Genus-level analysis - Faith's phylogenetic diversity (Supplementary Figure E9, panel b) ----
# Corresponds to the sum of the branch lengths of phylogenetic tree.
# Obtaining a DF with transposed data from the ASV table
AbsRfPatientGenusASVtabTranspDF <- as.data.frame(otu_table(AbsRfPatientGenusPs)) %>% t()
# Phylogenetic diversity will be measured from a rooted tree
AbsRfPatientGenusTree <- phy_tree(AbsRfPatientGenusPs)
# Using picante::pd to measure phylogenetic diversity
PhyloDivSputumSwabGenusDF <- pd(AbsRfPatientGenusASVtabTranspDF, AbsRfPatientGenusTree, include.root = T)

PhyloDivSputumSwabGenusDF %<>% dplyr::rename(FaithPhylogenDiv = PD, SpeciesRichness = SR)
PhyloDivSputumSwabGenusDF$SampleID <- rownames(PhyloDivSputumSwabGenusDF)

AbsRfPatientGenusMetaDivPhyloDivDF <- left_join(AbsRfPatientGenusMetaDivDF, 
                                           PhyloDivSputumSwabGenusDF, by = "SampleID")

# Plotting Faith's phylogenetic diversity for Sputum (LRT) and Swab (URT) (Supplementary Figure E9, panel b)
PhyloDivSputumSwabGenusBoxplot <- AbsRfPatientGenusMetaDivPhyloDivDF %>%
  ggplot(aes(x = TTTPrePost5gr, 
             y = FaithPhylogenDiv,
             fill = TTTPrePost5gr)) +
  geom_boxplot(width = .7,
               size = .3,
               alpha = .4,
               position = position_dodge(.8),
               outlier.shape = NA,
               show.legend = TRUE) +
  geom_jitter(position = position_jitterdodge(.8),
              shape = 21,
              colour = "black",
              size = 1,
              alpha = .7) +
  scale_fill_manual(values = TTTphaseCol5) +
  stat_compare_means(comparisons = CompTTT1) +
  stat_compare_means(label.y = 310) +
  facet_wrap(~Type, ncol = 2) +
  labs(fill = "Treatment phase",
       y = "Faith's phylogenetic diversity") +
  theme_classic() +
  theme(strip.background = element_blank()) +
  rotate_x_text(angle = 60, hjust = 1)

ggsave(filename = "PhyloDivSputumSwabGenusBoxplot.pdf",
       width = 15, height = 15, dpi = 200, units = "cm", device='pdf')

# Kruskal-Wallis test for Sputum samples
AbsRfPatientGenusMetaDivPhyloDivDF %>%
  filter(OriginType3 == "Sputum") %>% 
  droplevels() %>% 
  kruskal_test(data = .,
               FaithPhylogenDiv ~ TTTPrePost5gr) %>%
  adjust_pvalue(method = "BH") %>%
  add_significance("p.adj")
# Dunn's test for pairwise multiple comparisons
AbsRfPatientGenusMetaDivPhyloDivDF %>%
  filter(OriginType3 == "Sputum") %>% 
  droplevels() %>%
  dunn_test(data = .,
            FaithPhylogenDiv ~ TTTPrePost5gr,
            p.adjust.method = "BH")

# Kruskal-Wallis test for OPS samples
AbsRfPatientGenusMetaDivPhyloDivDF %>%
  filter(OriginType3 == "Swab") %>% 
  droplevels() %>% 
  kruskal_test(data = .,
               FaithPhylogenDiv ~ TTTPrePost5gr) %>%
  adjust_pvalue(method = "BH") %>%
  add_significance("p.adj")
# Dunn's test for pairwise multiple comparisons
AbsRfPatientGenusMetaDivPhyloDivDF %>%
  filter(OriginType3 == "Swab") %>% 
  droplevels() %>%
  dunn_test(data = .,
            FaithPhylogenDiv ~ TTTPrePost5gr,
            p.adjust.method = "BH")

# Genus-level analysis - Distance between OPS samples and the centroid of Sputum samples (Supplementary Figure E13, panel a) ----
HelRfPatientGenusPs <- HelRfPatientPs %>% tax_glom(., "Genus")

# Computing the distance matrix 
TTTPatientsGenus_UniFrac <- phyloseq::distance(HelRfPatientGenusPs, method = "unifrac")
# Data frame with Secondary data
PatientsGenusSecDataDF <- sample_data(HelRfPatientGenusPs) %>% data.frame()
# Adding a "TTTphase_Type" column
PatientsGenusSecData_TTT_TypeDF <- PatientsGenusSecDataDF %>% 
  mutate(TTTPrePost5grCopy = TTTPrePost5gr, TypeCopy = Type) %>% 
  unite(TTTphase_Type, c(TTTPrePost5grCopy, TypeCopy), sep = "_") %>% 
  convert(fct(TTTphase_Type))
# Computing distance between each sample and each centroid
PatientsGenus_centr_distDF <- usedist::dist_to_centroids(TTTPatientsGenus_UniFrac, 
                                                    PatientsGenusSecData_TTT_TypeDF$TTTphase_Type)
# Adding TTT group variable in data frame with distances to centroids
PatientsGenus_centr_distDF %<>% dplyr::rename(SampleID = Item)
PatientsGenus_centr_distDF2 <- left_join(PatientsGenus_centr_distDF, 
                                    PatientsGenusSecData_TTT_TypeDF, 
                                    by = "SampleID") %>% 
  select(SampleID, 
         CentroidGroup, 
         CentroidDistance, 
         TTTphase_Type, 
         TTTPrePost5gr, 
         Type)

# Data frame with Swab sample distances to Sputum centroid for each TTT phase
SwabSamples_to_SputumCentroidsGenusDF <- PatientsGenus_centr_distDF2 %>%
  filter(Type == "Swab") %>%
  filter(CentroidGroup %in% c("PreAZT_Sputum", "StartAZT_Sputum", "EndAZT_Sputum", 
                              "PostAZT_1mo_Sputum", "PostAZT>=3mo_Sputum")) %>% 
  mutate(CentroidGroup = factor(CentroidGroup, 
                                levels = c("PreAZT_Sputum", "StartAZT_Sputum", "EndAZT_Sputum", 
                                           "PostAZT_1mo_Sputum", "PostAZT>=3mo_Sputum"),
                                labels = c("PreAZT", "StartAZT", "EndAZT", 
                                           "PostAZT_1mo", "PostAZT>=3mo"))) %>% 
  arrange(CentroidGroup)

# Plotting Swab sample distances to Sputum centroid for each TTT phase (Supplementary Figure E13, panel a)
SwabSamples_to_SputumCentroidsGenusBoxplot <- SwabSamples_to_SputumCentroidsGenusDF %>%
  ggplot(aes(x = CentroidGroup, 
             y = CentroidDistance,
             fill = CentroidGroup)) +
  geom_boxplot(width = .7,
               size = .3,
               alpha = .4,
               position = position_dodge(.8),
               outlier.shape = NA,
               show.legend = TRUE) +
  geom_jitter(position = position_jitterdodge(.8),
              shape = 21,
              colour = "black",
              size = 1,
              alpha = .7) +
  scale_fill_manual(values = TTTphaseCol5) +
  stat_compare_means(comparisons = CompTTT2) +
  stat_compare_means(label.y = .87) +
  labs(fill = "Treatment phase",
       y = "OPS to sputum distance") +
  theme_classic() +
  theme(strip.background = element_blank()) +
  rotate_x_text(angle = 60, hjust = 1)

ggsave(filename = "SwabSamples_to_SputumCentroidsGenusBoxplot.pdf",
       width = 12, height = 12, dpi = 200, units = "cm", device='pdf')

# Kruskal-Wallis test
SwabSamples_to_SputumCentroidsGenusDF %>%
  kruskal_test(data = .,
               CentroidDistance ~ CentroidGroup) %>%
  adjust_pvalue(method = "BH") %>%
  add_significance("p.adj")
# Dunn's test for pairwise multiple comparisons
SwabSamples_to_SputumCentroidsGenusDF %>%
  dunn_test(data = .,
            CentroidDistance ~ CentroidGroup,
            p.adjust.method = "BH")

# Genus-level analysis - Microbiota dynamics in Sputum and OPS (Supplementary Figure E13, panel b) ----
# Analyses on all patients with available pairs of sputum samples independent of AR gene carriage
# 1A) Determining the proportion of taxa gained/conserved/lost between StartATZ and EndAZT
# Patients with a pair of sputum samples available
PatientsGenusStartEndSputum <- meta(AbsRfPatientGenusPs) %>% 
  filter(Type == "Sputum") %>%
  filter(TTTPrePost5gr == "StartAZT" | TTTPrePost7gr == "EndAZT") %>% 
  group_by(PatientID) %>% 
  dplyr::summarize(NumberSamplesAvailable = n()) %>% 
  filter(NumberSamplesAvailable == 2) %>% 
  select(PatientID)

PatientsGenusStartEndSputum$PatientID %<>% droplevels()

AbsRfGenusSputumStartEndPs <- subset_samples(AbsRfPatientGenusPs,
                                        PatientID %in% PatientsGenusStartEndSputum$PatientID &
                                          Type == "Sputum" &
                                          TTTPrePost5gr %in% c("StartAZT",
                                                               "EndAZT"))
# Keeping only the taxa represented after filtering, then counting them and listing their names
AbsRfGenusSputumStartEndPs %<>% prune_taxa(taxa_sums(.) > 0, .)
SputumStartEndGenusNumber <- AbsRfGenusSputumStartEndPs %>% ntaxa()
SputumStartEndGenusNames <- AbsRfGenusSputumStartEndPs %>% taxa_names()
# 1B) Inspecting taxa present in sputum samples at StartAZT
AbsRfGenusSputumStartAZTPs <- subset_samples(AbsRfGenusSputumStartEndPs,
                                             TTTPrePost5gr == "StartAZT")
# Keeping only the taxa represented after filtering, then counting them and listing their names
AbsRfGenusSputumStartAZTPs %<>% prune_taxa(taxa_sums(.) > 0, .)
SputumStartGenusNumber <- AbsRfGenusSputumStartAZTPs %>% ntaxa()
SputumStartGenusNames <- AbsRfGenusSputumStartAZTPs %>% taxa_names()
# 1C) Inspecting taxa present in sputum samples at EndAZT
AbsRfGenusSputumEndAZTPs <- subset_samples(AbsRfGenusSputumStartEndPs,
                                           TTTPrePost5gr == "EndAZT")
# Keeping only the taxa represented after filtering, then counting them and listing their names
AbsRfGenusSputumEndAZTPs %<>% prune_taxa(taxa_sums(.) > 0, .)
SputumEndGenusNumber <- AbsRfGenusSputumEndAZTPs %>% ntaxa()
SputumEndGenusNames <- AbsRfGenusSputumEndAZTPs %>% taxa_names()

# Preparing a list of vectors with taxa names for visualisation by Venn diagram
SputumStartEndAZTGenusList <- list(StartAZT = SputumStartGenusNames,
                              EndAZT = SputumEndGenusNames)
# Calculating the coordinates required for plotting a Venn diagram
SputumStartEndAZTGenusVenn <- ggVennDiagram::Venn(SputumStartEndAZTGenusList)
SputumStartEndAZTGenusVennData <- ggVennDiagram::process_data(SputumStartEndAZTGenusVenn)
# Obtaining the Venn diagram (For Supplementary Figure E13, panel b)
SputumStartEndAZTGenusVennDiagram <- ggVennDiagram(SputumStartEndAZTGenusList,
                                              set_color = "white", # Hiding original label
                                              label =  "both", # Showing "count" and "percent"
                                              label_size = 10,
                                              label_alpha = .6,
                                              edge_lty = "blank") +
  # to customize the label layer
  geom_sf_text(aes(label = name),
               fontface = "plain",
               size = 9,
               data = venn_setlabel(SputumStartEndAZTGenusVennData)) +
  # to customize the edge layer
  geom_sf(size = .1,
          color = "#ff7d1a",
          alpha = .6,
          data = venn_setedge(SputumStartEndAZTGenusVennData)) +
  scale_fill_gradient(low = "#d2ebf9", high = "#E69F00") +
  theme_void() +
  theme(legend.position = "none")

ggsave("SputumStartEndAZTGenusVennDiagram.pdf",
       width = 12, height = 10, dpi = 200, units = "cm", limitsize = FALSE)

# 2A) Determining the proportion of taxa gained/conserved/lost between PostAZT_1mo and PostAZT>=3mo
# Patients with a pair of sputum samples available
PatientsPost1moPost4moGenusSputum <- meta(AbsRfPatientGenusPs) %>% 
  filter(Type == "Sputum") %>%
  filter(TTTPrePost7gr == "PostAZT_1mo" | TTTPrePost7gr == "PostAZT_4mo") %>% 
  group_by(PatientID) %>% 
  dplyr::summarize(NumberSamplesAvailable = n()) %>% 
  filter(NumberSamplesAvailable >= 2) %>% 
  select(PatientID)

PatientsPost1moPost4moGenusSputum$PatientID %<>% droplevels()

AbsRfGenusSputumPost1moPost4moPs <- subset_samples(AbsRfPatientGenusPs,
                                              PatientID %in% PatientsPost1moPost4moGenusSputum$PatientID &
                                                Type == "Sputum" &
                                                TTTPrePost7gr %in% c("PostAZT_1mo",
                                                                     "PostAZT_4mo"))

# Keeping only the taxa represented after filtering, then counting them and listing their names
AbsRfGenusSputumPost1moPost4moPs %<>% prune_taxa(taxa_sums(.) > 0, .)
SputumPost1moPost4moGenusNumber <- AbsRfGenusSputumPost1moPost4moPs %>% ntaxa()
SputumPost1moPost4moGenusNames <- AbsRfGenusSputumPost1moPost4moPs %>% taxa_names()
# 2B) Inspecting taxa present in sputum samples at PostAZT_1mo
AbsRfGenusSputumPost1moPs <- subset_samples(AbsRfGenusSputumPost1moPost4moPs,
                                       TTTPrePost7gr == "PostAZT_1mo")
# Keeping only the taxa represented after filtering, then counting them and listing their names
AbsRfGenusSputumPost1moPs %<>% prune_taxa(taxa_sums(.) > 0, .)
SputumPost1moGenusNumber <- AbsRfGenusSputumPost1moPs %>% ntaxa()
SputumPost1moGenusNames <- AbsRfGenusSputumPost1moPs %>% taxa_names()
# 2C) Inspecting taxa present in sputum samples at PostAZT_4mo
AbsRfGenusSputumPost4moPs <- subset_samples(AbsRfGenusSputumPost1moPost4moPs,
                                       TTTPrePost7gr == "PostAZT_4mo")
# Keeping only the taxa represented after filtering, then counting them and listing their names
AbsRfGenusSputumPost4moPs %<>% prune_taxa(taxa_sums(.) > 0, .)
SputumPost4moGenusNumber <- AbsRfGenusSputumPost4moPs %>% ntaxa()
SputumPost4moGenusNames <- AbsRfGenusSputumPost4moPs %>% taxa_names()

# Preparing a list of vectors with taxa names for visualisation by Venn diagram
SputumPost1moPost4moGenusList <- list(PostAZT_1mo = SputumPost1moGenusNames,
                                 PostAZT_4mo = SputumPost4moGenusNames)
# Calculating the coordinates required for plotting a Venn diagram
SputumPost1moPost4moGenusVenn <- Venn(SputumPost1moPost4moGenusList)
SputumPost1moPost4moGenusVennData <- process_data(SputumPost1moPost4moGenusVenn)
# Obtaining the Venn diagram (For Supplementary Figure E13, panel b)
SputumPost1moPost4moGenusVennDiagram <- ggVennDiagram(SputumPost1moPost4moGenusList,
                                                 set_color = "white", # Hiding original label
                                                 label =  "both", # Showing "count" and "percent"
                                                 label_size = 10,
                                                 label_alpha = .6,
                                                 edge_lty = "blank") +
  # to customize the label layer
  geom_sf_text(aes(label = name),
               fontface = "plain",
               size = 9,
               data = venn_setlabel(SputumPost1moPost4moGenusVennData)) +
  # to customize the edge layer
  geom_sf(size = .1,
          color = "#ff7d1a",
          alpha = .6,
          data = venn_setedge(SputumPost1moPost4moGenusVennData)) +
  scale_fill_gradient(low = "#00CC96", high = "#CC79A7") +
  theme_void() +
  theme(legend.position = "none")

ggsave("SputumPost1moPost4moGenusVennDiagram.pdf",
       width = 12, height = 10, dpi = 200, units = "cm", limitsize = FALSE)

# 3A) Determining the proportion of taxa gained/conserved/lost between PreAZT_4mo and PreAZT_1mo
# in patients receiving placebo first.
# Keeping only patients with a pair of sputum samples available
PatientsGenusPlaceboFirst <- meta(AbsRfPatientGenusPs) %>% 
  filter(Type == "Sputum") %>%
  filter(TTTPrePost7gr == "PreAZT_4mo" | TTTPrePost7gr == "PreAZT_1mo") %>% 
  group_by(PatientID) %>% 
  dplyr::summarize(NumberSamplesAvailable = n()) %>% 
  filter(NumberSamplesAvailable == 2) %>% 
  select(PatientID)

PatientsGenusPlaceboFirst$PatientID %<>% droplevels()

AbsRfGenusSputumPlaceboFirst4Mo1MoPs <- subset_samples(AbsRfPatientGenusPs,
                                                  PatientID %in% PatientsGenusPlaceboFirst$PatientID &
                                                    Type == "Sputum" &
                                                    TTTPrePost7gr %in% c("PreAZT_1mo",
                                                                         "PreAZT_4mo"))
# Keeping only the taxa represented after filtering, then counting them and listing their names
AbsRfGenusSputumPlaceboFirst4Mo1MoPs %<>% prune_taxa(taxa_sums(.) > 0, .)
SputumPlaceboFirst4Mo1MoGenusNumber <- AbsRfGenusSputumPlaceboFirst4Mo1MoPs %>% ntaxa()
SputumPlaceboFirst4Mo1MoGenusNames <- AbsRfGenusSputumPlaceboFirst4Mo1MoPs %>% taxa_names()
# 3B) Inspecting taxa present at PreAZT_4mo in the sputum of patients receiving placebo first
AbsRfGenusSputumPlaceboFirstPre4MoPs <- subset_samples(AbsRfGenusSputumPlaceboFirst4Mo1MoPs,
                                                       TTTPrePost7gr == "PreAZT_4mo")
# Keeping only the taxa represented after filtering, then counting them and listing their names
AbsRfGenusSputumPlaceboFirstPre4MoPs %<>% prune_taxa(taxa_sums(.) > 0, .)
SputumPlaceboFirst4MoGenusNumber <- AbsRfGenusSputumPlaceboFirstPre4MoPs %>% ntaxa()
SputumPlaceboFirst4MoGenusNames <- AbsRfGenusSputumPlaceboFirstPre4MoPs %>% taxa_names()
# 3C) Inspecting taxa present at PreAZT_1mo in the sputum of patients receiving placebo first
AbsRfGenusSputumPlaceboFirstPre1MoPs <- subset_samples(AbsRfGenusSputumPlaceboFirst4Mo1MoPs,
                                                       TTTPrePost7gr == "PreAZT_1mo")
# Keeping only the taxa represented after filtering, then counting them and listing their names
AbsRfGenusSputumPlaceboFirstPre1MoPs %<>% prune_taxa(taxa_sums(.) > 0, .)
SputumPlaceboFirst1MoGenusNumber <- AbsRfGenusSputumPlaceboFirstPre1MoPs %>% ntaxa()
SputumPlaceboFirst1MoGenusNames <- AbsRfGenusSputumPlaceboFirstPre1MoPs %>% taxa_names()

# Preparing a list of vectors with taxa names for visualisation by Venn diagram
SputumPlaceboFirst4Mo1MoGenusList <- list("PreTTT_4mo" = SputumPlaceboFirst4MoGenusNames,
                                     "PreTTT_1moEndAZT" = SputumPlaceboFirst1MoGenusNames)
# Calculating the coordinates required for plotting a Venn diagram
SputumPlaceboFirst4Mo1MoGenusVenn <- Venn(SputumPlaceboFirst4Mo1MoGenusList)
SputumPlaceboFirst4Mo1MoGenusVennData <- process_data(SputumPlaceboFirst4Mo1MoGenusVenn)
# Obtaining the Venn diagram (For Supplementary Figure E13, panel b)
SputumPlaceboFirst4Mo1MoGenusVennDiagram <- ggVennDiagram(SputumPlaceboFirst4Mo1MoGenusList,
                                                     set_color = "white", # Hiding original label
                                                     label =  "both", # Showing "count" and "percent"
                                                     label_size = 10,
                                                     label_alpha = .6,
                                                     edge_lty = "blank") +
  # to customize the label layer
  geom_sf_text(aes(label = name),
               fontface = "plain",
               size = 9,
               data = venn_setlabel(SputumPlaceboFirst4Mo1MoGenusVennData)) +
  # to customize the edge layer
  geom_sf(size = .1,
          color = "#f2f2f2",
          alpha = .6,
          data = venn_setedge(SputumPlaceboFirst4Mo1MoGenusVennData)) +
  scale_fill_gradient(low = "#f0f0f0", high = "#999999") +
  theme_void() +
  theme(legend.position = "none")

ggsave("SputumPlaceboFirst4Mo1MoGenusVennDiagram.pdf",
       width = 12, height = 10, dpi = 200, units = "cm", limitsize = FALSE)

# Analyses on all patients with available pairs of OPS samples independent of AR gene carriage
# 1A) Determining the proportion of taxa gained/conserved/lost between StartATZ and EndAZT
# Patients with a pair of OPS samples available
PatientsGenusStartEndOPS <- meta(AbsRfPatientGenusPs) %>% 
  filter(Type == "Swab") %>%
  filter(TTTPrePost5gr == "StartAZT" | TTTPrePost7gr == "EndAZT") %>% 
  group_by(PatientID) %>% 
  dplyr::summarize(NumberSamplesAvailable = n()) %>% 
  filter(NumberSamplesAvailable == 2) %>% 
  select(PatientID)

PatientsGenusStartEndOPS$PatientID %<>% droplevels()

AbsRfGenusOPSStartEndPs <- subset_samples(AbsRfPatientGenusPs,
                                          PatientID %in% PatientsGenusStartEndOPS$PatientID &
                                            Type == "Swab" &
                                            TTTPrePost5gr %in% c("StartAZT",
                                                                 "EndAZT"))
# Keeping only the taxa represented after filtering, then counting them and listing their names
AbsRfGenusOPSStartEndPs %<>% prune_taxa(taxa_sums(.) > 0, .)
OPSStartEndGenusNumber <- AbsRfGenusOPSStartEndPs %>% ntaxa()
OPSStartEndGenusNames <- AbsRfGenusOPSStartEndPs %>% taxa_names()
# 1B) Inspecting taxa present in OPS samples at StartAZT
AbsRfGenusOPSStartAZTPs <- subset_samples(AbsRfGenusOPSStartEndPs,
                                          TTTPrePost5gr == "StartAZT")
# Keeping only the taxa represented after filtering, then counting them and listing their names
AbsRfGenusOPSStartAZTPs %<>% prune_taxa(taxa_sums(.) > 0, .)
OPSStartGenusNumber <- AbsRfGenusOPSStartAZTPs %>% ntaxa()
OPSStartGenusNames <- AbsRfGenusOPSStartAZTPs %>% taxa_names()
# 1C) Inspecting taxa present in OPS samples at EndAZT
AbsRfGenusOPSEndAZTPs <- subset_samples(AbsRfGenusOPSStartEndPs,
                                        TTTPrePost5gr == "EndAZT")
# Keeping only the taxa represented after filtering, then counting them and listing their names
AbsRfGenusOPSEndAZTPs %<>% prune_taxa(taxa_sums(.) > 0, .)
OPSEndGenusNumber <- AbsRfGenusOPSEndAZTPs %>% ntaxa()
OPSEndGenusNames <- AbsRfGenusOPSEndAZTPs %>% taxa_names()

# Preparing a list of vectors with taxa names for visualisation by Venn diagram
OPSStartEndAZTGenusList <- list(StartAZT = OPSStartGenusNames,
                                EndAZT = OPSEndGenusNames)
# Calculating the coordinates required for plotting a Venn diagram
OPSStartEndAZTGenusVenn <- ggVennDiagram::Venn(OPSStartEndAZTGenusList)
OPSStartEndAZTGenusVennData <- ggVennDiagram::process_data(OPSStartEndAZTGenusVenn)
# Obtaining the Venn diagram (For Supplementary Figure E13, panel b)
OPSStartEndAZTGenusVennDiagram <- ggVennDiagram(OPSStartEndAZTGenusList,
                                                set_color = "white", # Hiding original label
                                                label =  "both", # Showing "count" and "percent"
                                                label_size = 10,
                                                label_alpha = .6,
                                                edge_lty = "blank") +
  # to customize the label layer
  geom_sf_text(aes(label = name),
               fontface = "plain",
               size = 9,
               data = venn_setlabel(OPSStartEndAZTGenusVennData)) +
  # to customize the edge layer
  geom_sf(size = .1,
          color = "#ff7d1a",
          alpha = .6,
          data = venn_setedge(OPSStartEndAZTGenusVennData)) +
  scale_fill_gradient(low = "#d2ebf9", high = "#E69F00") +
  theme_void() +
  theme(legend.position = "none")

ggsave("OPSStartEndAZTGenusVennDiagram.pdf",
       width = 12, height = 10, dpi = 200, units = "cm", limitsize = FALSE)

# 2A) Determining the proportion of taxa gained/conserved/lost between PostAZT_1mo and PostAZT>=3mo
# Patients with a pair of OPS samples available
PatientsPost1moPost4moGenusOPS <- meta(AbsRfPatientGenusPs) %>% 
  filter(Type == "Swab") %>%
  filter(TTTPrePost7gr == "PostAZT_1mo" | TTTPrePost7gr == "PostAZT_4mo") %>% 
  group_by(PatientID) %>% 
  dplyr::summarize(NumberSamplesAvailable = n()) %>% 
  filter(NumberSamplesAvailable >= 2) %>% 
  select(PatientID)

PatientsPost1moPost4moGenusOPS$PatientID %<>% droplevels()

AbsRfGenusOPSPost1moPost4moPs <- subset_samples(AbsRfPatientGenusPs,
                                                PatientID %in% PatientsPost1moPost4moGenusOPS$PatientID &
                                                  Type == "Swab" &
                                                  TTTPrePost7gr %in% c("PostAZT_1mo",
                                                                       "PostAZT_4mo"))

# Keeping only the taxa represented after filtering, then counting them and listing their names
AbsRfGenusOPSPost1moPost4moPs %<>% prune_taxa(taxa_sums(.) > 0, .)
OPSPost1moPost4moGenusNumber <- AbsRfGenusOPSPost1moPost4moPs %>% ntaxa()
OPSPost1moPost4moGenusNames <- AbsRfGenusOPSPost1moPost4moPs %>% taxa_names()
# 2B) Inspecting taxa present in OPS samples at PostAZT_1mo
AbsRfGenusOPSPost1moPs <- subset_samples(AbsRfGenusOPSPost1moPost4moPs,
                                         TTTPrePost7gr == "PostAZT_1mo")
# Keeping only the taxa represented after filtering, then counting them and listing their names
AbsRfGenusOPSPost1moPs %<>% prune_taxa(taxa_sums(.) > 0, .)
OPSPost1moGenusNumber <- AbsRfGenusOPSPost1moPs %>% ntaxa()
OPSPost1moGenusNames <- AbsRfGenusOPSPost1moPs %>% taxa_names()
# 2C) Inspecting taxa present in OPS samples at PostAZT_4mo
AbsRfGenusOPSPost4moPs <- subset_samples(AbsRfGenusOPSPost1moPost4moPs,
                                         TTTPrePost7gr == "PostAZT_4mo")
# Keeping only the taxa represented after filtering, then counting them and listing their names
AbsRfGenusOPSPost4moPs %<>% prune_taxa(taxa_sums(.) > 0, .)
OPSPost4moGenusNumber <- AbsRfGenusOPSPost4moPs %>% ntaxa()
OPSPost4moGenusNames <- AbsRfGenusOPSPost4moPs %>% taxa_names()

# Preparing a list of vectors with taxa names for visualisation by Venn diagram
OPSPost1moPost4moGenusList <- list(PostAZT_1mo = OPSPost1moGenusNames,
                                   PostAZT_4mo = OPSPost4moGenusNames)
# Calculating the coordinates required for plotting a Venn diagram
OPSPost1moPost4moGenusVenn <- Venn(OPSPost1moPost4moGenusList)
OPSPost1moPost4moGenusVennData <- process_data(OPSPost1moPost4moGenusVenn)
# Obtaining the Venn diagram (For Supplementary Figure E13, panel b)
OPSPost1moPost4moGenusVennDiagram <- ggVennDiagram(OPSPost1moPost4moGenusList,
                                                   set_color = "white", # Hiding original label
                                                   label =  "both", # Showing "count" and "percent"
                                                   label_size = 10,
                                                   label_alpha = .6,
                                                   edge_lty = "blank") +
  # to customize the label layer
  geom_sf_text(aes(label = name),
               fontface = "plain",
               size = 9,
               data = venn_setlabel(OPSPost1moPost4moGenusVennData)) +
  # to customize the edge layer
  geom_sf(size = .1,
          color = "#ff7d1a",
          alpha = .6,
          data = venn_setedge(OPSPost1moPost4moGenusVennData)) +
  scale_fill_gradient(low = "#00CC96", high = "#CC79A7") +
  theme_void() +
  theme(legend.position = "none")

ggsave("OPSPost1moPost4moGenusVennDiagram.pdf",
       width = 12, height = 10, dpi = 200, units = "cm", limitsize = FALSE)

# 3A) Determining the proportion of taxa gained/conserved/lost between PreAZT_4mo and PreAZT_1mo
# in patients receiving placebo first.
# Keeping only patients with a pair of OPS samples available
PatientsGenusPlaceboFirstOPS <- meta(AbsRfPatientGenusPs) %>% 
  filter(Type == "Swab") %>%
  filter(TTTPrePost7gr == "PreAZT_4mo" | TTTPrePost7gr == "PreAZT_1mo") %>% 
  group_by(PatientID) %>% 
  dplyr::summarize(NumberSamplesAvailable = n()) %>% 
  filter(NumberSamplesAvailable == 2) %>% 
  select(PatientID)

PatientsGenusPlaceboFirstOPS$PatientID %<>% droplevels()

AbsRfGenusOPSPlaceboFirst4Mo1MoPs <- subset_samples(AbsRfPatientGenusPs,
                                                    PatientID %in% PatientsGenusPlaceboFirstOPS$PatientID &
                                                      Type == "Swab" &
                                                      TTTPrePost7gr %in% c("PreAZT_1mo",
                                                                           "PreAZT_4mo"))
# Keeping only the taxa represented after filtering, then counting them and listing their names
AbsRfGenusOPSPlaceboFirst4Mo1MoPs %<>% prune_taxa(taxa_sums(.) > 0, .)
OPSPlaceboFirst4Mo1MoGenusNumber <- AbsRfGenusOPSPlaceboFirst4Mo1MoPs %>% ntaxa()
OPSPlaceboFirst4Mo1MoGenusNames <- AbsRfGenusOPSPlaceboFirst4Mo1MoPs %>% taxa_names()
# 3B) Inspecting taxa present at PreAZT_4mo in OPS samples of patients receiving placebo first
AbsRfGenusOPSPlaceboFirstPre4MoPs <- subset_samples(AbsRfGenusOPSPlaceboFirst4Mo1MoPs,
                                                    TTTPrePost7gr == "PreAZT_4mo")
# Keeping only the taxa represented after filtering, then counting them and listing their names
AbsRfGenusOPSPlaceboFirstPre4MoPs %<>% prune_taxa(taxa_sums(.) > 0, .)
OPSPlaceboFirst4MoGenusNumber <- AbsRfGenusOPSPlaceboFirstPre4MoPs %>% ntaxa()
OPSPlaceboFirst4MoGenusNames <- AbsRfGenusOPSPlaceboFirstPre4MoPs %>% taxa_names()
# 3C) Inspecting taxa present at PreAZT_1mo in OPS samples of patients receiving placebo first
AbsRfGenusOPSPlaceboFirstPre1MoPs <- subset_samples(AbsRfGenusOPSPlaceboFirst4Mo1MoPs,
                                                    TTTPrePost7gr == "PreAZT_1mo")
# Keeping only the taxa represented after filtering, then counting them and listing their names
AbsRfGenusOPSPlaceboFirstPre1MoPs %<>% prune_taxa(taxa_sums(.) > 0, .)
OPSPlaceboFirst1MoGenusNumber <- AbsRfGenusOPSPlaceboFirstPre1MoPs %>% ntaxa()
OPSPlaceboFirst1MoGenusNames <- AbsRfGenusOPSPlaceboFirstPre1MoPs %>% taxa_names()

# Preparing a list of vectors with taxa names for visualisation by Venn diagram
OPSPlaceboFirst4Mo1MoGenusList <- list("PreTTT_4mo" = OPSPlaceboFirst4MoGenusNames,
                                       "PreTTT_1moEndAZT" = OPSPlaceboFirst1MoGenusNames)
# Calculating the coordinates required for plotting a Venn diagram
OPSPlaceboFirst4Mo1MoGenusVenn <- Venn(OPSPlaceboFirst4Mo1MoGenusList)
OPSPlaceboFirst4Mo1MoGenusVennData <- process_data(OPSPlaceboFirst4Mo1MoGenusVenn)
# Obtaining the Venn diagram (For Supplementary Figure E13, panel b)
OPSPlaceboFirst4Mo1MoGenusVennDiagram <- ggVennDiagram(OPSPlaceboFirst4Mo1MoGenusList,
                                                       set_color = "white", # Hiding original label
                                                       label =  "both", # Showing "count" and "percent"
                                                       label_size = 10,
                                                       label_alpha = .6,
                                                       edge_lty = "blank") +
  # to customize the label layer
  geom_sf_text(aes(label = name),
               fontface = "plain",
               size = 9,
               data = venn_setlabel(OPSPlaceboFirst4Mo1MoGenusVennData)) +
  # to customize the edge layer
  geom_sf(size = .1,
          color = "#f2f2f2",
          alpha = .6,
          data = venn_setedge(OPSPlaceboFirst4Mo1MoGenusVennData)) +
  scale_fill_gradient(low = "#f0f0f0", high = "#999999") +
  theme_void() +
  theme(legend.position = "none")

ggsave("OPSPlaceboFirst4Mo1MoGenusVennDiagram.pdf",
       width = 12, height = 10, dpi = 200, units = "cm", limitsize = FALSE)

# Genus-level analysis - Comparison of OPS and sputum samples ---- 
# 1A) Inspecting taxa present in OPS and sputum samples at 4 or 1 mo before treatment (patients starting with placebo)
# Patients with OPS and sputum samples available at 4 mo before treatment
Placebo4moPreAZTSputumSwabGenusPat <- meta(AbsRfPatientGenusPs) %>%
  filter(TTTPrePost7gr == "PreAZT_4mo") %>% 
  select(PatientID, 
         TTTPrePost7gr, 
         OriginType3) %>% 
  group_by(PatientID) %>%
  filter(n() == 2) %>% 
  distinct(PatientID) %>% 
  pull(PatientID) %>% 
  droplevels()

AbsRfPreAZT4moSputumSwabGenusPs <- subset_samples(AbsRfPatientGenusPs,
                                                  PatientID %in% Placebo4moPreAZTSputumSwabGenusPat &
                                                    TTTPrePost7gr == "PreAZT_4mo")

# Keeping only the taxa represented after filtering
AbsRfPreAZT4moSputumSwabGenusPs %<>% prune_taxa(taxa_sums(.) > 0, .)

# Patients with OPS and sputum samples available at 1 mo before treatment (end placebo)
Placebo1moPreAZTSputumSwabGenusPat <- meta(AbsRfPatientGenusPs) %>%
  filter(TTTPrePost7gr == "PreAZT_1mo") %>% 
  select(PatientID, TTTPrePost7gr, OriginType3) %>% 
  group_by(PatientID) %>%
  filter(n() == 2) %>% 
  distinct(PatientID) %>% 
  pull(PatientID) %>% 
  droplevels()

AbsRfPreAZT1moSputumSwabGenusPs <- subset_samples(AbsRfPatientGenusPs,
                                                  PatientID %in% Placebo1moPreAZTSputumSwabGenusPat &
                                                    TTTPrePost7gr == "PreAZT_1mo")

# Keeping only the taxa represented after filtering
AbsRfPreAZT1moSputumSwabGenusPs %<>% prune_taxa(taxa_sums(.) > 0, .)

# Patients with OPS and sputum samples available at 4 or 1 mo before treatment
AbsRfPreAZT4mo1moSputumSwabGenusPs <- merge_phyloseq(otu_table(AbsRfPreAZT4moSputumSwabGenusPs),
                                                     otu_table(AbsRfPreAZT1moSputumSwabGenusPs),
                                                     tax_table(AbsRfPreAZT4moSputumSwabGenusPs),
                                                     tax_table(AbsRfPreAZT1moSputumSwabGenusPs),
                                                     sample_data(AbsRfPreAZT4moSputumSwabGenusPs),
                                                     sample_data(AbsRfPreAZT1moSputumSwabGenusPs))

# 1B) Inspecting taxa present in OPS and sputum samples at 4 mo or 1 mo before treatment
SputumSwabPreAZT4mo1moGenusNumber <- AbsRfPreAZT4mo1moSputumSwabGenusPs %>% ntaxa()
SputumSwabPreAZT4mo1moGenusNames <- AbsRfPreAZT4mo1moSputumSwabGenusPs %>% taxa_names()
# 1C) Inspecting taxa present in OPS samples at 4 mo or 1 mo before treatment
AbsRfPreAZT4mo1moSwabGenusPs <- subset_samples(AbsRfPreAZT4mo1moSputumSwabGenusPs,
                                               OriginType3 == "Swab")
# Keeping only the taxa represented after filtering, then counting them and listing their names
AbsRfPreAZT4mo1moSwabGenusPs %<>% prune_taxa(taxa_sums(.) > 0, .)
SwabPreAZT4mo1moGenusNumber <- AbsRfPreAZT4mo1moSwabGenusPs %>% ntaxa()
SwabPreAZT4mo1moGenusNames <- AbsRfPreAZT4mo1moSwabGenusPs %>% taxa_names()
# 1D) Inspecting taxa present in sputum samples at 4 mo before or 1 mo before treatment
AbsRfPreAZT4mo1moSputumGenusPs <- subset_samples(AbsRfPreAZT4mo1moSputumSwabGenusPs,
                                                 OriginType3 == "Sputum")
# Keeping only the taxa represented after filtering, then counting them and listing their names
AbsRfPreAZT4mo1moSputumGenusPs %<>% prune_taxa(taxa_sums(.) > 0, .)
SputumPreAZT4mo1moGenusNumber <- AbsRfPreAZT4mo1moSputumGenusPs %>% ntaxa()
SputumPreAZT4mo1moGenusNames <- AbsRfPreAZT4mo1moSputumGenusPs %>% taxa_names()

# Preparing a list of vectors with taxa names for visualisation by Venn diagram
SputumSwabPreAZT4mo1moGenusList <- list(OPS = SwabPreAZT4mo1moGenusNames,
                                        Sputum = SputumPreAZT4mo1moGenusNames)
# Calculating the coordinates required for plotting
SputumSwabPreAZT4mo1moGenusVenn <- Venn(SputumSwabPreAZT4mo1moGenusList)
SputumSwabPreAZT4mo1moGenusVennData <- process_data(SputumSwabPreAZT4mo1moGenusVenn)
# Obtaining the Venn diagram
SputumSwabPreAZT4mo1moGenusVennDiagram <- ggVennDiagram(SputumSwabPreAZT4mo1moGenusList,
                                                        set_color = "white", # Hiding original label
                                                        label =  "both", # Showing "count" and "percent"
                                                        label_size = 10,
                                                        label_alpha = .6,
                                                        edge_lty = "blank") +
  # to customize the label layer
  geom_sf_text(aes(label = name),
               fontface = "plain",
               size = 9,
               data = venn_setlabel(SputumSwabPreAZT4mo1moGenusVennData)) +
  # to customize the edge layer
  geom_sf(size = .1,
          color = "#f2f2f2",
          alpha = .6,
          data = venn_setedge(SputumSwabPreAZT4mo1moGenusVennData)) +
  scale_fill_gradient(low = "#f0f0f0", high = "#999999") +
  theme_void() +
  theme(legend.position = "none")

ggsave("SputumSwabPreAZT4mo1moGenusVennDiagram.pdf",
       width = 12, height = 10, dpi = 200, units = "cm", limitsize = FALSE)

# 2A) Inspecting taxa present in OPS and sputum samples at start of treatment
# Patients with OPS and sputum samples available
StartAZTSputumSwabIncrStableARGenusPat <- meta(AbsRfPatientGenusPs) %>%
  filter(TTTPrePost7gr == "StartAZT") %>% 
  select(PatientID, 
         TTTPrePost5gr, 
         OriginType3) %>% 
  group_by(PatientID) %>%
  arrange(PatientID) %>% 
  filter(n() == 2) %>% 
  distinct(PatientID) %>% 
  pull(PatientID) %>% 
  droplevels()

StartAZTSputumSwabIncrStableARGenusPs <- subset_samples(AbsRfPatientGenusPs,
                                                        PatientID %in% StartAZTSputumSwabIncrStableARGenusPat &
                                                          TTTPrePost5gr == "StartAZT")

# Keeping only the taxa represented after filtering
StartAZTSputumSwabIncrStableARGenusPs %<>% prune_taxa(taxa_sums(.) > 0, .)

# Proportion of unique or shared taxa in OPS and sputum samples at the start of treatment
# 2B) Inspecting taxa present in OPS and sputum samples
StartAZTSputumSwabIncrStableGenusNumber <- StartAZTSputumSwabIncrStableARGenusPs %>% ntaxa()
StartAZTSputumSwabIncrStableGenusNames <- StartAZTSputumSwabIncrStableARGenusPs %>% taxa_names()
# 2C) Inspecting taxa present in OPS at the start of treatment
AbsRfStartAZTSwabIncrStableARGenusPs <- subset_samples(StartAZTSputumSwabIncrStableARGenusPs,
                                                       OriginType3 == "Swab")
# Keeping only the taxa represented after filtering, then counting them and listing their names
AbsRfStartAZTSwabIncrStableARGenusPs %<>% prune_taxa(taxa_sums(.) > 0, .)
StartAZTSwabIncrStableARGenusNumber <- AbsRfStartAZTSwabIncrStableARGenusPs %>% ntaxa()
StartAZTSwabIncrStableARGenusNames <- AbsRfStartAZTSwabIncrStableARGenusPs %>% taxa_names()
# 2D) Inspecting taxa present in sputum samples at the start of treatment
AbsRfStartAZTSputumIncrStableARGenusPs <- subset_samples(StartAZTSputumSwabIncrStableARGenusPs,
                                                         OriginType3 == "Sputum")
# Keeping only the taxa represented after filtering, then counting them and listing their names
AbsRfStartAZTSputumIncrStableARGenusPs %<>% prune_taxa(taxa_sums(.) > 0, .)
StartAZTSputumIncrStableARPGenusNumber <- AbsRfStartAZTSputumIncrStableARGenusPs %>% ntaxa()
StartAZTSputumIncrStableARPGenusNames <- AbsRfStartAZTSputumIncrStableARGenusPs %>% taxa_names()

# Preparing a list of vectors with taxa names for visualisation by Venn diagram
SputumSwabStartAZTIncrStableARGenusList <- list(OPS = StartAZTSwabIncrStableARGenusNames,
                                                Sputum = StartAZTSputumIncrStableARPGenusNames)
# Calculating the coordinates required for plotting
SputumSwabStartAZTIncrStableARGenusVenn <- Venn(SputumSwabStartAZTIncrStableARGenusList)
SputumSwabStartAZTIncrStableARGenusVennData <- process_data(SputumSwabStartAZTIncrStableARGenusVenn)
# Obtaining the Venn diagram
SputumSwabStartAZTIncrStableARGenusVennDiagram <- ggVennDiagram(SputumSwabStartAZTIncrStableARGenusList,
                                                                set_color = "white", # Hiding original label
                                                                label =  "both", # Showing "count" and "percent"
                                                                label_size = 10,
                                                                label_alpha = .6,
                                                                edge_lty = "blank") +
  # to customize the label layer
  geom_sf_text(aes(label = name),
               fontface = "plain",
               size = 9,
               data = venn_setlabel(SputumSwabStartAZTIncrStableARGenusVennData)) +
  # to customize the edge layer
  geom_sf(size = .1,
          color = "#EFCA81",
          alpha = .6,
          data = venn_setedge(SputumSwabStartAZTIncrStableARGenusVennData)) +
  scale_fill_gradient(low = "#F7DFAD", high = "#E69F00") +
  theme_void() +
  theme(legend.position = "none")

ggsave("SputumSwabStartAZTIncrStableARGenusVennDiagram.pdf",
       width = 12, height = 10, dpi = 200, units = "cm", limitsize = FALSE)

# 3A) Inspecting taxa present in OPS and sputum samples at the end of treatment
# Patients with OPS and sputum samples available
EndAZTSputumSwabIncrStableARGenusPat <- meta(AbsRfPatientGenusPs) %>%
  filter(TTTPrePost5gr == "EndAZT") %>% 
  select(PatientID, TTTPrePost5gr, OriginType3) %>% 
  group_by(PatientID) %>%
  arrange(PatientID) %>%
  filter(n() == 2) %>% 
  distinct(PatientID) %>% 
  pull(PatientID) %>% 
  droplevels()

EndAZTSputumSwabIncrStableARGenusPs <- subset_samples(AbsRfPatientGenusPs,
                                                      PatientID %in% EndAZTSputumSwabIncrStableARGenusPat &
                                                        TTTPrePost5gr == "EndAZT")

# Keeping only the taxa represented after filtering
EndAZTSputumSwabIncrStableARGenusPs %<>% prune_taxa(taxa_sums(.) > 0, .)

# Proportion of unique or shared taxa in OPS and sputum samples at the end of treatment
# 3B) Inspecting taxa present in OPS and sputum samples
EndAZTSputumSwabIncrStableGenusNumber <- EndAZTSputumSwabIncrStableARGenusPs %>% ntaxa()
EndAZTSputumSwabIncrStableGenusNames <- EndAZTSputumSwabIncrStableARGenusPs %>% taxa_names()
# 3C) Inspecting taxa present in OPS at the end of treatment
AbsRfEndAZTSwabIncrStableARGenusPs <- subset_samples(EndAZTSputumSwabIncrStableARGenusPs,
                                                     OriginType3 == "Swab")
# Keeping only the taxa represented after filtering, then counting them and listing their names
AbsRfEndAZTSwabIncrStableARGenusPs %<>% prune_taxa(taxa_sums(.) > 0, .)
EndAZTSwabIncrStableARGenusNumber <- AbsRfEndAZTSwabIncrStableARGenusPs %>% ntaxa()
EndAZTSwabIncrStableARGenusNames <- AbsRfEndAZTSwabIncrStableARGenusPs %>% taxa_names()
# 3D) Inspecting taxa present in sputum samples at the end of treatment
AbsRfEndAZTSputumIncrStableARGenusPs <- subset_samples(EndAZTSputumSwabIncrStableARGenusPs,
                                                       OriginType3 == "Sputum")
# Keeping only the taxa represented after filtering, then counting them and listing their names
AbsRfEndAZTSputumIncrStableARGenusPs %<>% prune_taxa(taxa_sums(.) > 0, .)
EndAZTSputumIncrStableARPGenusNumber <- AbsRfEndAZTSputumIncrStableARGenusPs %>% ntaxa()
EndAZTSputumIncrStableARPGenusNames <- AbsRfEndAZTSputumIncrStableARGenusPs %>% taxa_names()

# Preparing a list of vectors with taxa names for visualisation by Venn diagram
SputumSwabEndAZTIncrStableARGenusList <- list(OPS = EndAZTSwabIncrStableARGenusNames,
                                              Sputum = EndAZTSputumIncrStableARPGenusNames)
# Calculating the coordinates required for plotting
SputumSwabEndAZTIncrStableARGenusVenn <- Venn(SputumSwabEndAZTIncrStableARGenusList)
SputumSwabEndAZTIncrStableARGenusVennData <- process_data(SputumSwabEndAZTIncrStableARGenusVenn)
# Obtaining the Venn diagram
SputumSwabEndAZTIncrStableARGenusVennDiagram <- ggVennDiagram(SputumSwabEndAZTIncrStableARGenusList,
                                                              set_color = "white", # For hiding original label
                                                              label =  "both", # For showing "count" and "percent"
                                                              label_size = 10, # For region label size
                                                              label_alpha = .6, # For adding transparency to region labels
                                                              edge_lty = "blank") +
  # to customize the label layer
  geom_sf_text(aes(label = name),
               fontface = "plain",
               size = 9,
               data = venn_setlabel(SputumSwabEndAZTIncrStableARGenusVennData)) +
  # to customize the edge layer
  geom_sf(size = .1,
          color = "#a5d6f3",
          alpha = .6,
          data = venn_setedge(SputumSwabEndAZTIncrStableARGenusVennData)) +
  scale_fill_gradient(low = "#d2ebf9", high = "#34a4e5") +
  theme_void() +
  theme(legend.position = "none")

ggsave("SputumSwabEndAZTIncrStableARGenusVennDiagram.pdf",
       width = 12, height = 10, dpi = 200, units = "cm", limitsize = FALSE)

# 4A) Inspecting taxa present in OPS and sputum samples at 1 mo post-treatment
# Patients with OPS and sputum samples available
PostAZT1moSputumSwabIncrStableARGenusPat <- meta(AbsRfPatientGenusPs) %>%
  filter(TTTPrePost7gr == "PostAZT_1mo") %>% 
  select(PatientID, TTTPrePost5gr, OriginType3) %>% 
  group_by(PatientID) %>%
  arrange(PatientID) %>% 
  filter(n() == 2) %>% 
  distinct(PatientID) %>% 
  pull(PatientID) %>% 
  droplevels()

PostAZT1moSputumSwabIncrStableARGenusPs <- subset_samples(AbsRfPatientGenusPs,
                                                          PatientID %in% PostAZT1moSputumSwabIncrStableARGenusPat &
                                                            TTTPrePost5gr == "PostAZT_1mo")

# Keeping only the taxa represented after filtering
PostAZT1moSputumSwabIncrStableARGenusPs %<>% prune_taxa(taxa_sums(.) > 0, .)

# Proportion of unique or shared taxa in OPS and sputum samples at 1 mo post-treatment
# 4B) Inspecting taxa present in OPS and sputum samples
PostAZT1moSputumSwabIncrStableGenusNumber <- PostAZT1moSputumSwabIncrStableARGenusPs %>% ntaxa()
PostAZT1moSputumSwabIncrStableGenusNames <- PostAZT1moSputumSwabIncrStableARGenusPs %>% taxa_names()
# 4C) Inspecting taxa present in OPS at 1 mo post-treatment
AbsRfPostAZT1moSwabIncrStableARGenusPs <- subset_samples(PostAZT1moSputumSwabIncrStableARGenusPs,
                                                         OriginType3 == "Swab")
# Keeping only the taxa represented after filtering, then counting them and listing their names
AbsRfPostAZT1moSwabIncrStableARGenusPs %<>% prune_taxa(taxa_sums(.) > 0, .)
PostAZT1moSwabIncrStableARGenusNumber <- AbsRfPostAZT1moSwabIncrStableARGenusPs %>% ntaxa()
PostAZT1moSwabIncrStableARGenusNames <- AbsRfPostAZT1moSwabIncrStableARGenusPs %>% taxa_names()
# 4D) Inspecting taxa present in sputum samples at 1 mo post-treatment
AbsRfPostAZT1moSputumIncrStableARGenusPs <- subset_samples(PostAZT1moSputumSwabIncrStableARGenusPs,
                                                           OriginType3 == "Sputum")
# Keeping only the taxa represented after filtering, then counting them and listing their names
AbsRfPostAZT1moSputumIncrStableARGenusPs %<>% prune_taxa(taxa_sums(.) > 0, .)
PostAZT1moSputumIncrStableARPGenusNumber <- AbsRfPostAZT1moSputumIncrStableARGenusPs %>% ntaxa()
PostAZT1moSputumIncrStableARPGenusNames <- AbsRfPostAZT1moSputumIncrStableARGenusPs %>% taxa_names()

# Preparing a list of vectors with taxa names for visualisation by Venn diagram
SputumSwabPostAZT1moIncrStableARGenusList <- list(OPS = PostAZT1moSwabIncrStableARGenusNames,
                                                  Sputum = PostAZT1moSputumIncrStableARPGenusNames)

# Calculating the coordinates required for plotting
SputumSwabPostAZT1moIncrStableARGenusVenn <- Venn(SputumSwabPostAZT1moIncrStableARGenusList)
SputumSwabPostAZT1moIncrStableARGenusVennData <- process_data(SputumSwabPostAZT1moIncrStableARGenusVenn)
# Obtaining the Venn diagram
SputumSwabPostAZT1moIncrStableARGenusVennDiagram <- ggVennDiagram(SputumSwabPostAZT1moIncrStableARGenusList,
                                                                  set_color = "white", # Hiding original label
                                                                  label =  "both", # Showing "count" and "percent"
                                                                  label_size = 10,
                                                                  label_alpha = .6,
                                                                  edge_lty = "blank") +
  # to customize the label layer
  geom_sf_text(aes(label = name),
               fontface = "plain",
               size = 9,
               data = venn_setlabel(SputumSwabPostAZT1moIncrStableARGenusVennData)) +
  # to customize the edge layer
  geom_sf(size = .1,
          color = "#99ffe4",
          alpha = .6,
          data = venn_setedge(SputumSwabPostAZT1moIncrStableARGenusVennData)) +
  scale_fill_gradient(low = "#ccfff1", high = "#00cc96") +
  theme_void() +
  theme(legend.position = "none")

ggsave("SputumSwabPostAZT1moIncrStableARGenusVennDiagram.pdf",
       width = 12, height = 10, dpi = 200, units = "cm", limitsize = FALSE)

# 5A) Inspecting taxa present in OPS and sputum samples at 4 mo post-treatment
# Patients with OPS and sputum samples available
PostAZT4moSputumSwabIncrStableARGenusPat <- meta(AbsRfPatientGenusPs) %>%
  filter(TTTPrePost7gr == "PostAZT_4mo") %>%
  group_by(PatientID) %>%
  arrange(PatientID) %>% 
  filter(n() == 2) %>% 
  distinct(PatientID) %>% 
  pull(PatientID) %>% 
  droplevels()

PostAZT4moSputumSwabIncrStableARGenusPs <- subset_samples(AbsRfPatientGenusPs,
                                                          PatientID %in% PostAZT4moSputumSwabIncrStableARGenusPat &
                                                            TTTPrePost7gr == "PostAZT_4mo")

# Keeping only the taxa represented after filtering
PostAZT4moSputumSwabIncrStableARGenusPs %<>% prune_taxa(taxa_sums(.) > 0, .)

# Proportion of unique or shared taxa in OPS and sputum samples at 4 mo post-treatment
# 5B) Inspecting taxa present in OPS and sputum samples
PostAZT4moSputumSwabIncrStableGenusNumber <- PostAZT4moSputumSwabIncrStableARGenusPs %>% ntaxa()
PostAZT4moSputumSwabIncrStableGenusNames <- PostAZT4moSputumSwabIncrStableARGenusPs %>% taxa_names()
# 5C) Inspecting taxa present in OPS at 4 mo post-treatment
AbsRfPostAZT4moSwabIncrStableARGenusPs <- subset_samples(PostAZT4moSputumSwabIncrStableARGenusPs,
                                                         OriginType3 == "Swab")
# Keeping only the taxa represented after filtering, then counting them and listing their names
AbsRfPostAZT4moSwabIncrStableARGenusPs %<>% prune_taxa(taxa_sums(.) > 0, .)
PostAZT4moSwabIncrStableARGenusNumber <- AbsRfPostAZT4moSwabIncrStableARGenusPs %>% ntaxa()
PostAZT4moSwabIncrStableARGenusNames <- AbsRfPostAZT4moSwabIncrStableARGenusPs %>% taxa_names()
# 5D) Inspecting taxa present in sputum samples at 4 mo post-treatment
AbsRfPostAZT4moSputumIncrStableARGenusPs <- subset_samples(PostAZT4moSputumSwabIncrStableARGenusPs,
                                                           OriginType3 == "Sputum")
# Keeping only the taxa represented after filtering, then counting them and listing their names
AbsRfPostAZT4moSputumIncrStableARGenusPs %<>% prune_taxa(taxa_sums(.) > 0, .)
PostAZT4moSputumIncrStableARGenusPNumber <- AbsRfPostAZT4moSputumIncrStableARGenusPs %>% ntaxa()
PostAZT4moSputumIncrStableARPGenusNames <- AbsRfPostAZT4moSputumIncrStableARGenusPs %>% taxa_names()

# Preparing a list of vectors with taxa names for visualisation by Venn diagram
SputumSwabPostAZT4moIncrStableARGenusList <- list(OPS = PostAZT4moSwabIncrStableARGenusNames,
                                                  Sputum = PostAZT4moSputumIncrStableARPGenusNames)

# Calculating the coordinates required for plotting
SputumSwabPostAZT4moIncrStableARGenusVenn <- Venn(SputumSwabPostAZT4moIncrStableARGenusList)
SputumSwabPostAZT4moIncrStableARGenusVennData <- process_data(SputumSwabPostAZT4moIncrStableARGenusVenn)
# Obtaining the Venn diagram
SputumSwabPostAZT4moIncrStableARGenusVennDiagram <- ggVennDiagram(SputumSwabPostAZT4moIncrStableARGenusList,
                                                                  set_color = "white", # Hiding original label
                                                                  label =  "both", # Showing "count" and "percent"
                                                                  label_size = 10,
                                                                  label_alpha = .6,
                                                                  edge_lty = "blank") +
  # to customize the label layer
  geom_sf_text(aes(label = name),
               fontface = "plain",
               size = 9,
               data = venn_setlabel(SputumSwabPostAZT4moIncrStableARGenusVennData)) +
  # to customize the edge layer
  geom_sf(size = .1,
          color = "#eac8da",
          alpha = .6,
          data = venn_setedge(SputumSwabPostAZT4moIncrStableARGenusVennData)) +
  scale_fill_gradient(low = "#f1dae7", high = "#cc79a7") +
  theme_void() +
  theme(legend.position = "none")

ggsave("SputumSwabPostAZT3moIncrStableARGenusVennDiagram.pdf",
       width = 12, height = 10, dpi = 200, units = "cm", limitsize = FALSE)

# Genus-level analysis - Assessing in LRT the dynamics of bacteria cleared from LRT during treatment ----
# Sputum: Names of all ASVs present at StartAZT but not EndAZT
SputumStartOnlyGenusNames <- setdiff(SputumStartGenusNames, SputumEndGenusNames)

# Names of patients with a pair of sputum samples available at StartAZT and EndAZT
# and selection of all samples obtained from these patients
SputumStartEndPatientNames <- meta(AbsRfSputumStartAZTPs)$PatientID

AbsRfSputumPatSputumStartEndAllGenusPs <- subset_samples(AbsRfPatientGenusPs,
                                                         OriginType == "Patient_Sputum" &
                                                           PatientID %in% SputumStartEndPatientNames)
# Removing genera with zero count
AbsRfSputumPatSputumStartEndAllGenusPs2 <- prune_taxa(taxa_sums(AbsRfSputumPatSputumStartEndAllGenusPs) > 0, 
                                                      AbsRfSputumPatSputumStartEndAllGenusPs)
# Merging taxa cleared during treatment
AbsRfSputumAllGenusMergeClearedGenusPs2 <- microbiome::merge_taxa2(AbsRfSputumPatSputumStartEndAllGenusPs2, 
                                                                   taxa = SputumStartOnlyGenusNames,
                                                                   name = "Cleared during treatment")

# Working with relative abundance
RelRfSputumAllGenusMergeClearedGenusPs2 <- microbiome::transform(AbsRfSputumAllGenusMergeClearedGenusPs2, 
                                                                 "compositional")
# Data frame for subsequent plotting
RelRfSputumAllGenusMergeClearedGenusDF <- psmelt(RelRfSputumAllGenusMergeClearedGenusPs2)
# Keeping only taxa removed during treatment
RelRfSputumAllGenusMergeClearedGenusDF2 <- RelRfSputumAllGenusMergeClearedGenusDF %>%
  dplyr::rename(TaxaGroup = OTU, RelAbundance = Abundance) %>% 
  filter(TaxaGroup == "Cleared during treatment") %>% 
  select(SampleID, PatientID, TTTPrePost5gr, TTTPrePost7gr, RelAbundance, TaxaGroup) %>% 
  convert(fct(SampleID, PatientID, TTTPrePost5gr, TTTPrePost7gr, TaxaGroup))

SputumMergeClearedGenus_boxplot <- RelRfSputumAllGenusMergeClearedGenusDF2 %>%
  ggplot(aes(x = TTTPrePost7gr,
             y = RelAbundance,
             fill = TTTPrePost7gr)) +
  geom_boxplot(width = .5,
               size = .3,
               alpha = .4,
               position = position_dodge(.7),
               outlier.shape = NA,
               show.legend = FALSE) +
  geom_jitter(shape = 21,
              colour = "black",
              size = 1.5,
              alpha = .6,
              width = .2) +
  stat_compare_means(comparisons = CompTTT7) +
  stat_compare_means(label.y = .05) +
  ylim(0, .05) +
  scale_fill_manual(values = TTTphaseCol7) +
  labs(x = "", y = "Relative abundance") +
  theme_classic() +
  theme(strip.background = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1))

ggsave("SputumMergeClearedGenus_boxplot.pdf",
       width = 10, height = 12, dpi = 200, units = "cm")

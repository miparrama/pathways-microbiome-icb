### Script to correlate TopoScore Signature and the functional dysbiosis. 
# Figure S2F

# Setup -------------------------------------------------------------------

knitr::opts_chunk$set(warning = FALSE, message = FALSE) 
source('scripts/helpers.R')
library (janitor)
library (readxl)

# Load metadata:
drup.metadata = read_csv ("output/data/drup-gut-sampledata.csv") %>% 
  as.data.frame()

rownames (drup.metadata) <- drup.metadata$Sample

metaphlan.files = list.files("../workflows/biobakery/results-vJan21/metaphlan/read_stats/", full.names = T)
names(metaphlan.files) =  list.files("../workflows/biobakery/results-vJan21/metaphlan/read_stats/") %>% 
  basename() %>% 
  str_remove("-reads.txt")

metaphlan.files = metaphlan.files[names(metaphlan.files) != "NA"]

pseq.reads = make_phyloseq4(table_fn = metaphlan.files, 
                            sampledata = drup.metadata) %>% 
  subset_samples(Response != "MISSING") %>% 
  subset_samples(Response != "FOLLOWS") %>% 
  subset_samples(Response != "NE") %>% 
  subset_samples(msStatus != "MSS") %>% 
  prune_samples(sample_sums(.)>0, .) %>% 
  prune_taxa(taxa_sums(.)>0, .) %>% 
  subset_samples(patientId != "DRUP01290005") %>% 
  subset_taxa(Kingdom == "Bacteria") 

tse_reads = pseq.reads %>% 
  convertFromPhyloseq() %>% 
  mia::transformAssay(assay.type = "counts",  method = "relabundance") %>% 
  mia::transformAssay( assay.type = "relabundance", method = "clr", pseudocount = T)

species_score <- read_xlsx("resources/toposcore_species.xlsx",
                           skip = 1,sheet = 2) %>% 
  mutate (MGS = str_replace_all(MGS, " ", ""))

rowData(tse_reads)$MGS <- rowData(tse_reads)$Species %>% str_replace_all("_", "")
rowData(tse_reads)$SIG <- rowtibble(tse_reads) %>%
  left_join(species_score) %>% 
  pull(SIG)
rowtibble(tse_reads) %>%
  filter (!is.na (SIG)) %>% 
  distinct(Species)

species_score %>% filter (!MGS %in% rowData(tse_reads)$MGS ) #Only Akkermansia not found. 
species_score %>% filter (MGS %in% rowData(tse_reads)$MGS )

# Function to compute S-score
compute_S_score <- function(tse,
                            species_df,
                            assay_name = "counts",
                            agglomerate_rank = "Species") {
  
  #Akkermansia_muciniphila strain to classify gray zone samples. 
  fuso = tse[rowData(tse)$SGB == "SGB9226",] %>% 
    meltSE()
  
  if (dim(fuso)[1] == 0) {
    fuso = tse[rowData(tse)$Species == "Akkermansia_muciniphila",] %>% 
      meltSE()
  }
  
  tse = agglomerateByRank(tse, rank= agglomerate_rank)
  
  # Subset mgs by signature
  sig1_mgs <- species_df %>% filter(SIG == "SIG1") %>% pull(MGS)
  sig2_mgs <- species_df %>% filter(SIG == "SIG2") %>% pull(MGS)
  
  #Extract species:
  sig1_species <- rowData(tse) %>%
    as.data.frame() %>%
    filter (MGS %in% sig1_mgs) %>% 
    pull (agglomerate_rank)
  
  sig2_species <- rowData(tse) %>%
    as.data.frame() %>%
    filter (MGS %in% sig2_mgs) %>% 
    pull (agglomerate_rank)
  
  
  message(paste0(length(sig1_species), "/37 species of SIG1 found"))
  message(paste0(length(sig2_species), "/45 species of SIG2 found"))
  
  # Extract assay data
  abundance_data <- assay(tse, assay_name)
  
  # Subset abundance data for SIG1 and SIG2
  sig1_data <- abundance_data[rownames(abundance_data) %in% sig1_species, , drop = FALSE]
  sig2_data <- abundance_data[rownames(abundance_data) %in% sig2_species, , drop = FALSE]
  
  # Compute proportions of present species
  sig1_present <- colSums(sig1_data > 0) / 37 # SIG1 proportion of present
  sig2_present <- colSums(sig2_data > 0) / 45 # SIG2 proportion of present
  
  # Calculate S-score
  S <- ((sig2_present - sig1_present) + 1) / 2
  
  # Compute the class based on fuso abundance. 
  fuso$S = S
  fuso = fuso %>% 
    mutate( Sig = case_when (
      S > 0.7911 ~ "Sig2", 
      S < 0.5351 ~ "Sig1",
      .default = "gray_zone")) %>% 
    mutate(Sig = case_when (
      (Sig == "gray_zone") & (counts == 0) ~ "Sig1", 
      (Sig == "gray_zone") & (counts > 15) ~ "Sig1",
      (Sig == "gray_zone") & (counts < 15) & (counts > 0 ) ~ "Sig2",
      .default = Sig)) %>% 
    column_to_rownames("SampleID")
  
  sig = fuso$Sig
  names(sig) <- rownames (fuso)
  
  result = list ("s.score" = S, "toposcore" = sig)
  
  return(result)
}

# Correlation -------------------------------------------------------------

# Compute s-score. 
colData(tse_reads)$s_score = compute_S_score(tse_reads,
                                             species_score)$s.score
colData(tse_reads)$toposcore = compute_S_score(tse_reads,
                                               species_score)$toposcore

drup.dysbiosis <- read_csv("output/data/drup-functional-dysbiosis.csv")

S.score_dysbiosis_plot <- drup.dysbiosis %>% 
  select (Sample, dysbiosis) %>% 
  inner_join(coltibble(tse_reads)) %>% 
  ggplot (aes (dysbiosis, s_score)) +
  geom_smooth(method = "lm", color = "black", size = 2) +
  geom_point(size = 3, alpha = 0.5) +
  stat_cor(label.x.npc = "center", size = 5) +
  theme_bw(base_size = 14) + 
  labs (x = "Functional dysbiosis", 
        y = "S score (Derosa et al, 2024)") +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())

S.score_dysbiosis_plot

# Save plot: 
ggsave ("output/figures/toposcore_dysbiosis.pdf",S.score_dysbiosis_plot, height = 4, width = 5.5 )
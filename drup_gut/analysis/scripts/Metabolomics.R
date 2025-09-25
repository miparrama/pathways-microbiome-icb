# Metabolomics analysis. 
# Figure 2E

source ("scripts/helpers.R") 
library (janitor)

# Setup -------------------------------------------------------------------

tse.metabo.fil <- read_rds ("resources/metabolomics.se.fil.rds")

drup.functional.dysbiosis <- read_csv ("output/data/drup-functional-dysbiosis.csv")

# Filter data
colData(tse.metabo.fil) <- tse.metabo.fil %>% 
  coltibble() %>% 
  mutate (dysbiosis_level = ifelse (dysbiosis > 0, "High", "Low")) %>% 
  DataFrame()

tse.metabo.dysbiosis <- tse.metabo.fil[,colData(tse.metabo.fil)$dysbiosis_level %in% c("High", "Low")]

rowData(tse.metabo.dysbiosis)$prevalence <- rowMeans(assay(tse.metabo.dysbiosis) > rowData(tse.metabo.dysbiosis)$lod)

# Filter based on prevalence. 
tse.metabo.dysbiosis <- tse.metabo.dysbiosis[ rowData(tse.metabo.dysbiosis)$prevalence > 0.7, ]

colnames (tse.metabo.dysbiosis) <- colData(tse.metabo.dysbiosis)$Sample

# Filter uncharacterized compounds:
tse.metabo.dysbiosis <- tse.metabo.dysbiosis[rowData(tse.metabo.dysbiosis)$annotation_level %in% c("1", "2a", "2b", "3"),]

# Filter based on drugs. 
# Filter drugs: 
synthetic_drugs = c("Atenolol", "Saccharin", "4-Acetamidophenol", "Acesulfame", "Atenolol acid",
                    "Caprolactam", "Cyclamic acid", "DEET", "Dimetridazole", "Edoxaban",
                    "Famotidine", "iopromide", "Irbesartan", "Levetiracetam", "Losartan",
                    "Methoxyphenamine", "Metoprolol", "Oxycodone", "Pilocarpine", "Pregabalin",
                    "Quinine", "Telmisartan", "Trimethoprim", "Valsartan", "Venlafaxine", "Allopurinol")
# Filter PEG wich is usually a contamination source. 
tse.metabo.dysbiosis <- tse.metabo.dysbiosis[!(rowData(tse.metabo.dysbiosis)$name %in% synthetic_drugs) &
                                               !str_detect(rowData(tse.metabo.dysbiosis)$name, "PEG"),] 

# PCA metabolomics --------------------------------------------------------

core.data_mt <- tse.metabo.dysbiosis %>% 
  assay("log10") %>%
  t()

data_sd <- apply(core.data_mt,1,'sd')
core.data_mt <- core.data_mt[data_sd!=0,]

pca <- prcomp (core.data_mt, scale = TRUE, center = TRUE)
sum_pca <- summary(pca)

pca_df <-  pca$x %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "Sample") %>%
  left_join(., colData(tse.metabo.dysbiosis) %>% 
              as.data.frame(),
            by = "Sample")

dysbiosis_continious_pca <- pca_df %>% 
  mutate (dysbiosis_level = as.factor(dysbiosis_level)) %>% 
  mutate (dysbiosis_quantile = case_when(
    dysbiosis < quantile (dysbiosis, 1/3) ~ "Low", 
    dysbiosis >= quantile (dysbiosis, 2/3) ~ "High", 
    .default = "Medium"
  )) %>% 
  ggplot( aes (x = PC1, y = PC2)) +
  geom_point(size = 3.5, aes ( shape = time)) + 
  geom_point(size = 3, aes (color = dysbiosis, shape = time)) + 
  xlab(sprintf("PC1: (%.2f %%)",sum_pca$importance[2,1]*100 ))+
  ylab(sprintf("PC2: (%.2f %%)",sum_pca$importance[2,2]*100 ))+
  stat_ellipse(data = . %>% filter (dysbiosis_quantile == "High"),
               aes(group = dysbiosis_level), color = "red", type = "t", level = 0.95, linetype = "dashed", alpha = 0.65) +
  stat_ellipse(data = . %>% filter (dysbiosis_quantile == "Low"),
               aes(group = dysbiosis_level), color = "blue", type = "t", level = 0.95, linetype = "dashed", alpha = 0.65) +
  labs (title = "Untargeted Metabolomics", 
        color = "Functional dysbiosis", 
        shape = "Timepoint")+
  scale_colour_gradient2(high = "red",
                         mid = "white",
                         low = "blue") +
  annotate("text", x = 20, y = -25, size = 5, label = "Permanova p<0.0001")+
  theme_classic( base_size = 14 ) +
  theme(legend.position = "bottom")

dysbiosis_continious_pca

# Save plot: 
ggsave ("output/figures/metabolomics_dysbiosis_pca.pdf",dysbiosis_continious_pca, height = 6, width = 7 )


# Permanova ----------------------------------------------------------------

library (vegan)

otu_table <- assays(tse.metabo.dysbiosis)$log10  
otu_table <- t(otu_table) 

dist_matrix <- vegdist(otu_table, method = "euclidean")

set.seed(5)
adonis_result <- adonis2(dist_matrix ~ dysbiosis + ATB_use + tumor_class + Response, data = coltibble(tse.metabo.dysbiosis), permutations = 1000)
adonis_result



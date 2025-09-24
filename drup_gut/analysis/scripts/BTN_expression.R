# BTN expression in MSI cohort. 
# Figure 3H

source('scripts/helpers.R')
library (janitor)
library (patchwork)
library (ggpubr)
library (patchwork)
library (broom)
library (limma)
library (edgeR)

# Load metadata
msi.metadata = read_csv("resources/dMMR-MSIH-ICB-all.csv")
feature.counts = data.table::fread("resources/feature-counts-symbol.csv.gz")
clinical_data <- read_csv ("resources/clinical_data_elaborate_20250314.csv") %>% 
  mutate (patientId = str_replace_all(SubjectKey, "-", ""))

metadata_msi_rna <- msi.metadata %>% 
  filter (msStatus == "MSI") %>% 
  filter (sampleId %in% colnames (feature.counts)) %>% 
  select (patientId, cohort, sampleId, purity, tumor_groups, lymphnode) %>% 
  dplyr::left_join(clinical_data %>% select (patientId, CB, BOR_proxy)) %>% 
  filter (BOR_proxy != "NE") %>% 
  mutate (Response = ifelse (BOR_proxy %in% c("CR", "PR"), "Response", "Non-response"), 
          CB = ifelse (CB == 1, "CB", "NCB"))

fc.counts.symbol <- feature.counts %>% 
  column_to_rownames(var = "symbol")

fc.counts.symbol <- fc.counts.symbol[,metadata_msi_rna$sampleId]

d <- DGEList(fc.counts.symbol)
d <- calcNormFactors(d)

#Filter low expressed genes:
keep <- filterByExpr(d,
                     design = model.matrix( ~  CB, data = metadata_msi_rna),
                     min.count = 10, min.total.count = 15, large.n = 10, min.prop = 0.7
)
d <- d[keep,]
mm <- model.matrix( ~ 0 + CB + tumor_groups + lymphnode + purity , data =  metadata_msi_rna %>% mutate (tumor_groups = make.names(tumor_groups)))

colnames(mm)[1:2] <-c("CB", "NCB")

y <- voom(d, mm, plot = T)

#Fit of the model
fit <- lmFit(y, mm)

contr <- makeContrasts( 
  "Response" = CB - NCB,
  levels = mm
)

tmp <- contrasts.fit(fit, contr)
tmp <- eBayes(tmp)

top_table <- topTable(tmp, coef = "Response", sort.by = "P", n = Inf, adjust.method = "fdr" )

genes <- c("BTN2A1", "BTN3A1", "BTN3A2", "BTN3A3")

dea <- topTable(tmp, coef = "Response", sort.by = "P", n = Inf, adjust.method = "fdr" )[genes,] %>% 
  rownames_to_column(var = "symbol")

expression_data = y$E %>% 
  as.data.frame() %>% 
  rownames_to_column(var= "symbol") %>% 
  filter (symbol %in% genes) %>% 
  pivot_longer(., -symbol, values_to = "log2cpm", names_to = "sampleId") %>% 
  dplyr::left_join(., metadata_msi_rna) %>% 
  left_join(dea %>% select (symbol, P.Value))

btn_cb_boxplot <- expression_data %>% 
  ggplot (aes (x = symbol, y = log2cpm)) +
  geom_boxplot(aes (fill = CB),outliers = F) +
  geom_text(data = . %>% distinct(symbol, .keep_all = T),
            aes (label = sprintf("p=%.3f",  P.Value),
                 y = 9),
            size = 4.5) +
  geom_quasirandom(aes (group = CB), dodge.width=0.75, alpha = 0.5, size = 1) +
  theme_classic(base_size = 14) +
  scale_fill_manual(values = c("NCB" = "darkred", "CB" = "darkgreen")) +
  labs ( x = NULL, 
         y = "Gene expression (log2cpm)", 
         fill = NULL) +
  theme(axis.text.x = element_text(color = "black", size = 12))

btn_cb_boxplot
ggsave("output/figures/msi_btn_cb_boxplot.pdf", btn_cb_boxplot, height = 3.5, width = 6.5)

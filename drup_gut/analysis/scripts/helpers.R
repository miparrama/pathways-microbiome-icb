# Helper function for metagenomics analysis

# Load all libraries
library(ggsci)
library(remotes)
library(DivNet)
library(RColorBrewer)
library(lubridate)

library(ggthemes)
library(viridis)
library(scater)
library(MMUPHin)
library(mia)
library(miaViz)
library(tidyHeatmap)
library(ggsci)

library (broom)
library(ggbeeswarm)
library(ggpubr)
library(patchwork)
library(MicrobiomeStat)
library(vegan)
library(mixOmics)
library(ggrepel)
library(ggpubr)
library(NetCoMi)
library(phyloseq)
library (janitor)
library(microbiome)
library(tidyverse)


# Import metacyc Names
metacyc.names = vroom::vroom("../../resources/map_metacyc-pwy_name.txt.gz", col_names = c("ID", "Name"))

# Make phyloseq object from metaphlan output
make_phyloseq = function(table_fn, sampledata) {
  
  # Import metaphlan results (rel w readstats)
  metaphlan = table_fn %>% 
    map_dfr(.id = "sampleId", .f = function(x) {suppressMessages(read_tsv(x, skip = 4, col_types = cols(
      `#clade_name` = col_character(),
      clade_taxid = col_character(),
      relative_abundance = col_double(),
      coverage = col_double(),
      estimated_number_of_reads_from_the_clade = col_double()
    )))}) %>% 
    rename(clade_name = `#clade_name`,
           reads = estimated_number_of_reads_from_the_clade) %>% 
    filter(grepl("s__", clade_name)) %>% 
    select(-clade_taxid, -coverage)  %>% 
    separate(clade_name, into = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"), sep = "\\|") %>% 
    mutate(Kingdom = str_remove(Kingdom, "k__"),
           Phylum = str_remove(Phylum, "p__"),
           Class = str_remove(Class, "c__"),
           Order = str_remove(Order, "o__"),
           Family = str_remove(Family, "f__"),
           Genus = str_remove(Genus, "g__"), 
           Species = str_remove(Species, "s__")) 
  
  # Get OTU info (use mp2 read stats)
  otu.table.stats = metaphlan %>% 
    select(-relative_abundance) %>% 
    pivot_wider(names_from = sampleId, values_from = reads, values_fill = 0) %>% 
    select(-Kingdom:-Genus) %>% 
    column_to_rownames("Species") %>% 
    otu_table(taxa_are_rows = T)
  
  # Get taxa info
  taxa.table = metaphlan %>% 
    distinct(Species, .keep_all = T) %>% 
    select(Kingdom:Species) %>% 
    mutate(sp = Species) %>% 
    column_to_rownames("sp") %>% 
    mutate(Phylum = if_else(Phylum == "Eukaryota_unclassified", "Euryarchaeota", Phylum)) %>%
    mutate(Phylum = if_else(Phylum == "Viruses_unclassified", "Viruses", Phylum)) %>% 
    as.matrix() %>% 
    tax_table()

  # Get sample data
  sampledata = sampledata %>% 
    rownames_to_column("SAMPLE") %>% 
    filter(SAMPLE %in% colnames(otu.table.stats)) %>% 
    column_to_rownames("SAMPLE") %>% 
    sample_data()
  
  tree.data = ape::read.tree("../older-runs/biobakery_3_1_1/metaphlan/unifrac/mpa_v30_CHOCOPhlAn_201901_species_tree.nwk")
  tree.data$tip.label <- gsub(".+\\|s__", "", tree.data$tip.label)
  
  
  # Make phyloseq object (w/ estimated read counts)
  pseq.stats = phyloseq(otu.table.stats, taxa.table, sampledata, tree.data)
  return(pseq.stats)
  
}

make_phyloseq4 = function(table_fn, sampledata, tree = F, relab = F) {
  
  # Import metaphlan results (rel w readstats)
  if(relab){
    metaphlan = table_fn %>% 
      map_dfr(.id = "sampleId", .f = function(x) {suppressMessages(data.table::fread(x, skip = 4))}) %>% 
      rename(clade_name = `#clade_name`,
             reads = relative_abundance) %>% 
      filter(grepl("t__", clade_name)) %>% 
      select(-additional_species, -NCBI_tax_id)  %>% 
      separate(clade_name, into = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species", "SGB"), sep = "\\|") %>% 
      mutate(Kingdom = str_remove(Kingdom, "k__"),
             Phylum = str_remove(Phylum, "p__"),
             Class = str_remove(Class, "c__"),
             Order = str_remove(Order, "o__"),
             Family = str_remove(Family, "f__"),
             Genus = str_remove(Genus, "g__"), 
             Species = str_remove(Species, "s__"),
             SGB = str_remove(SGB, "t__")) %>% 
      mutate(SGB = str_remove(SGB, "_group"))
    
    # Get OTU info (use mp2 read stats)
    otu.table.stats = metaphlan %>% 
      pivot_wider(names_from = sampleId, values_from = reads, values_fill = 0) %>% 
      select(-Kingdom:-Species) %>% 
      column_to_rownames("SGB") %>% 
      otu_table(taxa_are_rows = T)
    
  } else{
    metaphlan = table_fn %>% 
      map_dfr(.id = "sampleId", .f = function(x) {suppressMessages(data.table::fread(x, skip = 5))}) %>% 
      rename(clade_name = `#clade_name`,
             reads = estimated_number_of_reads_from_the_clade) %>% 
      filter(grepl("t__", clade_name)) %>% 
      select(-coverage, -clade_taxid)  %>% 
      separate(clade_name, into = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species", "SGB"), sep = "\\|") %>% 
      mutate(Kingdom = str_remove(Kingdom, "k__"),
             Phylum = str_remove(Phylum, "p__"),
             Class = str_remove(Class, "c__"),
             Order = str_remove(Order, "o__"),
             Family = str_remove(Family, "f__"),
             Genus = str_remove(Genus, "g__"), 
             Species = str_remove(Species, "s__"),
             SGB = str_remove(SGB, "t__")) %>% 
      mutate(SGB = str_remove(SGB, "_group"))
    
    # Get OTU info (use mp2 read stats)
    otu.table.stats = metaphlan %>% 
      select(-relative_abundance) %>% 
      pivot_wider(names_from = sampleId, values_from = reads, values_fill = 0) %>% 
      select(-Kingdom:-Species) %>% 
      column_to_rownames("SGB") %>% 
      otu_table(taxa_are_rows = T)
  }
  
  # Get taxa info
  taxa.table = metaphlan %>% 
    distinct(SGB, .keep_all = T) %>% 
    select(Kingdom:SGB) %>% 
    mutate(sp = SGB) %>% 
    column_to_rownames("sp") %>% 
    mutate(Phylum = if_else(Phylum == "Eukaryota_unclassified", "Euryarchaeota", Phylum)) %>%
    mutate(Phylum = if_else(Phylum == "Viruses_unclassified", "Viruses", Phylum)) %>% 
    as.matrix() %>% 
    tax_table()
  
  # Get sample data
  sampledata = sampledata %>% 
    rownames_to_column("SAMPLE") %>% 
    filter(SAMPLE %in% colnames(otu.table.stats)) %>% 
    column_to_rownames("SAMPLE") %>% 
    sample_data()

  # Make phyloseq object (w/ estimated read counts)
  pseq.stats = phyloseq(otu.table.stats, taxa.table, sampledata)
  
  # Get phylogenetic tree
  if(tree){
    mpa_tree = ape::read.tree("http://cmprod1.cibio.unitn.it/biobakery4/metaphlan_databases/mpa_vOct22_CHOCOPhlAnSGB_202212.nwk")
    mpa_tree$tip.label = paste0("SGB", mpa_tree$tip.label)
    pseq.tree = merge_phyloseq(pseq.stats, mpa_tree)
    return(pseq.tree)
  } else{
    pseq.stats
  }
  
}

# Compute GMHI
# https://github.com/jaeyunsung/GMHI_2020/blob/master/GMHI.R
mh_species = read_tsv("https://raw.githubusercontent.com/jaeyunsung/GMHI_2020/master/MH_species.txt", col_names = "Species") %>% 
  mutate(Species = str_remove(Species, "s__"))

mn_species = read_tsv("https://raw.githubusercontent.com/jaeyunsung/GMHI_2020/master/MN_species.txt", col_names = "Species") %>% 
  mutate(Species = str_remove(Species, "s__"))

compute_ghmi = function(pseq){
  
  alpha <- function(x){sum((log(x[x>0]))*(x[x>0]))*(-1)}
  
  #tax.names = taxa_names(pseq) 
  #tax.names.keep = tax.names[!grepl("unclassified", tax.names, fixed = T)]
  
  species_profile = pseq %>% 
    #subset_taxa(Species %in% tax.names.keep) %>% 
    subset_taxa(Kingdom != "Viruses") %>% 
    transform(transform = "compositional")
  
  # Extracting Health-prevalent species present in metagenome
  # Extracting Health-scarce species present in metagenome
  MH_species_metagenome <- species_profile %>% 
    subset_taxa(Species %in% mh_species$Species) %>% 
    abundances()
  
  MN_species_metagenome <- species_profile %>% 
    subset_taxa(Species %in% mn_species$Species) %>% 
    abundances()
  
  # Diversity among Health-prevalent species
  # Diversity among Health-scarce species
  MH_shannon <- apply((MH_species_metagenome), 2, alpha) 
  MN_shannon <- apply((MN_species_metagenome), 2, alpha) 
  
  # Richness of Health-prevalent species
  # Richness of Health-scarce species
  R_MH <- apply(MH_species_metagenome, 2, function(i) (sum(i > 0))) 
  R_MN <- apply(MN_species_metagenome, 2, function(i) (sum(i > 0)))
  
  # Median RMH from 1% of the top-ranked samples (see Methods)
  # Median RMN from 1% of the bottom-ranked samples (see Methods)
  MH_prime <- 7
  MN_prime <- 31
  
  # Collective abundance of Health-prevalent species
  # Collective abundance of Health-scarce species
  psi_MH <- ((R_MH/MH_prime)*MH_shannon) 
  psi_MN <- ((R_MN/MN_prime)*MN_shannon)
  
  GMHI <- data.frame(log10((psi_MH+0.00001)/(psi_MN+0.00001))) # 0.00001 added to avoid having the denominator as 0
  colnames(GMHI) <- c("GMHI")
  GMHI.res = GMHI %>% 
    as.data.frame() %>% 
    rownames_to_column("Sample")
  
  return(GMHI.res)
}

# Function to compute B/F ratio from metaphlan3 input 
bf_ratio <- function(phylo){
  
  # Check to make sure input is a phyloseq object
  if(class(phylo)[1] != "phyloseq"){
    message("Error, input is not a phyloseq object!")
    return(NULL)
  }
  
  # Collapse on phyla
  phyla <- phyloseq::tax_glom(physeq = phylo, taxrank = "Phylum")
  
  # Find relative abundances
  phyla_rel <- phyloseq::transform_sample_counts(phyla, function(x) { x/sum(x) } )
  
  # Keep B/F taxa
  tax_table(phyla_rel)
  phyla_rel_bact <- otu_table(subset_taxa(phyla_rel, Phylum == "Bacteroidetes"))
  phyla_rel_firm <- suppressWarnings(otu_table(subset_taxa(phyla_rel, Phylum == "Firmicutes")))
  
  # OTU
  bf_ratio <- log2(phyla_rel_bact /  phyla_rel_firm)
  
  # Add to sample metadata
  phyloseq::sample_data(phylo)$log2_bf_ratio <- as.numeric(bf_ratio)
  
  # Return phyllseq object
  return(phylo)
}

# Get ancom results table
get_ancom = function(result, namesto = "Feature"){
  t1 = result$res$W %>% 
    as.data.frame() %>% 
    rownames_to_column("Taxa") %>% 
    pivot_longer(cols = -Taxa, names_to = namesto, values_to = "beta")
  
  t2 = result$res$p_val %>% 
    as.data.frame() %>% 
    rownames_to_column("Taxa") %>% 
    pivot_longer(cols = -Taxa, names_to = namesto, values_to = "p_val") %>% 
    select(-Taxa, -Feature)
  
  t3 = result$res$q_val %>% 
    as.data.frame() %>% 
    rownames_to_column("Taxa") %>% 
    pivot_longer(cols = -Taxa, names_to = namesto, values_to = "q_val") %>% 
    select(-Taxa, -Feature)
  
  t4 = result$res$se %>% 
    as.data.frame() %>% 
    rownames_to_column("Taxa") %>% 
    pivot_longer(cols = -Taxa, names_to = namesto, values_to = "se") %>% 
    select(-Taxa, -Feature)
  
  ancom.merged = bind_cols(t1, t2, t3, t4)
  return(ancom.merged)
}

calc_rclr <- function(pseq.in, ...){
  #https://github.com/microbiome/mia/blob/master/R/transformCounts.R
  pseq.pruned = prune_samples(sample_sums(pseq.in)>0, pseq.in)
  pseq.pruned = prune_taxa(taxa_sums(pseq.in)>0, pseq.in)
  mat = t(abundances(pseq.pruned))
  # Performs logarithmic transform
  log_mat <- log(mat)
  # If there are zeros, they are converted into infinite values. 
  # They are converted to NAs.
  log_mat[is.infinite(log_mat)] <- NA
  # Calculates means for every sample, does not take NAs into account
  mean_log_mat <- DelayedMatrixStats::colMeans2(log_mat, na.rm = TRUE)
  # Calculates exponential values from means, i.e., geometric means
  geometric_means_of_samples <- exp(mean_log_mat)
  # Divides all values by their sample-wide geometric means
  values_divided_by_geom_mean <- t(mat)/geometric_means_of_samples
  # Does logarithmic transform and transposes the table back to its original form
  return_mat <- t(log(values_divided_by_geom_mean))
  # If there were zeros, there are infinite values after logarithmic transform. 
  # They are converted to zero.
  return_mat[is.infinite(return_mat)] <- 0
  pseq.clr.res <- pseq.pruned
  otu_table(pseq.clr.res) <- otu_table(t(return_mat), taxa_are_rows = T)
  return(pseq.clr.res)
}

plot_nmds_stars = function(pseq.input, nmds, group = "CB", size = 2){
  
  scores = as.data.frame(nmds$points)
  group.vec = meta(pseq.input)[,group]
  
  gg <- data.frame(Group = group.vec, x=scores$MDS1, y=scores$MDS2)
  centroids <- aggregate(cbind(x,y) ~ Group, data = gg, median)
  gg <- merge(gg,centroids,by="Group",suffixes=c("",".centroid")) %>% 
    mutate(Group = fct_relevel(Group, "Responder"))
  
  gg.plot <- ggplot(gg) + 
    geom_point(aes(x=x,y=y,color=Group), size=size, pch = 1) +
    geom_point(data=centroids, aes(x=x, y=y, color=Group), size=size + 1) +
    geom_segment(aes(x=x.centroid, y=y.centroid, xend=x, yend=y, color=Group)) +
    stat_ellipse(mapping = aes(x = x, y = y, color = Group), type = "norm", linetype = 2, size = 0.2) +
    xlab("NMDS1") +
    ylab("NMDS2")
  return(gg.plot)
}


plot_pcoa_stars = function(pseq.input, pcoa, group = "CB", size = 2){
  
  scores = as.data.frame(pcoa$vectors)
  group.vec = meta(pseq.input)[,group]
  pc1_variance = (pcoa$values$Relative_eig[1] * 100) %>% round(2)
  pc2_variance = (pcoa$values$Relative_eig[2] * 100) %>% round(2)
  
  gg <- data.frame(Group = group.vec, x=scores$Axis.1, y=scores$Axis.2)
  centroids <- aggregate(cbind(x,y) ~ Group, data = gg, median)
  gg <- merge(gg,centroids,by="Group",suffixes=c("",".centroid")) %>% 
    mutate(Group = fct_relevel(Group, "Responder"))
  
  library(ggside)
  gg.plot <- ggplot(gg) + 
    geom_point(aes(x=x,y=y,color=Group), size=size, pch = 1) +
    geom_point(data=centroids, aes(x=x, y=y, color=Group), size=size + 1) +
    geom_segment(aes(x=x.centroid, y=y.centroid, xend=x, yend=y, color=Group)) +
    stat_ellipse(mapping = aes(x = x, y = y, color = Group), type = "norm", linetype = 2, size = 0.2) +
    xlab(paste0("PC1 (", pc1_variance, "%)")) +
    ylab(paste0("PC2 (", pc2_variance, "%)")) 
    return(gg.plot)
}


# Compute distance to healthy centoid
betadisper2 <- function(d, group, type = c("median","centroid"), bias.adjust=FALSE, sqrt.dist = FALSE, add = FALSE, healthy = NULL) {
  ## inline function for double centring. We used .C("dblcen", ...,
  ## PACKAGE = "stats") which does not dublicate its argument, but
  ## it was removed from R in r60360 | ripley | 2012-08-22 07:59:00
  ## UTC (Wed, 22 Aug 2012) "more conversion to .Call, clean up".
  
  ordimedian <- function(ord, groups, display = "sites", label = FALSE, ...) {
    ## Sum of distances from the statistic
    medfun <-
      function(x, ord) sum(sqrt(rowSums(sweep(ord, 2, x)^2)),
                           na.rm = TRUE)
    ## derivative of medfun (if NULL, optim will use numerical
    ## differentiation)
    dmedfun <- function(x, ord) {
      up <- -sweep(ord, 2, x)
      dn <- sqrt(rowSums(sweep(ord, 2, x)^2))
      colSums(sweep(up, 1, dn, "/"))
    }
    #dmedfun <- NULL
    pts <- scores(ord, display = display, ...)
    inds <- names(table(groups))
    medians <- matrix(NA, nrow = length(inds), ncol = ncol(pts))
    rownames(medians) <- inds
    colnames(medians) <- colnames(pts)
    for (i in inds) {
      X <- pts[groups == i, , drop = FALSE]
      if (NROW(X) > 0)
        medians[i, ] <- optim(apply(X, 2, median, na.rm = TRUE),
                              fn = medfun, gr = dmedfun,
                              ord = X, method = "BFGS")$par
      if(label)
        ordiArgAbsorber(medians[i,1], medians[i,2], label = i,
                        FUN = text, ...)
    }
    invisible(medians)
  }
  dblcen <- function(x, na.rm = TRUE) {
    cnt <- colMeans(x, na.rm = na.rm)
    x <- sweep(x, 2L, cnt, check.margin = FALSE)
    cnt <- rowMeans(x, na.rm = na.rm)
    sweep(x, 1L, cnt, check.margin = FALSE)
  }
  ## inline function for spatial medians
  spatialMed <- function(vectors, group, pos) {
    axes <- seq_len(NCOL(vectors))
    spMedPos <- ordimedian(vectors, group, choices = axes[pos])
    spMedNeg <- ordimedian(vectors, group, choices = axes[!pos])
    cbind(spMedPos, spMedNeg)
  }
  ## inline function for centroids
  centroidFUN <- function(vec, group) {
    cent <- apply(vec, 2,
                  function(x, group) tapply(x, INDEX = group, FUN = mean),
                  group = group)
    if(!is.matrix(cent)) { ## if only 1 group, cent is vector
      cent <- matrix(cent, nrow = 1,
                     dimnames = list(as.character(levels(group)),
                                     paste0("Dim", seq_len(NCOL(vec)))))
    }
    cent
  }
  ## inline function for distance computation
  Resids <- function(x, c) {
    if(is.matrix(c))
      d <- x - c
    else
      d <- sweep(x, 2, c)
    rowSums(d^2)
  }
  ## Tolerance for zero Eigenvalues
  TOL <- sqrt(.Machine$double.eps)
  ## uses code from stats:::cmdscale by R Core Development Team
  if(!inherits(d, "dist"))
    stop("distances 'd' must be a 'dist' object")
  ## Someone really tried to analyse correlation like object in range -1..+1
  if (any(d < -TOL, na.rm = TRUE))
    stop("dissimilarities 'd' must be non-negative")
  ## adjust to avoid negative eigenvalues (if they disturb you)
  if (sqrt.dist)
    d <- sqrt(d)
  if (is.logical(add) && isTRUE(add))
    add <- "lingoes"
  if (is.character(add)) {
    add <- match.arg(add, c("lingoes", "cailliez"))
    if (add == "lingoes") {
      ac <- addLingoes(as.matrix(d))
      d <- sqrt(d^2 + 2 * ac)
    }
    else if (add == "cailliez") {
      ac <- addCailliez(as.matrix(d))
      d <- d + ac
    }
  }
  if(missing(type))
    type <- "median"
  type <- match.arg(type)
  ## checks for groups - need to be a factor for later
  group <- if(!is.factor(group)) {
    as.factor(group)
  } else { ## if already a factor, drop empty levels
    droplevels(group, exclude = NA) # need exclude = NA under Rdevel r71113
  }
  n <- attr(d, "Size")
  x <- matrix(0, ncol = n, nrow = n)
  x[row(x) > col(x)] <- d^2
  ## site labels
  labs <- attr(d, "Labels")
  ## remove NAs in group
  if(any(gr.na <- is.na(group))) {
    group <- group[!gr.na]
    x <- x[!gr.na, !gr.na]
    ## update n otherwise C call crashes
    n <- n - sum(gr.na)
    ## update labels
    labs <- labs[!gr.na]
    message("missing observations due to 'group' removed")
  }
  ## remove NA's in d
  if(any(x.na <- apply(x, 1, function(x) any(is.na(x))))) {
    x <- x[!x.na, !x.na]
    group <- group[!x.na]
    ## update n otherwise C call crashes
    n <- n - sum(x.na)
    ## update labels
    labs <- labs[!x.na]
    message("missing observations due to 'd' removed")
  }
  x <- x + t(x)
  x <- dblcen(x)
  e <- eigen(-x/2, symmetric = TRUE)
  vectors <- e$vectors
  eig <- e$values
  ## Remove zero eigenvalues
  eig <- eig[(want <- abs(eig) > max(TOL, TOL * eig[1L]))]
  ## scale Eigenvectors
  vectors <- vectors[, want, drop = FALSE] %*% diag(sqrt(abs(eig)),
                                                    nrow = length(eig))
  ## store which are the positive eigenvalues
  pos <- eig > 0
  ## group centroids in PCoA space
  centroids <-
    switch(type,
           centroid = centroidFUN(vectors, group),
           median = spatialMed(vectors, group, pos)
    )
  
  # Filter
  centroids = centroids[row.names(centroids) == healthy, , drop = F]
  group = rep("dummy.group", length(group))
  group = as.factor(group)
  
  ## for each of the groups, calculate distance to centroid for
  ## observation in the group
  ## Uses in-line Resids function as we want LAD residuals for
  ## median method, and LSQ residuals for centroid method
  dist.pos <- Resids(vectors[, pos, drop=FALSE],
                     centroids[group, pos, drop=FALSE])
  dist.neg <- 0
  if(any(!pos))
    dist.neg <- Resids(vectors[, !pos, drop=FALSE],
                       centroids[group, !pos, drop=FALSE])
  
  ## zij are the distances of each point to its group centroid
  if (any(dist.neg > dist.pos)) {
    ## Negative squared distances give complex valued distances:
    ## take only the real part (which is zero). Github issue #306.
    warning("some squared distances are negative and changed to zero")
    zij <- Re(sqrt(as.complex(dist.pos - dist.neg)))
  } else {
    zij <- sqrt(dist.pos - dist.neg)
  }
  if (bias.adjust) {
    n.group <- as.vector(table(group))
    zij <- zij*sqrt(n.group[group]/(n.group[group]-1))
  }
  ## pre-compute group mean distance to centroid/median for `print` method
  grp.zij <- tapply(zij, group, "mean")
  ## add in correct labels
  if (any(want))
    colnames(vectors) <- names(eig) <-
    paste("PCoA", seq_along(eig), sep = "")
  if(is.matrix(centroids))
    colnames(centroids) <- names(eig)
  else
    names(centroids) <- names(eig)
  rownames(vectors) <- names(zij) <- labs
  retval <- list(eig = eig, vectors = vectors, distances = zij,
                 group = group, centroids = centroids,
                 group.distances = grp.zij, call = match.call())
  class(retval) <- "betadisper"
  attr(retval, "method") <- attr(d, "method")
  attr(retval, "type") <- type
  attr(retval, "bias.adjust") <- bias.adjust
  retval
}

# Get Adonis output
get_adonis = function(input, group = NULL, method = "euclidean") {
  meta.data = meta(input) 
  
  if(method %in% c("unifrac", "wunifrac", "dpcoa", "jsd")){
    dist.data = distance(input, method = method) 
  } else {
    
    dist.data = vegdist(t(abundances(input)), method = method) 
  }
 
  form.data = as.formula(paste0("dist.data ~ ", group))
  res = adonis2(form.data, data = meta.data, permutations = 999, by = "margin") %>% 
    as_tibble(rownames = "term")
  return(res)
}

# Parse ANCOM-1 output
process_ancom = function(ancom, alpha = 0.10, cutoff = 0.7){

  res2 = ancom$res
  q_val = ancom$q_data
  beta_val = ancom$beta_data
  
  # Only consider the effect sizes with the corresponding q-value less than alpha
  beta_val = beta_val * (q_val < alpha) 
  
  # Choose the maximum of beta's as the effect size
  beta_pos = apply(abs(beta_val), 2, which.max) 
  beta_max = vapply(seq_along(beta_pos), function(i) beta_val[beta_pos[i], i],
                    FUN.VALUE = double(1))
  
  # Choose the maximum of beta's as the effect size
  qval_pos = apply(q_val, 2, which.min) 
  qval_min = vapply(seq_along(qval_pos), function(i) q_val[qval_pos[i], i],
                    FUN.VALUE = double(1))
  
  # Number of taxa except structural zeros
  n_taxa = ifelse(is.null(ancom$zero_ind), 
                  nrow(ancom$res), 
                  sum(apply(ancom$zero_ind, 1, sum) == 0))
  
  # Cutoff values for declaring differentially abundant taxa
  cut_off_ntaxa = cutoff * (n_taxa - 1)
  message(paste0("Cutoff: ", cutoff))
  message(paste0("Cutoff ntaxa: ", cut_off_ntaxa))
  
  if(cutoff == 0.6) {
    df_fig_w = res2 %>%
      dplyr::mutate(beta = beta_max,
                    direct = case_when(
                      detected_0.6 == TRUE & beta > 0 ~ "Positive",
                      detected_0.6 == TRUE & beta <= 0 ~ "Negative",
                      TRUE ~ "Not Significant"
                    )) %>%
      dplyr::arrange(W)
  } else if (cutoff == 0.7) {
    df_fig_w = res2 %>%
      dplyr::mutate(beta = beta_max,
                    direct = case_when(
                      detected_0.7 == TRUE & beta > 0 ~ "Positive",
                      detected_0.7 == TRUE & beta <= 0 ~ "Negative",
                      TRUE ~ "Not Significant"
                    )) %>%
      dplyr::arrange(W)
  } else if (cutoff == 0.8) {
    df_fig_w = res2 %>%
      dplyr::mutate(beta = beta_max,
                    direct = case_when(
                      detected_0.8 == TRUE & beta > 0 ~ "Positive",
                      detected_0.8 == TRUE & beta <= 0 ~ "Negative",
                      TRUE ~ "Not Significant"
                    )) %>%
      dplyr::arrange(W)
  } else if (cutoff == 0.9) {
    df_fig_w = res2 %>%
      dplyr::mutate(beta = beta_max,
                    direct = case_when(
                      detected_0.9 == TRUE & beta > 0 ~ "Positive",
                      detected_0.9 == TRUE & beta <= 0 ~ "Negative",
                      TRUE ~ "Not Significant"
                    )) %>%
      dplyr::arrange(W)
  } else {
    message("Error. Incorrect cutoff")
  }
  
  df_fig_w$taxon_id = factor(df_fig_w$taxon_id, levels = df_fig_w$taxon_id)
  df_fig_w$W = replace(df_fig_w$W, is.infinite(df_fig_w$W), n_taxa - 1)
  df_fig_w$direct = factor(df_fig_w$direct, 
                           levels = c("Negative", "Positive", "Not Significant"))
  #df_fig_w$qval = qval_min

  p_w = df_fig_w %>%
    ggplot(aes(x = taxon_id, y = W, color = direct)) +
    geom_point(size = 2, alpha = 0.6) +
    labs(x = "Taxon", y = "W") +
    scale_color_discrete(name = NULL) + 
    geom_hline(yintercept = cut_off_ntaxa, linetype = "dotted", 
               color = "blue", size = 1.5) +
    geom_text(aes(x = 2, y = cut_off_ntaxa + 0.5, label = "W[0.7]"), 
              size = 5, vjust = -0.5, hjust = 0, color = "orange", parse = TRUE) +
    theme_bw() +
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          panel.grid.major = element_blank())
  print(p_w)
  
  return(df_fig_w %>% arrange(beta))
}



# Import humann3 pathway abundances
import_metacyc = function(table, sampledata, taxonomy = T, type = "metacyc", exclude = T) {
  pathways = data.table::fread(table)
  colnames(pathways)[1] = "metacyc"
  
  # Remove taxa
  if(!taxonomy){
    pathways = pathways %>% 
      filter(!(grepl("\\|", metacyc))) 
  }
  
  # Separate names
  pathways = pathways %>% 
    separate(metacyc, into = c("ID", "Name"), sep = ":")
  
  if(exclude) {
    pathways = pathways %>% 
      filter(!(ID %in% c("UNMAPPED", "UNINTEGRATED")))
  }
  
  # Fix colnames
  if(type == "metacyc") {
    colnames(pathways) = str_remove(colnames(pathways), "_Abundance")
  } else {
    colnames(pathways) = str_remove(colnames(pathways), "_Abundance-RPKs")
  }
  
  # Fix samples not found
  #sampledata.fil = sampledata %>% 
  #  filter(run %in% colnames(pathways))
  
  # Convert to phyloseq object
  pathways.otu = pathways %>% 
    select(-any_of(c("Name","NA"))) %>% 
    column_to_rownames("ID") %>% 
    phyloseq::otu_table(taxa_are_rows = T) 
  
  pathways.tax = pathways %>% 
    select(-any_of(c("Name","NA"))) %>%  
    column_to_rownames("ID") %>% 
    rownames_to_column("Species") %>% 
    select(Species) %>% 
    mutate(id = Species) %>% 
    column_to_rownames("id") %>% 
    as.matrix() %>% 
    phyloseq::tax_table()
  
  # Return pathways as phyloseq object
  pseq = phyloseq(pathways.otu, pathways.tax, sample_data(sampledata))
  return(pseq)
  
}


# Import humann3 pathway abundances
import_humann = function(table, sampledata, taxonomy = T) {
  gene.table = data.table::fread(table)
  colnames(gene.table)[1] = "gene"
  
  # Remove taxa
  if(!taxonomy){
    gene.table = gene.table %>% 
      filter(!(grepl("\\|", gene)))
  }
  
  # Fix colnames
  colnames(gene.table) = str_remove(colnames(gene.table), "_Abundance-RPKs")
  
  # Convert to phyloseq object
  gene.table.otu = gene.table %>% 
    column_to_rownames("gene") %>% 
    otu_table(taxa_are_rows = T) 
  
  gene.table.tax = gene.table %>% 
    column_to_rownames("gene") %>% 
    rownames_to_column("Species") %>% 
    select(Species) %>% 
    mutate(id = Species) %>% 
    column_to_rownames("id") %>% 
    as.matrix() %>% 
    tax_table()
  
  # Return pathways as phyloseq object
  pseq = phyloseq(gene.table.otu, gene.table.tax, sample_data(sampledata))
  return(pseq)
}

# Predict metabolitc compostion with melonpann
predict_metabolites = function(table, sampledata) {

  melonnpan.input = data.table::fread(file = table) 
  colnames(melonnpan.input)[1] = "gene_families"
  
  melonnpan.input = melonnpan.input %>% 
    dplyr::filter(!(grepl("|", gene_families, fixed = T))) %>% 
    dplyr::filter(gene_families != "UNMAPPED") %>% 
    dplyr::filter(gene_families %in% colnames(melonnpan::melonnpan.training.data)) %>% 
    column_to_rownames("gene_families") %>% 
    t() %>% 
    as.data.frame()  
  
  # Predict with melonpann
  melonpan.metabolites = melonnpan::melonnpan.predict(melonnpan.input, output = "data/melonpann.txt")
  
  # Convert to phyloseq object
  melon.otu = melonpan.metabolites$pred %>% 
    as.data.frame() %>% 
    dplyr::mutate(ID = str_remove(ID, "_Abundance-RPKs")) %>% 
    column_to_rownames("ID") %>% 
    t() %>% 
    otu_table(taxa_are_rows = T) 
  
  melon.tax = melonpan.metabolites$pred %>% 
    as.data.frame() %>% 
    dplyr::mutate(ID = str_remove(ID, "_Abundance-RPKs")) %>% 
    column_to_rownames("ID") %>% 
    t() %>% 
    as.data.frame() %>% 
    rownames_to_column("Species") %>% 
    dplyr::select(Species) %>% 
    dplyr::mutate(id = Species) %>% 
    column_to_rownames("id") %>% 
    as.matrix() %>% 
    tax_table()
  
  melon.pseq = phyloseq(melon.otu, melon.tax, sample_data(sampledata))
  gc()
  return(melon.pseq)
}

# Import pathseq counts files
import_pathseq = function (filepath, filenames, metadata, level = "genus", count_feature = "unambiguous", minimum = 2) {
  message("Importing PathSeq files...")
  #names(filepath) = filenames
  pathseq.data = filepath %>% 
    map_dfr(.id = "sampleid", .f = function(x){
      tbl = data.table::fread(x)
      if(nrow(tbl) == 0){
        return(NULL)
      } else { 
        return(tbl)
      }
    })
  
  # Only keep intersecting samples between metadata and counts
  message("Finding intersecting samples...")
  pathseq.data = pathseq.data %>% 
    filter(sampleid %in% row.names(metadata))
  
  message("Converting to phyloseq object...")
  otu.table = pathseq.data %>% 
    filter(type == level) %>% 
    filter(unambiguous >= minimum) %>% 
    dplyr::select(tax_id, count_feature, sampleid) %>% 
    pivot_wider(names_from = sampleid, values_from = count_feature) %>% 
    column_to_rownames("tax_id") %>% 
    replace(., is.na(.), 0) %>% as.matrix()
  
  TAX = NKIatlas:::TAX
  message(paste0("Detected ", nrow(otu.table), " taxa"))
  message(paste0("Detected ", ncol(otu.table), " samples"))
  message(paste0("Missing ", paste(setdiff(row.names(otu.table), row.names(TAX)), collapse = ","), " taxa ID's."))
  SAMPLE = metadata %>% sample_data()
  OTU = otu_table(otu.table, taxa_are_rows = TRUE)
  pseq = phyloseq(OTU, TAX, SAMPLE)
  pseq.fil = prune_taxa(taxa_sums(pseq) > 0, pseq)
  return(pseq.fil)
}


# Function to import Bracken counts into phyloseq object
import_bracken = function (filepath, filenames, metadata) {
  message("Importing Kraken2/Bracken files...")
  kraken2.import = data.frame(sampleid = filenames, filelocation = filepath) %>% 
    mutate(file_contents = map(filelocation, function(x) read.delim(x))) %>% 
    unnest(cols = c(file_contents)) %>% dplyr::select(-filelocation)
  
  # Only keep intersecting samples between metadata and counts
  message("Finding intersecting samples...")
  kraken2.import = kraken2.import %>% 
    filter(sampleid %in% row.names(metadata))
  
  # Convert to phyloseq object
  message("Converting to phyloseq object...")
  otu.table = kraken2.import %>% 
    dplyr::select(taxonomy_id, new_est_reads, sampleid) %>% 
    pivot_wider(names_from = sampleid, values_from = new_est_reads) %>% column_to_rownames("taxonomy_id") %>% 
    replace(., is.na(.), 0) %>% as.matrix()
  TAX = NKIatlas:::TAX
  message(paste0("Detected ", nrow(otu.table), " taxa"))
  message(paste0("Detected ", ncol(otu.table), " samples"))
  message(paste0("Missing ", paste(setdiff(row.names(otu.table), 
                                           row.names(TAX)), collapse = ","), " taxa ID's."))
  SAMPLE = metadata %>% sample_data()
  OTU = otu_table(otu.table, taxa_are_rows = TRUE)
  pseq = phyloseq(OTU, TAX, SAMPLE)
  pseq.fil = prune_taxa(taxa_sums(pseq) > 0, pseq)
  return(pseq.fil)
}


# Function to null-out the Subject-OTU entries not found in Pathseq
get_intersection = function(kraken2, pathseq, pcutoff = 1, kcutoff = 10){
  
  # For each sample, filter kraken2 based on the microbes found in pathseq
  kraken2.abun = abundances(kraken2) %>% as.data.frame() %>% rownames_to_column("Id")
  pathseq.abun = abundances(pathseq) %>% as.data.frame() %>% rownames_to_column("Id")
  
  # Iterate over samples in Kraken2 table
  pb <- progress_estimated(length(sample_names(kraken2)))
  filtered = sample_names(kraken2) %>% 
    map(.f = function(x){
      
      pb$tick()$print()
      to.fil = select(pathseq.abun, "Id" , x) %>% 
        filter(.data[[x]] >= pcutoff) %>% 
        pull(Id)
      
      select(kraken2.abun, "Id", x) %>% 
        filter(Id %in% to.fil) %>% 
        filter(.data[[x]] >= kcutoff) 
    }) 
  
  # Merge all tables together
  joined = filtered %>% 
    purrr::reduce(full_join, by = "Id") %>% 
    replace(is.na(.), 0)
  
  # Convert to phyloseq object
  mat = joined %>% 
    column_to_rownames("Id") %>% 
    as.matrix()
  
  # Add back to phylsoeq object
  results = kraken2
  otu_table(results) = otu_table(mat, taxa_are_rows = T)
  
  # Prune samples with no reads
  results = prune_samples(sample_sums(results)>0, results)
  results = prune_taxa(taxa_sums(results)>0, results)
  
  return(results)
}

# Function to get metacyc table split by taxa
get_protein = function(input, taxa = T) {
  
  # Import table
  input = data.table::fread(input)
  colnames(input)[1] = "pathway"
  
  # Fix colnames
  colnames(input) = str_remove(colnames(input), "_Abundance-RPKs")
  
  # Separate names
  if(taxa) {
    res = input %>% 
      filter(grepl("g__", pathway)) %>% 
      separate(pathway, into = c("Name", "Taxa"), sep = "\\|") %>% 
      separate(Taxa, into = c("Genus", "Species"), sep = "\\.") %>% 
      mutate(Genus = str_remove(Genus, "g__")) %>% 
      mutate(Species = str_remove(Species, "s__")) %>% 
      as_tibble()
  } else {
    res = input %>% 
      filter(!grepl("g__", pathway)) %>% 
      filter(!grepl("\\|", pathway)) %>% 
      as_tibble()
  }


  return(res)
}

# Function to get metacyc table split by taxa
get_metacyc_taxa = function(metacyc) {
  
  # Import metacyc table
  metacyc = data.table::fread(metacyc)
  colnames(metacyc)[1] = "pathway"
  
  # Fix colnames
  colnames(metacyc) = str_remove(colnames(metacyc), "_Abundance")
  
  # Separate names
  metacyc = metacyc %>% 
    filter(grepl("g__", pathway)) %>% 
    separate(pathway, into = c("ID", "Name"), sep = ": ") %>% 
    separate(Name, into = c("Name", "Taxa"), sep = "\\|") %>% 
    separate(Taxa, into = c("Genus", "Species"), sep = "\\.") %>% 
    mutate(Genus = str_remove(Genus, "g__")) %>% 
    mutate(Species = str_remove(Species, "s__")) %>% 
    as_tibble()
  return(metacyc)
}

# Run maaslin2
run_maaslin = function(pseq.in, 
                       transform = "LOG", 
                       analysis_method = "LM", 
                       normalization = "TSS", 
                       fixed_effects = NULL, 
                       random_effects = NULL, 
                       max_significance = 0.25, 
                       min_prevalence = 0.10, 
                       min_variance = 0.0, 
                       correction = "BH",
                       min_abundance = 0.0,
                       reference = c("primaryTumorLocation,CUP", "biopsySite,CNS")) {
  
  # get temporary directory
  tmp.dir = tempdir()
  
  # Format counts table
  pseq.in %>% 
    #aggregate_taxa(level = "Species") %>% 
    abundances() %>% 
    t() %>% 
    as.data.frame() %>% 
    rownames_to_column("ID") %>% 
    write_tsv(paste0(tmp.dir, "/maaslin2.counts.tsv"))
  
  # Format metadata
  meta(pseq.in) %>% 
    rownames_to_column("ID") %>% 
    write_tsv(paste0(tmp.dir, "/maaslin2.metadata.tsv"))
  
  # Run command
  library(Maaslin2)
  fit_data <- Maaslin2( input_data = paste0(tmp.dir, "/maaslin2.counts.tsv"), 
                        input_metadata = paste0(tmp.dir, "/maaslin2.metadata.tsv"), 
                        output = paste0(tmp.dir, "/outdir"), 
                        transform = transform,
                        normalization = normalization,
                        fixed_effects = fixed_effects,
                        random_effects = random_effects,
                        analysis_method = analysis_method,
                        standardize = T,
                        cores = 2, 
                        correction = correction,
                        min_variance = min_variance,
                        min_abundance = min_abundance,
                        min_prevalence = min_prevalence,
                        max_significance = max_significance,
                        plot_heatmap = F, 
                        plot_scatter = F,
                        reference = reference)
  
  
  # Import results
  results = read_tsv( paste0(tmp.dir, "/outdir/all_results.tsv"))
  sig_results = read_tsv( paste0(tmp.dir, "/outdir/significant_results.tsv"))
  #normalized = read_tsv( paste0(tmp.dir, "/outdir/features/filtered_data_norm_transformed.tsv"))
  #model = readRDS( paste0(tmp.dir, "/outdir/fits/models.rds"))
  #residuals = readRDS( paste0(tmp.dir, "/outdir/fits/residuals.rds"))
  #fitted = readRDS( paste0(tmp.dir, "/outdir/fits/fitted.rds"))
  #ranef = readRDS( paste0(tmp.dir, "/outdir/fits/ranef.rds"))
  #scatter_plots = readRDS( paste0(tmp.dir, "/outdir/figures/scatter_plots.rds"))
  
  # Return
  return(list(results = results, sig_results = sig_results))
  
}


# Filter based on variance
filter_variance = function(pseq.in, var = 0.20) {
  
  vars.list = abundances(pseq.in) %>% 
    genefilter::rowVars() %>% 
    sort(decreasing = T) %>% 
    .[1:(var * ntaxa(pseq.in))] %>% 
    names()


  return(vars.list)

}


get_residuals2 = function(input.table, formula = c("primaryTumorLocation", "biopsySite")) {
  input.table.long = input.table %>% 
    pivot_longer(cols = -sampleId, names_to = "signatures", values_to = "y") %>% 
    inner_join(purple.metadata) %>% 
    as_tibble() 
  
  unique(input.table.long$signatures) %>% 
    map_dfr(.f = function(x) {
      input.table.long %>% 
        filter(signatures == x) %>% 
        tibble::column_to_rownames("sampleId") %>% 
        glm(as.formula(paste("y", paste(formula, collapse=" + "), sep=" ~ ")), data = ., family = "gaussian") %>% 
        broom::augment() %>% 
        rename(sampleId = `.rownames`,
               residual = `.resid`) %>% 
        mutate(diff = residual - mean(residual)) %>% 
        select(sampleId, diff) %>% 
        mutate(signatures = x)
    }) %>% 
    pivot_wider(id_cols = sampleId, names_from = signatures, values_from = diff) %>% 
    as.data.frame() 
}



# Map colors to names of brewer paletter. 
map_colors <- function(variable, palette_name = "Set3") {
  # Check if the input is a factor; if not, convert it to a factor
  if (!is.factor(variable)) {
    variable <- factor(variable)
  }
  
  # Get the unique levels of the factor
  levels <- levels(variable)
  num_levels <- length(levels)
  
  # Check if the palette has enough colors
  if (num_levels > brewer.pal.info[palette_name, "maxcolors"]) {
    message( sprintf("The selected palette (%s) does not have enough colors for %s levels.",
                     palette_name, 
                     num_levels)
    )
    message ("Extending the number of colors in the palette....")
    
    # Generate the color palette with the same number of colors as the number of levels
    color_palette <- colorRampPalette(brewer.pal(8, palette_name))(num_levels)
    
  } else {
    
    # Generate the color palette with the same number of colors as the number of levels
    color_palette <- brewer.pal(num_levels, palette_name)
  }
  
  # Create a named vector that maps each level to a color
  color_mapping <- setNames(color_palette, levels)
  
  return(color_mapping)
}


# Updated survminer::ggforest fixing a bug. 
ggforest2 <- function (model,
                       data = NULL,
                       main = "Hazard ratio",
                       cpositions = c(0.02, 0.22, 0.4),
                       fontsize = 0.7, refLabel = "reference", noDigits = 2) 
{
  conf.high <- conf.low <- estimate <- NULL
  stopifnot(inherits(model, "coxph"))
  data <- survminer:::.get_data(model, data = data)
  terms <- attr(model$terms, "dataClasses")[-1]
  coef <- as.data.frame(broom::tidy(model, conf.int = TRUE))
  gmodel <- broom::glance(model)
  allTerms <- lapply(seq_along(terms), function(i) {
    var <- names(terms)[i]
    if (terms[i] %in% c("factor", "character")) {
      adf <- as.data.frame(table(data[, var]))
      cbind(var = var, adf, pos = 1:nrow(adf))
    }
    else if (terms[i] == "numeric") {
      data.frame(var = var, Var1 = "", Freq = nrow(data), 
                 pos = 1)
    }
    else {
      vars <- grep(paste0("^", var, "*."), coef$term, 
                   value = TRUE)
      data.frame(var = vars, Var1 = "", Freq = nrow(data), 
                 pos = seq_along(vars))
    }
  })
  allTermsDF <- do.call(rbind, allTerms)
  colnames(allTermsDF) <- c("var", "level", "N", "pos")
  inds <- apply(allTermsDF[, 1:2], 1, paste0, collapse = "")
  rownames(coef) <- gsub(coef$term, pattern = "`", replacement = "")
  toShow <- cbind(allTermsDF, coef[inds, ])[, c("var", "level", 
                                                "N", "p.value", "estimate", "conf.low", "conf.high", 
                                                "pos")]
  toShowExp <- toShow[, 5:7]
  toShowExp[is.na(toShowExp)] <- 0
  toShowExp <- format(exp(toShowExp), digits = noDigits)
  toShowExpClean <- data.frame(toShow, pvalue = signif(toShow[, 
                                                              4], noDigits + 1), toShowExp)
  toShowExpClean$stars <- paste0(round(toShowExpClean$p.value, 
                                       noDigits + 1), " ", ifelse(toShowExpClean$p.value < 
                                                                    0.05, "*", ""), ifelse(toShowExpClean$p.value < 0.01, 
                                                                                           "*", ""), ifelse(toShowExpClean$p.value < 0.001, "*", 
                                                                                                            ""))
  toShowExpClean$ci <- paste0("(", toShowExpClean[, "conf.low.1"], 
                              " - ", toShowExpClean[, "conf.high.1"], ")")
  toShowExpClean$estimate.1[is.na(toShowExpClean$estimate)] <- refLabel
  toShowExpClean$stars[which(toShowExpClean$p.value < 0.001)] <- "<0.001 ***"
  toShowExpClean$stars[is.na(toShowExpClean$estimate)] <- ""
  toShowExpClean$ci[is.na(toShowExpClean$estimate)] <- ""
  toShowExpClean$estimate[is.na(toShowExpClean$estimate)] <- 0
  toShowExpClean$var <- as.character(toShowExpClean$var)
  toShowExpClean$var[duplicated(toShowExpClean$var)] <- ""
  toShowExpClean$N <- paste0("(N=", toShowExpClean$N, ")")
  toShowExpClean <- toShowExpClean[nrow(toShowExpClean):1, 
  ]
  ci_h <- toShowExpClean$conf.high[is.finite(as.numeric(toShowExpClean$conf.high.1))]
  ci_l <- toShowExpClean$conf.low[!near(as.numeric(toShowExpClean$conf.low.1), 
                                        0)]
  toShowExpClean$conf.high[is.infinite(as.numeric(toShowExpClean$conf.high.1))] <- max(ci_h, 
                                                                                       na.rm = TRUE)
  toShowExpClean$conf.low[near(as.numeric(toShowExpClean$conf.low.1), 
                               0)] <- min(ci_l, na.rm = TRUE)
  rangeb <- range(toShowExpClean$conf.low, toShowExpClean$conf.high, 
                  na.rm = TRUE)
  breaks <- axisTicks(rangeb/2, log = TRUE, nint = 7)
  rangeplot <- rangeb
  rangeplot[1] <- rangeplot[1] - diff(rangeb)
  rangeplot[2] <- rangeplot[2] + 0.15 * diff(rangeb)
  width <- diff(rangeplot)
  y_variable <- rangeplot[1] + cpositions[1] * width
  y_nlevel <- rangeplot[1] + cpositions[2] * width
  y_cistring <- rangeplot[1] + cpositions[3] * width
  y_stars <- rangeb[2]
  x_annotate <- seq_len(nrow(toShowExpClean))
  annot_size_mm <- fontsize * as.numeric(grid::convertX(unit(theme_get()$text$size, 
                                                             "pt"), "mm"))
  p <- ggplot(toShowExpClean, aes(seq_along(var), exp(estimate))) + 
    geom_rect(aes(xmin = seq_along(var) - 0.5, xmax = seq_along(var) + 
                    0.5, ymin = exp(rangeplot[1]), ymax = exp(rangeplot[2]), 
                  fill = ordered(seq_along(var)%%2 + 1))) + scale_fill_manual(values = c("#FFFFFF33", 
                                                                                         "#00000033"), guide = "none") + geom_point(pch = 15, 
                                                                                                                                    size = 4) + geom_errorbar(aes(ymin = exp(conf.low), 
                                                                                                                                                                  ymax = exp(conf.high)), width = 0.15) + geom_hline(yintercept = 1, 
                                                                                                                                                                                                                     linetype = 3) + coord_flip(ylim = exp(rangeplot)) + 
    ggtitle(main) + scale_y_log10(name = "", labels = sprintf("%g", 
                                                              breaks), expand = c(0.02, 0.02), breaks = breaks) + 
    theme_light() + theme(panel.grid.minor.y = element_blank(), 
                          panel.grid.minor.x = element_blank(), panel.grid.major.y = element_blank(), 
                          legend.position = "none", panel.border = element_blank(), 
                          axis.title.y = element_blank(), axis.text.y = element_blank(), 
                          axis.ticks.y = element_blank(), plot.title = element_text(hjust = 0.5)) + 
    xlab("") + annotate(geom = "text", x = x_annotate, y = exp(y_variable), 
                        label = toShowExpClean$var, fontface = "bold", hjust = 0, 
                        size = annot_size_mm) + annotate(geom = "text", x = x_annotate, 
                                                         y = exp(y_nlevel), hjust = 0, label = toShowExpClean$level, 
                                                         vjust = -0.1, size = annot_size_mm) + annotate(geom = "text", 
                                                                                                        x = x_annotate, y = exp(y_nlevel), label = toShowExpClean$N, 
                                                                                                        fontface = "italic", hjust = 0, vjust = ifelse(toShowExpClean$level == 
                                                                                                                                                         "", 0.5, 1.1), size = annot_size_mm) + annotate(geom = "text", 
                                                                                                                                                                                                         x = x_annotate, y = exp(y_cistring), label = toShowExpClean$estimate.1, 
                                                                                                                                                                                                         size = annot_size_mm, vjust = ifelse(toShowExpClean$estimate.1 == 
                                                                                                                                                                                                                                                "reference", 0.5, -0.1)) + annotate(geom = "text", 
                                                                                                                                                                                                                                                                                    x = x_annotate, y = exp(y_cistring), label = toShowExpClean$ci, 
                                                                                                                                                                                                                                                                                    size = annot_size_mm, vjust = 1.1, fontface = "italic") + 
    annotate(geom = "text", x = x_annotate, y = exp(y_stars), 
             label = toShowExpClean$stars, size = annot_size_mm, 
             hjust = -0.2, fontface = "italic") + annotate(geom = "text", 
                                                           x = 0.5, y = exp(y_variable), label = paste0("# Events: ", 
                                                                                                        gmodel$nevent, "; Global p-value (Log-Rank): ", 
                                                                                                        format.pval(gmodel$p.value.log, eps = ".001"), " \nAIC: ", 
                                                                                                        round(gmodel$AIC, 2), "; Concordance Index: ", round(gmodel$concordance, 
                                                                                                                                                             2)), size = annot_size_mm, hjust = 0, vjust = 1.2, 
                                                           fontface = "italic")
  gt <- ggplot_gtable(ggplot_build(p))
  gt$layout$clip[gt$layout$name == "panel"] <- "off"
  ggpubr::as_ggplot(gt)
}

# Functions to quickly access the 
rowtibble <- function (tse){
  row_data_tibble <- rowData(tse) %>% as_tibble()
  return(row_data_tibble)
}

coltibble <- function (tse){
  col_data_tibble <- colData(tse) %>% as_tibble()
  return(col_data_tibble)
}

rowdf <- function (tse){
  row_data_df <- rowData(tse) %>% as.data.frame()
  return(row_data_df)
}

coldf <- function (tse){
  col_data_df <- colData(tse) %>% as.data.frame()
  return(col_data_df)
}



check_lm <- function(lmodel){
  
  require (performance)
  
  checks_df = data.frame(
    normality = check_normality (lmodel) %>% capture.output() %>% paste(., collapse = "\n"),
    heteroscedasticity = check_heteroscedasticity(lmodel)  %>% capture.output() %>% paste(., collapse = "\n"),
    collinearity = check_collinearity(lmodel) %>% capture.output() %>% paste(., collapse = "\n"),
    outliers = check_outliers(lmodel)  %>% capture.output() %>% paste(., collapse = "\n"),
    residuals = check_residuals(lmodel)  %>% capture.output()  %>% paste(., collapse = "\n") 
  )
  
  return (checks_df)
  
}


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



















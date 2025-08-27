#!/usr/bin/bash

# Initialize nextflow variables to run gcp
export NXF_MODE=gcp
export GOOGLE_PROJECT_ID=<project_id>
export SERVICE_ACCOUNT_KEY_FILE=<service_account_key_file>
export GOOGLE_APPLICATION_CREDENTIALS=<google_application_credentials>
export SERVICE_ACCOUNT_EMAIL=<service_account_email>

nextflow -Dnxf.pool.type=sync run \
biobakery-nf/main.nf \
--outdir 'results' \
--csv 'SraRunTable.csv' \
--isSE \
-work-dir <WORK-DIRECTORY> \
--gcp_project $GOOGLE_PROJECT_ID \
--gcp_service_account $SERVICE_ACCOUNT_EMAIL \
--print_files false \
-resume

# - - - - - - - - - - - - - - - - -

# Merge tables together (metaphlan)
merge_metaphlan_tables.py \
results/metaphlan/taxa_profile/*.txt > \
results/metaphlan/metaphlan-profile.txt

# Merge tables together (metaphlan)
merge_metaphlan_tables.py \
results/metaphlan/read_stats/*.txt > \
results/metaphlan/metaphlan-readstats.txt

# Merge tables together (gene families)
humann_join_tables \
--input results/humann/gene_families \
--output results/humann/gene_families.txt

humann_join_tables \
--input results/humann/pathway_abundance \
--output results/humann/pathway_abundance.txt

humann_join_tables \
--input results/humann/pathway_coverage \
--output results/humann/pathway_coverage.txt


# - - - - - - - - - - - - - - - - - 
# Normalize data
humann_renorm_table \
--input results/humann/pathway_abundance.txt \
--output results/humann/pathway_abundance-cpm.txt \
--units 'cpm' \
--special 'n'

# Normalize data
humann_renorm_table \
--input results/humann/gene_families.txt \
--output results/humann/gene_families-relab.txt \
--units 'relab' \
--special 'n'

# - - - - - - - - - - - - - - - - - 
# Convert to other databases
mkdir -p results/humann/regroup

## EC
humann_regroup_table \
--input results/humann/gene_families.txt  \
--output results/humann/regroup/level4ec_families.txt \
--groups uniref90_level4ec

## KEGG
humann_regroup_table \
--input results/humann/gene_families.txt  \
--output results/humann/regroup/kegg_families.txt \
--groups uniref90_ko

## Pfam
humann_regroup_table \
--input results/humann/gene_families.txt  \
--output results/humann/regroup/pfam_families.txt \
--groups uniref90_pfam

## Eggnog
humann_regroup_table \
--input results/humann/gene_families.txt  \
--output results/humann/regroup/eggnog_families.txt \
--groups uniref90_eggnog

## RXN
humann_regroup_table \
--input results/humann/gene_families.txt  \
--output results/humann/regroup/rxn_families.txt \
--groups uniref90_rxn

## GO
humann_regroup_table \
--input results/humann/gene_families.txt  \
--output results/humann/regroup/infogo_families.txt \
--groups uniref90_go

# KEGG compounds
humann_split_table \
--input results/humann/regroup/level4ec_families.txt \
--output results/humann/regroup/ec_files

mkdir -p results/humann/kegg_compounds
for sample in $(ls results/humann/regroup/ec_files)
do
    bn=$(basename "$sample" _Abundance_RPKs.tsv)
    echo $bn

    mv results/humann/regroup/ec_files/$sample results/humann/regroup/ec_files/$bn
    humann \
    --input results/humann/regroup/ec_files/$bn \
    --output results/humann/kegg_compounds \
    --pathways-database ../../assets/compoundc \
    --remove-stratified-output \
    --threads 16
done

humann_join_tables \
--input results/humann/kegg_compounds \
--output results/humann/regroup/kegg_compounds.txt \
--file_name "_pathabundance.tsv"

# KEGG modules
humann_split_table \
--input results/humann/regroup/kegg_families.txt \
--output results/humann/regroup/kegg_files

mkdir -p results/humann/kegg_modules
for sample in $(ls results/humann/regroup/kegg_files)
do
    bn=$(basename "$sample" _Abundance_RPKs.tsv)
    echo $bn

    mv results/humann/regroup/kegg_files/$sample results/humann/regroup/kegg_files/$bn
    humann \
    --input results/humann/regroup/kegg_files/$bn \
    --output results/humann/kegg_modules \
    --pathways-database ../../assets/modulep \
    --remove-stratified-output \
    --threads 16
done

humann_join_tables \
--input results/humann/kegg_modules \
--output results/humann/regroup/kegg_modules.txt \
--file_name "_pathabundance.tsv"
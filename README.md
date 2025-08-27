## Microbial metabolic pathways guide response to immune checkpoint blockade therapy

This repository contains the code pertinent to our publication. Please refer to the Materials and Methods section of the article for details.
>Iris L. Mimpen, Thomas W. Battaglia, Miguel Parra Martinez, Catherine Toner-Bartelds, Laurien J. Zeverijn, Birgit S. Geurts, Karlijn Verkerk, Louisa R. Hoes, Allard W. J. van Renterghem, Michaël Noë, Ingrid Hofland, Annegien Broeks, Vincent van der Noort, Edwin C.A. Stigter, Can M.C. Gulersonmez, Boudewijn M.T. Burgering, Merel van Gogh, Marcel R. de Zoete, Hans Gelderblom, Krijn K. Dijkstra, Lodewyk F.A. Wessels, Emile E. Voest.
Microbial metabolic pathways guide response to immune checkpoint blockade therapy, Cancer Discovery (Accepted)

## Abstract
Studies have identified a link between specific microbiome-derived bacteria and immune checkpoint blockade (ICB) efficacy. However, these species lack consistency across studies and their immunomodulatory mechanisms remain elusive. To understand the influence of the microbiome on ICB response we studied its functional capacity. Using pan-cancer metagenomics data of ICB-treated patients, we showed that community-level metabolic pathways are stable across individuals, making them suitable to predict ICB response. We identified several microbial metabolic processes significantly associated with response, including the methylerythritol phosphate (MEP) pathway, which was associated with response and induced Vδ2 T cell-mediated anti-tumor responses in patient-derived tumor organoids. In contrast, riboflavin synthesis was associated with ICB resistance, and its intermediates induced mucosal-associated invariant T (MAIT) cell-mediated immune suppression. Moreover, gut metabolomics revealed that high riboflavin levels were linked to worse survival in patients with abundant intratumoral MAIT cells. Collectively, our results highlight the relevance of metabolite-mediated microbiome-immune cell crosstalk. 

## Analysis Overview and files. 
The analyses performed include:

- `drup-gut/`: Analysis of the DRUP-MSI metagenomic dataset. 
- `microbeICB`/: Analysis of the publicly available metagenomic datasets. 

## Pipeline
This analysis makes use of a Google Cloud (GCP) based Nextflow pipeline to process whole metagenome sequences with the Biobakery tools Humann and Metaphlan.
The pipeline can be found at: https://github.com/miparrama/biobakery-nf

## Data
Compositions of the data can be found:

## Contributors

- [Miguel Parra-Martinez](https://github.com/miparrama) - [@miparrama](https://github.com/miparrama)
- [Thomas W. Battaglia](https://github.com/twbattaglia) - [@twbattaglia](https://github.com/twbattaglia)
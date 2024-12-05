# Analysis for "Elevated rates of molecular evolution genome-wide in mutualist legumes and rhizobia" 

Tia Harrison 
November 20 2024

This repository stores data and code for our analysis in the paper:

Harrison TL, Stinchcombe JR, Frederickson ME. "Elevated rates of molecular evolution genome-wide in mutualist legumes and rhizobia"

## Abstract 

Rates of molecular evolution vary greatly among even closely related species. Although theory predicts that antagonistic interactions between species increase rates of molecular evolution, predictions for how mutualism affects evolutionary rates are mixed. We compared rates of molecular evolution between 1) mutualistic and non-mutualistic legumes, 2) an independent set of symbiotic rhizobia and their non-symbiotic close relatives, and 3) symbiotic and non-symbiotic clades within Ensifer, a diverse genus of bacteria with various lifestyles. We assembled transcriptomes de novo for 12 legume species and calculated dN/dS ratios at orthologous genes in all species to determine if genes in mutualistic plants evolve faster or slower than in their non-mutualistic relatives. We also calculated dN/dS ratios in genes known to be important for symbiosis. We found that mutualists have higher rates of molecular evolution genome-wide compared to non-mutualist legumes, but this pattern did not hold in symbiosis genes. We next calculated dN/dS ratios in 14 bacteria species across the proteobacteria phylogeny that differ in whether they associate mutualistically with plants, using published data. In most pairs, symbiotic rhizobia show higher dN/dS values compared to their non-symbiotic relatives. Within a bacterial genus with many well-characterized mutualist species (Ensifer), we calculated dN/dS ratios in symbiotic and non-symbiotic clades and found that symbiotic lineages have higher rates of molecular evolution genome-wide, but not at genes on the symbiotic plasmid pSymB. Our results suggest that although mutualism between legumes and rhizobia is associated with elevated rates of molecular evolution genome-wide, symbiosis genes may be evolutionarily stagnant.

## Description of folders 
- Plants
  - dN/dS output data from plant transcriptomes analyses done in PAML
  - R code available for analysis of dN/dS ratios from plant genomes
  - dN/dS ratios from symbiosis genes in plant genomes
  - invasion data for plant species in dataset
 
- Rhizobia
  - dN/dS output data from species pairs of rhizobia and non-symbiotic bacteria (PAML analysis) 
  - dN/dS output data from Ensifer genus analysis done in PAML
  - R code available for analysis of dN/dS ratios of rhizobia and Ensifer genomes
  - dN/dS ratios from symbiosis genes in bacteria genomes
  - dN/dS ratios for Ensifer chromosome and pSymB plasmid 


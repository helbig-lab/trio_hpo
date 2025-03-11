# Phenotypic analysis of Trio exomes

This repo contains all the code for doing the analysis on this paper, "Phenotypic analysis of 11,125 trio exomes in neurodevelopmental disorders". 

For each patient, phenotypes were standardized using HPO (Human Phenotype Ontology) terms along with the filtered _de novo_ variants. 

Different steps and the scripts:

- [HPO data manipulation](https://github.com/helbig-lab/trio_hpo/blob/main/scripts/harmonize_data_hpo.R). This script checks the input data in HPO format and reports an overview of the dataset. We also compute Information Content (IC) for each HPO terms based on the frequency in the cohort.

- [Computation of similarity score](https://github.com/helbig-lab/trio_hpo/blob/main/scripts/compute_similarity_score.R). This script guides us in computing the similarity score based on all different algorithms. The output of this step will be a matrix with pair-wise similarity score.

- Filtering Genes and permutation analysis [tbd](tbd). This script focus on filtering de novo variants and counting the no. of individuals with each genetic etiology. For each genetic etiology, we compute median similarity score for each genes. For the permutation analysis, this script generates a random distribution of similarity scores for each n. The output of this file generates p-value for each genetic etiology based on similiraity score. We also compute the p-values based on _de novo_ burden using [denovolyzeR](https://github.com/jamesware/denovolyzeR). 

-  Clustering [tbd](tbd). This workflow focuses on finding internal clusters within the genes based on the pairwise similarity scores. We also look at the phenotypes that are contributing the different clusters.

-  Twin index [tbd](tbd). This workflow computes the twin indices for the genes based on the neighborhood rank. To create the volcano and phenograms, the code is also available in this script. 

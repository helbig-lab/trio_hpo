######## Permutation 
library(memoise)
library(tidyverse)
library(dplyr)
library(readr)
library(data.table)
library(ontologyIndex)
library(doParallel)
library(dqrng)

# input files ####
new_hpo_ancs <- read.csv("../outputs/hpo_ancs.csv",stringsAsFactors = F)
new_hpo_name <- read.csv("../outputs/hpo_def.csv",stringsAsFactors = F)
base <- read.csv("../outputs/base_new_hpo.csv",stringsAsFactors = F)
prop <- read.csv("../outputs/prop_new_hpo.csv",stringsAsFactors = F)
IC <- read.csv("../outputs/IC.csv",stringsAsFactors = F)
pairwise_similarity_mat <- read.csv("../outputs/sim_score.csv",stringsAsFactors = F)

# genes ####

all_genes <- base %>% select(id,gene) %>% unique %>% dplyr::count(gene) %>%
  dplyr::rename(Freq = n) %>%
  filter(Freq > 1) %>% filter(!is.na(gene)) 

pat_pairs = all_genes$Freq %>% unique

perm = 100000
for( i in 1:length(pat_pairs)) {
sim_pat_draw = function(sim_score, num_pats,perm)  {
  pat_vect <- names(sim_score)
  r_100k = rep(NA_real_,perm)
  r_100k_id <- NULL
  for(n in seq_along(r_100k)){
    IDs = dqsample(pat_vect,num_pats,replace=F)
    sub_sim = sim_score[ (rownames(sim_score) %in% IDs), (names(sim_score) %in% IDs)]
    r_100k[n] = median(unlist(sub_sim),na.rm=T)
    r_100k_id[n] <- paste(IDs,collapse=";")
    }
  return(list(r_100k,r_100k_id))
}

v_name <- paste("n",pat_pairs[i],"_perm",sep = "")

assign(v_name, sim_pat_draw(pairwise_similarity_mat,pat_pairs[i],perm)
       ,envir=globalenv())
}

gene_df <- function(gene_name) {
  gene_mat_df <- pairwise_similarity_df %>% filter(pat1 %in% (base %>% filter(gene == gene_name) %>% pull(id) %>% unique),
                                    pat2 %in% (base %>% filter(gene == gene_name) %>% pull(id) %>% unique))
  gene_mat_df$gene <- gene_name
  return(gene_mat_df)
}


pairwise_gene <- do.call(rbind,lapply(1:nrow(all_genes), function(xx) gene_df(all_genes$gene[xx])))

write.csv(pairwise_gene,"../outputs/pairwise_gene.csv",row.names = F)

gene_count <- as.data.frame(matrix(ncol=5,nrow=length(all_genes$gene)))
names(gene_count) <- c("gene","n_pats","pairs","median_sim","p_median")
gene_count <- gene_count %>% mutate(gene = all_genes$gene)

#Finding the average, median, and mode similarity,for each gene in the cohort
for (i in 1:nrow(gene_count)) {
  name_x <- pairwise_gene %>% filter(gene == gene_count$gene[i])
  gene_count[i,c('n_pats')] <- c(name_x$pat1, name_x$pat2) %>% unique %>% length()
  gene_count[i,c('pairs')] <- nrow(name_x)
  gene_count[i,c('median_sim')] <- median(name_x$sim_score)
}


get_pval <- function(sample_size,sim_score_med,perm){
  #builds matrix with nx random draws with n random individuals, example SCN1A with 21 pairs
  
  v_name <- paste("n",sample_size,"_perm",sep = "")
  ax <- get(v_name)
  return( pnorm(sim_score_med,mean = mean(ax[[1]]),sd = sd(ax[[1]]),lower.tail = F) )
}

gene_count$p_median <- lapply(1:nrow(gene_count), function(xx) get_pval(gene_count$n_pats[xx],gene_count$median_sim[xx],perm)) %>% unlist()

write.csv(gene_count,"../outputs/gene_similarity_score.csv",row.names = F)


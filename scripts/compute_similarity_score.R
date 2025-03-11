######## Similarity methods
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



# similarity score ####

path <- new_hpo_ancs %>%
  tidyr::separate_rows(ancs, sep = ";") %>%
  filter(ancs != "") %>%
  mutate(Counter1 = ancs) %>%
  rename(Term = ID, HPO1 = ancs)


# pairwise similarity calculation ####

resnik_mod <- function(pat1, pat2){
  mica <- function(hpo1, hpo2){
    
    path1_unique <- path %>% dplyr::filter(Term == hpo1)
    # path1_unique <- path[which(path$Term == hpo1),]
    
    path2_unique <- path %>% dplyr::filter(Term == hpo2)
    # path2_unique <- path[which(path$Term == hpo2),]
    
    joint1 <- path1_unique %>% inner_join(path2_unique, by = 'Counter1')
    
    joint2 <- joint1 %>% left_join(IC, by = c('Counter1' = 'HPO'))
    
    mica_ic <- joint2$prop.IC %>% max
    
    
    return(mica_ic)
    
  }
  
memo_mica <- memoise(mica)

  hpo_pat1 <- base %>% dplyr::filter(id == pat1) %>% unique
  hpo_pat2 <- base %>% dplyr::filter(id == pat2) %>% unique
  
  #create data frame with HPO of pat1 in x, HPO of pat2 in y
  
  x_length <- length(unique(hpo_pat1$new_hpo_id))
  y_length <- length(unique(hpo_pat2$new_hpo_id))
  
  ic_matrix <- as.data.frame(matrix(ncol=x_length, nrow=y_length))
  names(ic_matrix) <- unique(hpo_pat1$new_hpo_id)
  rownames(ic_matrix) <- unique(hpo_pat2$new_hpo_id)
  
  for (i in 1:y_length){
    for(j in 1:x_length){
      ic_matrix[i,j] <- memo_mica(hpo_pat2$new_hpo_id[i], hpo_pat1$new_hpo_id[j])
    }
  }
  max_col <- apply(ic_matrix,2,max)
  max_row <- apply(ic_matrix,1,max)
  max_complete <- sum(max_col,max_row)/2
  
  return(max_complete)
}

resnik_min <- function(pat1, pat2){
  mica <- function(hpo1, hpo2){
    
    path1_unique <- path %>% dplyr::filter(Term == hpo1)
    # path1_unique <- path[which(path$Term == hpo1),]
    
    path2_unique <- path %>% dplyr::filter(Term == hpo2)
    # path2_unique <- path[which(path$Term == hpo2),]
    
    joint1 <- path1_unique %>% inner_join(path2_unique, by = 'Counter1')
    
    joint2 <- joint1 %>% left_join(IC, by = c('Counter1' = 'HPO'))
    
    mica_ic <- joint2$prop.IC %>% max
    
    
    return(mica_ic)
    
  }
  
  minimal_set_sg <- function(ancs,terms){
    redundant <- unlist(use.names = FALSE, lapply(terms, function(x) setdiff( ancs %>% filter(ID == x) %>% pull(ancs) %>% str_split(.,";") %>% unlist,
                                                                              x)))
    setdiff(terms, redundant)
  }
  
  
  memo_mica <- memoise(mica)
  
  hpo_pat1 <- base %>% dplyr::filter(id == pat1) %>% unique
  hpo_pat2 <- base %>% dplyr::filter(id == pat2) %>% unique
  
  hpo_pat1 <- minimal_set_sg(new_hpo_ancs,hpo_pat1$new_hpo_id)
  hpo_pat2 <- minimal_set_sg(new_hpo_ancs,hpo_pat2$new_hpo_id)
  
  #create data frame with HPO of pat1 in x, HPO of pat2 in y
  
  x_length <- length(unique(hpo_pat1$new_hpo_id))
  y_length <- length(unique(hpo_pat2$new_hpo_id))
  
  ic_matrix <- as.data.frame(matrix(ncol=x_length, nrow=y_length))
  names(ic_matrix) <- unique(hpo_pat1$new_hpo_id)
  rownames(ic_matrix) <- unique(hpo_pat2$new_hpo_id)
  
  for (i in 1:y_length){
    for(j in 1:x_length){
      ic_matrix[i,j] <- memo_mica(hpo_pat2$new_hpo_id[i], hpo_pat1$new_hpo_id[j])
    }
  }
  max_col <- apply(ic_matrix,2,max)
  max_row <- apply(ic_matrix,1,max)
  max_complete <- sum(max_col,max_row)/2
  
  return(max_complete)
}

cube <- function(pat1, pat2){
  
  hpo_pat1 <- prop %>% dplyr::filter(id == pat1) %>% unique %>% left_join(IC %>% select(HPO,prop.IC))
  hpo_pat2 <- prop %>% dplyr::filter(id == pat2) %>% unique %>% left_join(IC %>% select(HPO,prop.IC))
  
  return(hpo_pat1 %>% inner_join(hpo_pat2,by=c("HPO","definition","prop.IC")) %>% pull(prop.IC) %>% sum(.,na.rm = T))
}


# compute patient matrixx ####
patients=unique(base$id)
pairwise_similarity_df <- t(combn(patients,2)) %>% data.frame %>% rename("pat1"="X1","pat2"="X2")

## using the resnik mod algorithm as default
pairwise_similarity_df$sim_score <- unlist(lapply(1:nrow(pairwise_similarity_df), function(xx) resnik_mod(pairwise_similarity_df$pat1[xx],pairwise_similarity_df$pat2[xx])))

pairwise_similarity_mat <- pairwise_similarity_df %>% pivot_wider(names_from = "pat2",values_from = "sim_score") %>% data.frame
row.names(pairwise_similarity_mat) <- pairwise_similarity_mat$pat1

names(pairwise_similarity_mat)[1] <- unique(pairwise_similarity_df$pat1)[!unique(pairwise_similarity_df$pat1) %in% unique(pairwise_similarity_df$pat2)]

pairwise_similarity_mat <- rbind(pairwise_similarity_mat,pairwise_similarity_mat[1,])
row.names(pairwise_similarity_mat)[nrow(pairwise_similarity_mat)] <- unique(pairwise_similarity_df$pat2)[!unique(pairwise_similarity_df$pat2) %in% unique(pairwise_similarity_df$pat1)]

pairwise_similarity_mat[,unique(pairwise_similarity_df$pat1)[!unique(pairwise_similarity_df$pat1) %in% unique(pairwise_similarity_df$pat2)]] <- NA
pairwise_similarity_mat[unique(pairwise_similarity_df$pat2)[!unique(pairwise_similarity_df$pat2) %in% unique(pairwise_similarity_df$pat1)],] <- NA
diag(pairwise_similarity_mat) = NA_real_


write.csv(pairwise_similarity_mat,"../outputs/sim_score.csv",row.names = F)


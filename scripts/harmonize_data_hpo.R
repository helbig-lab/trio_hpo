# packages ####

library(tidyverse)
library(dplyr)
library(readr)
library(data.table)
library(ontologyIndex)
library(doParallel)
library(dqrng)

# OBO File ####
### get new HPO obo file
new_hpo <- get_OBO("../files/hp.obo") ### you can get this file from hpo.jax.org

## create an id and definition data frame
new_hpo_name <- data.frame(new_hpo$name)
new_hpo_name$id <- row.names(new_hpo_name)
row.names(new_hpo_name) <- NULL
names(new_hpo_name)[1] = "definition"

## exclude obsolete terms
new_hpo_name <- new_hpo_name %>% filter(!id %in% (new_hpo_name %>% filter(grepl("obsolete",definition)) %>% pull(2)))

# create ancestor data frame
ancestors <- do.call(rbind,lapply( 1:length(new_hpo$ancestors), function(xx) paste(new_hpo$ancestors[[xx]],collapse = ";")  )) %>% as.data.frame()
terms <- do.call(rbind,lapply( 1:length(new_hpo$ancestors), function(xx) names(new_hpo$ancestors[xx])  )) %>% as.data.frame()
names(terms) = "ID"
names(ancestors)= "ancs"

new_hpo_ancs <- cbind(terms,ancestors)
new_hpo_ancs$ID <- as.character(new_hpo_ancs$ID)
new_hpo_ancs$ancs <- as.character(new_hpo_ancs$ancs)

new_hpo_ancs <- new_hpo_ancs %>% filter(ID %in% new_hpo_name$id)

write.csv(new_hpo_ancs,"../outputs/hpo_ancs.csv",row.names = F)
write.csv(new_hpo_name,"../outputs/hpo_def.csv",row.names = F)


# Phenotypic data ####
bc_base <- read_csv("../files/clinical_data.csv")

bc_base <- bc_base %>% left_join(new_hpo_name, by =c("HPO"="id")) 

## create a new prop 
bc_prop <- bc_base %>% 
  left_join(new_hpo_ancs,by=c("new_hpo_id"="ID")) %>%
  select(famID, ancs) %>%
  separate_rows(ancs, sep=";") %>%
  rename(HPO = ancs) %>%
  filter(!is.na(HPO)) %>%
  left_join(new_hpo_name,by=c("HPO"="id")) %>% 
  unique

write.csv(bc_base,"../outputs/bc_base_new_hpo.csv",row.names = F)
write.csv(bc_prop,"../outputs/bc_prop_new_hpo.csv",row.names = F)

# Create Information Content ####

base_ct <- bc_base %>%
  count(new_hpo_id) %>%
  mutate(freq_base = n/length(unique(bc_base$famID))) %>%
  mutate(base.IC = -log2(freq_base)) %>%
  select(-n) %>% 
  rename(HPO=new_hpo_id)

prop_ct <- bc_prop %>%
  count(HPO) %>%
  mutate(freq_prop = n/length(unique(bc_prop$famID))) %>%
  mutate(prop.IC = -log2(freq_prop)) %>%
  select(-n)

IC <- base_ct %>%
  full_join(prop_ct) %>%
  left_join(new_hpo_name %>% unique, by =c("HPO"="id"))

IC[is.na(IC)] <- 0

write.csv(IC,"../outputs/bc2_IC_v2.csv",row.names = F)

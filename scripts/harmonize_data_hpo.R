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
base <- read_csv("../files/clinical_data.csv")

base <- base %>% left_join(new_hpo_name, by =c("new_hpo_id"="id")) 

## create a new prop 
prop <- base %>% 
  left_join(new_hpo_ancs,by=c("new_hpo_id"="ID")) %>%
  select(id, ancs) %>%
  separate_rows(ancs, sep=";") %>%
  rename(HPO = ancs) %>%
  filter(!is.na(HPO)) %>%
  left_join(new_hpo_name,by=c("HPO"="id")) %>% 
  unique

write.csv(base,"../outputs/base_new_hpo.csv",row.names = F)
write.csv(prop,"../outputs/prop_new_hpo.csv",row.names = F)

# Create Information Content ####

base_ct <- base %>%
  count(new_hpo_id) %>%
  mutate(freq_base = n/length(unique(base$id))) %>%
  mutate(base.IC = -log2(freq_base)) %>%
  select(-n) %>% 
  rename(HPO=new_hpo_id)

prop_ct <- prop %>%
  count(HPO) %>%
  mutate(freq_prop = n/length(unique(prop$id))) %>%
  mutate(prop.IC = -log2(freq_prop)) %>%
  select(-n)

IC <- base_ct %>%
  full_join(prop_ct) %>%
  left_join(new_hpo_name %>% unique, by =c("HPO"="id"))

IC[is.na(IC)] <- 0

write.csv(IC,"../outputs/IC.csv",row.names = F)

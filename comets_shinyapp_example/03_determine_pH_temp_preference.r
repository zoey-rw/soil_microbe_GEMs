library(tidyverse)
library(broom)
library(data.table)

# First run 01_merge_env_data.r

# Read in sample-level abundances with biome info
sample_abundance <- fread("./sample_abundance_data.csv")


# Function to find the peak pH/temperature value
loess_maximum <- function(x, y, tol = .Machine$double.eps^0.5){
    model <- loess(y ~ x,span = 0.8)
    yfit <- model$fitted
    inx <- which(abs(yfit - max(yfit)) < tol)[[1]]
    list(x = x[inx], y.fitted = yfit[inx], ix = inx)
}

# Get pH of maximum abundance trends per taxon
pH_max = sample_abundance %>% 
    group_by(taxonomy_id, taxon) %>% 
    filter(!is.na(soilInCaClpH)) %>% 
    summarize(pH_preference = loess_maximum(soilInCaClpH, percentage)[[1]]) %>% 
    rename("Species of interest" = taxon)

# Get temp of maximum abundance trends per taxon
temp_max = sample_abundance %>% 
    group_by(taxonomy_id, taxon) %>% 
    filter(!is.na(soilTemp)) %>% 
    summarize(temperature_preference = loess_maximum(soilTemp, percentage)[[1]])  %>% 
    rename("Species of interest" = taxon)

# Combine and save
organism_preference = full_join(temp_max, pH_max)
write_csv(organism_preference, "./intermediate_data/organism_pH_temp_preference.csv")


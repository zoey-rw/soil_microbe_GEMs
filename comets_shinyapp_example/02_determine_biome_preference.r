
# First run 01_merge_env_data.r

library(tidyverse)
library(data.table)
library(broom)

# For future_map
library(segmented)
library(furrr)
plan(strategy = "multisession")


#setwd("/projectnb/frpmars/soil_microbe_db/")

# Function to assign biome preference to soil organisms based on prevalence 
assign_biome_presence = function(taxon_df, prevalence_cutoff=.05) {
    
    # get counts of samples per biome 
    sample_counts = taxon_df %>% group_by(nlcdClass) %>% tally(name = "biome_count")
    
    # get counts of non-zero abundances per biome
    presence_counts = taxon_df %>% 
        mutate(present = ifelse(percentage > 0, 1, 0)) %>% 
        group_by(nlcdClass, present) %>% tally(name = "presence_count")
    
    # ensure all combinations are present
    presence_counts_expanded = presence_counts %>% 
        ungroup %>%  
        expand(nlcdClass, present)
    # set NA values to zero
    presence_counts_expanded = merge(presence_counts, presence_counts_expanded, all=T) %>% 
        replace_na(list(presence_count = 0))
    
    # Get prevalence (percent of biome samples with nonzero abundance)
    df_merged = merge(presence_counts_expanded, sample_counts, all=T) %>% 
        filter(present == 1) %>% 
        mutate(prevalence = presence_count/biome_count,
               present_in_biome = ifelse(prevalence > prevalence_cutoff, 1, 0))
    
    df_out = df_merged  %>% 
        dplyr::select(nlcdClass, biome_count, prevalence, present_in_biome)
    return(df_out)
}

# Read in sample-level abundances with biome info
sample_abundance <- fread("./sample_abundance_data.csv")

# Nest by taxon before running function
sample_abundance_nest = sample_abundance %>% 
    group_by(db_name, taxon, taxonomy_id, lineage) %>% nest()

#species_of_interest = unique(merged_df_nest$taxon) %>% sample(1000)
#data_nest_biome <- merged_df_nest %>% filter(taxon %in% species_of_interest) 

# Do this in chunks otherwise the map() function freezes
df_biome_pt1 <- sample_abundance_nest[1:5000,] %>% 
    mutate(biome_presence = future_map(data, assign_biome_presence)) 

df_biome_pt2 <- sample_abundance_nest[5001:10000,] %>% 
    mutate(biome_presence = future_map(data, assign_biome_presence)) 

df_biome_pt3 <- sample_abundance_nest[10001:nrow(sample_abundance_nest),] %>% 
    mutate(biome_presence = future_map(data, assign_biome_presence)) 

# Combine and unlist
df_biome = data.table::rbindlist(list(df_biome_pt1, df_biome_pt2, df_biome_pt3)) %>% 
    dplyr::select(-data) %>% 
    unnest(cols = c(biome_presence)) %>% 
    ungroup 

# Create separate column for each biome's "present in biome" (true/false = 1/0)
df_biome_wide = df_biome %>% 
    dplyr::select(-c(biome_count, prevalence)) %>% 
    pivot_wider(names_from = nlcdClass, values_from = present_in_biome)

#nrow(df_biome_wide) # 18576 taxa
write_csv(df_biome_wide, "./intermediate_data/organism_biome_preference.csv")




table(df_biome_wide$cultivatedCrops)
df_biome %>% group_by(nlcdClass, present_in_biome) %>% tally()

head(df_biome_pt1$biome_presence)

p3 = ggplot(merged_df %>% 
                dplyr::filter(taxon %in% species_of_interest_carbon[1:4]),
            aes(x = organicd13C, y = percentage, color=taxon)) +
    geom_point(alpha=.5, 
               position=position_jitter(width = .01, height=0), size=2, 
               show.legend = F) + 
    geom_smooth(show.legend = F) +
    facet_wrap(~taxon, scales = "free") + theme_bw(base_size = 22) + 
    scale_y_sqrt() + xlab("Soil organic carbon") + ylab("Microbial abundance")


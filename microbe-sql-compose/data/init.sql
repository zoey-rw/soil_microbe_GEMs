-- Table: organism_data_to_print
CREATE TABLE organism_data_to_print (
    "Row" SERIAL PRIMARY KEY,
    "Speciesofinterest" TEXT,
    "match_by" TEXT,
    "GEM_ID" TEXT,
    "Genomesource" TEXT,
    "Functional in COMETS" TEXT
);

-- Table: organism_data_to_subset
CREATE TABLE organism_data_to_subset (
    "ROW ID" SERIAL PRIMARY KEY,
    "GEM ID" TEXT,
    "Kingdom" TEXT,
    "Genus" TEXT,
    "Species.strain" TEXT,
    "Citation" TEXT,
    "Source" TEXT,
    "Filepath" TEXT,
    "Match criteria" TEXT,
    "Species of interest" TEXT,
    "Functional in COMETS?" TEXT,
    "taxonomy_id" INT,
    "lineage" TEXT,
    "db_name" TEXT,
    "cultivatedCrops" FLOAT,
    "deciduousForest" INT,
    "dwarfScrub" FLOAT,
    "emergentHerbaceousWetlands" FLOAT,
    "evergreenForest" INT,
    "grasslandHerbaceous" INT,
    "mixedForest" FLOAT,
    "pastureHay" FLOAT,
    "sedgeHerbaceous" FLOAT,
    "shrubScrub" INT,
    "woodyWetlands" FLOAT,
    "source" TEXT,
    "pH_preference" FLOAT,
    "temperature_preference" FLOAT,
    "genome_link" FLOAT,
    "GEM_ID" TEXT,
    "match_by" TEXT
);

-- Table: organism_taxonomy
CREATE TABLE organism_taxonomy (
    "Unnamed: 0" SERIAL PRIMARY KEY,
    "taxon" TEXT,
    "Kingdom" TEXT,
    "Phylum" TEXT,
    "Class" TEXT,
    "Order" TEXT,
    "Family" TEXT,
    "Genus" TEXT,
    "accession" FLOAT
);

-- Table: species_abundance_filt
CREATE TABLE species_abundance_filt (
    id SERIAL PRIMARY KEY,
    "name" TEXT,
    "taxonomy_id" INT,
    "percentage" FLOAT,
    "lineage" TEXT,
    "source" TEXT,
    "is_MAG" BOOLEAN,
    "taxid_lineage" TEXT,
    "genomicsSampleID" TEXT,
    "d15N" FLOAT,
    "organicd13C" FLOAT,
    "nitrogenPercent" FLOAT,
    "organicCPercent" FLOAT,
    "soilTemp" FLOAT,
    "soilMoisture" FLOAT,
    "soilInWaterpH" FLOAT,
    "soilInCaClpH" FLOAT,
    "latitude" FLOAT,
    "longitude" FLOAT,
    "elevation" FLOAT,
    "sampleTiming" TEXT,
    "nlcdClass" TEXT,
    "db_name" TEXT,
    "taxon" TEXT,
    "n_samples" INT
);

-- Load data into the tables
COPY organism_data_to_subset 
FROM '/docker-entrypoint-initdb.d/cleaned_organism_data_to_subset.csv' 
WITH (FORMAT csv, HEADER true, QUOTE '"', DELIMITER ',');
COPY organism_data_to_print FROM '/docker-entrypoint-initdb.d/organism_data_to_print.csv' WITH (FORMAT csv, HEADER true, QUOTE '"');
COPY organism_taxonomy FROM '/docker-entrypoint-initdb.d/cleaned_organism_taxonomy.csv' WITH CSV HEADER;
COPY species_abundance_filt (
    "name", 
    "taxonomy_id", 
    "percentage", 
    "lineage", 
    "source", 
    "is_MAG", 
    "taxid_lineage", 
    "genomicsSampleID", 
    "d15N", 
    "organicd13C", 
    "nitrogenPercent", 
    "organicCPercent", 
    "soilTemp", 
    "soilMoisture", 
    "soilInWaterpH", 
    "soilInCaClpH", 
    "latitude", 
    "longitude", 
    "elevation", 
    "sampleTiming", 
    "nlcdClass", 
    "db_name", 
    "taxon", 
    "n_samples"
)
FROM '/docker-entrypoint-initdb.d/cleaned_species_abundance_filt.csv'
WITH (FORMAT csv, HEADER true);


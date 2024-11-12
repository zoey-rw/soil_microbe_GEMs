# Use the official R image as a base
FROM rocker/shiny:4

# system libraries of general use
## install debian packages
RUN apt-get update -qq && apt-get -y --no-install-recommends install \
    libxml2-dev \
    libcairo2-dev \
    libsqlite3-dev \
    libmariadbd-dev \
    libpq-dev \
    libssh2-1-dev \
    unixodbc-dev \
    libcurl4-openssl-dev \
    libssl-dev \
    libharfbuzz-dev \
    libfribidi-dev \
    libfreetype6-dev \
    libpng-dev \
    libtiff5-dev \
    libjpeg-dev


## update system libraries
RUN apt-get update && \
    apt-get upgrade -y && \
    apt-get clean


# Install necessary R packages
COPY ./renv.lock /renv.lock

RUN Rscript -e 'install.packages("renv")'
RUN Rscript -e 'renv::restore()'


# Copy the app files into the image
COPY ./app.R /srv/shiny-server/app.R
COPY ./species_abundance_filt.csv /srv/shiny-server/species_abundance_filt.csv
COPY ./organism_data_to_subset.csv /srv/shiny-server/organism_data_to_subset.csv
COPY ./organism_data_to_print.csv /srv/shiny-server/organism_data_to_print.csv
COPY ./nlcd_key.csv /srv/shiny-server/nlcd_key.csv
COPY ./organism_taxonomy.csv /srv/shiny-server/organism_taxonomy.csv

# Expose the port the app runs on
EXPOSE 3838

# Run the Shiny app
CMD ["R", "-e", "shiny::runApp('/srv/shiny-server', host = '0.0.0.0', port = 3838)"]
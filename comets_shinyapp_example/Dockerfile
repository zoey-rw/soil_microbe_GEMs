# Use the official R image as a base
FROM rocker/shiny:4

# Set environment variables to reduce interactive prompts during installation
ENV DEBIAN_FRONTEND=noninteractive

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


# Copy the app files into the image. The R code reads from
# ./reference_data/ and ./intermediate_data/ — preserve that layout
# rather than the legacy flat copy that didn't match what app.R expects.
COPY ./app.R /srv/shiny-server/app.R
COPY ./reference_data /srv/shiny-server/reference_data
COPY ./intermediate_data /srv/shiny-server/intermediate_data

# Set permissions for the shiny-server directory
RUN chown -R shiny:shiny /srv/shiny-server

# Expose the Shiny app port
EXPOSE 3838

# Switch to the 'shiny' user to improve security
USER shiny

# Run the Shiny app
CMD ["R", "-e", "shiny::runApp('/srv/shiny-server', host = '0.0.0.0', port = 3838)"]
FROM rocker/rstudio:4.1

COPY . .

# Install R packages
RUN Rscript ./install_r_packages.r
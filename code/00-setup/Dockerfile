FROM rocker/rstudio:latest

COPY . .

RUN apt-get clean all && \
  apt-get update && \
  apt-get upgrade -y && \
  apt-get install -y libglpk40 \
  libz-dev

# Install R packages
RUN Rscript ./install_r_packages.r
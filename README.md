# cytokineFinder

## Local

- assumes docker desktop is installed, running, and you are logged in

### make docker image to set up env for data curation rmds
- if additional R packages required: add to install_r_packages.r

```bash
make build
```

### run rstudio instance

```bash
make run
```

> log into RStudio by navigating to [localhost](http://localhost:8787/) with user: rstudio and password: 123

### push to dockerhub

- modified DOCKERHUB_USERNAME to your dockerhub account name

```bash
make push
```

## HPC (Sockeye)

### clone repo to hpc project folder
- log into sockeye and go to project user folder

```bash
module load git
git clone https://github.com/CompBio-Lab/cytokineFinder.git
```

### pull docker image

```bash
cd multi_omics/data
make sockeye_pull
mkdir sif
mv cytokinefinder sif/
```

### download geo data
* the following examples is based on a certain GEO dataset (GSE92415)
* change to your user

```bash
cd sif
module load apptainer
apptainer shell cytokinefinder
R
setwd("/arc/project/st-singha53-1/jt1013/cytokineFinder/gse/")
source("gse71669.R")
quit()
exit
```
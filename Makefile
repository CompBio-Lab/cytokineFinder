DOCKERHUB_USERNAME=jeffreysjtang
IMAGE_VERSION=v0.1
IMAGE_NAME=cytokineFinder

# Docker 	
build:
	docker build -t $(DOCKERHUB_USERNAME)/$(IMAGE_NAME):$(IMAGE_VERSION) .

run:
	docker run --rm -it -p 8787:8787 -e PASSWORD=123 -v $(shell pwd):/home/rstudio/ $(DOCKERHUB_USERNAME)/$(IMAGE_NAME):$(IMAGE_VERSION)
	
push:
	docker push $(DOCKERHUB_USERNAME)/$(IMAGE_NAME):$(IMAGE_VERSION)

sockeye_pull:
	module load singularity; \
	singularity pull --name modc docker://$(DOCKERHUB_USERNAME)/$(IMAGE_NAME):$(IMAGE_VERSION)
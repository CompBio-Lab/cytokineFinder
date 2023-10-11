DOCKERHUB_USERNAME=jeffreysjtang
IMAGE_VERSION=v0.2
IMAGE_NAME=cytokinefinder

# Docker 	
build:
	docker build -t $(DOCKERHUB_USERNAME)/$(IMAGE_NAME):$(IMAGE_VERSION) .

# go to localhost:8787 for local instance of docker image. Username is rstudio, password is stated in the run command
run:
	docker run --rm -it -p 8787:8787 -e PASSWORD=123 -v $(shell pwd):/home/rstudio/ $(DOCKERHUB_USERNAME)/$(IMAGE_NAME):$(IMAGE_VERSION)
	
push:
	docker push $(DOCKERHUB_USERNAME)/$(IMAGE_NAME):$(IMAGE_VERSION)

sockeye_pull:
	module load singularity; \
	singularity pull --name cytokinefinder docker://$(DOCKERHUB_USERNAME)/$(IMAGE_NAME):$(IMAGE_VERSION)

remove:
	rm *[0-9]+*
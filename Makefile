REGISTRY ?= ghcr.io/philiplinden
IMAGE ?= cremons-et-al-2022
TAG ?= main

build:
	docker build . -t  $(REGISTRY)/$(IMAGE):$(TAG)

push: build
	docker push $(REGISTRY)/$(IMAGE):$(TAG)

pull:
	docker pull $(REGISTRY)/$(IMAGE):$(TAG)

run:
	docker run -v $(PWD):"/opt" $(REGISTRY)/$(IMAGE):$(TAG) repro/main.py

bacalhau: pull
	bacalhau docker run $(REGISTRY)/$(IMAGE):$(TAG) repro/main.py

IMAGE ?= cremons-et-al-2022
TAG ?= latest

build:
	docker build . -t $(IMAGE):$(TAG)

run:
	docker run -v $(PWD):"/opt" $(IMAGE):$(TAG) repro/main.py

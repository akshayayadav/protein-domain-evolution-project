# Protein Domain Evolution

## Goals

1. Creat a docker file to build a container which will not require any software installation to run.
2. Write a Snakemake file to generate a workflow and run the analysis in parallel.
3. Dockerize the project to make it reproducible.

## How to start

1. Launch an instance on Jetstream, using Ubuntu 18.04 Devel and Docker v1.22, with m1.xlarge (CPU: 24, Mem: 60 GB, Disk: 60 GB)
2. Fork this github repository with scripts (python scripts, snakefile, .sh) along with the dockerfile.
3. Create a new folder "data" in the directory, download two datasets (fasta and pfam) into ../data directory.

## Create a container

- 


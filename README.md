Protein Domain Evolution
===


## 1. Goals

If you are a total beginner to this, start here!

1. Creat a docker file to build a container which will not require any software installation to run.
2. Write a Snakemake file to generate a workflow and run the analysis in parallel.
3. Dockerize the project to make it reproducible.

## 2. How to start

1. Launch an instance on Jetstream, using Ubuntu 18.04 Devel and Docker v1.22, with m1.xlarge (CPU: 24, Mem: 60 GB, Disk: 60 GB)
2. Fork this github repository with scripts (python scripts, snakefile, .sh) along with the dockerfile.
3. Create a new folder "data" in the directory, download two datasets (fasta and pfam) into ../data directory.

## 3. Build a Docker container
### 3.1 Starting from Dockerfile
1. install make, perl, hmmer, pfamscan
```bash
FROM ubuntu:16.04
RUN apt-get update && \
    apt-get install -y wget build-essential make perl hmmer && \
    cd /root/ && \
    wget "http://ftp.ebi.ac.uk/pub/databases/Pfam/Tools/OldPfamScan/PfamScan1.5/PfamScan.tar.gz"
```
2. install python3 and libs to run the scripts
```
RUN apt-get update \
  && apt-get install -y python3-pip python3-dev
RUN pip3 install pandas
RUN pip3 install rpy2
RUN pip3 install scipy
RUN pip3 install sklearn
RUN pip3 install matplotlib
```
3. install snakemake as the workflow management system
```
RUN pip3 install snakemake
```
4. add scripts for data analysis inside the container
```
ADD scripts /usr/local/bin
# make the scripts executable
RUN chmod +x /usr/local/bin/* 
```
### 3.2 Making the snakemake workflow file

1. The analysis includes three steps. 
- assigning pfam protein domains to species fasta
- calculating the domain matrices (content, duplication, abundance, versatility) from domain assignments
- analyzing domain matrices for filtering out significantly evolving domains

![](https://i.imgur.com/NXAbrww.jpg)

2. The snakefile 
- 1st step, run pfamscan.pl 
```
rule pfamscan:
```
- 2nd step, run python scripts to get matrices
```
rule domain_content_matrix:
rule domain_duplication_matrix:
rule domain_abundance_matrix:
rule domain_versatility_matrix:
```
- 3rd step, run python scripts to analyze the matrices

```

```
### 3.3 Build the container using Docker

Once the dockerfile and snakefile are ready, build the container using virtual machine as:

```bash
docker build -t domainevolution -f Docker/Dockerfile .
```
Use ```docker images``` to check the built images

### 3.4 Create and run a writeable container layer over the built image

Since the data directory is not built into the container, we need to bind mount a volume with the data directory into the container. 

```
docker run -v /home/$USERNAME/protein-domain-evolution-project/data:/data domainevolution run_analysis.sh
```

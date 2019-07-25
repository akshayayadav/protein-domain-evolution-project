Protein Domain Evolution
===


## 1. Goals

If you are a total beginner to this, start here!

1. Creat a docker file to build a container which will not require any software installation to run.
2. Write a Snakemake file to generate a workflow and run the analysis in parallel.
3. Dockerize the project to make it reproducible.

## 2. How to start

1. Launch an instance on Jetstream, using Ubuntu 18.04 Devel and Docker v1.22, with m1.xlarge (CPU: 24, Mem: 60 GB, Disk: 60 GB)
   try 
   ```
   #check it matches the username, if not, change any $USERNAME into the username
   echo $USER 
   ```
   ssh to the VM using
   ```
   # get the username and IP address
   ssh $USER@xxx.xxx.xxx.xxx
   ```
2. Fork this github repository with scripts (python scripts, snakefile, .sh) along with the dockerfile.
```
git clone https://github.com/cyber-carpentry/Group5-protein-domain-evolution-project.git
```
3. Download the dataset (fasta, pfam, species.label; test data available) 
```https://drive.google.com/file/d/1yr3_NfQ6lpcGGN1tzJnb8RIBaEdAaLmK/view?usp=sharing```
copy into the project directory.
  $USER is the username showed in echo $USER in VM. 
```
scp -r <download dir>/test_data $USER@xxx.xxx.xxx.xx:/home/$USER/Group5-protein-domain-evolution-project/
```
then upzip it using 
```
unzip test_data.zip
```
## 3. Build a Docker container
### 3.1 Starting from Dockerfile (explanation)
1. Install make, perl #v5.22.1, hmmer, pfamscan
```bash
FROM ubuntu:16.04
RUN apt-get update && \
    apt-get install -y wget build-essential make perl hmmer && \
    cd /root/ && \
    wget "http://ftp.ebi.ac.uk/pub/databases/Pfam/Tools/OldPfamScan/PfamScan1.5/PfamScan.tar.gz"
```
2. Install python3 and libs to run the scripts
```
RUN apt-get update \
  && apt-get install -y python3-pip python3-dev  #Version:Python 3.5.2
RUN pip3 install pandas #Version:0.24.2
RUN pip3 install rpy2   #Version:3.0.5
RUN pip3 install scipy  #Version:1.3.0
RUN pip3 install sklearn  #Version:0.21.2
RUN pip3 install matplotlib #Version:3.0.3
```
3. Install snakemake as the workflow management system
```
RUN pip3 install snakemake #Version:5.5.4
```
4. Add scripts for data analysis inside the container
```
ADD scripts /usr/local/bin
# make the scripts executable
RUN chmod +x /usr/local/bin/* 
```
### 3.2 Making the snakemake workflow file (explanation)

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
rule domain_content_matrix_analysis:
rule domain_duplication_matrix_analysis:
rule domain_abundance_matrix_analysis:
rule domain_versatility_matrix_analysis:
```
### 3.3 Build the container using Docker (hands on)

Once the dockerfile and snakefile are ready, build the container using virtual machine as:

```bash
docker build -t domainevolution -f Docker/Dockerfile .
```
Use ```docker images``` to check the built images

### 3.4 Create and run a writeable container layer over the built image (hands on)

Since the data directory is not built into the container, we need to bind mount a volume with the data directory into the container. 

```
docker run -v /home/$USER/Group5-protein-domain-evolution-project/test_data:/data domainevolution run_analysis.sh -c 24
```
The number "24" gives the number of cores passed to snakemake to run the analysis.


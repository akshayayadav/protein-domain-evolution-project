Protein domain evolution analysis pipeline
===

Protein domains are independent sections of protein sequences that can have function distint functions. One of the major ways in which proteins can evovle is through domain insertion/delection/duplication. In this project we will attempt to build an analysis pipeline that will take in 2 groups of species proteomes and find differences in domain compositions between the two groups. The whole pipeline will be packaged inside a docker container which can executed  on any given data in any machine environment.

## Inputs
* Proteome fasta file for each species. The number of fasta files will depend of which two groups of species the user decides to compare and how many species the user the user wants in each group. It is recommended that there should be atleast 10 species per group to obtain significant results. These can be collected from following databases.
    * [Phytozome](https://phytozome.jgi.doe.gov/pz/portal.html) for plant species
    * [Uniprot](https://www.uniprot.org/proteomes/) for any species

* The HMM file containing Pfam domains from  [Pfam](https://pfam.xfam.org/) database which contains registry of all the domains found in all the organisms. The HMM file must be processed using hmmpress program to create a HMM database. For more details on how to use the hmmpress tool please see the HMMER [user manual](http://hmmer.org/documentation.html).

* A file with two columns `species` and `species_label`. The `species` column contains fasta file names of individual species append by string "pfamscan". The `species_label` column contains labels (0 or 1) classifying the species in different groups.


## Execute the container using user defined datasets
   - Pull the docker image from https://hub.docker.com/r/akshayayadav/protein-domain-evolution-project . 
   ```
   docker pull akshayayadav/protein-domain-evolution-project
   ```
   - Since the data directory is not built into the container, you need to bind mount a volume with the data directory into the container. 

   ```
   docker run -v <path to the data directory>:/data akshayayadav/protein-domain-evolution-project run_analysis.sh -c 10
   ```
   The number "10" gives the number of cores passed to snakemake to run the analysis.



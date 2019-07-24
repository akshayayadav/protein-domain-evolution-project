#!/bin/bash

mkdir /data/pfamscan_results
mkdir /data/matrix_results
snakemake -s /usr/local/bin/snakefile



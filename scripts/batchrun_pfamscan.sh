#!/bin/bash

mkdir /data/pfamscan_results

for fasta_file in /data/fasta/* ;do
	fasta_file_name=$(basename -- "$fasta_file")
	 pfam_scan.pl -fasta /data/fasta/$fasta_file_name -dir /data/pfam_database -outfile /data/pfamscan_results/$fasta_file_name.pfamscan -e_seq 1e-5 -e_dom 1e-5
done


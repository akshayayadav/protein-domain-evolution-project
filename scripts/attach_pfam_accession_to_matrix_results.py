#!/usr/bin/python3

import sys
import pandas as pd

def read_pfam_domainID_accession_mapping_file(domainID_accession_mapping_fileName):
	domainID_accession_mapping_matrix = pd.read_csv(domainID_accession_mapping_fileName, sep=",", header=None, names=["domainID","domainAcc"])
	return(domainID_accession_mapping_matrix)

def read_results_matrix(results_matrix_fileName):
	results_matrix = pd.read_csv(results_matrix_fileName, sep=",", header=0)
	return(results_matrix)

def join_results_matrix_domainID_accession_mapping_matrix(results_matrix, domainID_accession_mapping_matrix, joined_matrix_outfile):
	joined_matrix = results_matrix.join(domainID_accession_mapping_matrix.set_index('domainID'), on='domainID')
	joined_matrix.to_csv(joined_matrix_outfile, index=False)

######################################################################################################################################
results_matrix_fileName = sys.argv[1]
domainID_accession_mapping_fileName = sys.argv[2]
joined_matrix_outfile = sys.argv[3]

domainID_accession_mapping_matrix = read_pfam_domainID_accession_mapping_file(domainID_accession_mapping_fileName)
results_matrix = read_results_matrix(results_matrix_fileName)
join_results_matrix_domainID_accession_mapping_matrix(results_matrix, domainID_accession_mapping_matrix, joined_matrix_outfile)

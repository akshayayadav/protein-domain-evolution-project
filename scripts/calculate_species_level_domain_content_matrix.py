#!/usr/bin/env python3

import re
import sys
import os
import pandas as pd

def get_species_domain_dict(species_pfamscan_files_dirName):
	species_domain_dict = {}
	for species_pfamscan_fileName in os.listdir(species_pfamscan_files_dirName):
		species_pfamscan_file = open(species_pfamscan_files_dirName+"/"+species_pfamscan_fileName, "r")
		for line in species_pfamscan_file:
			line = line.rstrip()
			if(re.match(r'^\#', line) or line == ""):
				continue
			linearr = re.split(r'\s+', line)
			domain_name = linearr[6]
			if(species_pfamscan_fileName in species_domain_dict):
				species_domain_dict[species_pfamscan_fileName][domain_name]=1
			else:
				species_domain_dict[species_pfamscan_fileName]={}
				species_domain_dict[species_pfamscan_fileName][domain_name]=1

		species_pfamscan_file.close()
	return(species_domain_dict)

def get_domain_feature_vector(species_domain_dict):
	domain_feature_vector=list()
	domain_feature_dict = {}
	for species in species_domain_dict:
		for domain in species_domain_dict[species]:
			domain_feature_dict[domain] = 1
	
	domain_feature_vector = list(domain_feature_dict.keys())
	return(domain_feature_vector)

def get_domain_feature_matrix_for_species(species_domain_dict, domain_feature_vector):
	domain_feature_vector.insert(0, 'species')
	domain_feature_matrix = pd.DataFrame(columns = domain_feature_vector)
	domain_feature_vector.pop(0)
	count=0
	for species in species_domain_dict:
		print (species)
		feature_vector = list()
		feature_vector.append(species)
		for domain in domain_feature_vector:
			if(domain in species_domain_dict[species]):
				feature_vector.append(1)
			else:
				feature_vector.append(0)
		domain_feature_matrix.loc[count] = list(feature_vector)
		count+=1

	return(domain_feature_matrix)

def remove_common_domains(domain_feature_matrix):
	domain_feature_matrix = domain_feature_matrix.loc[:,domain_feature_matrix.apply(pd.Series.nunique) != 1]
	return(domain_feature_matrix)		
		
def print_domain_content_matrix(domain_feature_matrix):
	#domain_feature_matrix.to_csv("/data/matrix_results/domain_content.matrix", index=False)
	domain_feature_matrix.to_csv(sys.argv[2], index=False)
###########################################################################################################################
#species_pfamscan_files_dirName = "/data/pfamscan_results/"
species_pfamscan_files_dirName = sys.argv[1]

print ("Calculating domain content matrix......")

species_domain_dict = get_species_domain_dict(species_pfamscan_files_dirName)
domain_feature_vector = get_domain_feature_vector(species_domain_dict)
domain_feature_matrix = get_domain_feature_matrix_for_species(species_domain_dict, domain_feature_vector)
domain_feature_matrix = remove_common_domains(domain_feature_matrix)
print_domain_content_matrix(domain_feature_matrix)

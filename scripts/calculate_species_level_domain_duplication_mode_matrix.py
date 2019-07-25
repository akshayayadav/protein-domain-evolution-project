#!/usr/bin/env python3

from __future__ import division
import re
import os
import sys
import pandas as pd
import numpy as np

def get_species_seqid_domarr_dict(species_pfamscan_files_dirName):
	species_seqid_domarr_dict = {}
	domain_feature_vector = list()
	domain_feature_dict = {}
	for species_pfamscan_fileName in os.listdir(species_pfamscan_files_dirName):
		species_pfamscan_file = open(species_pfamscan_files_dirName+"/"+species_pfamscan_fileName, "r")
		for line in species_pfamscan_file:
			line = line.rstrip()
			if(re.match(r'^\#', line) or line == ""):
				continue
			linearr = re.split(r'\s+', line)
			seqid = linearr[0]
			domain_name = linearr[6]
			domain_feature_dict[domain_name] = 1
			if(species_pfamscan_fileName in species_seqid_domarr_dict):
				if(seqid in species_seqid_domarr_dict[species_pfamscan_fileName]):
					species_seqid_domarr_dict[species_pfamscan_fileName][seqid].append(domain_name)
				else:
					species_seqid_domarr_dict[species_pfamscan_fileName][seqid]=list()
					species_seqid_domarr_dict[species_pfamscan_fileName][seqid].append(domain_name)
			else:
				species_seqid_domarr_dict[species_pfamscan_fileName]={}
				species_seqid_domarr_dict[species_pfamscan_fileName][seqid]=list()
				species_seqid_domarr_dict[species_pfamscan_fileName][seqid].append(domain_name)

	
	domain_feature_vector = list(domain_feature_dict.keys())
	return([species_seqid_domarr_dict, domain_feature_vector])

def get_species_domDuplicationCountarr_dict(species_seqid_domarr_dict):
	species_domDuplicationCountarr_dict = {}
	for species in species_seqid_domarr_dict:
		species_domDuplicationCountarr_dict[species]={}
		for seq in species_seqid_domarr_dict[species]:
			domDuplicationCount_dict = {}
			for domain in species_seqid_domarr_dict[species][seq]:
				if(domain in domDuplicationCount_dict):
					domDuplicationCount_dict[domain]+=1
				else:
					domDuplicationCount_dict[domain]=1
			
		
			for domain in domDuplicationCount_dict:
				if(domain in species_domDuplicationCountarr_dict[species]):
					species_domDuplicationCountarr_dict[species][domain].append(domDuplicationCount_dict[domain])
				else:
					species_domDuplicationCountarr_dict[species][domain] = list()
					species_domDuplicationCountarr_dict[species][domain].append(domDuplicationCount_dict[domain])
	
	return(species_domDuplicationCountarr_dict)

def get_domain_duplication_feature_matrix_for_species(species_domDuplicationCountarr_dict, domain_feature_vector):
	domain_feature_vector.insert(0, 'species')
	domain_feature_matrix = pd.DataFrame(columns = domain_feature_vector)
	domain_feature_vector.pop(0)
	count=0
	
	for species in species_domDuplicationCountarr_dict:
		print (species)
		feature_vector = list()
		feature_vector.append(species)
		for domain in domain_feature_vector:
			if(domain in species_domDuplicationCountarr_dict[species]):
				domain_duplication_count = species_domDuplicationCountarr_dict[species][domain]
				#feature_vector.append(calculate_mean_duplication(domain_duplication_count))
				feature_vector.append(calculate_mode_duplication(domain_duplication_count))
			else:
				feature_vector.append(0)
		domain_feature_matrix.loc[count] = list(feature_vector)
		count+=1

	return(domain_feature_matrix)

def calculate_mean_duplication(count_vector): 
    return sum(count_vector) / len(count_vector)

def calculate_mode_duplication(count_vector):
	return max(set(count_vector), key=count_vector.count)

def remove_constantly_duplicated_domains(domain_feature_matrix):
	domain_feature_matrix = domain_feature_matrix.loc[:,domain_feature_matrix.apply(pd.Series.nunique) != 1]
	return(domain_feature_matrix)

def remove_nonduplicated_domains(domain_feature_matrix):
	first_col = domain_feature_matrix['species']
	domain_feature_matrix.drop(domain_feature_matrix.columns[0], axis=1, inplace=True)
	domain_feature_matrix = domain_feature_matrix.loc[:,domain_feature_matrix.apply(detect_nonduplicated_domains)>1]
	domain_feature_matrix.insert(loc=0, column='species', value=first_col)
	return(domain_feature_matrix)

def detect_nonduplicated_domains(domain_col):
	domain_col = pd.Series(domain_col)
	return np.sum(domain_col.unique())

def print_domain_feature_matrix(domain_feature_matrix):
	#domain_feature_matrix.to_csv("/data/matrix_results/domain_duplication_mode.matrix", index=False)
	domain_feature_matrix.to_csv(sys.argv[2], index=False)
	
################################################################################################################################
#species_pfamscan_files_dirName = "/data/pfamscan_results/"
species_pfamscan_files_dirName = sys.argv[1]

print ("Calculating domain duplication mode matrix......")

species_seqid_domarr_dict, domain_feature_vector = get_species_seqid_domarr_dict(species_pfamscan_files_dirName)
species_domDuplicationCountarr_dict = get_species_domDuplicationCountarr_dict(species_seqid_domarr_dict)
domain_feature_matrix = get_domain_duplication_feature_matrix_for_species(species_domDuplicationCountarr_dict, domain_feature_vector)
domain_feature_matrix = remove_constantly_duplicated_domains(domain_feature_matrix)
domain_feature_matrix = remove_nonduplicated_domains(domain_feature_matrix)
print_domain_feature_matrix(domain_feature_matrix)
#print domain_feature_matrix

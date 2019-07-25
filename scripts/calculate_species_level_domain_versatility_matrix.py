#!/usr/bin/env python3

from __future__ import division
import re
import os
import sys
import pandas as pd
import numpy as np
from math import log

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

def get_species_domvrstlty_dict(species_seqid_domarr_dict):
	species_domvrstlty_dict = {}
	for species in species_seqid_domarr_dict:
		species_domvrstlty_dict[species]={}
		for seq in species_seqid_domarr_dict[species]:
			dom_arr = species_seqid_domarr_dict[species][seq]
			dom1_dom2_dict = get_domain_adjacency_dict(dom_arr)

			for dom1 in dom1_dom2_dict:
				if not (dom1 in species_domvrstlty_dict[species]):
					species_domvrstlty_dict[species][dom1]={}
				for dom2 in dom1_dom2_dict[dom1]:
					species_domvrstlty_dict[species][dom1][dom2]=1
	
	return(species_domvrstlty_dict)
					
def get_domain_adjacency_dict(dom_arr):
	dom1_dom2_dict = {}
	if(len(dom_arr)<2):
		return(dom1_dom2_dict)
	for i in range(0, len(dom_arr)):
		if not (dom_arr[i] in dom1_dom2_dict):
			dom1_dom2_dict[dom_arr[i]]={}
		
		if(i==0):
			#if(dom_arr[i] != dom_arr[i+1]):
			dom1_dom2_dict[dom_arr[i]][dom_arr[i+1]]=1
		elif(i==len(dom_arr)-1):
			#if(dom_arr[i] != dom_arr[i-1]):		
			dom1_dom2_dict[dom_arr[i]][dom_arr[i-1]]=1
		else:
			#if(dom_arr[i] != dom_arr[i+1]):
			dom1_dom2_dict[dom_arr[i]][dom_arr[i+1]]=1
			#if(dom_arr[i] != dom_arr[i-1]):
			dom1_dom2_dict[dom_arr[i]][dom_arr[i-1]]=1

	return(dom1_dom2_dict)


def get_domain_versatility_feature_matrix_for_species(species_domvrstlty_dict, domain_feature_vector):
	domain_feature_vector.insert(0, 'species')
	domain_feature_matrix = pd.DataFrame(columns = domain_feature_vector)
	domain_feature_vector.pop(0)
	count=0
	
	for species in species_domvrstlty_dict:
		print (species)
		feature_vector = list()
		feature_vector.append(species)
		for domain in domain_feature_vector:
			if(domain in species_domvrstlty_dict[species]):
				if(len(species_domvrstlty_dict[species][domain])>0):
					feature_vector.append(1/len(species_domvrstlty_dict[species][domain]))
				else:
					feature_vector.append(0)
			else:
				feature_vector.append(0)
		domain_feature_matrix.loc[count] = list(feature_vector)
		count+=1

	return(domain_feature_matrix)

def remove_constantly_versatile_domains(domain_feature_matrix):
	domain_feature_matrix = domain_feature_matrix.loc[:,domain_feature_matrix.apply(pd.Series.nunique) != 1]
	return(domain_feature_matrix)

def print_domain_feature_matrix(domain_feature_matrix):
	#domain_feature_matrix.to_csv("/data/matrix_results/domain_versatility.matrix", index=False)
	domain_feature_matrix.to_csv(sys.argv[2], index=False)

#########################################################################################################################################
#species_pfamscan_files_dirName = "/data/pfamscan_results/"
species_pfamscan_files_dirName = sys.argv[1]

print ("Calculating domain versatility matrix......")

species_seqid_domarr_dict, domain_feature_vector = get_species_seqid_domarr_dict(species_pfamscan_files_dirName)
species_domvrstlty_dict = get_species_domvrstlty_dict(species_seqid_domarr_dict)
domain_feature_matrix = get_domain_versatility_feature_matrix_for_species(species_domvrstlty_dict, domain_feature_vector)
domain_feature_matrix = remove_constantly_versatile_domains(domain_feature_matrix)
print_domain_feature_matrix(domain_feature_matrix)


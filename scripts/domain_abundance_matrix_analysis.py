#!/usr/bin/python3
import pandas as pd
import numpy as np
import operator
import sys
from sklearn.feature_selection import mutual_info_classif
from sklearn.preprocessing import normalize
from scipy.stats import ranksums
import matplotlib.pyplot as plt
from rpy2.robjects.packages import importr
from rpy2.robjects.vectors import FloatVector

def read_matrix_file(matrix_fileName, species_label_fileName):
	domain_matrix = pd.read_csv(matrix_fileName, sep=",", header=0)
	species_table = pd.read_csv(species_label_fileName, sep=",", header=0)
	domain_matrix = domain_matrix.join(species_table.set_index('species'), on='species')
	return(domain_matrix)

def get_data_matrix(domain_matrix):
	X_domain_matrix = domain_matrix.iloc[:,1:len(domain_matrix.columns)-1]
	return(X_domain_matrix)

def get_species_labels(domain_matrix):
	y_labels = domain_matrix['species_label']
	return(y_labels)

def normalize_domain_matrix(X_domain_matrix):
	X_domain_matrix_norm = normalize(X_domain_matrix, axis = 0, norm='l2')
	X_domain_matrix_norm = pd.DataFrame(X_domain_matrix_norm, columns=X_domain_matrix.columns)
	return(X_domain_matrix_norm)

def calculate_domain_mutual_info_scores(X_domain_matrix, y_labels):
	domain_mutual_info_scores = mutual_info_classif(X_domain_matrix, y_labels,random_state=8)
	domain_mutual_info_scores = dict(zip(X_domain_matrix.columns, domain_mutual_info_scores))
	sorted_domain_mutual_info_scores = sorted(domain_mutual_info_scores.items(), key=operator.itemgetter(1))
	sorted_domain_mutual_info_scores = pd.DataFrame(sorted_domain_mutual_info_scores)
	sorted_domain_mutual_info_scores.columns = ['domainID', 'MI-score']
	return(sorted_domain_mutual_info_scores)

def domain_ranksum_test(domain_col, y_labels):
    grp1 = domain_col.loc[y_labels==1]
    grp2 = domain_col.loc[y_labels==0]
    ranksum_test = ranksums(grp1,grp2)
    pvalue = ranksum_test[1]
    return(pvalue)

def run_signifance_tests_on_domain_cols(X_domain_matrix, y_labels):
	domain_signifance_test_results = pd.DataFrame(X_domain_matrix.apply(domain_ranksum_test, y_labels=y_labels))
	domain_signifance_test_results = domain_signifance_test_results.reset_index()
	domain_signifance_test_results.columns = ['domainID', 'p-value']
	return(domain_signifance_test_results)

def join_MI_results_with_significance_results(sorted_domain_mutual_info_scores, domain_signifance_test_results):
	domain_results = sorted_domain_mutual_info_scores.join(domain_signifance_test_results.set_index('domainID'), on='domainID')
	return(domain_results)

def multiple_testing_correction(domain_results):
	stats = importr('stats')
	p_adjust = stats.p_adjust(FloatVector(np.array(domain_results['p-value'])), method = 'fdr')
	domain_results['adjusted-p-values']=np.array(p_adjust)
	return(domain_results)

def print_significant_domains(domain_results, results_outfile):
	sig_domains = domain_results.loc[domain_results['adjusted-p-values']<0.05]
	sig_domains.to_csv(results_outfile, index=False)

def export_plot(domain_results, plot_outfile):
	plt.rcParams["figure.figsize"] = (30,20)
	plt.plot(domain_results['MI-score'],-1*np.log10(domain_results['adjusted-p-values']), 'ro')
	plt.axhline(y=-1*np.log10(0.05), color='g', linestyle='-')
	plt.xlabel('Mutual-Information Score', fontsize=30)
	plt.ylabel('-log10(fdr-adjusted-pvalues)', fontsize=30)
	ax = plt.gca()
	ax.tick_params(axis = 'both', which = 'major', labelsize = 20)
	ax.tick_params(axis = 'both', which = 'minor', labelsize = 20)
	plt.grid(True)
	plt.title("Protein domain abundance results", fontsize=50)
	plt.savefig(plot_outfile)
################################################################################################################################################
matrix_fileName = sys.argv[1]
species_label_fileName = sys.argv[2]
results_outfile = sys.argv[3]
plot_outfile = sys.argv[4]

print ("Analysing domain abundance matrix......")

domain_matrix = read_matrix_file(matrix_fileName, species_label_fileName)
X_domain_matrix = get_data_matrix(domain_matrix)
y_labels = get_species_labels(domain_matrix)
X_domain_matrix = normalize_domain_matrix(X_domain_matrix)
sorted_domain_mutual_info_scores = calculate_domain_mutual_info_scores(X_domain_matrix, y_labels)
domain_signifance_test_results = run_signifance_tests_on_domain_cols(X_domain_matrix, y_labels)
domain_results = join_MI_results_with_significance_results(sorted_domain_mutual_info_scores, domain_signifance_test_results)
domain_results = multiple_testing_correction(domain_results)
print_significant_domains(domain_results, results_outfile)
export_plot(domain_results, plot_outfile)

import os
import random
import pandas as pd
import numpy as np
from random import sample
from scipy.stats import ttest_ind

from loguru import logger
logger.info('Import OK')

correlation_path = 'results/tpemi_proteomics/correlations/max_per_protein_correlations.csv'
feature_path = 'results/tpemi_proteomics/protein_pathways/KEGG_summary.csv'
output_folder = 'results/tpemi_proteomics/pathway_correlations/'

if not os.path.exists(output_folder):
    os.mkdir(output_folder)


def sample(data):
    return [random.choice(data) for _ in range(len(data))]

def bootstrap_t_test(treatment, control, nboot = 1000, direction = "less"):
    ones = np.vstack((np.ones(len(treatment)),treatment))
    treatment = ones.conj().transpose()
    zeros = np.vstack((np.zeros(len(control)), control))
    control = zeros.conj().transpose()
    Z = np.vstack((treatment, control))
    tstat = np.mean(treatment[:,1])-np.mean(control[:,1])
    tboot = np.zeros(nboot)
    for i in range(nboot):
        sboot = sample(Z)
        sboot = pd.DataFrame(np.array(sboot), columns=['treat', 'vals'])
        tboot[i] = np.mean(sboot['vals'][sboot['treat'] == 1]) - np.mean(sboot['vals'][sboot['treat'] == 0]) - tstat
    if direction == "greater":
        pvalue = np.sum(tboot>=tstat-0)/nboot
    elif direction == "less":
        pvalue = np.sum(tboot<=tstat-0)/nboot
    else:
        logger.info('Enter a valid arg for direction')


# Read in correlation data
raw_data = pd.read_csv(f'{correlation_path}')
raw_data.drop([col for col in raw_data.columns.tolist() if 'Unnamed: ' in col], axis=1, inplace=True)
protein_corr = pd.melt(
    raw_data,
    id_vars=['Proteins'],
    value_vars=[col for col in raw_data.columns if col != 'Proteins'],
    var_name='Protein_B',
    value_name='corr_val',
)
protein_corr['key'] = [tuple(sorted([prot_A, prot_B])) for prot_A, prot_B in protein_corr[['Proteins', 'Protein_B']].values]
protein_corr.drop_duplicates(subset='key', inplace=True)
protein_corr = protein_corr[protein_corr['Proteins'] != protein_corr['Protein_B']]

# Read in interaction data
feature_data = pd.read_csv(f'{feature_path}')
feature_data.drop([col for col in feature_data.columns.tolist() if 'Unnamed: ' in col], axis=1, inplace=True)


# for a given complex, assign proteins inside versus outside the complex
pathway_correlations = protein_corr.copy()
pathway_names = dict(feature_data[['ko_pathway_term', 'name']].values)

# assign proteins into pathways
for pathway_name, df in feature_data.groupby('ko_pathway_term'):
    # pathway_name = 'mmu04120'
    pathway_proteins = feature_data[feature_data['ko_pathway_term'] == pathway_name]['Proteins'].unique().tolist()

    pathway_correlations['pathway_A'] = [1 if protein in pathway_proteins else 0 for protein in pathway_correlations['Proteins']]
    pathway_correlations['pathway_B'] = [1 if protein in pathway_proteins else 0 for protein in pathway_correlations['Protein_B']]
    pathway_correlations[pathway_name] = pathway_correlations[['pathway_A', 'pathway_B']].sum(axis=1)

pathway_correlations.drop(['pathway_A', 'pathway_B'], axis=1, inplace=True)

correlation_tests = []
for pathway in pathway_names:
    pathway_corr = pathway_correlations[['corr_val', 'key', pathway]].copy()
    pathway_corr = {membership: df['corr_val'].tolist() for membership, df in pathway_corr.groupby(pathway)}

    stat, pval = ttest_ind(pathway_corr[2], pathway_corr[1], equal_var=False)
    bpval = bootstrap_t_test(pathway_corr[2], pathway_corr[1], nboot=100, direction="less")

    correlation_tests.append(pd.DataFrame([pathway, np.mean(pathway_corr[2]), np.mean(pathway_corr[1]), np.mean(pathway_corr[2]), stat, pval, bpval]).T)
correlation_tests = pd.concat(correlation_tests)
correlation_tests.columns = ['pathway', 'members', 'external', 'nomembers', 'ttest_stat', 'ttest_pval', 'bootstrap_pval']
correlation_tests['pathway_name'] = correlation_tests['pathway'].map(pathway_names)

correlation_tests.to_csv(f'{output_folder}pathway_correlations.csv')

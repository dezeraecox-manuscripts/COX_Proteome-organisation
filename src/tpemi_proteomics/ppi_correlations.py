
import os
import random
import pandas as pd
import numpy as np
from scipy.stats import ttest_ind

from loguru import logger
logger.info('Import OK')

correlation_path = 'results/tpemi_proteomics/correlations/max_per_protein_correlations.csv'
feature_path = 'results/tpemi_proteomics/protein_interactions/STRING_protein_interactions.xlsx'
output_folder = 'results/tpemi_proteomics/ppi_correlations/'

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

    logger.info(f'The p-value is {pvalue}')


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
protein_corr['key'] = [f'{prot_a}_{prot_b}' if prot_a < prot_b else f'{prot_b}_{prot_a}' for prot_a, prot_b in protein_corr[['Proteins', 'Protein_B']].values]
protein_corr = protein_corr[protein_corr['Proteins'] != protein_corr['Protein_B']]


# Read in interaction data
feature_data = pd.read_excel(f'{feature_path}', sheet_name=None)
feature_data.update({key: df.drop([col for col in df.columns.tolist() if 'Unnamed: ' in str(col)], axis=1) for key, df in feature_data.items()})

# assign proteins which do vs do not interact
interactions = feature_data['all_cys_ppis'].copy()
interactions = interactions[interactions['score'] >= 0.4][['Protein_A', 'Protein_B']].copy().drop_duplicates()
interactions['key'] = [f'{prot_a}_{prot_b}' if prot_a < prot_b else f'{prot_b}_{prot_a}' for prot_a, prot_b in interactions[['Protein_A', 'Protein_B']].values]
interactions['corr_type'] = 1

ppi_correlations = protein_corr.copy()
ppi_correlations = pd.merge(ppi_correlations, interactions[['key', 'corr_type']], on='key', how='left').fillna(0)

mean_corr = ppi_correlations.groupby(['Proteins', 'corr_type']).mean().reset_index()
mean_corr['corr_val'] = mean_corr['corr_val'].abs()
mean_corr = pd.pivot_table(mean_corr, index='Proteins', columns='corr_type', values='corr_val').dropna() # remove proteins with no PPIs

# Are proteins that interact more correlated than those that dont?
tstat, pval = ttest_ind(
    mean_corr[0].tolist(), 
    mean_corr[1].tolist())
test_results = pd.DataFrame([tstat, pval], index=['tstat', 'pval']).T
test_results['mean_0'] = mean_corr.mean()[0]
test_results['mean_1'] = mean_corr.mean()[1]
test_results['std_0'] = mean_corr.std()[0]
test_results['std_1'] = mean_corr.std()[1]
test_results['count'] = len(mean_corr)

mean_corr.to_csv(f'{output_folder}ppi_correlation_vals.csv')
test_results.to_csv(f'{output_folder}ppi_correlation_summary.csv')

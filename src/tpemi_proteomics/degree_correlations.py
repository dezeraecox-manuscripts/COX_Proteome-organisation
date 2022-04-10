import os
import pandas as pd
import numpy as np
from scipy.stats.stats import spearmanr

from loguru import logger
logger.info('Import OK')

input_intersections = 'results/tpemi_proteomics/intersections/protein_intersection_degree.csv'
input_path = 'results/tpemi_proteomics/peptide_summary/peptide_summary.xlsx'
output_folder = 'results/tpemi_proteomics/degree_correlations/'

if not os.path.exists(output_folder):
    os.mkdir(output_folder)


# Read in summary data - proccessed same as PPIs
treatments = ['Celastrol', 'Novobiocin', 'MG132', 'Ver155008', 'Staurosporine']
# Read in summary data
proteins = pd.read_csv(f'{input_intersections}')
proteins.drop([col for col in proteins.columns.tolist() if 'Unnamed: ' in str(col)], axis=1, inplace=True)
proteins.sort_values('Degree')
degree_dict = pd.melt(proteins, id_vars=['Proteins', 'Degree'], value_vars=treatments, value_name='changed', var_name='treatment')
degree_dict = {protein: str(degree) if degree > 1 else treatment for protein, degree, treatment in degree_dict[degree_dict['changed']][['Proteins', 'Degree', 'treatment']].values}

# Read in quantitative data from summary, which contains the mean and max per protein
quant_data = pd.read_excel(f'{input_path}', sheet_name='changed_protein_summary')
quant_data.drop([col for col in quant_data.columns.tolist() if 'Unnamed: ' in col], axis=1, inplace=True)

proteins_max = quant_data.copy().dropna()
proteins_max = pd.pivot_table(
    proteins_max, values='max_log2_pval_thresh', index='Proteins', columns='treatment').reset_index()
# Select only proteins with at least one sign change 
proteins_max_change = proteins_max.replace(0, np.nan).dropna(how='all', subset=treatments)

# Map quant data onto degree information
proteins_max_change['degree_type'] = proteins_max_change['Proteins'].map(degree_dict)
proteins_max_change['direction'] = ['up' if all(val > 0 for val in vals[~np.isnan(vals)]) else (
    'down' if all(val < 0 for val in vals[~np.isnan(vals)]) else 'mixed') for vals in proteins_max_change[treatments].values]


# Fill na's with 0 --> no change
proteins_max_change.fillna(0, inplace=True)
proteins_max_change['degree_type'] = [int(
    degree) if degree not in treatments else 1 for degree in proteins_max_change['degree_type']]

proteins_max_change.to_csv(f'{output_folder}proteins_max_change.csv')


# Simple spearmans R
for_corr = pd.melt(
    proteins_max_change,
    id_vars=['Proteins', 'degree_type', 'direction'],
    value_vars=treatments,
    value_name='ratio',
    var_name='treatment')
for_corr.dropna(inplace=True)
from scipy.stats import spearmanr

spearmans_vals = {treatment: list(spearmanr(df['degree_type'], df['ratio'])) for treatment, df in for_corr.groupby('treatment')}

# Add complete dataset
spearmans_vals['complete'] = list(spearmanr(for_corr['degree_type'],for_corr['ratio']))
spearmans_vals = pd.DataFrame(dict([ (k,pd.Series(v)) for k,v in spearmans_vals.items() ])).T
spearmans_vals.columns = ['spearmans_r', 'spearmans_pval']

spearmans_vals.to_csv(f'{output_folder}degree_spearmans_correlations.csv')
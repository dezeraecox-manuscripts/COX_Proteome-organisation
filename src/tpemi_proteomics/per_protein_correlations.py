import os
import pandas as pd
import numpy as np

from loguru import logger

logger.info('Import OK')

input_folder = 'results/tpemi_proteomics/peptide_summary/'
output_folder = 'results/tpemi_proteomics/correlations/'

if not os.path.exists(output_folder):
    os.mkdir(output_folder)

treatments = ['Celastrol', 'Novobiocin', 'MG132', 'Ver155008', 'Staurosporine']

# Read in raw data
proteins = pd.read_excel(
    f'{input_folder}peptide_summary.xlsx', sheet_name='max_proteins')
proteins.drop([col for col in proteins.columns.tolist() if 'Unnamed: ' in col], axis=1, inplace=True)

proteins = proteins.replace(0, np.nan).dropna(how='all', subset=treatments).fillna(0) # select only proteins quantified in all treatments with at least one change

# calculate cross-correlation of all proteins
protein_corr = proteins.set_index('Proteins').T.corr(method='spearman')

# Calculate cross-correlation of all treatments
treatment_corr = proteins.set_index('Proteins').corr(method='spearman')

protein_corr.to_csv(f'{output_folder}max_per_protein_correlations.csv')
treatment_corr.to_csv(f'{output_folder}treatment_max_protein_correlations.csv')

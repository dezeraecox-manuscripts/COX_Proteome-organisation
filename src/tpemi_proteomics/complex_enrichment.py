
import os
import pandas as pd
import numpy as np

from utilities.statistical_tests import apply_fisher_test

from loguru import logger
logger.info('Import OK')

input_folder = 'results/tpemi_proteomics/peptide_summary/'
feature_path = 'results/tpemi_proteomics/protein_complexes/ontology_summary.xlsx'
background_path = 'results/tpemi_proteomics/peptide_normalisation/normalised_summary.xlsx'

output_folder = 'results/tpemi_proteomics/complex_enrichment/'

if not os.path.exists(output_folder):
    os.mkdir(output_folder)

# ----Prepare peptide database----
# Read in raw data
proteins = pd.read_excel(
    f'{input_folder}peptide_summary.xlsx', sheet_name='max_proteins')
proteins.drop([col for col in proteins.columns.tolist()
              if 'Unnamed: ' in col], axis=1, inplace=True)
# collect proteins changed in at least one treatment
proteins = proteins.replace(0, np.nan).dropna(how='all', subset=['MG132', 'Celastrol', 'Ver155008', 'Novobiocin', 'Staurosporine'])
proteins = proteins['Proteins'].unique().tolist()
proteins = pd.DataFrame([proteins], index=['change_in_one']).T

# Read in complex data
feature_data = pd.read_excel(f'{feature_path}', sheet_name=None)
feature_data.update({key: df.drop([col for col in df.columns.tolist() if 'Unnamed: ' in str(col)], axis=1) for key, df in feature_data.items()})
feature_data.keys()
complexes = feature_data['GO_complexes'].copy()
complex_names = dict(feature_data['GO_summary'][['go_id', 'name']].values)

# Read in background data
background_genes = pd.read_excel(f'{background_path}', sheet_name=None)
background = background_genes['cys_peptides']['Proteins'].unique().tolist()

# Test enrichment of individual complex ontology terms within changed proteins compared to total background
enrichments = []
for complex_name in complexes.columns:
    enrichment = apply_fisher_test(proteins.rename(columns={'change_in_one': complex_name}), background, complexes[complex_name], complex_name)
    enrichments.append(enrichment['test_summary'])
enrichments = pd.concat(enrichments).reset_index().rename(columns={'index': 'GO_id'})
enrichments['complex_name'] = enrichments['GO_id'].map(complex_names)
enrichments[enrichments['pval'] < 0.05].sort_values('pval')

enrichments.to_csv(f'{output_folder}complex_fishers_enrichment.csv')
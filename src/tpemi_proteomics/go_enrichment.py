import os
import pandas as pd
import numpy as np

from utilities.statistical_tests import apply_enrichment

from loguru import logger
logger.info('Import OK')

input_path = 'results/tpemi_proteomics/peptide_summary/peptide_summary.xlsx'
background_path = 'results/tpemi_proteomics/peptide_normalisation/normalised_summary.xlsx'
output_folder = 'results/tpemi_proteomics/go_enrichment/'

if not os.path.exists(output_folder):
    os.mkdir(output_folder)


# -----------------Read in standard components-----------------
# cleaned peptide summary
peptides = pd.read_excel(f'{input_path}', sheet_name='summary')
peptides.drop([col for col in peptides.columns.tolist() if 'Unnamed: ' in col], axis=1, inplace=True)

# raw data as background
background_genes = pd.read_excel(f'{background_path}', sheet_name=None)
background = background_genes['cys_peptides']['Proteins'].unique().tolist()

# ----------------Perform Panther enrichment test----------------

proteins = {}
for treatment, df in peptides.groupby('treatment'):
    pval_changed = df[['Proteins', 'log2_thresh_pval_ratio']].replace(0, np.nan).dropna().copy()
    pval_protected = pval_changed[pval_changed['log2_thresh_pval_ratio'] > 0].copy()
    pval_exposed = pval_changed[pval_changed['log2_thresh_pval_ratio'] < 0].copy()
    proteins[f'{treatment}_changed'] = pval_changed['Proteins'].unique().tolist()
    proteins[f'{treatment}_protected'] = pval_protected['Proteins'].unique().tolist()
    proteins[f'{treatment}_exposed'] = pval_exposed['Proteins'].unique().tolist()
combinations = pd.DataFrame(dict([ (k,pd.Series(v)) for k,v in proteins.items() ]))

# perform enrichment test
enrichment = apply_enrichment(combinations, searches=None, obo_path='resources/bioinformatics_databases/PANTHERGOslim.obo', organism='10090', refOrganism='10090', enrichmentTestType='FISHER', correction='FDR', min_proteins=2, reference=background)

# Save all to excel
enrichment.to_csv(f'{output_folder}go_enrichment.csv')


# -----------Read in cytoscape-clustered degree summary-----------
input_path = 'results/tpemi_proteomics/plot_degree_ppis/cytoscape_node_summary.csv'
output_folder = 'results/tpemi_proteomics/go_enrichment_degree_clusters/'

if not os.path.exists(output_folder):
    os.mkdir(output_folder)

node_summary = pd.read_csv(f'{input_path}')

clusters = {}
for cluster, df in node_summary.groupby('glay_cluster'):
    if len(df) < 2:
        continue
    clusters[f'complete_cluster_{cluster}'] = df['Proteins'].tolist()
    for direction, dataframe in df.groupby('direction'):
        clusters[f'{direction}_cluster_{cluster}'] = dataframe['Proteins'].tolist()
cluster_combinations = pd.DataFrame(
    dict([(k, pd.Series(v)) for k, v in clusters.items()]))
    

# perform enrichment test
enrichment = apply_enrichment(
    cluster_combinations, 
    searches=None, 
    obo_path='resources/bioinformatics_databases/PANTHERGOslim.obo', 
    organism='10090', 
    refOrganism='10090', 
    enrichmentTestType='FISHER', 
    correction='FDR', 
    min_proteins=2, 
    reference=node_summary['Proteins'].tolist()
    )

# Save all to excel
enrichment.to_csv(f'{output_folder}go_enrichment.csv')

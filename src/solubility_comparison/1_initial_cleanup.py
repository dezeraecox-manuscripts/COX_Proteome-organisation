import os, re
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

from utilities.databases import  uniprot_summary, create_uniprot_xref

from loguru import logger
from GEN_Utils import FileHandling

logger.info('Import OK')

input_folder = 'data/solubility_comparison/'
output_folder = 'results/solubility_comparison/initial_cleanup/'

if not os.path.exists(output_folder):
    os.makedirs(output_folder)

# --------Read in raw datasets from XS--------
file_list = [filename for filename in os.listdir(f'{input_folder}xs_datasets/') if '.xlsx' in filename]
drugs = {'Hsp70 inhibition': 'Ver155008', 'Novobiocin': 'Novobiocin', 'Proteasome inhibition': 'MG132'}

xs_data = []
for filename, sheet_name in zip(file_list, drugs.keys()):
    # columns are merged as header here, so skip row 1 
    # and note that repeated column names are in order of Total proteome, pSup and pellet proteome
    raw_data = pd.read_excel(f'{input_folder}xs_datasets/{filename}', sheet_name=sheet_name, skiprows=1)

    total_data = raw_data[['Accession', f'Average Abundance Ratio ({drugs[sheet_name]}:Control)', 'Significant']].copy()
    total_data.columns = ['Accession', 'abundance_ratio', 'significant']
    total_data['type'] = 'total'

    pSup_data = raw_data[['Accession', [col for col in raw_data.columns if f'pSup (Control - {drugs[sheet_name]})' in col][0], 'Significant.1']].copy()
    pSup_data.columns = ['Accession', 'abundance_ratio', 'significant']
    pSup_data['type'] = 'pSup'

    pellet_data = raw_data[['Accession', f'Average Abundance Ratio ({drugs[sheet_name]}:Control).1', 'Significant.2']].copy()
    pellet_data.columns = ['Accession', 'abundance_ratio', 'significant']
    pellet_data['type'] = 'pellet'

    df = pd.concat([total_data, pSup_data, pellet_data])
    df['drug'] = drugs[sheet_name]

    xs_data.append(df)
    
xs_data = pd.concat(xs_data)
xs_data['data_source'] = 'xs'

# ----Convert gene IDs to accession number----
uniprot_ids = uniprot_summary(
                    tax_id='10090',
                    resource_folder='resources/bioinformatics_databases/',
                    genes=[],
                    reviewed=False)
uniprot_map = dict(uniprot_ids[['UniProtKB-ID', 'Entry']].values)
xs_data['Proteins'] = xs_data['Accession'].map(uniprot_map)


# --------Read in raw datasets from DC--------
dc_data = pd.read_excel(f'{input_folder}dc_datasets/peptide_summary.xlsx', sheet_name=None)
dc_data.update({key: df.drop([col for col in df.columns.tolist() if 'Unnamed: ' in col], axis=1) for key, df in dc_data.items()})

# Collect max change per protein as overall comparison to xs
max_data = dc_data['max_proteins'].copy()

# melt to match columns to xs
max_data = pd.melt(
    max_data,
    id_vars=['Proteins'],
    value_vars=['MG132', 'Ver155008', 'Novobiocin'],
    value_name='max_cys_ratio',
    var_name='treatment'
)


max_data.columns = ['Proteins', 'drug', 'abundance_ratio']
max_data['significant'] = [np.nan if val == 0 else '+' for val in max_data['abundance_ratio']]

max_data['data_source'] = 'dc'

# --------Generate compiled dataframe---------

clean_data = pd.concat([xs_data.drop('Accession', axis=1), max_data])

# ----------------Save to csv-----------------
clean_data.to_csv(f'{output_folder}compiled_summary.csv')


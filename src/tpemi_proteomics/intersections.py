import os
import pandas as pd
import numpy as np
from upsetplot import from_contents, plot
import matplotlib.pyplot as plt

from GEN_Utils import FileHandling
from loguru import logger

logger.info('Import OK')

input_folder = 'results/tpemi_proteomics/peptide_summary/'
output_folder = 'results/tpemi_proteomics/intersections/'


if not os.path.exists(output_folder):
    os.makedirs(output_folder)

# Read in raw data
cys_peptides = pd.read_excel(f'{input_folder}peptide_summary.xlsx', sheet_name=None)
peptides = cys_peptides['summary'].copy()
peptides.drop([col for col in peptides.columns.tolist() if 'Unnamed: ' in col], axis=1, inplace=True)

treatments = ['Celastrol', 'Novobiocin', 'MG132', 'Ver155008', 'Staurosporine']
proteins = pd.pivot_table(peptides.copy(), index=['Proteins'], columns='treatment', values='log2_thresh_pval_ratio')
proteins = proteins[treatments].dropna() # select only proteins quantified in all treatments
proteins_all_quant = pd.melt(proteins.reset_index(), id_vars='Proteins', value_vars=treatments, var_name='treatment', value_name='log2_thresh_pval_ratio')
proteins_all_quant_one_change = pd.melt(proteins.replace(0, np.nan).dropna(how='all').replace(np.nan, 0).reset_index(), id_vars='Proteins', value_vars=treatments, var_name='treatment', value_name='log2_thresh_pval_ratio')

peptide_quant = pd.pivot_table(peptides.copy(), index=[
                          'Proteins', 'Sequence'], columns='treatment', values='log2_thresh_pval_ratio')
# select only proteins quantified in all treatments
peptide_quant = peptide_quant[treatments].dropna()
peptides_all_quant = pd.melt(peptide_quant.reset_index(), id_vars=['Proteins', 'Sequence'], value_vars=treatments,var_name='treatment', value_name='log2_thresh_pval_ratio')
peptides_all_quant_one_change = pd.melt(peptide_quant.replace(0, np.nan).dropna(how='all').replace(np.nan, 0).reset_index(), id_vars=['Proteins', 'Sequence'], value_vars=treatments, var_name='treatment', value_name='log2_thresh_pval_ratio')
# -----------------------------------UpSet plotting-----------------------------------
# Generate UpSet figure for significant proteins i.e. any protein that had at least one significant change 
# turn into dictionary to feed to pyupset
changed_prots = {df_name: df.replace(0, np.nan).dropna(subset=['log2_thresh_pval_ratio'])[
    'Proteins'].unique().tolist() for df_name, df in proteins_all_quant_one_change.groupby('treatment')}
changed_prots = from_contents(changed_prots)

fig = plot(changed_prots, sort_by='degree', show_counts=True)
plt.savefig(f'{output_folder}protein_upset.svg')

# Generate UpSet figure for significant peptides
changed_peps = {df_name: df.replace(0, np.nan).dropna(subset=['log2_thresh_pval_ratio'])[
    'Sequence'].unique().tolist() for df_name, df in peptides_all_quant_one_change.groupby('treatment')}
changed_peps = from_contents(changed_peps)

fig = plot(changed_peps, sort_by='degree', show_counts=True)


#-------------------Generating degree count results-----
degree_peps = changed_peps.reset_index().set_index('id').sum(axis=1).reset_index()
degree_peps.columns = ['Sequence', 'Degree']
degree_peps['Proteins'] = degree_peps['Sequence'].map(dict(peptides[['Sequence', 'Proteins']].values))
peps_summary = pd.merge(changed_peps.reset_index().rename(columns={'id': 'Sequence'}), degree_peps, on=['Sequence'])

degree_prots = changed_prots.reset_index().set_index('id').sum(axis=1).reset_index()
degree_prots.columns = ['Proteins', 'Degree']
prots_summary = pd.merge(changed_prots.reset_index().rename(
    columns={'id': 'Proteins'}), degree_prots, on=['Proteins'])

# Save summary info
peps_summary.to_csv(f'{output_folder}peptide_intersection_degree.csv')
prots_summary.to_csv(f'{output_folder}protein_intersection_degree.csv')

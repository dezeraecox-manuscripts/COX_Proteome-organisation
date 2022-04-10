import os
import pandas as pd
import numpy as np

from loguru import logger
from GEN_Utils import FileHandling

logger.info('Import OK')

input_normalised = 'results/tpemi_proteomics/peptide_normalisation/normalised_summary.xlsx'
input_thresholded = 'results/tpemi_proteomics/baseline_thresholding/thresholded_summary.xlsx'
input_ncthresholded = 'results/tpemi_proteomics/baseline_thresholding/noncys_thresholded_summary.xlsx'
output_folder = 'results/tpemi_proteomics/peptide_summary/'

quant_cols = ['Celastrol', 'MG132', 'Novobiocin', 'Staurosporine', 'Ver155008',]

if not os.path.exists(output_folder):
    os.mkdir(output_folder)

# Read in raw data
normalised = pd.read_excel(input_normalised, sheet_name=None)
normalised.update({key: df.drop([col for col in df.columns.tolist() if 'Unnamed: ' in str(col)], axis=1) for key, df in normalised.items()})

thresholded = pd.read_excel(input_thresholded, sheet_name=None)
thresholded.update({key: df.drop([col for col in df.columns.tolist() if 'Unnamed: ' in str(col)], axis=1) for key, df in thresholded.items()})

nc_thresholded = pd.read_excel(input_ncthresholded, sheet_name=None)
nc_thresholded.update({key: df.drop([col for col in df.columns.tolist(
) if 'Unnamed: ' in str(col)], axis=1) for key, df in nc_thresholded.items()})

"""------------------------Generate summary df------------------------"""

# raw_abundance
raw_vals = normalised['raw'].copy()
raw_vals = raw_vals[raw_vals['Sequence'].str.contains('C')]
raw_vals = pd.melt(raw_vals, id_vars=['Sequence', 'Proteins'], value_vars=[col for col in [f'{col}_{x}' for col in quant_cols for x in range(1, 5)] if col in raw_vals.columns.tolist()], var_name='treatment', value_name='mean_raw_ratio')
raw_vals[['treatment', 'replicate']] = raw_vals['treatment'].str.split('_', expand=True)
raw_vals = raw_vals.groupby(['Sequence', 'Proteins', 'treatment']).mean().reset_index()
raw_vals['log2_mean_raw_ratio'] = np.log2(raw_vals['mean_raw_ratio'])
summary = raw_vals.copy()

# noncys vals
noncys_vals = normalised['noncys_peptides'].copy()
noncys_vals = noncys_vals.groupby(['Proteins']).mean()[quant_cols].reset_index()
noncys_vals = pd.melt(noncys_vals, id_vars=['Proteins'], value_vars=quant_cols, var_name='treatment', value_name='mean_noncys_ratio')
noncys_vals['log2_mean_noncys_ratio'] = np.log2(noncys_vals['mean_noncys_ratio'])

summary = pd.merge(summary, noncys_vals, on=['Proteins', 'treatment'], how='left')

# corrected ratio
cys_vals = normalised['cys_noncys_peptides'].copy()
cys_vals = cys_vals.groupby(['Sequence', 'Proteins']).mean()[quant_cols].reset_index()
cys_vals = pd.melt(cys_vals, id_vars=['Sequence', 'Proteins'], value_vars=quant_cols, var_name='treatment', value_name='mean_corrected_cys_ratio')
cys_vals['log2_mean_corrected_cys_ratio'] = np.log2(cys_vals['mean_corrected_cys_ratio'])

summary = pd.merge(summary, cys_vals, on=['Sequence', 'Proteins', 'treatment'], how='left')


# pval_smooth - cys
pval_vals = thresholded['pval_smooth'].copy().dropna() #remove peptides which were not quantified
pval_vals['TPE_thresholded'] = [val if thresh == 1 else 0 for  thresh, val in pval_vals[['TPE_thresholded', 'log2_mean_ratio']].values]
pval_vals = pval_vals[['Sequence', 'Proteins', 'treatment', 'log2_mean_ratio', 'TPE_thresholded']]
pval_vals.columns = ['Sequence', 'Proteins', 'treatment', 'log2_pval_ratio', 'log2_thresh_pval_ratio']

summary = pd.merge(summary, pval_vals, on=['Sequence', 'Proteins', 'treatment'], how='left')

# pval_smooth - noncys
pval_vals = nc_thresholded['pval_smooth'].copy()
pval_vals['Proteins'] = pval_vals['Sequence'] # fix issue with taking mean of nc per protein
pval_vals.dropna(inplace=True)  # remove peptides which were not quantified
pval_vals['TPE_thresholded'] = [val if thresh == 1 else 0 for  thresh, val in pval_vals[['TPE_thresholded', 'log2_mean_ratio']].values]
pval_vals = pval_vals[['Proteins', 'treatment', 'log2_mean_ratio', 'TPE_thresholded']]
pval_vals.columns = ['Proteins', 'treatment', 'log2_noncys_pval_ratio', 'log2_thresh_noncys_pval_ratio']

summary = pd.merge(summary, pval_vals, on=['Proteins', 'treatment'], how='left')

# significant_changes
sign_vals = thresholded['significant_changes'].copy().dropna() #remove peptides which were not quantified
sign_vals['TPE_thresholded'] = [val if thresh == 1 else 0 for  thresh, val in sign_vals[['TPE_thresholded', 'log2_mean_ratio']].values]
sign_vals = sign_vals[['Sequence', 'Proteins', 'treatment', 'log2_mean_ratio', 'TPE_thresholded']]
sign_vals.columns = ['Sequence', 'Proteins', 'treatment', 'log2_significant_ratio', 'log2_thresh_significant_ratio']

summary = pd.merge(summary, sign_vals, on=['Sequence', 'Proteins', 'treatment'], how='left')

"""------------------------Assign common changes------------------------"""

# collect peptides identified in all treatments of interest for measure of interest - currently using pval smoothed dataset
# currently using ------"""pval smoothed"""------ dataset, alternative is log2_thresh_significant_ratio
treatments = ['MG132', 'Celastrol', 'Ver155008', 'Novobiocin', 'Staurosporine']
per_peptide = pd.pivot_table(summary.copy(), index=['Sequence', 'Proteins'], columns='treatment', values='log2_thresh_pval_ratio')

# Select only peptides which change in at least one treatment --> no change will now be represented by nan
changed_peptides = per_peptide[treatments].dropna() # select only peptides quantified in all treatments
changed_peptides = changed_peptides.replace(0, np.nan).dropna(how='all', subset=treatments)
len(changed_peptides.reset_index().groupby('Proteins').count()['Sequence'].value_counts())

# Calculate number of treatments protein was changed in using binary method
binary_peptides = changed_peptides.copy()
binary_peptides[treatments] = np.where(binary_peptides[treatments].isnull(), 0, 1)
binary_peptides['num_changed'] = binary_peptides.sum(axis=1)


# Produce summary measure per protein --> explore mean, max and simple binary measures
mean_proteins = per_peptide.reset_index().copy().groupby(
    'Proteins').mean().dropna()  # select only proteins quantified in all treatments.dropna()  # select only proteins quantified in all treatments

binary_proteins = mean_proteins.copy()
binary_proteins[treatments] = np.where(binary_proteins[treatments].replace(0, np.nan).isnull(), 0, 1)

max_proteins = []
for protein, df in per_peptide.reset_index().groupby('Proteins'):
    if len(df) == 1:
        max_proteins.append(df.reset_index()[['Proteins']+treatments])
    else:
        treatment_vals = []
        for treatment in treatments:
            protein_vals = df[treatment].dropna().sort_values().tolist()
            if len(protein_vals) == 0:
                treatment_vals.append(np.nan)
            else:
                max_val = protein_vals[-1] if abs(protein_vals[-1]) > abs(protein_vals[0]) else protein_vals[0]
                treatment_vals.append(max_val)
        max_vals = pd.DataFrame([[protein] + treatment_vals])
        max_vals.columns = ['Proteins'] + treatments
                
        max_proteins.append(max_vals)

# select only proteins quantified in all treatments
max_proteins = pd.concat(max_proteins).dropna()

# Melt and summarise per protein changes
protein_summary = pd.melt(mean_proteins.copy().reset_index(), id_vars='Proteins', value_vars=treatments, var_name='treatment', value_name='mean_log2_pval_thresh')
protein_summary = pd.merge(protein_summary, pd.melt(max_proteins.copy().reset_index(), id_vars='Proteins', value_vars=treatments, var_name='treatment', value_name='max_log2_pval_thresh'), on=['Proteins', 'treatment'], how='outer')
protein_summary = pd.merge(protein_summary, pd.melt(binary_proteins.copy().reset_index(), id_vars='Proteins', value_vars=treatments, var_name='treatment', value_name='binary_pval_thresh'), on=['Proteins', 'treatment'], how='outer')

# Determine number of treatments with changes in common
binary_proteins['num_changed'] = binary_proteins.sum(axis=1)


FileHandling.df_to_excel(
    output_path=f'{output_folder}peptide_summary.xlsx',
    sheetnames=['summary', 'per_peptide', 'changed_peptides', 'binary_peptides', 'mean_proteins', 'max_proteins', 'binary_proteins', 'changed_protein_summary', ],
    data_frames=[summary, per_peptide, changed_peptides, binary_peptides, mean_proteins, max_proteins, binary_proteins, protein_summary, ]
)

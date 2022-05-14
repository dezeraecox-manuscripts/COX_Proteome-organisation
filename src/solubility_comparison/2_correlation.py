import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats

from loguru import logger

logger.info('Import OK')

input_path = 'results/solubility_comparison/initial_cleanup/compiled_summary.csv'
output_folder = 'results/solubility_comparison/correlation/'

if not os.path.exists(output_folder):
    os.mkdir(output_folder)

def r2(x, y, round_pos=False):
    r = stats.pearsonr(x, y)[0] ** 2
    return round(r, round_pos) if round_pos else r


import matplotlib
font = {'family' : 'normal',
'weight' : 'normal',
'size'   : 12 }
matplotlib.rc('font', **font)
plt.rcParams['svg.fonttype'] = 'none'

# Read in compiled summary
compiled = pd.read_csv(input_path)
compiled.drop([col for col in compiled.columns.tolist() if 'Unnamed: ' in col], axis=1, inplace=True)
compiled['type'] = compiled['type'].fillna('pval')


# Determine scale of protein changes as proportion of total number ID proteins
scale_of_changes = compiled.groupby(['data_source', 'type', 'drug']).count().reset_index()
scale_of_changes['proportion_changes'] = scale_of_changes['significant'] / scale_of_changes['abundance_ratio'] * 100


# pivot df based on protein ID
abundance_correlation = pd.pivot_table(
    data=compiled[compiled['type'].isin(['pval', 'pellet'])].dropna(subset=['significant']).copy(),
    index=['Proteins', 'drug'],
    columns=['data_source'],
    values=['abundance_ratio']
)
abundance_correlation.columns = ['_'.join(col).strip() for col in abundance_correlation.columns.values]
abundance_correlation = abundance_correlation.reset_index().dropna()

for drug, df in abundance_correlation.groupby('drug'):
    logger.info(df.corr())

    fig, ax = plt.subplots(figsize=(5, 5))
    sns.regplot(
        x=df['abundance_ratio_dc'],
        y=np.log2(df['abundance_ratio_xs']),
        color='black'
        )
    plt.annotate(r2(df['abundance_ratio_dc'],
                 np.log2(df['abundance_ratio_xs']), 4),
                 xy=(0.8, 0.9),
                 xycoords='axes fraction')
    plt.title(drug)
    plt.ylabel("Sui dataset\n"r"$log_2$(Abundance Ratio)")
    plt.xlabel(r"TPE-MI dataset\nCorrected cys Ratio)")
    plt.savefig(f'{output_folder}regression_{drug}.svg')
    plt.show()

# For each data_type/treatment, determine correlation overall
treatment_correlation = pd.pivot_table(
    data=compiled[(compiled['type'].isin(['pval', 'total'])) & compiled['drug'].isin(['Ver155008', 'MG132', 'Novobiocin'])].dropna(subset=['significant']).copy(),
    index=['Proteins'],
    columns=['data_source', 'drug'],
    values=['abundance_ratio']
)
treatment_correlation.columns = ['_'.join(col).strip() for col in treatment_correlation.columns.values]
treatment_correlation.dropna(how='all').corr()

sns.pairplot(
    np.log2(treatment_correlation)
)


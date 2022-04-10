import os
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib_venn

from loguru import logger

logger.info('Import OK')

input_path = 'results/solubility_comparison/initial_cleanup/compiled_summary.csv'
output_folder = 'results/solubility_comparison/intersections/'

if not os.path.exists(output_folder):
    os.makedirs(output_folder)

import matplotlib
font = {'family' : 'normal',
'weight' : 'normal',
'size'   : 18 }
matplotlib.rc('font', **font)
plt.rcParams['svg.fonttype'] = 'none'

# Read in raw data
compiled = pd.read_csv(f'{input_path}')
compiled.drop([col for col in compiled.columns.tolist() if 'Unnamed: ' in col], axis=1, inplace=True)
compiled['type'] = compiled['type'].fillna('pval')

# Plot Venn for all proteins quantified
for drug in ['Ver155008', 'MG132', 'Novobiocin']:
    dc_set = set(compiled[(compiled['drug'] == drug) & (
        compiled['data_source'] == 'dc') & (compiled['type'] == 'pval')]['Proteins'])
    xs_set = set(compiled[(compiled['drug'] == drug) & (
        compiled['data_source'] == 'xs') & (compiled['type'] == 'pellet')]['Proteins'])
    v = matplotlib_venn.venn2(subsets=[dc_set, xs_set], set_labels=[
                          f'dc_{drug}', f'xs_{drug}'])
    v.get_patch_by_id('10').set_color('#d1d1d1')
    v.get_patch_by_id('11').set_color('black')
    v.get_patch_by_id('01').set_color('#d1d1d1')
    plt.savefig(f'{output_folder}dcVxs_{drug}_all.svg')
    plt.show()


# Plot Venn for significant changes

sig_proteins = compiled.dropna(subset=['significant']).copy()

for drug in ['Ver155008', 'MG132', 'Novobiocin']:
    dc_set = set(sig_proteins[(sig_proteins['drug'] == drug) & (sig_proteins['data_source'] == 'dc') & (sig_proteins['type'] == 'pval')]['Proteins'])
    xs_set = set(sig_proteins[(sig_proteins['drug'] == drug) & (sig_proteins['data_source'] == 'xs') & (sig_proteins['type'] == 'pellet')]['Proteins'])
    v = matplotlib_venn.venn2(subsets=[dc_set, xs_set], set_labels=[f'dc_{drug}', f'xs_{drug}'], )
    v.get_patch_by_id('10').set_color('#d1d1d1')
    v.get_patch_by_id('11').set_color('black')
    v.get_patch_by_id('01').set_color('#d1d1d1')
    plt.savefig(f'{output_folder}dcVxs_{drug}_changed.svg')
    plt.show()

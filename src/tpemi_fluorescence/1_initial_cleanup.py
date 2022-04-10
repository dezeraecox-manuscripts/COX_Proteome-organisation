import os
import pandas as pd

from loguru import logger

logger.info('Import OK')

input_folder = 'data/tpemi_fluorescence/'
output_folder = 'results/tpemi_fluorescence/initial_cleanup/'

if not os.path.exists(output_folder):
    os.makedirs(output_folder)

# Read in EXPA data, layout
exp_A = pd.read_excel(f'{input_folder}EXPA.xls')
exp_A.drop([col for col in exp_A.columns.tolist() if 'Unnamed: ' in col], axis=1, inplace=True)
exp_A.columns = ['tube_name', 'median_fluorescence']

layout_A = pd.read_excel(f'{input_folder}EXPA_layout.xlsx')
layout_A = pd.melt(layout_A, value_name='well', var_name='treatment')
layout_A['tube_name'] = layout_A['well'].str.replace('_', '')
layout_A.dropna(subset=['tube_name'], inplace=True)

# map treatment to tube name, select treatments of interest
exp_A['treatment'] = exp_A['tube_name'].map(
    dict(layout_A[['tube_name', 'treatment']].values))
exp_A.dropna(subset=['treatment'], inplace=True)
exp_A = exp_A[exp_A['treatment'].isin(['Control', 'Ver155008', 'MG132', 'Celastrol'])]

# normalise to control
control_mean = exp_A[exp_A['treatment'] == 'Control'].mean()['median_fluorescence']
exp_A['norm_fluorescence'] = exp_A['median_fluorescence'] / control_mean

# Repeat for EXPB
exp_B = pd.read_excel(f'{input_folder}EXPB.xls')
exp_B.drop([col for col in exp_B.columns.tolist()
            if 'Unnamed: ' in col], axis=1, inplace=True)
exp_B.drop('PLATE NAME', axis=1, inplace=True)
exp_B.columns = ['tube_name', 'median_fluorescence']

layout_B = pd.read_excel(f'{input_folder}EXPB_layout.xlsx', header=None)
layout_B = pd.melt(layout_B, id_vars=[0], value_name='well', var_name='column').drop('column', axis=1)
layout_B.columns = ['treatment', 'tube_name']
layout_B['tube_name'] = layout_B['tube_name'].str.lstrip('1_').str.replace('_', '0')
layout_B.dropna(subset=['tube_name'], inplace=True)

# map treatment to tube name
exp_B['treatment'] = exp_B['tube_name'].map(
    dict(layout_B[['tube_name', 'treatment']].values))
exp_B.dropna(subset=['treatment'], inplace=True)

# normalise to control
mQ_mean = exp_B[exp_B['treatment'] == 'mQ'].mean()['median_fluorescence']
DMSO_mean = exp_B[exp_B['treatment'] == 'DMSO'].mean()['median_fluorescence']
controls = {
    'mQ': mQ_mean, 
    'Novo': mQ_mean, 
    'DMSO': DMSO_mean, 
    'Stauro': DMSO_mean
}

exp_B['norm_fluorescence'] = [val / controls[treatment] for treatment, val in exp_B[['treatment', 'median_fluorescence']].values]

exp_B['treatment'] = exp_B['treatment'].map(
    {'mQ': 'Control', 'DMSO': 'Control', 'Novo': 'Novobiocin', 'Stauro': 'Staurosporine'})

# Compile, save
summary = pd.concat([exp_A, exp_B])
summary.to_csv(f'{output_folder}compiled.csv')

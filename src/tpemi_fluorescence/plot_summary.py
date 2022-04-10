import os
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

from loguru import logger

logger.info('Import OK')

input_path = 'results/tpemi_fluorescence/initial_cleanup/compiled.csv'
output_folder = 'results/tpemi_fluorescence/plot_summary/'

if not os.path.exists(output_folder):
    os.makedirs(output_folder)

import matplotlib
font = {'family' : 'normal',
'weight' : 'normal',
'size'   : 12 }
matplotlib.rc('font', **font)
plt.rcParams['svg.fonttype'] = 'none'

# Read in summarised data
compiled = pd.read_csv(input_path)
compiled.drop([col for col in compiled.columns.tolist()
              if 'Unnamed: ' in col], axis=1, inplace=True)

# Generate summary plot
xorder = ['Control', 'MG132', 'Ver155008',
          'Staurosporine', 'Celastrol', 'Novobiocin']
fig, ax = plt.subplots(figsize=(5, 5))
sns.boxplot(
    data=compiled,
    x='treatment',
    y='norm_fluorescence',
    color='white',
    fliersize=0,
    order=xorder
    )
plt.setp(ax.artists, edgecolor='darkgrey', facecolor='w')
plt.setp(ax.lines, color='darkgrey')
sns.stripplot(
    data=compiled,
    x='treatment',
    y='norm_fluorescence',
    dodge=True,
    color='black',
    order=xorder)
ax.axhline(1, color='lightgrey', linestyle='--')
plt.ylabel(r'Relative TPE-MI Fluorescence (A.U.)')
plt.xlabel('')
ax.set_xticklabels(ax.get_xticklabels(), rotation=45)
plt.ylim(0.95, 1.65)
plt.savefig(f'{output_folder}boxplot_fluorescence.svg')
plt.show()

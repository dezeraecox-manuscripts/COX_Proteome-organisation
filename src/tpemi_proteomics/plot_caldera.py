
import os
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

from loguru import logger
logger.info("Import ok")

import matplotlib
font = {'family' : 'normal',
'weight' : 'normal',
'size'   : 22 }
matplotlib.rc('font', **font)
plt.rcParams['svg.fonttype'] = 'none'

input_treated = 'results/tpemi_proteomics/peptide_summary/peptide_summary.xlsx'
input_controls = 'results/tpemi_proteomics/baseline/significance/cys_scaled_summary.xlsx'
thresholds_cys = 'results/tpemi_proteomics/baseline_thresholding/thresholded_summary.xlsx'
thresholds_noncys = 'results/tpemi_proteomics/baseline_thresholding/noncys_thresholded_summary.xlsx'

output_folder = 'results/tpemi_proteomics/plot_caldera/'

if not os.path.exists(output_folder):
    os.makedirs(output_folder)    


# Read in summary data
peptides = pd.read_excel(f'{input_treated}', sheet_name='summary')
peptides.drop([col for col in peptides.columns.tolist() if 'Unnamed: ' in str(col)], axis=1, inplace=True)

cys_thresholds = pd.read_excel(thresholds_cys, sheet_name='CIs')
cys_thresholds = cys_thresholds[(cys_thresholds['data_type'] == 'pval_smooth') & (cys_thresholds['sample'] == 'TPE')][['lower_CI', 'upper_CI']].values.flatten()

noncys_thresholds = pd.read_excel(thresholds_noncys, sheet_name='CIs')
noncys_thresholds = noncys_thresholds[(noncys_thresholds['data_type'] == 'pval_smooth') & (noncys_thresholds['sample'] == 'TPE')][['lower_CI', 'upper_CI']].values.flatten()

# Prepare control data
cys_controls = 'results/tpemi_proteomics/baseline/significance/cys_scaled_summary.xlsx'
cys_controls = pd.read_excel(input_controls, sheet_name='pval_smooth')
cys_controls.drop([col for col in cys_controls.columns.tolist()
                  if 'Unnamed: ' in col], axis=1, inplace=True)
cys_controls = pd.melt(
    cys_controls,
    id_vars=['Sequence', 'Proteins'],
    value_vars=['TPE', 'NOTPE'],
    value_name='log2_pval_ratio',
    var_name='treatment'
)
noncys_controls = 'results/tpemi_proteomics/baseline/significance/noncys_scaled_summary.xlsx'
noncys_controls = pd.read_excel(input_controls, sheet_name='pval_smooth')
noncys_controls.drop([col for col in noncys_controls.columns.tolist()
                      if 'Unnamed: ' in col], axis=1, inplace=True)
noncys_controls = pd.melt(
    noncys_controls,
    id_vars=['Sequence', 'Proteins'],
    value_vars=['TPE', 'NOTPE'],
    value_name='log2_noncys_pval_ratio',
    var_name='treatment'
)

control_data = pd.merge(
    cys_controls,
    noncys_controls.groupby(['Proteins', 'treatment']).mean().reset_index(),
    on=['Proteins', 'treatment'],
    how='left'
)

# plot caldera for each stress
def caldera_plot(df, x='log2_mean_noncys_ratio', y='log2_pval_ratio', xthresholds=(-1, 1), ythresholds=(-1, 1), palette=False, xlabel=False, ylabel=False, title=False, color_scheme='reactivity', save_fig=False):
    """Generates caldera plot for changes in cys peptides, according to thresholds provided by TPE control samples

    Parameters
    ----------
    df : DataFrame
        melted dataframe containing only treatment of interest, x and y columns
    x : str, optional
        column label of x values, by default 'log2_mean_noncys_ratio'
    y : str, optional
        column label of y values, by default 'log2_pval_ratio'
    thresholds : tuple, optional
        lower and upper y thresholds according to control values, by default (-1, 1)
    xlabel : bool, optional
        If provided, alter x axis label, by default False
    ylabel : bool, optional
        If provided, alter y axis label, by default False
    title : bool, optional
        If provided, add title to plot, by default False
    save_fig : bool, optional
        If True, save fig to current output_folder. Title must also be supplied, by default False
    """
    # assign thresholded colours
    if not palette:
        palette = {0: '#b8b8b8', 1: '#404040', 2: '#000000'}
    df['color'] = 0
    if color_scheme == 'quadrants':
        df['color'] = [1 if (
            (xval < xthresholds[0]) | (yval < ythresholds[0]) |
            (xval > xthresholds[1]) | (yval > ythresholds[1])) else color for xval, yval, color in df[[x, y, 'color']].values]
        df['color'] = [2 if (
            (xval < xthresholds[0]) & (yval < ythresholds[0]) |
            (xval > xthresholds[1]) & (yval > ythresholds[1]) |
            (xval < xthresholds[0]) & (yval > ythresholds[1]) |
            (xval > xthresholds[1]) & (yval < ythresholds[0])) else color for xval, yval, color in df[[x, y, 'color']].values]
    elif color_scheme == 'reactivity':
        df['color'] = [1 if (yval < ythresholds[0]) else color for yval, color in df[[y, 'color']].values]
        df['color'] = [2 if (yval > ythresholds[1]) else color for yval, color in df[[ y, 'color']].values]
    print(len(df[df['color'] != 0]['Proteins'].unique()))

    # plot scatterplot
    fig, ax = plt.subplots(figsize=(6, 6))

    sns.scatterplot(data=df, x=x, y=y, hue='color', palette=palette,
                    linewidth=0.25, alpha=0.6)

    for thresh in xthresholds:
        ax.axvline(thresh, color='red', linestyle='--', linewidth=0.5)
    for thresh in ythresholds:
        ax.axhline(thresh, color='red', linestyle='--', linewidth=0.5)

    plt.legend('')
    plt.ylim(-8, 8)
    plt.xlim(-4, 4)
    if xlabel:
        plt.xlabel(xlabel)    
    if ylabel:
        plt.ylabel(ylabel)
    if title:
        plt.title(title)
    if save_fig:
        plt.savefig(f'{output_folder}{title}.png')
        plt.savefig(f'{output_folder}{title}.svg')


# Visualise Caldera plot --> using pval smoothed values without thresholding (to show distribution of values within threshold)
   
for treatment_name in ['TPE', 'NOTPE']:
    df = control_data[control_data['treatment'] == treatment_name].copy()
    caldera_plot(
        df,
        x='log2_noncys_pval_ratio',
        y='log2_pval_ratio',
        xthresholds=noncys_thresholds,
        ythresholds=cys_thresholds,
        xlabel='Log$_2$ non-cys *ratio',
        ylabel='Log$_2$ corrected cys *ratio',
        title=treatment_name,
        palette={0: '#404040', 1: '#404040', 2: '#404040'},
        save_fig=True)

for treatment_name in ['MG132', 'Ver155008', 'Novobiocin', 'Celastrol', 'Staurosporine']:
    df = peptides[peptides['treatment'] == treatment_name].copy()
    caldera_plot(
        df,
        x='log2_noncys_pval_ratio',
        y='log2_pval_ratio',
        xthresholds=noncys_thresholds,
        ythresholds=cys_thresholds,
        xlabel='Log$_2$ non-cys *ratio',
        ylabel='Log$_2$ corrected cys *ratio',
        title=treatment_name,
        palette={0: '#5d5d5d', 1: '#a10000', 2: '#063478'},
        save_fig=True)

# ---------Plot distribution on TPE-MI control dataset for Supp info---------
after = control_data[control_data['treatment'] == 'TPE'].copy().dropna()
after.columns

input_before = 'results/tpemi_proteomics/baseline/peptide_normalisation/normalised_summary.xlsx'
# read in raw data
before_cys = pd.read_excel(f'{input_before}', sheet_name='cys_noncys_peptides')
before_cys.drop(
    [col for col in before_cys.columns.tolist() if 'Unnamed: ' in col], axis=1, inplace=True)

before_noncys = pd.read_excel(f'{input_before}', sheet_name='noncys_peptides')
before_noncys.drop(
    [col for col in before_noncys.columns.tolist() if 'Unnamed: ' in col], axis=1, inplace=True)


# Cys
fig, ax = plt.subplots()
sns.distplot(
    np.log2(before_cys['TPE']),
    bins=np.arange(-2.5, 2.5, 0.25),
    label='Before scaling',
    color='darkgrey',
    # hist=False,
    kde_kws={'bw':0.25}
)
sns.distplot(
    after['log2_pval_ratio'],
    bins=np.arange(-2.5, 2.5, 0.25),
    label='After scaling',
    color='black',
    # hist=False,
    kde_kws={'bw':0.25}
)
plt.xlabel('Corrected cys ratio')
plt.ylabel('Density')
for threshold in cys_thresholds:
    ax.axvline(threshold, color='firebrick', linestyle='--')
plt.legend(bbox_to_anchor=(1.0, 1.0))
plt.savefig(f'{output_folder}distribution_TPE_cys.svg')

# NonCys
fig, ax = plt.subplots()
sns.distplot(
    np.log2(before_noncys['TPE']),
    bins=np.arange(-2.5, 2.5, 0.25),
    label='Before scaling',
    color='darkgrey',
    # hist=False,
    kde_kws={'bw': 0.25}
)
sns.distplot(
    after['log2_noncys_pval_ratio'],
    bins=np.arange(-2.5, 2.5, 0.25),
    label='After scaling',
    color='black',
    # hist=False,
    kde_kws={'bw': 0.25}
)
plt.xlabel('Corrected cys ratio')
plt.ylabel('Density')
for threshold in noncys_thresholds:
    ax.axvline(threshold, color='firebrick', linestyle='--')
plt.legend(bbox_to_anchor=(1.0, 1.0))
plt.savefig(f'{output_folder}distribution_TPE_noncys.svg')

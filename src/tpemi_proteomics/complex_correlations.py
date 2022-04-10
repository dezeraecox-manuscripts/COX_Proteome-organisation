import os
import pandas as pd
from scipy.stats import ttest_ind

from loguru import logger
logger.info('Import OK')

correlation_path = 'results/tpemi_proteomics/correlations/max_per_protein_correlations.csv'
feature_path = 'results/tpemi_proteomics/protein_complexes/ontology_summary.xlsx'
output_folder = 'results/tpemi_proteomics/complex_correlations/'

if not os.path.exists(output_folder):
    os.mkdir(output_folder)


# Read in correlation data
raw_data = pd.read_csv(f'{correlation_path}')
raw_data.drop([col for col in raw_data.columns.tolist() if 'Unnamed: ' in col], axis=1, inplace=True)
protein_corr = pd.melt(
    raw_data,
    id_vars=['Proteins'],
    value_vars=[col for col in raw_data.columns if col != 'Proteins'],
    var_name='Protein_B',
    value_name='corr_val',
).rename(columns={'Proteins': 'Protein_A'})
protein_corr['key'] = [tuple(sorted([prot_A, prot_B])) for prot_A, prot_B in protein_corr[['Protein_A', 'Protein_B']].values]
protein_corr = protein_corr[protein_corr['Protein_A'] != protein_corr['Protein_B']]

# Read in interaction data
feature_data = pd.read_excel(f'{feature_path}', sheet_name='GO_complexes')
feature_data.drop([col for col in feature_data.columns.tolist() if 'Unnamed: ' in col], axis=1, inplace=True)
complex_names = pd.read_excel(f'{feature_path}', sheet_name='GO_summary')
complex_names = dict(complex_names[['go_id', 'name']].values)

# for a given complex, assign proteins inside versus outside the complex
complex_correlations = []
for complex_name in feature_data.columns:
    complex_proteins = feature_data[complex_name].dropna().tolist()
    complex_corrs = protein_corr[protein_corr['Protein_A'].isin(complex_proteins)].copy()
    complex_corrs['corr_type'] = [1 if protein_b in complex_proteins else 0 for protein_b in complex_corrs['Protein_B']]
    complex_corrs = pd.pivot(complex_corrs.groupby(['Protein_A', 'corr_type']).mean(
    ).reset_index(), columns='corr_type', index='Protein_A').reset_index()
    complex_corrs['complex_id'] = complex_name
    complex_correlations.append(complex_corrs)
complex_correlations = pd.concat(complex_correlations)
complex_correlations.columns = ['Proteins', 'complex_id', '0', '1']

# for each complex, test the difference in average correlation of protein with proteins within the complex vs without
mean_comparisons_tests = []
for complex_name, df in complex_correlations.groupby('complex_id'):
    if len(df) < 3:
        # Only process complexes with at least 3 members quantified
        continue 
    stat, pval = ttest_ind(df['0'].abs(), df['1'].abs())
    summary = df.copy()
    summary['tstat'] = stat
    summary['pval'] = pval
    mean_comparisons_tests.append(summary)

mean_comparisons_tests = pd.concat(
    mean_comparisons_tests).reset_index(drop=True)
mean_comparisons_tests['complex_name'] = mean_comparisons_tests['complex_id'].map(complex_names)

mean_comparisons_tests['sign'] = ['***' if pval < 0.001 else ('**' if pval < 0.01 else ('*' if pval < 0.05 else 'ns')) for pval in mean_comparisons_tests['pval']]


complex_correlations.to_csv(f'{output_folder}protein_mean_comparisons.csv')
mean_comparisons_tests.to_csv(f'{output_folder}protein_mean_comparisons_tests.csv')


# Generate stats summary
summary = mean_comparisons_tests.copy()
summary['0'] = summary['0'].abs()
summary['1'] = summary['1'].abs()
summary = summary.groupby(
    ['complex_id', 'complex_name', 'tstat', 'pval', 'sign']).agg(['mean', 'std', 'count']).reset_index()
summary.columns = [f'{i}_{j}' if j !=
                   '' else f'{i}' for i, j in summary.columns]
summary.to_csv(
    f'{output_folder}complex_correlation_summary.csv')

# compile graph-pad friendly table
corr_summary = mean_comparisons_tests.copy()
corr_summary = {f'{cluster}_{corr_type}': df['value'].tolist() for (cluster, corr_type), df in pd.melt(
    corr_summary, id_vars=['complex_id'], value_vars=['0', '1'], var_name='corr_type').groupby(['complex_id', 'corr_type'])}
corr_summary = pd.DataFrame(
    dict([(k, pd.Series(v)) for k, v in corr_summary.items()]))
corr_summary[corr_summary.columns.tolist(
)] = corr_summary[corr_summary.columns.tolist()].abs()
corr_summary.to_csv(
    f'{output_folder}complex_corr_vals.csv')

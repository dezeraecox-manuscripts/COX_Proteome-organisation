import os
import pandas as pd
from scipy.stats import ttest_ind

from loguru import logger
logger.info('Import OK')

correlation_path = 'results/tpemi_proteomics/correlations/max_per_protein_correlations.csv'
input_clusters = 'results/tpemi_proteomics/plot_degree_ppis/cytoscape_node_summary.csv'

output_folder = 'results/tpemi_proteomics/cluster_correlations/'

if not os.path.exists(output_folder):
    os.mkdir(output_folder)


# Read in correlation data
raw_data = pd.read_csv(f'{correlation_path}')
raw_data.drop([col for col in raw_data.columns.tolist()
              if 'Unnamed: ' in col], axis=1, inplace=True)
protein_corr = pd.melt(
    raw_data,
    id_vars=['Proteins'],
    value_vars=[col for col in raw_data.columns if col != 'Proteins'],
    var_name='Protein_B',
    value_name='corr_val',
).rename(columns={'Proteins': 'Protein_A'})
protein_corr['key'] = [tuple(sorted([prot_A, prot_B]))
                       for prot_A, prot_B in protein_corr[['Protein_A', 'Protein_B']].values]
protein_corr = protein_corr[protein_corr['Protein_A']
                            != protein_corr['Protein_B']].copy()
                            
# Read in cluster info
ppi_nodes = pd.read_csv(input_clusters)
ppi_nodes.drop([col for col in ppi_nodes.columns.tolist()
                if 'Unnamed: ' in col], axis=1, inplace=True)
# select clustered proteins, all else assigned to 'unclustered' 0
ppi_nodes['cluster'] = [cluster if cluster in [1.0, 2.0, 3.0,
                                               4.0, 6.0] else 0 for cluster in ppi_nodes['glay_cluster']]
ppi_clusters = dict(ppi_nodes[['Proteins', 'cluster']].values)
# organise according to text structure
cluster_map = {0: 0, 1: 4, 2: 2, 3: 1, 4: 5, 6: 3}
protein_corr['cluster_A'] = protein_corr['Protein_A'].map(
    ppi_clusters).map(cluster_map)
protein_corr['cluster_B'] = protein_corr['Protein_B'].map(
    ppi_clusters).map(cluster_map)
ppi_clusters = dict(protein_corr[['Protein_A', 'cluster_A']].values)


# Test average correlation for all members of a cluster with proteins inside cluster, versus outside cluster
protein_corr['corr_type'] = [1 if clust_a == clust_b else 0 for clust_a, clust_b in protein_corr[['cluster_A', 'cluster_B']].values]

mean_comparisons = []
for protein, df in protein_corr.groupby('Protein_A'):
    mean_corr = df.groupby('corr_type').mean()['corr_val'].reset_index().set_index('corr_type').T
    mean_corr['Proteins'] = protein
    mean_corr['cluster'] = ppi_clusters[protein]
    mean_comparisons.append(mean_corr)
mean_comparisons = pd.concat(mean_comparisons).reset_index(drop=True)

# for clusters, test the difference in average correlation of protein with proteins within the cluster vs without
mean_comparisons_tests = []
for cluster_name, df in mean_comparisons.groupby('cluster'):
    stat, pval = ttest_ind(df[0].abs(), df[1].abs())
    summary = df.copy()
    summary['tstat'] = stat
    summary['pval'] = pval
    mean_comparisons_tests.append(summary)

summary = mean_comparisons.copy()
stat, pval = ttest_ind(summary[0].abs().tolist(), summary[1].abs().tolist())
summary['cluster'] = 'all'
summary['tstat'] = stat
summary['pval'] = pval
mean_comparisons_tests.append(summary)

mean_comparisons_tests = pd.concat(
    mean_comparisons_tests).reset_index(drop=True)

mean_comparisons_tests[mean_comparisons_tests['pval']
                       < 0.05].groupby('cluster').mean()
mean_comparisons_tests['sign'] = ['***' if pval < 0.001 else ('**' if pval < 0.01 else ('*' if pval < 0.05 else 'ns')) for pval in mean_comparisons_tests['pval']]

mean_comparisons.to_csv(f'{output_folder}protein_mean_comparisons.csv')
mean_comparisons_tests.to_csv(
    f'{output_folder}protein_mean_comparisons_tests.csv')

# Generate stats summary
summary = mean_comparisons_tests.copy()
summary[0] = summary[0].abs()
summary[1] = summary[1].abs()
summary = summary.groupby(
    ['cluster', 'tstat', 'pval', 'sign']).agg(['mean', 'std', 'count']).reset_index()
summary.columns = [f'{i}_{j}' if j !=
                   '' else f'{i}' for i, j in summary.columns]
summary.to_csv(
    f'{output_folder}cluster_correlation_summary.csv')

# compile graph-pad friendly table
corr_summary = mean_comparisons_tests.copy()
corr_summary = {f'{cluster}_{corr_type}': df['value'].tolist() for (cluster, corr_type), df in pd.melt(
    corr_summary, id_vars=['cluster'], value_vars=[0, 1]).groupby(['cluster', 'corr_type'])}
corr_summary = pd.DataFrame(
    dict([(k, pd.Series(v)) for k, v in corr_summary.items()]))
corr_summary[corr_summary.columns.tolist()] = corr_summary[corr_summary.columns.tolist()].abs()
corr_summary.to_csv(
    f'{output_folder}cluster_corr_vals.csv')


import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

from loguru import logger
logger.info('Import OK')

input_proteins = 'results/tpemi_proteomics/peptide_summary/peptide_summary.xlsx'
input_clusters = 'results/tpemi_proteomics/plot_degree_ppis/cytoscape_node_summary.csv'
input_complexes = 'results/tpemi_proteomics/complex_correlations/protein_mean_comparisons_tests.csv'
input_interactions = 'results/tpemi_proteomics/ppi_correlations/ppi_correlation_vals.csv'
output_folder = 'results/tpemi_proteomics/plot_correlations/'


if not os.path.exists(output_folder):
    os.mkdir(output_folder)

import matplotlib
font = {'family' : 'normal',
'weight' : 'normal',
'size'   : 12 }
matplotlib.rc('font', **font)
plt.rcParams['svg.fonttype'] = 'none'

# Read in summary data - proccessed same as PPIs
treatments = ['Celastrol', 'Novobiocin', 'MG132', 'Ver155008', 'Staurosporine']
# Read in summary data
proteins = pd.read_excel(f'{input_proteins}', sheet_name='max_proteins')
proteins.drop([col for col in proteins.columns.tolist()
              if 'Unnamed: ' in str(col)], axis=1, inplace=True)
proteins = proteins.replace(0, np.nan).dropna(how='all', subset=treatments)
proteins['degree'] = proteins[treatments].count(axis=1)

for_plotting = pd.melt(proteins, id_vars=['Proteins', 'degree'], var_name='treatment', value_vars=treatments, value_name='max_cys_ratio')

# Read in cluster info
ppi_clusters = pd.read_csv(input_clusters)
ppi_clusters.drop([col for col in ppi_clusters.columns.tolist() if 'Unnamed: ' in col], axis=1, inplace=True)
ppi_clusters['cluster'] = [cluster if cluster in [1.0, 2.0, 3.0, 4.0, 6.0] else 0 for cluster in ppi_clusters['glay_cluster']] # select clustered proteins, all else assigned to 'unclustered' 0
ppi_clusters = dict(ppi_clusters[['Proteins', 'cluster']].values)

# Read in per-complex correlations
complex_correlations = pd.read_csv(input_complexes)
complex_correlations.drop([col for col in complex_correlations.columns.tolist() if 'Unnamed: ' in col], axis=1, inplace=True)

# Read in average interaction correlations
interactions_correlations = pd.read_csv(input_interactions)
interactions_correlations.drop([col for col in interactions_correlations.columns.tolist(
) if 'Unnamed: ' in col], axis=1, inplace=True)
interactions_correlations = pd.melt(interactions_correlations, id_vars=['Proteins'], value_name='corr_mag', value_vars=['0.0', '1.0'], var_name='interaction')
# ---------Visualise linear regression ~ correlation---------
palette = {
    'Celastrol': '#0055d4',
    'Novobiocin': '#d40000',
    'MG132': '#ff6600',
    'Ver155008': '#ffcc00',
    'Staurosporine': '#217821'
}


# fig, ax = plt.subplots(figsize=(5, 8))
fig = sns.lmplot(
    data=for_plotting,
    x='degree',
    y='max_cys_ratio',
    hue='treatment',
    scatter=False,
    palette=palette,
)
sns.regplot(
    data=for_plotting,
    x='degree',
    y='max_cys_ratio',
    color='black',
    scatter=False,
    line_kws={'linestyle': '--'}
)
sns.stripplot(
    data=for_plotting,
    x='degree',
    y='max_cys_ratio',
    hue='treatment',
    palette=palette,
    alpha=0.5
)
plt.legend(bbox_to_anchor=(1.0, 1.0))
plt.xlim(-0.5, 5.5)
plt.ylabel('Corrected cys ratio')
plt.xlabel('Intersection degree')
plt.savefig(f'{output_folder}per_treatment_regplot.svg')
plt.show()

# -------------Plot cluster-sorted correlation heatmap-------------
corr_df = proteins.fillna(0).copy()
corr_df['cluster'] = corr_df['Proteins'].map(ppi_clusters)

treatment_corr = corr_df[treatments].corr()
# Draw the heatmap
g = sns.clustermap(treatment_corr, center=0, cmap="vlag",
                #    row_colors=network_colors, col_colors=network_colors,
                   dendrogram_ratio=(.1, .2),
                   cbar_pos=(.02, .32, .03, .2),
                   linewidths=.75, figsize=(12, 13))

g.ax_row_dendrogram.remove()

# ----------------repeat for clustered correlations----------------
protein_corr = corr_df.set_index(['Proteins', 'cluster']).T.corr()
clusters = corr_df['cluster'].unique().tolist()

# Create a categorical palette to identify the clusters
cluster_pal = sns.husl_palette(len(clusters), s=.45)
cluster_lut = dict(zip(clusters, cluster_pal))

# Convert the palette to vectors that will be drawn on the side of the matrix
clusters = protein_corr.columns.get_level_values("cluster")
clusters_colors = pd.Series(
    clusters, index=protein_corr.columns).map(cluster_lut)

g = sns.clustermap(protein_corr, center=0, cmap="vlag",
                   row_colors=clusters_colors, col_colors=clusters_colors,
                   dendrogram_ratio=(.1, .2),
                   cbar_pos=(.02, .32, .03, .2),
                   linewidths=.05, figsize=(50, 50))

g.ax_row_dendrogram.remove()

# --------------Plot degree-sorted correlation heatmap--------------
degree_df = proteins.fillna(0).copy()
protein_corr = degree_df.set_index(['Proteins', 'degree']).T.corr().abs()
clusters = degree_df['degree'].unique().tolist()

# Create a categorical palette to identify the clusters
cluster_pal = sns.husl_palette(len(clusters), s=.45)
cluster_lut = dict(zip(clusters, cluster_pal))

# Convert the palette to vectors that will be drawn on the side of the matrix
clusters = protein_corr.columns.get_level_values("degree")
clusters_colors = pd.Series(
    clusters, index=protein_corr.columns).map(cluster_lut)

g = sns.clustermap(protein_corr, center=0, cmap="vlag",
                   row_colors=clusters_colors, col_colors=clusters_colors,
                   dendrogram_ratio=(.1, .2),
                   cbar_pos=(.02, .32, .03, .2),
                   linewidths=.05, figsize=(50, 50))

g.ax_row_dendrogram.remove()


# -------Plot boxplots for mean per-protein complex correlations-------
complex_correlations['corr_mag_0'] = complex_correlations['0'].abs()
complex_correlations['corr_mag_1'] = complex_correlations['1'].abs()
for (go_id, complex_name), df in complex_correlations.groupby(['complex_id', 'complex_name']):

    fig, ax = plt.subplots(figsize=(3, 6))
    sns.boxplot(
        data=pd.melt(df, id_vars=['Proteins'], value_vars=[
            'corr_mag_0', 'corr_mag_1'], var_name='complex_membership', value_name='corr_val'),
        x='complex_membership',
        y='corr_val',
        color='darkgrey',
        fliersize=0
    )
    plt.setp(ax.artists, edgecolor='darkgrey', facecolor='w')
    plt.setp(ax.lines, color='darkgrey')
    sns.stripplot(
        data=pd.melt(df, id_vars=['Proteins'], value_vars=[
            'corr_mag_0', 'corr_mag_1'], var_name='complex_membership', value_name='corr_val'),
        x='complex_membership',
        y='corr_val',
        color='black',
        dodge=True,
    )
    plt.xlabel(f'{go_id}: {complex_name}')
    plt.ylabel(r'Mean absolute $r_{s}$ per protein')
    plt.ylim(-0.05, 1.05)
    plt.savefig(f'{output_folder}complex_{go_id.split(":")[1]}.svg')

palette = {'corr_mag_0': 'lightgrey', 'corr_mag_1': 'black',}
for_plotting = pd.melt(complex_correlations, id_vars=['Proteins', 'complex_id', 'complex_name'], value_vars=['corr_mag_0', 'corr_mag_1'], var_name='complex_membership', value_name='corr_val')
complex_names = dict(for_plotting[['complex_id', 'complex_name']].values)

fig, ax = plt.subplots(figsize=(16, 3))
sns.boxplot(
    data=for_plotting,
    x='complex_id',
    y='corr_val',
    hue='complex_membership',
    palette=palette,
    fliersize=0
)
# plt.setp(ax.artists, facecolor='w')
sns.stripplot(
    data=for_plotting,
    x='complex_id',
    y='corr_val',
    hue='complex_membership',
    palette=palette,
    dodge=True,
    s=5
)

ax.axhline(0, linestyle='--', color='grey')

ax_t = ax.secondary_xaxis('top')
ax_t.set_xticklabels(ax.get_xticklabels(), rotation=90)
ax_t.set_xticks(range(len(ax.get_xticklabels())))

ax.set_xticks(range(len(ax.get_xticklabels())))
ax.set_xticklabels([complex_names[x.get_text()] for x in ax.get_xticklabels()], rotation=90)
plt.xlabel('')

plt.ylim(-0.05, 1.05)
plt.ylabel(r'Mean absolute $r_{s}$ per protein')
plt.legend(bbox_to_anchor=(1.0, 1.0))
plt.savefig(f'{output_folder}complex_summary.svg')
plt.show()

# -----Plot summary of corerelation between proteins with or without PPI-----

fig, ax = plt.subplots(figsize=(3, 6))
sns.boxplot(data=interactions_correlations, x='interaction', y='corr_mag', fliersize=0, color='darkgrey',)
plt.setp(ax.artists, edgecolor='darkgrey', facecolor='w')
plt.setp(ax.lines, color='darkgrey')
sns.stripplot(data=interactions_correlations, x='interaction', y='corr_mag', alpha=0.3, color='black')
plt.xlabel('Interaction')
plt.ylabel(r'Mean absolute $r_{s}$ per protein')
plt.ylim(-0.05, 1.05)
plt.savefig(f'{output_folder}interactions.svg')
plt.show()

# -------Plot boxplots for mean per-protein complex correlations-------
cluster_correlations = pd.read_csv(
    'results/tpemi_proteomics/cluster_correlations/protein_mean_comparisons_tests.csv')
cluster_correlations['corr_mag_0'] = cluster_correlations['0'].abs()
cluster_correlations['corr_mag_1'] = cluster_correlations['1'].abs()

fig, ax = plt.subplots(figsize=(3, 6))
sns.boxplot(
    data=pd.melt(cluster_correlations, id_vars=['Proteins'], value_vars=[
        'corr_mag_0', 'corr_mag_1'], var_name='cluster_membership', value_name='corr_val'),
    x='cluster_membership',
    y='corr_val',
    color='darkgrey',
    fliersize=0,
)
plt.setp(ax.artists, edgecolor='darkgrey', facecolor='w')
plt.setp(ax.lines, color='darkgrey')
sns.stripplot(
    data=pd.melt(cluster_correlations, id_vars=['Proteins'], value_vars=[
        'corr_mag_0', 'corr_mag_1'], var_name='cluster_membership', value_name='corr_val'),
    x='cluster_membership',
    y='corr_val',
    color='black',
    dodge=True,
    alpha=0.3
)
plt.xlabel('Cluster')
plt.ylabel(r'Mean absolute $r_{s}$ per protein')
plt.ylim(-0.05, 1.05)
plt.savefig(f'{output_folder}cluster_correlation.svg')
plt.show()

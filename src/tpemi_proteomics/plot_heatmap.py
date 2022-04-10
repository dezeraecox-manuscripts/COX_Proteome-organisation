
import os, functools
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

from GEN_Utils import FileHandling
from loguru import logger

logger.info("Import ok")

input_proteins = 'results/tpemi_proteomics/peptide_summary/peptide_summary.xlsx'
input_clusters = 'results/tpemi_proteomics/plot_degree_ppis/cytoscape_node_summary.csv'
output_folder = 'results/tpemi_proteomics/plot_heatmap/'

if not os.path.exists(output_folder):
    os.makedirs(output_folder)    

font = {'family': 'normal',
        'weight': 'normal',
        'size': 12}
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

# Read in cluster info
ppi_nodes = pd.read_csv(input_clusters)
ppi_nodes.drop([col for col in ppi_nodes.columns.tolist()
                  if 'Unnamed: ' in col], axis=1, inplace=True)
# select clustered proteins, all else assigned to 'unclustered' 0
ppi_nodes['cluster'] = [cluster if cluster in [1.0, 2.0, 3.0,
                                               4.0, 6.0] else 0 for cluster in ppi_nodes['glay_cluster']]
ppi_clusters = dict(ppi_nodes[['Proteins', 'cluster']].values)
cluster_map = {0: 0, 1: 4, 2: 2, 3: 1, 4: 5, 6: 3} # organise according to text structure
proteins['cluster'] = proteins['Proteins'].map(ppi_clusters).map(cluster_map)
ppi_nodes['display_name'] = ppi_nodes['display_name'].str.upper()
display_names = dict(ppi_nodes[['Proteins', 'display_name']].values)


treatments = ['Celastrol', 'Novobiocin', 'MG132', 'Ver155008', 'Staurosporine']
for_plotting = proteins.copy()
for_plotting['sort'] = for_plotting[treatments].sum(axis=1)


# -------------Visualise clustered heatmap for degree > 2-------------
df = for_plotting[for_plotting['degree'] > 2].copy().set_index(['Proteins', 'cluster', 'degree'])[treatments].fillna(0).T
# Create a categorical palette to identify the degrees
# Convert the palette to vectors that will be drawn on the side of the matrix
degree_lut = dict(
    zip([1, 2, 3, 4, 5], ['#e3d7f4', '#c6afe9', '#aa87de', '#8d5fd3', '#442178', ]))
degree_labels = pd.Series(
    df.columns.get_level_values("degree"), index=df.columns).map(degree_lut)

# Repeat for glay clusters
cluster_lut = dict(
    zip([0, 1, 2, 3, 4, 5], ['#ffffff', '#cccccc', '#b3b3b3', '#808080', '#666666', '#4d4d4d']))

cluster_labels = pd.Series(
    df.columns.get_level_values("cluster"), index=df.columns).map(cluster_lut)


font = {'family': 'normal',
        'weight': 'normal',
        'size': 12}
matplotlib.rc('font', **font)
plt.rcParams['svg.fonttype'] = 'none'

for label_name, labels in zip(['cluster', 'degree'], [cluster_labels, degree_labels]):
    g = sns.clustermap(
        df,
        center=0,
        cmap="RdBu",
        vmin=-1.5, vmax=1.5,
        col_colors=cluster_labels,
        dendrogram_ratio=(.01, .02),
        cbar_pos=(.02, .32, .03, .2),
        linewidths=.05,
        figsize=(40, 5)
    )

    g.ax_col_dendrogram.remove()
    plt.setp(g.ax_heatmap.xaxis.set_ticklabels([display_names[protein.get_text(
    ).split('-')[0]] for protein in g.ax_heatmap.xaxis.get_ticklabels()]))
    plt.savefig(f'{output_folder}heatmap_agglomerative_clustered_{label_name}.svg')

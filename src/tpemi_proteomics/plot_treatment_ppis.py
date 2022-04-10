""" NOTE: current versions of the py4cytoscape require cytoscape to be open and running."""

import os
import subprocess
import time
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from py4cytoscape import cytoscape_system
import seaborn as sns
import matplotlib
import py4cytoscape as cyto

from loguru import logger

logger.info('Import OK')

input_treated = 'results/tpemi_proteomics/peptide_summary/peptide_summary.xlsx'
output_folder = 'results/tpemi_proteomics/plot_treatment_ppis/'

if not os.path.exists(output_folder):
    os.makedirs(output_folder)

font = {'family': 'normal',
        'weight': 'normal',
        'size': 14}
matplotlib.rc('font', **font)
plt.rcParams['svg.fonttype'] = 'none'


def cytoscape_STRING_import(proteins, species='Mus Musculus', confidence=0.4, add_interactors=0):
    proteins = ','.join(proteins)
    # Collect network from STRING
    string_cmd = f'string protein query species="{species}" query="{proteins}" cutoff="{confidence}" limit="{add_interactors}"'
    cyto.commands_run(string_cmd)

    # Get table of mapped nodes
    nodes = cyto.get_table_columns()

    return nodes


def cytoscape_map_table(df, df_key='Proteins', node_key='query term'):
    # Map onto node table
    cyto.load_table_data(df, data_key_column=df_key,
                         table_key_column=node_key)
    # nodes var is not automatically updated with new info
    nodes = cyto.get_table_columns()
    return nodes


def cytoscape_create_style(style_name, color_col, map_vals=[-1, 0, 1], map_colors=['#05008a', '#adadad', '#a80000']):

    # Create new style with some mappings
    defaults = {'NODE_SHAPE': "circle", 'EDGE_TRANSPARENCY': 70,
                'NODE_LABEL_POSITION': "W,E,c,0.00,0.00", "NODE_BORDER_WIDTH": 2, "NODE_FILL_COLOR": '#ffffff', "NODE_SIZE": 10}
    node_labels = cyto.map_visual_property(
        visual_prop='node label', table_column='display name', mapping_type='p')  # 'p' means 'passthrough' mapping
    node_color = cyto.map_visual_property(
        visual_prop='node fill color', table_column=color_col, mapping_type='c',
        table_column_values=map_vals, visual_prop_values=map_colors)  # 'p' means 'passthrough' mapping, col_vals are the extremes to map, vis_prop are the colors
    cyto.create_visual_style(style_name, defaults, [node_labels, node_color])


# def cytoscape_create_style(style_name, color_col, map_vals=[-1, 0, 1], map_colors=['#05008a', '#adadad', '#a80000']):

#     # Create new style with some mappings
#     defaults = {'NODE_SHAPE': "circle", 'EDGE_TRANSPARENCY': 70,
#                 'NODE_LABEL_POSITION': "W,E,c,0.00,0.00", "NODE_BORDER_WIDTH": 2, "NODE_FILL_COLOR": '#ffffff', "NODE_SIZE": 20}
#     node_labels = cyto.map_visual_property(
#         visual_prop='node label', table_column='display name', mapping_type='p')  # 'p' means 'passthrough' mapping
#     node_color = cyto.map_visual_property(
#         visual_prop='node fill color', table_column=color_col, mapping_type='c',
#         table_column_values=map_vals, visual_prop_values=map_colors)  # 'p' means 'passthrough' mapping, col_vals are the extremes to map, vis_prop are the colors
#     node_size = cyto.map_visual_property(
#         visual_prop='node size', table_column='degree', mapping_type='c',
#         table_column_values=[0, 25], visual_prop_values=[10, 250])  # 'p' means 'passthrough' mapping, col_vals are the extremes to map, vis_prop are the colors
#     cyto.create_visual_style(style_name, defaults, [
#                              node_labels, node_color, node_size])


# Check cytoscape has been started
# start_cyto = subprocess.Popen(
#     r"C:/Program Files/Cytoscape_v3.9.0/cytoscape.exe", shell=True)
cyto.cytoscape_ping()

# Read in summary data
protein_summary = pd.read_excel(f'{input_treated}', sheet_name='summary')
protein_summary.drop([col for col in protein_summary.columns.tolist()
              if 'Unnamed: ' in str(col)], axis=1, inplace=True)
protein_summary = pd.pivot_table(
    protein_summary,
    index=['Proteins', 'Sequence'],
    columns='treatment',
    values='log2_thresh_pval_ratio')   
protein_summary.reset_index(inplace=True)
# Prepare per-protein max (but NOT removing proteins not quantified in all treatments)
treatments = ['MG132', 'Celastrol', 'Ver155008', 'Novobiocin', 'Staurosporine']
max_proteins = []
for protein, df in protein_summary.groupby('Proteins'):
    if len(df) == 1:
        max_proteins.append(df[['Proteins']+treatments])
    else:
        treatment_vals = []
        for treatment in treatments:
            protein_vals = df[treatment].dropna().sort_values().tolist()
            if len(protein_vals) == 0:
                treatment_vals.append(np.nan)
            else:
                max_val = protein_vals[-1] if abs(protein_vals[-1]) > abs(
                    protein_vals[0]) else protein_vals[0]
                treatment_vals.append(max_val)
        max_vals = pd.DataFrame([[protein] + treatment_vals])
        max_vals.columns = ['Proteins'] + treatments

        max_proteins.append(max_vals)
max_proteins = pd.concat(max_proteins)
protein_summary = pd.melt(
    max_proteins,
    id_vars=['Proteins'],
    value_vars=treatments,
    value_name='max_log2_pval_thresh',
    var_name='treatment'
)

network_ids = {}
tables = {}
for treatment, df in protein_summary.groupby('treatment'):

    treatment_df = df[['Proteins', 'max_log2_pval_thresh']].copy()
    # Remove proteins not outside control thresholds
    treatment_df['max_log2_pval_thresh'] = treatment_df['max_log2_pval_thresh'].replace(0, np.nan)
    treatment_df.dropna(inplace=True)
    treatment_df['sort_col'] = treatment_df['max_log2_pval_thresh'] * -1
    # Get STRING map for proteins of interest
    proteins = treatment_df['Proteins'].unique().tolist()

    cytoscape_STRING_import(proteins, species='Mus Musculus',
                            confidence=0.4, add_interactors=0)

    cytoscape_map_table(treatment_df, df_key='Proteins', node_key='query term')

    edges = cyto.get_table_columns('edge')
    interactions = edges[['name']].copy()
    interactions[['protein_a', 'protein_b']
                 ] = interactions['name'].str.split(' \(pp\) ', expand=True)
    interaction_degree = pd.melt(interactions, id_vars='name', value_vars=['protein_a', 'protein_b'], value_name='protein', var_name='degree')
    interaction_degree = interaction_degree.groupby(
        'protein').count()['degree'].reset_index()

    tables[treatment] = interaction_degree

    cytoscape_map_table(interaction_degree, df_key='protein', node_key='display name')

    if len(network_ids) < 1:
        # Create a new style for range of vals in peptides
        cytoscape_create_style(style_name='log2_pVal_color', color_col='max_log2_pval_thresh', map_vals=[-1, 0, 1], map_colors=['#a10000', '#d1d1d1', '#063478', ])
        cyto.style_mappings.set_node_size_mapping(
            'degree', table_column_values=[0, 25], sizes=[10, 100], mapping_type='c', style_name='log2_pVal_color')
    cyto.set_visual_style('log2_pVal_color')

    logger.info(f'Network created for {treatment}')
    # cyto.notebook_export_show_image()

    cyto.networks.rename_network(treatment)
    network_ids[treatment] = cyto.get_network_suid()


for network in network_ids:
    cyto.set_current_network(f'{network}')
    cyto.export_image(f'{network}',
                    type='SVG', network=f'{network}')
    # copy image file to Notebook directory
    cyto.sandbox_get_from(
        f'{network}.svg', f'{output_folder}{network}.svg')


# Save current version, overwriting previous versions of the session!
# copy session file from Notebook directory to workstation
cyto.save_session('treatment_STRING_layout')
cyto.sandbox_get_from('treatment_STRING_layout.cys',
                      f'{output_folder}treatment_STRING_layout.cys')


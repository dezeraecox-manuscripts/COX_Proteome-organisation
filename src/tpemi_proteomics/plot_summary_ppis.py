""" NOTE: current versions of the py4cytoscape require cytoscape to be open and running."""

import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import py4cytoscape as cyto

from loguru import logger

logger.info('Import OK')

input_intersections = 'results/tpemi_proteomics/intersections/protein_intersection_degree.csv'
input_path = 'results/tpemi_proteomics/peptide_summary/peptide_summary.xlsx'
output_folder = 'results/tpemi_proteomics/plot_degree_ppis/'

if not os.path.exists(output_folder):
    os.makedirs(output_folder)    

# Check cytoscape has been started
# start_cyto = subprocess.Popen(
#     r"C:/Program Files/Cytoscape_v3.9.0/cytoscape.exe", shell=True)
cyto.cytoscape_ping()

font = {'family' : 'normal',
'weight' : 'normal',
'size'   : 14 }
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
    defaults = {'NODE_SHAPE': "circle", 'NODE_SIZE': 10,
                'EDGE_TRANSPARENCY': 70, 'NODE_LABEL_POSITION': "W,E,c,0.00,0.00"}
    node_labels = cyto.map_visual_property(
        visual_prop='node label', table_column='display name', mapping_type='p')  # 'p' means 'passthrough' mapping
    node_color = cyto.map_visual_property(
        visual_prop='node fill color', table_column=color_col, mapping_type='c',
        table_column_values=map_vals, visual_prop_values=map_colors)  # 'p' means 'passthrough' mapping, col_vals are the extremes to map, vis_prop are the colors
    cyto.create_visual_style(style_name, defaults, [node_labels, node_color])


# -------------------------------Prepare quantitative data-------------------------------
treatments = ['Celastrol', 'Novobiocin', 'MG132', 'Ver155008', 'Staurosporine']
# Read in summary data
proteins = pd.read_csv(f'{input_intersections}')
proteins.drop([col for col in proteins.columns.tolist() if 'Unnamed: ' in str(col)], axis=1, inplace=True)
proteins.sort_values('Degree')
degree_dict = pd.melt(proteins, id_vars=['Proteins', 'Degree'], value_vars=treatments, value_name='changed', var_name='treatment')
degree_dict = {protein: str(degree) if degree > 1 else treatment for protein, degree, treatment in degree_dict[degree_dict['changed']][['Proteins', 'Degree', 'treatment']].values}

# Read in quantitative data from summary, which contains the mean and max per protein
quant_data = pd.read_excel(f'{input_path}', sheet_name='changed_protein_summary')
quant_data.drop([col for col in quant_data.columns.tolist() if 'Unnamed: ' in col], axis=1, inplace=True)

proteins_max = quant_data.copy().dropna()
proteins_max = pd.pivot_table(
    proteins_max, values='max_log2_pval_thresh', index='Proteins', columns='treatment').reset_index()
# Select only proteins with at least one sign change 
proteins_max_change = proteins_max.replace(0, np.nan).dropna(how='all', subset=treatments)

# Map quant data onto degree information
proteins_max_change['degree_type'] = proteins_max_change['Proteins'].map(degree_dict)
proteins_max_change['direction'] = ['up' if all(val > 0 for val in vals[~np.isnan(vals)]) else (
    'down' if all(val < 0 for val in vals[~np.isnan(vals)]) else 'mixed') for vals in proteins_max_change[treatments].values]

# ---------------------------Prepare STRING PPI map--------------------------- 

# Get STRING map for proteins of interest
protein_ids = proteins_max_change['Proteins'].unique().tolist()

cytoscape_STRING_import(protein_ids, species='Mus Musculus',
                        confidence=0.4, add_interactors=0)

# replace NaNs with 0.0 to enable mapping
cytoscape_map_table(proteins_max_change.replace(np.nan, 0.0), df_key='Proteins', node_key='query term')


# Create a new style for range of vals in peptides
# Create new style with some mappings
defaults = {'NODE_SHAPE': "circle", 'EDGE_TRANSPARENCY': 70,
            'NODE_LABEL_POSITION': "W,E,c,0.00,0.00", "NODE_BORDER_WIDTH": 2, "NODE_FILL_COLOR": '#e3d7f4', "NODE_SIZE": 20}
node_labels = cyto.map_visual_property(
    visual_prop='node label', table_column='display name', mapping_type='p')  # 'p' means 'passthrough' mapping
node_color = cyto.map_visual_property(
    visual_prop='node fill color', 
    table_column='degree_type', 
    mapping_type='d',
    table_column_values=['2', '3', '4', '5'], 
    visual_prop_values=['#c6afe9', '#aa87de', '#8d5fd3', '#442178']
    )
node_border = cyto.map_visual_property(
    visual_prop='node border color',
    table_column='degree_type',
    mapping_type='d',
    table_column_values=treatments,
    visual_prop_values=['#0055d4ff', '#d40000', '#ff6600', '#ffcc00', '#217821']
)

cyto.create_visual_style('degree_type', defaults, [node_labels, node_color, node_border])
cyto.style_mappings.set_node_border_width_mapping('degree_type', table_column_values=treatments, widths=[10]*len(treatments), style_name='degree_type', mapping_type='d')
cyto.style_mappings.set_node_size_mapping('degree_type', table_column_values=['2', '3', '4', '5'], sizes=[30, 40, 50, 60], style_name='degree_type', mapping_type='d')


cyto.set_visual_style('degree_type')
# cyto.notebook_export_show_image()


cyto.networks.rename_network('degree')

# NOTE: "Create new clustered network" toggle should be activated manually in cytoscape from within the ClusterMaker ->Community (Glay) menu, as there appears to be no option to programmatically enable this. 
cyto.commands.commands_run(
    'cluster glay clusterAttribute="glay" ')
# cyto.notebook_export_show_image()

# Get cluster data from table and remap onto node table (this should be automatic but appears broken)
nodes = cyto.get_table_columns('node')
glay_clusters = cyto.commands.commands_run(
    'node get attribute columnList="query term,glay"')[0]
glay_clusters = pd.DataFrame(glay_clusters.split('Node suid: ')).drop(0)
glay_clusters[['suid', 'discard', 'Proteins', 'glay_cluster']] = glay_clusters[0].str.split(':', expand=True)
glay_clusters.drop('discard', axis=1, inplace=True)
glay_clusters['Proteins'] = glay_clusters['Proteins'].str.replace(',glay', '')
glay_clusters['glay_cluster'] = glay_clusters['glay_cluster'].str.replace('}', '').str.replace(',', '').astype(int)
proteins_max_change['glay_cluster'] = proteins_max_change['Proteins'].map(dict(glay_clusters[['Proteins', 'glay_cluster']].values))
proteins_max_change['display_name'] = proteins_max_change['Proteins'].map(dict(nodes[['query term', 'display name']].values))
proteins_max_change.to_csv(f'{output_folder}cytoscape_node_summary.csv')

# At this point, switch to the cytoscape window and apply 'Layout -> yFiles Organic'
# note this is not available via API due to proprietary issues.....
# This also means that network will be reordered every time the layout is applied


# Select orphan node
edges = cyto.get_table_columns('edge')
edges[['Protein A', 'Protein B']] = edges['shared name'].str.split(r' \(pp\) ', expand=True)
edges.to_csv(f'{output_folder}cytoscape_edge_summary.csv')


# select orphan nodes (single cluster, non-connected nodes) and apply stacked layout to produce linear layout
orphans = nodes[~nodes['display name'].isin((edges['Protein A'].tolist() + edges['Protein B'].tolist()))]
cyto.select_nodes(orphans['display name'].tolist(), by_col='display name')

# Save session file
# copy session file from Notebook directory to workstation
cyto.save_session('degree_layout')
cyto.sandbox_get_from('degree_layout.cys', f'{output_folder}degree_layout.cys')
# Save image
cyto.export_image('degree', type='SVG', network='degree')
# copy image file to Notebook directory
cyto.sandbox_get_from('degree.svg', f'{output_folder}degree.svg')
# create copies of the network to be coloured according to changes in each treatment

for treatment in treatments:
    cyto.clone_network('degree--clustered')
    import time
    time.sleep(10)
    cyto.networks.rename_network(f'{treatment}')

# Create new style with some mappings
cyto.copy_visual_style('degree_type', 'treatment_max')
cyto.set_node_color_default('#ffffff', style_name='treatment_max')
cyto.set_node_border_width_default(2, style_name='treatment_max')
cyto.set_node_border_color_default('#000000', style_name='treatment_max')
cyto.style_mappings.delete_style_mapping(style_name='treatment_max', visual_prop='node border color')
cyto.style_mappings.delete_style_mapping('treatment_max', 'node border width')

for treatment in treatments:
    cyto.set_current_network(treatment)
    time.sleep(10)
    cyto.set_visual_style('treatment_max')
    # set color mapping - based on proteins_max_change[treatments].describe()
    # more than 50% proteins fall well within the -1 -> 1 range so map within this for now
    cyto.style_mappings.set_node_color_mapping(treatment, table_column_values=[-1, 0, 1], colors=['#a10000', '#ffffff', '#063478'], style_name='treatment_max', mapping_type='c')
    logger.info(treatment)
    # Save image
    cyto.export_image(f'degree_{treatment}', type='SVG', network=f'{treatment}')
    # copy image file to Notebook directory
    cyto.sandbox_get_from(f'degree_{treatment}.svg', f'{output_folder}degree_{treatment}.svg')


# Save current version, overwriting previous versions of the session!
# copy session file from Notebook directory to workstation
# cyto.save_session('degree_STRING_layout')
# cyto.sandbox_get_from('degree_STRING_layout.cys', f'{output_folder}degree_STRING_layout.cys')

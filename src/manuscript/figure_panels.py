import os
from shutil import copyfile

from loguru import logger
logger.info('Import OK')

input_folder = 'results/'
output_folder = 'results/manuscript/figures/'

if not os.path.exists(output_folder):
    os.makedirs(output_folder)


def collect_panels(panel_paths, output_path):

    if not os.path.exists(output_path):
        os.mkdir(output_path)

    for position, panel in panel_paths.items():
        old_path = f'{input_folder}{panel}'
        __, file_extension = os.path.splitext(old_path)
        new_path = f'{output_path}panel_{position}{file_extension}'
        new_path
        copyfile(old_path, new_path)


def collect_supp(supp_paths, output_path):

    if not os.path.exists(output_path):
        os.mkdir(output_path)
    for file_name, file_path in supp_paths.items():
        old_path = f'{input_folder}{file_path}'
        new_path = f'{output_path}{file_name}.{file_path.split(".")[-1]}'
        copyfile(old_path, new_path)


"""Collect individual svg versions of figure panels and deposit here for editing. ONLY cosmetic edits were performed manually once inside this folder."""


# ---------------------------------Figure 1---------------------------------
output_path = f'{output_folder}F1_fluorescence/'
panel_paths = {
    'B': 'tpemi_fluorescence/plot_summary/boxplot_fluorescence.svg', 
}

collect_panels(panel_paths, output_path)

# ---------------------------------Figure 2---------------------------------
output_path = f'{output_folder}F2_case-study-MG132/'
panel_paths = {
    'A': 'tpemi_proteomics/plot_caldera/MG132.svg',
    'B': 'tpemi_proteomics/plot_treatment_ppis/MG132.svg',
    'C': 'tpemi_proteomics/plot_go_enrichment/MG132_bubblechart.svg',
}

collect_panels(panel_paths, output_path)


# ---------------------------------Figure 3---------------------------------
output_path = f'{output_folder}F3_modifier-summary/'
panel_paths = {
    'A': 'tpemi_proteomics/intersections/protein_upset.svg', 
    'B': 'tpemi_proteomics/plot_degree_ppis/degree.svg', 
    'C1': 'tpemi_proteomics/plot_heatmap/heatmap_agglomerative_clustered_degree.svg', 
    'C2': 'tpemi_proteomics/plot_heatmap/heatmap_agglomerative_clustered_cluster.svg'
}

collect_panels(panel_paths, output_path)


# ---------------------------------Figure 4---------------------------------
output_path = f'{output_folder}F4_complex-correlations/'
panel_paths = {
    'A': 'tpemi_proteomics/plot_correlations/cluster_correlation.svg',
    'B': 'tpemi_proteomics/plot_correlations/interactions.svg', 
    'C': 'tpemi_proteomics/plot_correlations/complex_0032991.svg', 
    'D': 'tpemi_proteomics/plot_correlations/complex_summary.svg', 
    # 'E1': 'tpemi_proteomics/pymol_complexes/proteasome/proteasome_6MSB_MG132.png', 
    # 'E2': 'tpemi_proteomics/pymol_complexes/proteasome/proteasome_6MSB_Ver155008.png',
    # 'E3': 'tpemi_proteomics/pymol_complexes/proteasome/proteasome_6MSB_Novobiocin.png',
}

collect_panels(panel_paths, output_path)

# --------------------------------Supp Figure 1--------------------------------
output_path = f'{output_folder}S1_summary-caldera/'
panel_paths = {
    'A': 'tpemi_proteomics/plot_caldera/Ver155008.svg',
    'B': 'tpemi_proteomics/plot_caldera/Novobiocin.svg',
    'C': 'tpemi_proteomics/plot_caldera/Staurosporine.svg',
    'D': 'tpemi_proteomics/plot_caldera/Celastrol.svg',
}

collect_panels(panel_paths, output_path)

# --------------------------------Supp Figure 2--------------------------------
output_path = f'{output_folder}S2_summary-go/'
panel_paths = {
    'A': 'tpemi_proteomics/plot_go_enrichment/degree_Bio Process_bubblechart.svg',
    'B': 'tpemi_proteomics/plot_go_enrichment/degree_Mol Function_bubblechart.svg',
    'C': 'tpemi_proteomics/plot_go_enrichment/degree_Cell Component_bubblechart.svg',
}

collect_panels(panel_paths, output_path)


# --------------------------------Supp Figure 3--------------------------------
output_path = f'{output_folder}S3_constellations/'
panel_paths = {
    'A': 'tpemi_proteomics/plot_degree_ppis/degree_MG132.svg',
    'B': 'tpemi_proteomics/plot_degree_ppis/degree_Ver155008.svg',
    'C': 'tpemi_proteomics/plot_degree_ppis/degree_Novobiocin.svg',
    'D': 'tpemi_proteomics/plot_degree_ppis/degree_Staurosporine.svg',
    'E': 'tpemi_proteomics/plot_degree_ppis/degree_Celastrol.svg',
}

collect_panels(panel_paths, output_path)

# --------------------------------Supp Figure 4--------------------------------
output_path = f'{output_folder}S4_sui-correlations/'
panel_paths = {
    'A1': 'solubility_comparison/intersections/dcVxs_MG132_all.svg',
    'A2': 'solubility_comparison/intersections/dcVxs_MG132_changed.svg',
    'A3': 'solubility_comparison/correlation/regression_MG132.svg',
    'B1': 'solubility_comparison/intersections/dcVxs_Novobiocin_all.svg',
    'B2': 'solubility_comparison/intersections/dcVxs_Novobiocin_changed.svg',
    'B3': 'solubility_comparison/correlation/regression_Novobiocin.svg',
    'C1': 'solubility_comparison/intersections/dcVxs_Ver155008_all.svg',
    'C2': 'solubility_comparison/intersections/dcVxs_Ver155008_changed.svg',
}


collect_panels(panel_paths, output_path)

# --------------------------------Supp Figure 5--------------------------------
output_path = f'{output_folder}S5_treatment-correlations/'
panel_paths = {
    'A': 'tpemi_proteomics/plot_correlations/per_treatment_regplot.svg',
}

collect_panels(panel_paths, output_path)

# --------------------------------Supp Figure 6--------------------------------
output_path = f'{output_folder}S6_pval-data-thresholding/'
panel_paths = {
    'A': 'tpemi_proteomics/plot_caldera/TPE.svg',
    'B': 'tpemi_proteomics/plot_caldera/distribution_TPE_cys.svg', 
    'C': 'tpemi_proteomics/plot_caldera/distribution_TPE_noncys.svg', 
}

collect_panels(panel_paths, output_path)


"""Collect supp info and deposit here for editing. ONLY cosmetic edits were be performed manually once inside this folder."""

# -----------------------Supplementary datasets-----------------------

# Supp Table 1: Statistical information summary
output_path = f'{output_folder}ST1_Statistics-summary/'
supp_paths = {
    'Fig_4A': 'tpemi_proteomics/cluster_correlations/cluster_correlation_summary.csv',
    'Fig_4B': 'tpemi_proteomics/ppi_correlations/ppi_correlation_summary.csv',
    'Fig_4C': 'tpemi_proteomics/complex_correlations/complex_correlation_summary.csv',
    'Fig_4D': 'tpemi_proteomics/complex_correlations/complex_correlation_summary.csv',
    }

collect_supp(supp_paths, output_path)

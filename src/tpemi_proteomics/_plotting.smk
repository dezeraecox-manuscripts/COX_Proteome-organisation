
rule proteomics_plot_caldera:
    input: 
        'results/tpemi_proteomics/peptide_summary/peptide_summary.xlsx',
        'results/tpemi_proteomics/baseline/significance/cys_scaled_summary.xlsx',
        'results/tpemi_proteomics/baseline_thresholding/thresholded_summary.xlsx',
        'results/tpemi_proteomics/baseline_thresholding/noncys_thresholded_summary.xlsx',
        'results/tpemi_proteomics/baseline/peptide_normalisation/normalised_summary.xlsx'
    output: 
        'results/tpemi_proteomics/plot_caldera/distribution_TPE_cys.svg',
        'results/tpemi_proteomics/plot_caldera/distribution_TPE_noncys.svg',
        expand("results/tpemi_proteomics/plot_caldera/{treatment}.svg", treatment=['MG132', 'Celastrol', 'Ver155008', 'Novobiocin', 'Staurosporine', 'TPE', 'NOTPE'])
    script:
        "plot_caldera.py"



rule proteomics_plot_correlations:
    input: 
        'results/tpemi_proteomics/peptide_summary/peptide_summary.xlsx',
        'results/tpemi_proteomics/plot_degree_ppis/cytoscape_node_summary.csv',
        'results/tpemi_proteomics/complex_correlations/protein_mean_comparisons_tests.csv',
        'results/tpemi_proteomics/ppi_correlations/ppi_correlation_vals.csv',
        'results/tpemi_proteomics/cluster_correlations/protein_mean_comparisons_tests.csv'
    output: 
        'results/tpemi_proteomics/plot_correlations/per_treatment_regplot.svg',
        'results/tpemi_proteomics/plot_correlations/complex_summary.svg',
        'results/tpemi_proteomics/plot_correlations/interactions.svg',
        'results/tpemi_proteomics/plot_correlations/cluster_correlation.svg',
    script:
        "plot_correlations.py"


rule proteomics_plot_go_enrichment:
    input: 
        'results/tpemi_proteomics/go_enrichment/go_enrichment.csv',
        'results/tpemi_proteomics/go_enrichment_degree_clusters/go_enrichment.csv'
    output: 
        expand('results/tpemi_proteomics/plot_go_enrichment/degree_{search_type}_bubblechart.svg', search_type=['Mol Function', 'Bio Process', 'Cell Component']),
        expand('results/tpemi_proteomics/plot_go_enrichment/{treatment}_bubblechart.svg', treatment=['MG132', 'Ver155008', 'Novobiocin', 'Staurosporine', 'Celastrol']),
    script:
        "plot_go_enrichment.py"


rule proteomics_plot_heatmap:
    input: 
        'results/tpemi_proteomics/peptide_summary/peptide_summary.xlsx',
        'results/tpemi_proteomics/plot_degree_ppis/cytoscape_node_summary.csv'
    output: 
        'results/tpemi_proteomics/plot_heatmap/heatmap_agglomerative_clustered_cluster.svg', 
        'results/tpemi_proteomics/plot_heatmap/heatmap_agglomerative_clustered_degree.svg', 
    script:
        "plot_heatmap.py"


rule proteomics_plot_summary_ppis:
    input: 
        'results/tpemi_proteomics/intersections/protein_intersection_degree.csv',
        'results/tpemi_proteomics/peptide_summary/peptide_summary.xlsx'
    output: 
        # 'results/tpemi_proteomics/plot_degree_ppis/degree_STRING_layout.cys',
        'results/tpemi_proteomics/plot_degree_ppis/degree_layout.cys',
        'results/tpemi_proteomics/plot_degree_ppis/cytoscape_node_summary.csv'

    script:
        "plot_summary_ppis.py"


rule proteomics_plot_treatment_ppis:
    input: 
        'results/tpemi_proteomics/peptide_summary/peptide_summary.xlsx'
    output: 
        'results/tpemi_proteomics/plot_treatment_ppis/treatment_STRING_layout.cys'
    script:
        "plot_treatment_ppis.py"
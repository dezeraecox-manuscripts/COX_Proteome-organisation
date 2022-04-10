rule proteomics_cluster_correlations:
    input: 
        'results/tpemi_proteomics/correlations/max_per_protein_correlations.csv',
        'results/tpemi_proteomics/plot_degree_ppis/cytoscape_node_summary.csv'
    output: 
        'results/tpemi_proteomics/cluster_correlations/protein_mean_comparisons.csv',
        'results/tpemi_proteomics/cluster_correlations/protein_mean_comparisons_tests.csv',
        'results/tpemi_proteomics/cluster_correlations/cluster_correlation_summary.csv',
        'results/tpemi_proteomics/cluster_correlations/cluster_corr_vals.csv'
        
    script:
        "cluster_correlations.py"


rule proteomics_complex_correlations:
    input: 
        'results/tpemi_proteomics/correlations/max_per_protein_correlations.csv',
        'results/tpemi_proteomics/protein_complexes/ontology_summary.xlsx'
    output: 
        'results/tpemi_proteomics/complex_correlations/protein_mean_comparisons.csv',
        'results/tpemi_proteomics/complex_correlations/protein_mean_comparisons_tests.csv',
        'results/tpemi_proteomics/complex_correlations/complex_correlation_summary.csv',
        'results/tpemi_proteomics/complex_correlations/complex_corr_vals.csv'
    script:
        "complex_correlations.py"


rule proteomics_complex_enrichment:
    input: 
        'results/tpemi_proteomics/peptide_summary/peptide_summary.xlsx',
        'results/tpemi_proteomics/protein_complexes/ontology_summary.xlsx',
        'results/tpemi_proteomics/peptide_normalisation/normalised_summary.xlsx'
    output: 
        'results/tpemi_proteomics/complex_enrichment/complex_fishers_enrichment.csv'
    script:
        "complex_enrichment.py"


rule proteomics_degree_correlations:
    input: 
        'results/tpemi_proteomics/intersections/protein_intersection_degree.csv',
        'results/tpemi_proteomics/peptide_summary/peptide_summary.xlsx'
    output: 
        'results/tpemi_proteomics/degree_correlations/proteins_max_change.csv',
        'results/tpemi_proteomics/degree_correlations/degree_spearmans_correlations.csv'
    script:
        "degree_correlations.py"


rule proteomics_go_enrichment:
    input: 
        'results/tpemi_proteomics/peptide_summary/peptide_summary.xlsx',
        'results/tpemi_proteomics/peptide_normalisation/normalised_summary.xlsx',
        'results/tpemi_proteomics/plot_degree_ppis/cytoscape_node_summary.csv'
    output: 
        'results/tpemi_proteomics/go_enrichment/go_enrichment.csv',
        'results/tpemi_proteomics/go_enrichment_degree_clusters/go_enrichment.csv'
    script:
        "go_enrichment.py"


rule proteomics_intersections:
    input: 
        'results/tpemi_proteomics/peptide_summary/peptide_summary.xlsx'
    output: 
        'results/tpemi_proteomics/intersections/protein_upset.svg',
        'results/tpemi_proteomics/intersections/peptide_intersection_degree.csv',
        'results/tpemi_proteomics/intersections/protein_intersection_degree.csv'

    script:
        "intersections.py"


rule proteomics_pathway_correlations:
    input: 
        'results/tpemi_proteomics/correlations/max_per_protein_correlations.csv',
        'results/tpemi_proteomics/protein_pathways/KEGG_summary.csv'
    output: 
        'results/tpemi_proteomics/pathway_correlations/pathway_correlations.csv'
    script:
        "pathway_correlations.py"


rule proteomics_ppi_correlations:
    input: 
        'results/tpemi_proteomics/correlations/max_per_protein_correlations.csv',
        'results/tpemi_proteomics/protein_interactions/STRING_protein_interactions.xlsx'
    output: 
        'results/tpemi_proteomics/ppi_correlations/ppi_correlation_vals.csv',
        'results/tpemi_proteomics/ppi_correlations/ppi_correlation_summary.csv'
    script:
        "ppi_correlations.py"


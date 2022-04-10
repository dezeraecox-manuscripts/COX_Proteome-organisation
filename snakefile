import os

rule all:
    input:
        [
            # --------------------------------Resources--------------------------------
            'resources/bioinformatics_databases/10090_idmapping.tab',
            'resources/bioinformatics_databases/10090_idmapping.dat.gz',
            'resources/bioinformatics_databases/10090.tab',
            'resources/bioinformatics_databases/PANTHERGOslim.obo',

            'results/tpemi_fluorescence/plot_summary/boxplot_fluorescence.svg', 
            'results/tpemi_proteomics/plot_caldera/MG132.svg',
            'results/tpemi_proteomics/plot_treatment_ppis/MG132.svg',
            'results/tpemi_proteomics/plot_go_enrichment/MG132_bubblechart.svg',
            'results/tpemi_proteomics/intersections/protein_upset.svg', 
            'results/tpemi_proteomics/plot_degree_ppis/degree.svg', 
            'results/tpemi_proteomics/plot_heatmap/heatmap_agglomerative_clustered_degree.svg', 
            'results/tpemi_proteomics/plot_heatmap/heatmap_agglomerative_clustered_cluster.svg',
            'results/tpemi_proteomics/plot_correlations/cluster_correlation.svg',
            'results/tpemi_proteomics/plot_correlations/interactions.svg', 
            'results/tpemi_proteomics/plot_correlations/complex_0032991.svg', 
            'results/tpemi_proteomics/plot_correlations/complex_summary.svg', 
            'results/tpemi_proteomics/plot_caldera/Ver155008.svg',
            'results/tpemi_proteomics/plot_caldera/Novobiocin.svg',
            'results/tpemi_proteomics/plot_caldera/Staurosporine.svg',
            'results/tpemi_proteomics/plot_caldera/Celastrol.svg',
            'results/tpemi_proteomics/plot_go_enrichment/degree_Bio Process_bubblechart.svg',
            'results/tpemi_proteomics/plot_go_enrichment/degree_Mol Function_bubblechart.svg',
            'results/tpemi_proteomics/plot_go_enrichment/degree_Cell Component_bubblechart.svg',
            'results/tpemi_proteomics/plot_degree_ppis/degree_MG132.svg',
            'results/tpemi_proteomics/plot_degree_ppis/degree_Ver155008.svg',
            'results/tpemi_proteomics/plot_degree_ppis/degree_Novobiocin.svg',
            'results/tpemi_proteomics/plot_degree_ppis/degree_Staurosporine.svg',
            'results/tpemi_proteomics/plot_degree_ppis/degree_Celastrol.svg',
            'results/solubility_comparison/intersections/dcVxs_MG132_all.svg',
            'results/solubility_comparison/intersections/dcVxs_MG132_changed.svg',
            'results/solubility_comparison/correlation/regression_MG132.svg',
            'results/solubility_comparison/intersections/dcVxs_Novobiocin_all.svg',
            'results/solubility_comparison/intersections/dcVxs_Novobiocin_changed.svg',
            'results/solubility_comparison/correlation/regression_Novobiocin.svg',
            'results/solubility_comparison/intersections/dcVxs_Ver155008_all.svg',
            'results/solubility_comparison/intersections/dcVxs_Ver155008_changed.svg',
            'results/tpemi_proteomics/plot_correlations/per_treatment_regplot.svg',
            'results/tpemi_proteomics/plot_caldera/TPE.svg',
            'results/tpemi_proteomics/plot_caldera/distribution_TPE_cys.svg', 
            'results/tpemi_proteomics/plot_caldera/distribution_TPE_noncys.svg', 
            'results/tpemi_proteomics/cluster_correlations/cluster_correlation_summary.csv',
            'results/tpemi_proteomics/ppi_correlations/ppi_correlation_summary.csv',
            'results/tpemi_proteomics/complex_correlations/complex_correlation_summary.csv',
            'results/tpemi_proteomics/complex_correlations/complex_correlation_summary.csv',
        ]
    shell:
        'echo "Initiating all rules"'


include: "src/tpemi_fluorescence/_workflow.smk"

include: "src/tpemi_proteomics/_preprocessing.smk"
include: "src/tpemi_proteomics/_analysis.smk"
include: "src/tpemi_proteomics/_plotting.smk"

include: "src/solubility_comparison/_workflow.smk"



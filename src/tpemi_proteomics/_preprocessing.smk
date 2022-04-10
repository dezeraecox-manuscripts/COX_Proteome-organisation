
rule proteomics_resources:
    output: 
        'resources/bioinformatics_databases/10090_idmapping.tab',
        'resources/bioinformatics_databases/10090_idmapping.dat.gz',
        'resources/bioinformatics_databases/10090.tab',
        'resources/bioinformatics_databases/PANTHERGOslim.obo'
    script:
        "utilities/databases.py"


rule proteomics_initial_cleanup:
    input: 
        'data/tpemi_proteomics/baseline/peptides.txt', 'data/tpemi_proteomics/baseline/proteinGroups.txt', 'data/tpemi_proteomics/treatments/peptides.txt', 'data/tpemi_proteomics/treatments/proteinGroups.txt', 
    output: 
        'results/tpemi_proteomics/baseline/initial_cleanup/cleaned_peptides.csv', 
        'results/tpemi_proteomics/baseline/initial_cleanup/cleaned_proteins.csv', 
        'results/tpemi_proteomics/initial_cleanup/cleaned_peptides.csv', 
        'results/tpemi_proteomics/initial_cleanup/cleaned_proteins.csv'
    script:
        "1_initial_cleanup.py"


rule proteomics_peptide_normalisation:
    input: 
        'results/tpemi_proteomics/baseline/initial_cleanup/cleaned_peptides.csv', 
        'results/tpemi_proteomics/baseline/initial_cleanup/cleaned_proteins.csv', 
        'results/tpemi_proteomics/initial_cleanup/cleaned_peptides.csv', 
        'results/tpemi_proteomics/initial_cleanup/cleaned_peptides.csv', 
    output: 
        'results/tpemi_proteomics/baseline/peptide_normalisation/normalised_summary.xlsx',
        'results/tpemi_proteomics/peptide_normalisation/normalised_summary.xlsx'
    script:
        "2_peptide_normalisation.py"


rule proteomics_significance_and_scaling:
    input: 
        'results/tpemi_proteomics/baseline/peptide_normalisation/normalised_summary.xlsx',
        'results/tpemi_proteomics/peptide_normalisation/normalised_summary.xlsx'
    output: 
        'results/tpemi_proteomics/baseline/significance/cys_scaled_summary.xlsx',
        'results/tpemi_proteomics/baseline/significance/noncys_scaled_summary.xlsx',
        'results/tpemi_proteomics/significance/cys_scaled_summary.xlsx',
        'results/tpemi_proteomics/significance/noncys_scaled_summary.xlsx'
    script:
        "3_significance_and_scaling.py"


rule proteomics_baseline_thresholding:
    input: 
        'results/tpemi_proteomics/baseline/significance/cys_scaled_summary.xlsx',
        'results/tpemi_proteomics/baseline/significance/noncys_scaled_summary.xlsx',
        'results/tpemi_proteomics/significance/cys_scaled_summary.xlsx',
        'results/tpemi_proteomics/significance/noncys_scaled_summary.xlsx'        
    output: 
        'results/tpemi_proteomics/baseline_thresholding/noncys_thresholded_summary.xlsx',
        'results/tpemi_proteomics/baseline_thresholding/thresholded_summary.xlsx'
    script:
        "4_baseline_thresholding.py"


rule proteomics_summary:
    input: 
        'results/tpemi_proteomics/peptide_normalisation/normalised_summary.xlsx',
        'results/tpemi_proteomics/baseline_thresholding/thresholded_summary.xlsx',
        'results/tpemi_proteomics/baseline_thresholding/noncys_thresholded_summary.xlsx'
    output: 
        'results/tpemi_proteomics/peptide_summary/peptide_summary.xlsx'
    script:
        "5_summary.py"


rule proteomics_per_protein_correlations:
    input: 
        'results/tpemi_proteomics/peptide_summary/peptide_summary.xlsx'
    output: 
        'results/tpemi_proteomics/correlations/max_per_protein_correlations.csv',
        'results/tpemi_proteomics/correlations/treatment_max_protein_correlations.csv',
    script:
        "per_protein_correlations.py"


rule proteomics_protein_complexes:
    input: 
        'results/tpemi_proteomics/significance/cys_scaled_summary.xlsx',
        'resources/bioinformatics_databases/PANTHERGOslim.obo'
    output: 
        'results/tpemi_proteomics/protein_complexes/ontology_summary.xlsx'
    script:
        "protein_complexes.py"


rule proteomics_protein_interactions:
    input: 
        'results/tpemi_proteomics/peptide_normalisation/normalised_summary.xlsx'
    output: 
        'results/tpemi_proteomics/protein_interactions/STRING_protein_interactions.xlsx'
    script:
        "protein_interactions.py"


rule proteomics_protein_pathways:
    output: 
        'results/tpemi_proteomics/protein_pathways/KEGG_pathways.xlsx',
        'results/tpemi_proteomics/protein_pathways/KEGG_summary.csv',
    script:
        "protein_pathways.py"


# rule proteomics_0_raw_data:
#     input: 
        
#     output: 
        
#     script:
#         "0_raw_data.py"



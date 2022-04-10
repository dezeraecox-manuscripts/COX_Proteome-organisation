
rule comparison_0_raw_data:
    input: 
        'results/tpemi_proteomics/peptide_summary/peptide_summary.xlsx'
    output: 
        'data/solubility_comparison/dc_datasets/peptide_summary.xlsx',
        'data/solubility_comparison/xs_datasets/pnas.1912897117.sd04.xlsx',
        'data/solubility_comparison/xs_datasets/pnas.1912897117.sd05.xlsx',
        'data/solubility_comparison/xs_datasets/pnas.1912897117.sd06.xlsx',
        
    script:
        "0_raw_data.py"


rule comparison_1_initial_cleanup:
    input: 
        'data/solubility_comparison/xs_datasets/pnas.1912897117.sd04.xlsx',
        'data/solubility_comparison/xs_datasets/pnas.1912897117.sd05.xlsx',
        'data/solubility_comparison/xs_datasets/pnas.1912897117.sd06.xlsx',
        'data/solubility_comparison/dc_datasets/peptide_summary.xlsx',
        'resources/bioinformatics_databases/10090.tab',
        'resources/bioinformatics_databases/10090_idmapping.tab',
    output: 
        'results/solubility_comparison/initial_cleanup/compiled_summary.csv'
        
    script:
        "1_initial_cleanup.py"


rule comparison_2_correlation:
    input: 
        'results/solubility_comparison/initial_cleanup/compiled_summary.csv'
        
    output: 
        expand('results/solubility_comparison/correlation/regression_{drug}.svg', drug=['Ver155008', 'Novobiocin', 'MG132'])
    script:
        "2_correlation.py"


rule comparison_3_intersections:
    input: 
        'results/solubility_comparison/initial_cleanup/compiled_summary.csv'
        
    output: 
        expand('results/solubility_comparison/intersections/dcVxs_{drug}_all.svg', drug=['Ver155008', 'Novobiocin', 'MG132']),
        expand('results/solubility_comparison/intersections/dcVxs_{drug}_changed.svg', drug=['Ver155008', 'Novobiocin', 'MG132'])
    script:
        "3_intersections.py"


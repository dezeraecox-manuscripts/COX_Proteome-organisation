# rule fluorescence_0_raw_data:
#     input: 
        
#     output: 
        
#     script:
#         "0_raw_data.py"


rule fluorescence_initial_cleanup:
    input: 
        [f'data/tpemi_fluorescence/{filename}' for filename in os.listdir('data/tpemi_fluorescence/')]
    output: 
        "results/tpemi_fluorescence/initial_cleanup/compiled.csv"
    script:
        "1_initial_cleanup.py"


rule fluorescence_comparison:
    input: 
        "results/tpemi_fluorescence/initial_cleanup/compiled.csv"
    output: 
        [
            "results/tpemi_fluorescence/comparison/anova.csv",
            "results/tpemi_fluorescence/comparison/multi_comparison.csv",
            "results/tpemi_fluorescence/comparison/ttest_result.csv",
        ]
    script:
        "2_comparison.py"


rule fluorescence_plot_summary:
    input: 
        'results/tpemi_fluorescence/initial_cleanup/compiled.csv'
    output: 
        'results/tpemi_fluorescence/plot_summary/boxplot_fluorescence.svg'
    script:
        "plot_summary.py"

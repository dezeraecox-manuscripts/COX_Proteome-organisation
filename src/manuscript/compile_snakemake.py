import os, re
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

from loguru import logger

logger.info('Import OK')

input_folder = 'src/'
output_folder = 'results/manuscript/compile_snakemake/'

if not os.path.exists(output_folder):
    os.makedirs(output_folder)


# list of python scripts
paths = [os.path.join(root, filename).replace('\\', '/') for root, dirs, files in os.walk(input_folder, topdown=False) for filename in files if 'utilities' not in root]

paths = [path for path in paths if 'snakemake' not in path]

# Read each script and generate example rules
components = [
    "rule all:\n",
    "    input:\n",
    "        ['results/one.csv', \n",
    "        'results/two.csv',  \n",
    "        'results/three.svg', \n",
    "        'results/four.png']\n",
    "    shell:\n",
    "        'echo \"Initiating all rules\"'\n\n\n"
]
for filepath in paths:
    script_name = filepath.split('/')[-1].replace('.py', '')
    if 'fluorescence' in filepath:
        components.extend([f'rule fluorescence_{script_name}:\n'])
    elif 'proteomics' in filepath:
        components.extend([f'rule proteomics_{script_name}:\n'])
    elif 'comparison' in filepath:
        components.extend([f'rule comparison_{script_name}:\n'])
    elif 'manuscript' in filepath:
        components.extend([f'rule manuscript_{script_name}:\n'])
    components.extend(['    input: \n        ', '\n    output: \n        ', '\n    script:\n        ', f'"{filepath}"\n\n\n'])


# Write to  template file with some basic snakemake boilerplate
with open(f'{output_folder}snakefile', 'w') as f:
    for line in components:
        f.write(line)


# -----------CURATE SNAKEMAKE FILE WITH INPUTS AND OUTPUTS FOR EACH STAGE-----------
# Move compiled.txt into root directory and relabel as 'snakemake'
# Run using ```snakemake -p --cores 1``` at the command line

## NOTE: The compile file must be edited manually to ensure accurate curation of the inputs and outputs for each script!
## NOTE: A initial rule "all" is added at the top of the script - it should contain as input all the output files, and prints a success message:
# rule all:
#     input:
#         ["results/one.csv", 
#         "results/two.csv",  
#         "results/three.svg", 
#         "results/four.png"]
#     shell:
#         'echo "Initiating all rules"'
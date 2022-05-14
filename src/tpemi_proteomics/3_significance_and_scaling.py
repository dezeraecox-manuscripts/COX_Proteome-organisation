import os
import numpy as np
import pandas as pd

from GEN_Utils import FileHandling

# Snakemake requires a truncated version of the relative import - to use directly, use complete path
# from src.proteomics.utilities import smoothing_and_scaling as sms
from utilities import smoothing_and_scaling as sms

from loguru import logger
logger.info('Import OK')


def main(input_path, output_folder, sheet_name='cys_noncys_peptides'):
    
    # read in raw data
    raw_data = pd.read_excel(f'{input_path}', sheet_name=sheet_name)
    compiled = raw_data.copy().drop([col for col in raw_data.columns.tolist() if 'Unnamed: ' in col], axis=1)
    compiled.drop(['cys_rank', 'Gene names', 'Protein names', 'Unique (Groups)', 'Unique (Proteins)'], inplace=True, axis=1)

    raw = compiled.copy()
    compiled['sample_name'] = 'AllStresses'

    # define useful data bits
    if sheet_name == 'noncys_peptides':
        # if provessing the non-cys peptides, we want to smooth the mean per protein
        compiled = compiled.groupby(
            ['Proteins', 'sample_name', 'replicate']).mean().reset_index()
        compiled['Sequence'] = compiled['Proteins']
    protein_map = dict(raw.reset_index()[['Sequence', 'Proteins']].values)
    info_cols = ['Sequence', 'Proteins', 'sample_name', 'replicate']
    sample_cols = list(dict.fromkeys([col for col in compiled.columns.tolist() if col not in info_cols]))

    """--------------SIGNIFICANCE FILTERING-----------------------"""
    # Complete standard one-sample t-test to identify sig. diff from 1
    ttest_results = ttest_results = sms.one_sample_ttest(compiled, sample_cols)
    compiled_pvals = pd.pivot_table(ttest_results, index=['Sequence'], columns=['channel'], values=['p-val'])
    compiled_pvals.columns = [col for _, col in compiled_pvals.columns.tolist()]
    
    significant_changes = compiled.groupby(['Sequence']).mean().dropna()[sample_cols].copy()
    significant_changes = significant_changes.where(compiled_pvals < 0.05).fillna(1)

    significant_changes['Proteins'] = significant_changes.reset_index()['Sequence'].map(protein_map).tolist()
    significant_changes[sample_cols] = np.log2(significant_changes[sample_cols])

    """ ------------------PVAL SCALING------------------"""
    log_compiled = compiled.copy()
    log_compiled[sample_cols] = np.log2(log_compiled[sample_cols])
    pval_smooth = sms.pval_smoothing(log_compiled, sample_cols, center=0)
    # Add back protein info, remove unnecessary columns
    pval_smooth.reset_index(inplace=True)
    pval_smooth['Proteins'] = pval_smooth['Sequence'].map(protein_map)

    # save original and processed data to excel
    data_dict = {'raw': raw, 'pval_smooth': pval_smooth, 'significant_changes': significant_changes}

    FileHandling.df_to_excel(output_path=f'{output_folder}scaled_summary.xlsx', sheetnames=list(data_dict.keys()), data_frames=list(data_dict.values()))


if __name__ == "__main__":

    # -----------------------process stress datasets-----------------------

    input_path = 'results/tpemi_proteomics/peptide_normalisation/normalised_summary.xlsx'
    output_folder = 'results/tpemi_proteomics/significance/'

    filter_cols = ['']
    info_cols = ['Sequence', 'Proteins', 'replicate', 'sample_name']

    if not os.path.exists(output_folder):
        os.makedirs(output_folder)

    main(input_path,
         output_folder=f'{output_folder}cys_', sheet_name='cys_noncys_peptides')
    main(input_path,
         output_folder=f'{output_folder}noncys_', sheet_name='noncys_peptides')

    # -----------------------process baseline datasets-----------------------

    input_path = 'results/tpemi_proteomics/baseline/peptide_normalisation/normalised_summary.xlsx'
    output_folder = 'results/tpemi_proteomics/baseline/significance/'

    filter_cols = ['']
    info_cols = ['Sequence', 'Proteins', 'replicate', 'sample_name']

    if not os.path.exists(output_folder):
        os.makedirs(output_folder)

    main(input_path,
        output_folder=f'{output_folder}cys_', sheet_name='cys_noncys_peptides')
    main(input_path,
         output_folder=f'{output_folder}noncys_', sheet_name='noncys_peptides')

    

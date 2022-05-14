import pandas as pd
import numpy as np
import os

from loguru import logger
logger.info("Import OK")


def mq_peptides(input_path, sample_names=None, pep_cols=[], quant_col='Reporter intensity corrected'):

    peptides = pd.read_table(input_path, sep='\t')
    logger.info(
        f'Imported peptides from {input_path}. {peptides.shape[0]} entries found.')
    cleaned_dfs = {}
    standard_cols = ['Sequence', 'Proteins', 'Gene names',
                     'Protein names', 'Unique (Groups)', 'Unique (Proteins)', ] + pep_cols

    # filter out non-unique, reverse and contaminant peptides
    filtered_pep = peptides[(peptides['Reverse'] != '+') & (
        peptides['Potential contaminant'] != '+') & (peptides['Unique (Proteins)'] == 'yes')]
    # add cys rank
    filtered_pep['cys_rank'] = [
        1 if 'C' in pep else 0 for pep in filtered_pep['Sequence']]

    if sample_names is None:
        logger.info(f'Sample names not set. Collecting all samples.')
        logger.debug(f'Columns found: {peptides.columns.tolist()}')
        sample_names = [x.replace('Experiment ', '').split(
            '_')[:-1] for x in peptides.columns.tolist() if 'Experiment ' in x]
        sample_names = list({('_').join(x) for x in sample_names})
        logger.info(f'Samples detected: {sample_names}')

    for sample in sample_names:
        sample_df = filtered_pep[standard_cols +
                                 [x for x in filtered_pep.columns.tolist() if f' {sample}_' in x]]
        sample_cols = [col for col in sample_df.columns.tolist(
        ) if f'{quant_col} {sample}' in col]
        # filter those without any values for value in variable:
        sample_df = sample_df.replace(0, np.nan).dropna(
            axis=0, how='all', subset=sample_cols)
        sample_df = sample_df[standard_cols + sample_cols]
        cleaned_dfs[sample] = sample_df
    logger.info(f'Successfully cleaned peptide dataframe.')

    return cleaned_dfs, sample_names


def mq_proteins(input_path, sample_names=None, prot_cols=[], quant_col='Reporter intensity corrected'):

    logger.info(f'Collecting proteins')
    proteins = pd.read_table(input_path, sep='\t')
    logger.info(
        f'Imported proteins from {input_path}. {proteins.shape[0]} entries found.')

    # remove contaminant and reverse proteins
    proteins = proteins[(proteins['Reverse'] != '+') &
                        (proteins['Potential contaminant'] != '+')]
    logger.info(
        f'Removed contaminant and reverse proteins: {proteins.shape[0]} entries remain.')

    cleaned_prot_dfs = {}
    standard_cols = ['Protein IDs', 'Gene names',
                     'Protein names', 'Number of proteins'] + prot_cols

    for sample in sample_names:
        sample_cols = standard_cols + \
            [x for x in proteins.columns.tolist() if f' {sample}_' in x]
        sample_df = proteins[sample_cols]
        ratio_cols = [col for col in sample_df.columns.tolist(
        ) if f'{quant_col} {sample}' in col]
        logger.debug(f'Ratio cols: {ratio_cols}')
        #collect columns of interest
        sample_vals = proteins[sample_cols + ratio_cols]
        #collect only proteins with at least one quantification in that sample
        sample_reps = sample_df.replace(0, np.nan).dropna(
            axis=0, how='all', subset=ratio_cols)
        logger.debug(f'Sample reps: {sample_reps.head(10)}')
        # collect only proteins which are master proteins
        master_proteins = sample_reps[sample_reps['Number of proteins'] == 1]
        logger.debug(f'Master proteins: {master_proteins.head(10)}')

        cleaned_prot_dfs[sample] = master_proteins

    logger.info(f'Successfully cleaned proteins dataframe.')

    return cleaned_prot_dfs


def mq_cleaner(input_folder, output_path, sample_names=None, proteins_file='proteinGroups.txt', peptides_file='peptides.txt', prot_cols=[], pep_cols=[], quant_col='Reporter intensity corrected', quant_name='abundance'):

    cleaned_dfs = {}

    cleaned_peptides, sample_names = mq_peptides(input_path=f'{input_folder}{peptides_file}', sample_names=sample_names, pep_cols=pep_cols, quant_col=quant_col)
    cleaned_proteins = mq_proteins(input_path=f'{input_folder}{proteins_file}', sample_names=sample_names, prot_cols=prot_cols, quant_col=quant_col)

    logger.info(f'Sorting cleaned data per sample...')
    ## Collecting specific results for each set of samples for further processing
    for sample in sample_names:
        #collect peptide dataframe, rename relevant columns
        pep_dataframe = cleaned_peptides[sample].copy()
        pep_dataframe.columns = [col.strip(quant_col) if quant_col in col else col for col in pep_dataframe]

        prot_dataframe = cleaned_proteins[sample].copy()
        prot_dataframe.columns = [col.strip(quant_col) if quant_col in col else col for col in prot_dataframe]

    logger.info(f'Proteins and peptides successfully cleaned. Dataframes save to {output_folder}.')

    cleaned_peptides = compiler(cleaned_peptides, quant_col, quant_name)
    cleaned_proteins = compiler(cleaned_proteins, quant_col, quant_name)


    cleaned_peptides.to_csv(f'{output_folder}cleaned_peptides.csv')
    cleaned_proteins.to_csv(f'{output_folder}cleaned_proteins.csv')
    
    return cleaned_peptides, cleaned_proteins


def compiler(cleaned_data, quant_col, quant_name='abundance'):
    compiled_dfs = []
    for sample in cleaned_data.keys():
        df = cleaned_data[sample].copy()
        
        compiled_df = pd.melt(
            df,
            id_vars=[col for col in df.columns.tolist() if quant_col not in col],
            value_vars=[col for col in df.columns.tolist() if quant_col in col],
            var_name='sample_name',
            value_name=quant_name,)
        compiled_df['replicate'] = compiled_df['sample_name'].str.split('_').str[-1]
        compiled_df['sample_name'] = sample
        compiled_dfs.append(compiled_df)
    return pd.concat(compiled_dfs)


if __name__ == "__main__":


    ## Cleanup raw MQ data to generate peptide and protein information
    input_folder = 'data/tpemi_proteomics/treatments/'
    output_folder = 'results/tpemi_proteomics/initial_cleanup/'

    if not os.path.exists(output_folder):
        os.makedirs(output_folder)

    mq_cleaner(input_folder, output_folder, proteins_file='proteinGroups.txt', peptides_file='peptides.txt', prot_cols=[], pep_cols=[], quant_col='Ratio H/L normalized', quant_name='ratio')


    ## Repeat cleanup raw MQ data to generate peptide and protein information for baseline thresholds
    input_folder = 'data/tpemi_proteomics/baseline/'
    output_folder = 'results/tpemi_proteomics/baseline/initial_cleanup/'

    if not os.path.exists(output_folder):
        os.makedirs(output_folder)

    mq_cleaner(input_folder, output_folder, proteins_file='proteinGroups.txt', peptides_file='peptides.txt', prot_cols=[], pep_cols=[], quant_col='Ratio H/L normalized', quant_name='ratio')
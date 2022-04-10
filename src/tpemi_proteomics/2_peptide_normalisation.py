import os
import numpy as np
import pandas as pd
from GEN_Utils import FileHandling
from collections import defaultdict

from loguru import logger

logger.info('Import OK')


def cys_noncys_filter(peptides):
    """ Assigns cys_rank, then separates cys and noncys peptides.
    Collects only peptides from proteins that have at least one cys and one non-cys peptide"""

    # Add cys rank to identify cys peptides
    peptides['cys_rank'] = [
        1 if 'C' in sequence else 0 for sequence in peptides['Sequence']]

    # separate cys and noncys peptides
    cys_peptides = peptides[peptides['cys_rank'] == 1]
    noncys_peptides = peptides[peptides['cys_rank'] == 0]

    # Collect only proteins with at least one cys, one noncys
    common_prots = set(cys_peptides['Proteins']).intersection(
        set(noncys_peptides['Proteins']))
    cys_peptides = cys_peptides[cys_peptides['Proteins'].isin(common_prots)]
    noncys_peptides = noncys_peptides[noncys_peptides['Proteins'].isin(
        common_prots)]

    return cys_peptides, noncys_peptides


def noncys_ci_calculator(noncys_peptides, sample_cols):
    """Calculates relevant info per protein from noncys peptides, including overall mean, SD
    and the per-channel mean (used for cys/noncys calculation)"""

    # for each protein, collect all noncys values and determine mean + SD
    noncys_cis = {}
    for protein, df in noncys_peptides.groupby(['Proteins']):
        values = df[sample_cols].values.flatten()
        # Remove NaN values for calculations
        values = values[~np.isnan(values)]
        noncys_cis[protein] = {'df': df,
                               'noncys_means': df[sample_cols].mean(),
                               'num_peptides': df.shape[0],
                               'values': values,
                               'pop_mean': np.mean(values),
                               'pop_stdev': np.std(values),
                               'num_vals': len(values)}
    return noncys_cis


def cys_ratio(cys_peptides, noncys_cis, sample_cols):
    """Calculates cys/noncys ratio per cys peptide, using the noncys per protein mean"""
    # Generate cys/noncys ratio
    norm_cys = []
    for protein, df in cys_peptides.groupby('Proteins'):
        data = df.copy()
        noncys_vals = noncys_cis[protein]['noncys_means']
        data[sample_cols] = data[sample_cols] / noncys_vals
        norm_cys.append(data)

    cys_noncys_peptides = pd.concat(norm_cys)

    return cys_noncys_peptides


def med_normalise(peptides, control_plex):
    """Calculates correction factor from median of column 1, then normalises all other channels to this factor"""

    medians = peptides.median()
    control_factor = medians[control_plex] / medians

    # normalise to control channel
    peptides = peptides * control_factor

    return peptides


def main(sample, info_cols=[], med_norm=True):

    info_cols = ['Sequence', 'Proteins', 'Gene names',
            'Protein names', 'cys_rank', 'replicate'] + info_cols
    if type(sample) == str:
        raw_data = pd.read_excel(f'{input_folder}{sample}_Compiled.xlsx', sheet_name=None)

        peptides = raw_data['Peptides'].copy().set_index([col for col in raw_data['Peptides'].columns.tolist() if col in info_cols])
        # remove outlier samples as identified in initial cleanup/batch effects
        peptides = peptides.copy().drop(filter_cols, axis=1)
        peptides.drop([col for col in peptides.columns.tolist()
                    if 'Unnamed: ' in col], axis=1, inplace=True)
    elif type(sample) == pd.DataFrame:
        peptides = sample.copy()
    
    else:
        logger.info(f'Datatype not recognised. You provided a sample of type {type(sample)}.')
        return

    peptides.drop([col for col in peptides.columns.tolist()
                if 'Unique ' in col], axis=1, inplace=True)

    replicates = list(dict.fromkeys(([col.split(
        '_')[1] for col in peptides.columns.tolist()])))

    # collecting and sorting individual replicates, calculate cy/noncys
    replicate_normalised = defaultdict(dict)
    for replicate in replicates:
        # collect relevant channels
        replicate_cols = [col for col in peptides.columns.tolist() if f'_{replicate}' in col]
        rep_peptides = peptides[replicate_cols]
        replicate_normalised[replicate]['raw'] = rep_peptides.copy()

        if pooled_plex:
            pooled_col = f'Reporter intensity corrected {pooled_plex} {sample_name}_{replicate}'
            # drop control channel
            rep_peptides = rep_peptides.drop(pooled_col, axis=1)
            replicate_cols.remove(pooled_col)

        # filter out cols seen to be outliers via PCA
        if filter_cols:
            for col in filter_cols:
                if col in rep_peptides.columns.tolist():
                    # drop control channel
                    rep_peptides = rep_peptides.drop(col, axis=1)

        # collect column names for quantitative processes, remove replicate number as this will be added in a column
        sample_cols = [col for col in rep_peptides.columns.tolist() if col not in info_cols]
        rep_peptides.rename(columns={old: new for old, new in zip(sample_cols, [col.split('_')[0] for col in sample_cols])}, inplace=True)
        sample_cols = [col.split('_')[0] for col in sample_cols]

        # median-normalise across all samples
        if med_norm:
            rep_peptides = med_normalise(rep_peptides, control_plex)
            replicate_normalised[replicate]['med_norm'] = rep_peptides.copy()

        # separate cys and noncys, generate ratio
        cys_peptides, noncys_peptides = cys_noncys_filter(
            rep_peptides.reset_index())
        replicate_normalised[replicate]['cys_peptides'] = cys_peptides
        replicate_normalised[replicate]['noncys_peptides'] = noncys_peptides

        # Calculate relevant noncys info
        noncys_cis = noncys_ci_calculator(noncys_peptides, sample_cols)

        # Calculate cys/noncys ratio
        cys_noncys_peptides = cys_ratio(cys_peptides, noncys_cis, sample_cols)
        replicate_normalised[replicate]['cys_noncys_peptides'] = cys_noncys_peptides

        # Calculate noncys/noncys ratio
        noncys_noncys_peptides = cys_ratio(
            noncys_peptides, noncys_cis, sample_cols)
        replicate_normalised[replicate]['noncys_noncys_peptides'] = noncys_noncys_peptides


        # Normalise to control sample such that this is 1 (as there should be no difference in control)
        if control_plex:
            norm_cys_noncys_peptides = cys_noncys_peptides.copy()
            norm_cys_noncys_peptides[sample_cols] = (
                norm_cys_noncys_peptides[sample_cols].T / norm_cys_noncys_peptides[control_plex].T).T
            replicate_normalised[replicate]['control_norm_cys_ratio'] = norm_cys_noncys_peptides


        # remove any peptides quantified in less than quant threshold in this replicate
        dfs = []
        for key, df in replicate_normalised[replicate].items():
            if key in ['cys_peptides', 'noncys_peptides', 'cys_noncys_peptides', 'noncys_noncys_peptides', 'control_norm_cys_ratio']:
                replicate_normalised[replicate][key] = df.dropna(
            subset=sample_cols, thresh=quant_threshold)

        # Add column with replicate info
        for key, df in replicate_normalised[replicate].items():
            df['replicate'] = replicate
            if type(sample) == str:
                df['sample_name'] = sample
            replicate_normalised[replicate][key].update(df)

    # compile replicates
    compiled = {}
    for df_type in replicate_normalised[replicate].keys():
        compiled_df = pd.concat([replicate_normalised[replicate][df_type] for replicate in replicates])
        compiled[df_type] = compiled_df

    # save to excel
    if type(sample) != str:
        sample_name = ''
    else:
        sample_name = f'{sample}_'
    FileHandling.df_to_excel(
        output_path=f'{output_folder}{sample_name}normalised_summary.xlsx', 
        sheetnames=list(compiled.keys()),
        data_frames=list(compiled.values()))

    return compiled


#---------------------------------------------------------------------------------------------

if __name__ == "__main__":

    # -----------------------process stress datasets-----------------------
    input_folder = 'results/tpemi_proteomics/initial_cleanup/'
    output_folder = 'results/tpemi_proteomics/peptide_normalisation/'
    sample_names = ['Novobiocin', 'Ver155008', 'MG132', 'Celastrol', 'Staurosporine']
    filter_cols = []
    pooled_plex = False
    control_plex = False
    quant_threshold = 1

    if not os.path.exists(output_folder):
        os.makedirs(output_folder)

    # Read in raw data
    raw_data = pd.read_csv(f'{input_folder}cleaned_peptides.csv')
    raw_data.drop([col for col in raw_data.columns.tolist() if 'Unnamed: ' in col], axis=1, inplace=True)

    # Select only samples of interest
    raw_data = raw_data[raw_data['sample_name'].isin(sample_names)]

    # format to mimic TMT dataset
    pivot_cols = ['sample_name', 'replicate', 'ratio']
    pivot_data = pd.pivot_table(
        raw_data,
        index=[col for col in raw_data.columns.tolist() if col not in pivot_cols],
        columns=['sample_name', 'replicate'],
        values='ratio'
        )
    pivot_data.columns = ['_'.join([str(name) for name in col]) for col in pivot_data.columns.values]

    summary = main(pivot_data, med_norm=False) #used med norm values from maxquant

    # -----------------------process baseline datasets-----------------------
    input_folder = 'results/tpemi_proteomics/baseline/initial_cleanup/'
    output_folder = 'results/tpemi_proteomics/baseline/peptide_normalisation/'
    sample_names = ['TPE', 'NOTPE']
    filter_cols = []
    pooled_plex = False
    control_plex = False
    quant_threshold = 1

    if not os.path.exists(output_folder):
        os.makedirs(output_folder)

    # Read in raw data
    raw_data = pd.read_csv(f'{input_folder}cleaned_peptides.csv')
    raw_data.drop([col for col in raw_data.columns.tolist() if 'Unnamed: ' in col], axis=1, inplace=True)

    # Select only samples of interest
    raw_data = raw_data[raw_data['sample_name'].isin(sample_names)]

    # format to mimic TMT dataset
    pivot_cols = ['sample_name', 'replicate', 'ratio']
    pivot_data = pd.pivot_table(
        raw_data,
        index=[col for col in raw_data.columns.tolist() if col not in pivot_cols],
        columns=['sample_name', 'replicate'],
        values='ratio'
        )
    pivot_data.columns = ['_'.join([str(name) for name in col]) for col in pivot_data.columns.values]

    summary = main(pivot_data, med_norm=False) #used med norm values from maxquant
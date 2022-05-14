import os
import pandas as pd
import numpy as np

from GEN_Utils import FileHandling

from loguru import logger
logger.info("Import OK")

input_path = 'results/tpemi_proteomics/significance/cys_scaled_summary.xlsx'
thresholds_path = 'results/tpemi_proteomics/baseline/significance/cys_scaled_summary.xlsx'
output_folder = 'results/tpemi_proteomics/baseline_thresholding/'

if not os.path.exists(output_folder):
    os.makedirs(output_folder)

def distribution_z_threshold(data, z_val=1.96):
    """
    Calculate the lower and upper values encompassing a given proportion of the population).

    Common vals:
    ============
    0.5: 38.2%
    1.0: 68.2%
    1.5: 86.6%
    1.96: 95%
    2.0: 95.4%
    2.5: 98.8%
    2.58: 99%
    3.0: 99.8%

    For more explanation: https://upload.wikimedia.org/wikipedia/commons/thumb/2/25/The_Normal_Distribution.svg/1024px-The_Normal_Distribution.svg.png
    """
    range_val = z_val * np.std(data)
    return np.mean(data) - range_val, np.mean(data) + range_val

def main(input_path, thresholds_path, output_folder, z_val=1.96):

    raw_data = pd.read_excel(f'{input_path}', sheet_name=None)
    raw_data.update({key: df.drop([col for col in df.columns.tolist() if 'Unnamed: ' in col], axis=1) for key, df in raw_data.items()})

    # ---------------------Calculate CIs, apply as threshold---------------------
    # Calculate C.I. thresholds for the control TPE-positive and TPE-negative samples 
    if input_path == thresholds_path:
        CIs = []
        for data_type, df in raw_data.items():
            if data_type != 'raw':
                for sample in ['TPE', 'NOTPE']:
                    lower, upper = distribution_z_threshold(df[f'{sample}'].dropna().tolist(), z_val=z_val)
                    CIs.append([data_type, sample, lower, upper])
        CIs = pd.DataFrame(CIs)
        CIs.columns = ['data_type', 'sample', 'lower_CI', 'upper_CI']
    else:
        ci_data = pd.read_excel(f'{thresholds_path}', sheet_name=None)
        ci_data.update({key: df.drop([col for col in df.columns.tolist() if 'Unnamed: ' in col], axis=1) for key, df in ci_data.items()})
        CIs = []
        for data_type, df in ci_data.items():
            if data_type != 'raw':
                for sample in ['TPE', 'NOTPE']:
                    lower, upper = distribution_z_threshold(df[f'{sample}'].dropna().tolist(), z_val=z_val)
                    CIs.append([data_type, sample, lower, upper])
        CIs = pd.DataFrame(CIs)
        CIs.columns = ['data_type', 'sample', 'lower_CI', 'upper_CI']


    # Add thresholds to relevant datasets
    sample_cols = ['Celastrol', 'MG132', 'Novobiocin', 'Staurosporine', 'Ver155008',]
    thresholds = {}
    for data_type, df in raw_data.items():
        if data_type != 'raw':
            thresholded = df.set_index(['Sequence', 'Proteins'])[sample_cols].copy()
            thresholded = pd.melt(thresholded.reset_index(), id_vars=['Sequence', 'Proteins'], value_vars=sample_cols, var_name='treatment', value_name='log2_mean_ratio')
            for sample in ['TPE', 'NOTPE']:
                (lower_ci, upper_ci) = CIs[(CIs['sample'] == sample) & (CIs['data_type'] == data_type)][['lower_CI', 'upper_CI']].values[0]
                thresholded[f'{sample}_lower'] = lower_ci
                thresholded[f'{sample}_upper'] = upper_ci

                # assign thresholded changes as those outside threshold limits
                thresholded[f'{sample}_thresholded'] = [1 if lower_ci > value or value > upper_ci else 0 for value in thresholded['log2_mean_ratio'] ]

            thresholds[data_type] = thresholded

    thresholds['CIs'] = CIs

    # ---------------------------Save summaries---------------------------
    FileHandling.df_to_excel(
        output_path=f'{output_folder}thresholded_summary.xlsx',
        sheetnames=list(thresholds.keys()),
        data_frames=list(thresholds.values())
    )

if __name__ == '__main__':


    output_folder = 'results/tpemi_proteomics/baseline_thresholding/'

    if not os.path.exists(output_folder):
        os.makedirs(output_folder)

    # process cys smoothed values
    input_path = 'results/tpemi_proteomics/significance/cys_scaled_summary.xlsx'
    thresholds_path = 'results/tpemi_proteomics/baseline/significance/cys_scaled_summary.xlsx'
    main(input_path, thresholds_path, output_folder=f'{output_folder}', z_val=1.96)

    # process noncys smoothed values
    input_path = 'results/tpemi_proteomics/significance/noncys_scaled_summary.xlsx'
    thresholds_path = 'results/tpemi_proteomics/baseline/significance/noncys_scaled_summary.xlsx'
    main(input_path, thresholds_path, output_folder=f'{output_folder}noncys_', z_val=1.96)

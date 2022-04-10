import os
import shutil
from shutil import copyfile
from loguru import logger
from contextlib import closing
import urllib.request as request


logger.info('Import OK')


def download_resources(filename, url, output_folder):
    """
    Worker function to download and save file from URL.
    
    inputs
    ======
    filename: (str) name of output file (including extension)
    url: (str) complete location of file to be downloaded
    output_path: (str) relative or complete path to directory where folder will be saved.
    returns:
    ======
    None
    """
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)

    try:
        with closing(request.urlopen(url)) as r:
            with open(f'{output_folder}{filename}', 'wb') as f:
                shutil.copyfileobj(r, f)
        logger.info(f'Downloaded {filename}')
    except:
        logger.info(f'Downloaded failed for {filename}.')



if __name__ == "__main__":

    # # ----------Collect XS proteomics datasets----------

    url_root = 'https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7007570/bin/'
    files = [
        'pnas.1912897117.sd04.xlsx',
        'pnas.1912897117.sd05.xlsx',
        'pnas.1912897117.sd06.xlsx'
    ]
    output_folder = 'data/solubility_comparison/xs_datasets/'

    if not os.path.exists(output_folder):
        os.makedirs(output_folder)

    # Download files from repository - note this currently does not work without authentication due to journal access policies
    for filename in files:
        download_resources(
            filename=f'{filename.split("/")[-1]}',
            url=f'{url_root}{filename}', 
            output_folder=output_folder) 


    # ------------Collect DC summary results------------

    input_path = 'results/tpemi_proteomics/peptide_summary/peptide_summary.xlsx'
    output_folder = 'data/solubility_comparison/dc_datasets/'

    if not os.path.exists(output_folder):
        os.makedirs(output_folder)

    copyfile(input_path, f'{output_folder}peptide_summary.xlsx')
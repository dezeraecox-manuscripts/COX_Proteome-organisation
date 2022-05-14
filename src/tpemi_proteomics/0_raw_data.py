import os
import zipfile
import shutil
from contextlib import closing
import urllib.request as request
from src.tpemi_proteomics.utilities import databases

from loguru import logger
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

    # Download static bioinformatics databases of interest
    resource_folder='resources/bioinformatics_databases/'
    databases.download_databases(tax_ids={'MOUSE': '10090'}, resource_folder=resource_folder) 
    
    # Download preprocessed raw data from repository
    url = 'https://zenodo.org/record/6439170/files/tpemi_proteomics.zip?download=1'
    folder_name = 'tpemi_proteomics'
    output_folder = 'data/'

    if not os.path.exists(output_folder):
        os.makedirs(output_folder)
    
    download_resources(filename=f'{folder_name}.zip', url=url, resource_folder=output_folder) 
    with zipfile.ZipFile(f'{output_folder}{folder_name}.zip', 'r') as zip_ref:
        zip_ref.extractall(f'{output_folder}')
    

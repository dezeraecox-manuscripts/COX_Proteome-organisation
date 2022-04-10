import os, re
import zipfile
from shutil import copyfile
import shutil
from loguru import logger
import json
from contextlib import closing
import urllib.request as request
import requests

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


    output_folder = 'experimental_data/'

    if not os.path.exists(output_folder):
        os.makedirs(output_folder)


        copyfile(f'', f'{output_folder}')

import os, re
import shutil
import pandas as pd
import numpy as np

from loguru import logger
from GEN_Utils import FileHandling

logger.info('Import OK')

output_folder = 'raw_data/'

if not os.path.exists(output_folder):
    os.makedirs(output_folder)



# Zenodo dataset

# Bioinformatics resources
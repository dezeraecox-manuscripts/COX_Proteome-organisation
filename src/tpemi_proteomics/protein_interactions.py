import os
import pandas as pd
from loguru import logger
from GEN_Utils import FileHandling

from utilities.databases import network_interactions
from utilities.databases import create_uniprot_xref

logger.info('Import OK')

input_path = 'results/tpemi_proteomics/peptide_normalisation/normalised_summary.xlsx'
output_folder = 'results/tpemi_proteomics/protein_interactions/'
resource_folder = 'resources/bioinformatics_databases/'

confidence_threshold = 0.4

if not os.path.exists(output_folder):
    os.makedirs(output_folder)


# As collecting all PPIs for all proteins can be tricky, and potentially not useful in this case, restricted to all direct PPI between proteins which were quantified here i.e. no first-shell interactions

# -----Read in peptides data for all proteins identified-----

raw_data = pd.read_excel(f'{input_path}', sheet_name='cys_peptides')
proteins = raw_data['Proteins'].unique().tolist()

# Map genes to STRING ids
string_map = create_uniprot_xref(resource_folder, tax_id='10090', gene_ids=[], id_type='STRING')
string_map = dict(string_map[['UniProtKB-AC', 'ID']].values)
string_ids = [string_map[gene] for gene in proteins if gene in list(string_map.keys()) ]

# Colect all network interactions - remember limit is 2000 genes?
interactions = network_interactions(genes=string_ids, tax_id='10090', id_type='string')

# Add back original gene IDs based on mapping
inv_map = {v: k for k, v in string_map.items()}
interactions['Protein_A'] = interactions['stringId_A'].map(inv_map)
interactions['Protein_B'] = interactions['stringId_B'].map(inv_map)

FileHandling.df_to_excel(
    output_path=f'{output_folder}STRING_protein_interactions.xlsx',
    sheetnames=['all_cys_ppis'],
    data_frames=[interactions]
)
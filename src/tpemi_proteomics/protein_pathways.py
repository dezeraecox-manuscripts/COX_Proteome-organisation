import os
import pandas as pd
import numpy as np
import KEGGutils as kg
import requests

from utilities.databases import create_uniprot_xref, uniprot_summary
from utilities.databases import create_uniprot_xref, uniprot_summary

from loguru import logger
from GEN_Utils import FileHandling

logger.info('Import OK')

resource_folder='resources/bioinformatics_databases/'
output_folder = 'results/tpemi_proteomics/protein_pathways/'

if not os.path.exists(output_folder):
    os.mkdir(output_folder)


def find_pathway_genes(pathway_ids, tax_id='10090', resource_folder='resources/bioinformatics_databases/'):
    """Collects and converts all kegg id's for genes associated with a list of kegg terms, outputting mapped UniProt Accessions.

    Parameters
    ----------
    pathway_ids : list
        KEGG pathway ids of the form e.g. 'mmu00010'
    tax_id : str, optional
        taxonomy id as specified by UniProt, by default '10090' for Mus Musculus
    resource_folder : str, optional
        directory where bioinformatics databases are stored. If not already downloaded, download standard databases for tax id of interest using database_collection scripts, by default 'resources/bioinformatics_databases/'

    Returns
    -------
    df
        mapped KEGG pathway term to KEGG genes, mapped to UniProt Accessions where protein is reviewed else NaN
    """
    # Generate orthology database
    uniprot = uniprot_summary(tax_id=tax_id, resource_folder=resource_folder, genes=[], reviewed=True)
    keggo = create_uniprot_xref(input_path=resource_folder, tax_id=tax_id, gene_ids=[], id_type=None) # note some non-unique keys here...
    keggo = keggo[keggo['ID_type'] == 'KEGG']
    uniprot['KEGG_id'] = uniprot['Entry'].map(dict(keggo[['UniProtKB-AC', 'ID']].values))
    kegg_id_map = dict(uniprot[['KEGG_id', 'Entry']].values)

    results = []
    for x, pathway_term in enumerate(pathway_ids):
        pathway_term
        if x % 10 == 0:
            logger.info(f'Processing pathway number {x}')
        try:
            pathway = kg.KEGGpathway(pathway_id = pathway_term)
        except:
            logger.info(f'{pathway_term} not found')

        try:
            genes = [pathway.genes[key]['gene'].strip() for key in pathway.genes.keys()]
            kegg_gene_ids = [kegg_id_map.get(gene, np.nan) for gene in genes]
            genes = pd.DataFrame([genes, kegg_gene_ids], index=['kegg_ids', 'Proteins']).T
            genes['ko_pathway_term'] = pathway_term
            results.append(genes)
        except:
            logger.info(f'Genes for {pathway_term} not found')



    logger.info('Gene collection complete')

    return pd.concat(results)


# Collect all KEGG pathways
response = requests.get("http://rest.kegg.jp/list/pathway/mmu")
pathways = pd.DataFrame([line.split('\t') for line in response.text.split('\n')])
pathways.columns = ['kegg_id', 'name']
pathways['name'] = pathways['name'].str.strip(' - Mus musculus (mouse)')
pathways['kegg_id'] = pathways['kegg_id'].str.strip('path:')
pathways = pathways.dropna()
pathways_dict = dict(pathways[['kegg_id', 'name']].values)

# Find genes associated with each of the KEGG pathway terms
genes = find_pathway_genes(pathway_ids=pathways['kegg_id'].tolist(), tax_id='10090', resource_folder='resources/bioinformatics_databases/')

pathway_summary = genes.dropna().drop_duplicates()
pathway_summary['name'] = pathway_summary['ko_pathway_term'].map(pathways_dict)

pathway_summary.to_csv(f'{output_folder}KEGG_summary.csv')

# Save to excel
FileHandling.df_to_excel(
    output_path=f'{output_folder}KEGG_pathways.xlsx',
    sheetnames=['pathways', 'pathway_summary'],
    data_frames=[pathways, pathway_summary]
)
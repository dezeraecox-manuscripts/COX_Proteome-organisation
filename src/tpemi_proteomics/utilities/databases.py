import gzip
import json
import os
import shutil
import tarfile
import time
import urllib.request as request
from collections import OrderedDict
from contextlib import closing

import numpy as np
import obonet
import pandas as pd
import requests
from GEN_Utils import FileHandling
from loguru import logger
import time
import sys

logger.info('Import OK')

# Certificate expiry error with updated python packages - temporary fix
# import urllib.request
# import ssl

# ssl._create_default_https_context = ssl._create_unverified_context
# response = urllib.request.urlopen('https://www.python.org')
# print(response.read().decode('utf-8'))



class ProgressBarPrinter:
    def __init__(self, width, step, stream, fname):
        self.width = width
        self.block_progress = 0
        self.current_progress = 0
        self.start_time = time.time()
        self.step = step
        self.stream = stream

        # Prints the progress bar's layout
        print(f"{fname}: [{'-' * self.width}] 0%% <0.00 seconds>",
              flush=True, end='', file=self.stream)
        print("\b" * (self.width + len("] 0%% <0.00 seconds>")),
              end='', file=self.stream)

    def update(self, progress):
        # Parsing input
        if progress is None:
            progress = self.current_progress + self.step
        if not isinstance(progress, float):
            raise TypeError("ProgressBar: input must be float or None")

        # Keep the progress bar under 99% until end() has been called
        self.current_progress = min(progress, 0.99)
        self.print_bar(self.current_progress)

    def print_bar(self, progress):
        block = int(round(self.width * progress)) - self.block_progress
        self.block_progress += block
        bar = ('#' * block) + ('-' * (self.width - self.block_progress))
        progress = int(progress * 100)
        elapsed_time = round(time.time() - self.start_time, 2)
        text = f"{bar}] {progress}% <{elapsed_time} seconds>"
        print(text + ("\b" * (len(text) - block)),
              flush=True, end='', file=self.stream)

    def end(self):
        self.print_bar(1.0)
        print(flush=True, file=self.stream)


def ProgressBar(width=50, step=0.1, stream=sys.stdout):
    """Decorator, prints a progress bar when a decored function yields it's
    current progress.

    When you want the progress bar to be updated you should yield the progress
    of your function between 0 and 1. The general calcul for this is:
    (current_iteration + 1) / total_iterations.

    When yielding None, the progress bar goes up by `current progress + step`.
    This is usefull to show some feedback instead of a dead terminal when it's
    not possible to calculate the progress of a function.

    Limitation: It uses yield statements as callbacks for the decorator. That
    means you can't yield your result, meaning this progress bar doesn't
    work if your function is intended to be a generator.
    """
    def decorator(func):
        def wrapper(*args, **kwargs):
            pb = ProgressBarPrinter(width, step, stream, func.__name__)
            progress_generator = func(*args, **kwargs)
            try:
                while True:
                    progress = next(progress_generator)
                    pb.update(progress)
            except StopIteration as result:
                pb.end()
                return result.value
        return wrapper
    return decorator


def download_resources(filename, url, resource_folder):
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
    if not os.path.exists(resource_folder):
        os.makedirs(resource_folder)

    try:
        with closing(request.urlopen(url)) as r:
            with open(f'{resource_folder}{filename}', 'wb') as f:
                shutil.copyfileobj(r, f)
        logger.info(f'Downloaded {filename}')
    except:
        logger.info(f'Downloaded failed for {filename}.')


def chunks(lst, n):
    """Yield successive n-sized chunks from lst."""
    for i in range(0, len(lst), n):
        yield lst[i:i + n]


def uniprot_download(genes, output_folder):

    if not os._exists(output_folder):
        os.makedirs(output_folder)

    for gene in genes:
        url = f'https://www.uniprot.org/uniprot/{gene}.xml'
        response = requests.get(url)
        with open(f'{output_folder}{gene}.xml', 'w') as f:
            f.write(response.text)


def uniprot_features(genes, resource_folder='resources/bioinformatics_databases/'):

    def collect_feature(feature_dict, feature_name):
        try:
            return feature_dict[feature_name]
        except:
            return np.nan

    if not os.path.exists(f'{resource_folder}uniprot_download/'):
        os.makedirs(f'{resource_folder}uniprot_download/')

    protein_features = {}
    genes_unprocessed = []
    for gene in genes:
        try:
            with open(f'{resource_folder}uniprot_download/{gene}.xml') as fd:
                gene_dict = xmltodict.parse(fd.read())
        except FileNotFoundError:
            uniprot_download(genes, f'{resource_folder}uniprot_download/')
            with open(f'{resource_folder}uniprot_download/{gene}.xml') as fd:
                gene_dict = xmltodict.parse(fd.read())
        try:
            feature_names = ['location', '@type', '@description', '@evidence']
            new_features = []
            for feature in gene_dict['uniprot']['entry']['feature']:
                new_dict = {}
                for feature_name in feature_names:
                    new_dict[feature_name] = collect_feature(feature, feature_name)
                new_features.append(new_dict)
            features = pd.DataFrame(new_features)
            features = pd.merge(features[['@type', '@description', '@evidence']], features['location'].apply(pd.Series), left_index=True, right_index=True)
            for location_col in ['begin', 'end', 'position']:
                if location_col in features.columns.tolist():
                    features[f'{location_col}_key'] = [list(pos.keys())[0] if type(pos) == OrderedDict else np.nan for pos in features[location_col] ]
                    features = features[~features[f'{location_col}_key'].astype(str).str.contains('@status')]
                    features[location_col] = [int(pos['@position']) if type(pos) == OrderedDict else np.nan for pos in features[location_col] ]
            if len(features) < 1:
                raise Exception(f"No features found for {gene}.")
            features['entry'] = gene
            cols = ['@type', '@description', '@evidence', 'begin', 'end', 'position', 'entry']
            protein_features[gene] = features[[col for col in features.columns.tolist() if col in cols]].rename(columns={'@type': 'feature_type', '@description': 'feature_description', '@evidence': 'evidence'})
        except:
            logger.info(f'{gene} not processed.')
            genes_unprocessed.append(gene)
    
    return protein_features, genes_unprocessed


def uniprot_function(genes, resource_folder='resources/bioinformatics_databases/'):

    protein_function = {}
    genes_unprocessed = []
    for gene in genes:
        try:
            with open(f'{resource_folder}uniprot_download/{gene}.xml') as fd:
                gene_dict = xmltodict.parse(fd.read())
        except FileNotFoundError:
            uniprot_download(genes, f'{resource_folder}uniprot_download/')
            with open(f'{resource_folder}uniprot_download/{gene}.xml') as fd:
                gene_dict = xmltodict.parse(fd.read())
        try:
            protein_function[gene] = {words['@id'] : words['#text'] for words in gene_dict['uniprot']['entry']['keyword']}
        except:
            logger.info(f'{gene} not processed.')
            genes_unprocessed.append(gene)

    return protein_function


def uniprot_sequences(genes, resource_folder='resources/bioinformatics_databases/'):

    protein_sequence = {}
    genes_unprocessed = []
    for gene in genes:
        try:
            with open(f'{resource_folder}uniprot_download/{gene}.xml') as fd:
                gene_dict = xmltodict.parse(fd.read())
        except:
            uniprot_download(genes, f'{resource_folder}uniprot_download/')
            with open(f'{resource_folder}uniprot_download/{gene}.xml') as fd:
                gene_dict = xmltodict.parse(fd.read())
        try:
            protein_sequence[gene] = gene_dict['uniprot']['entry']['sequence']['#text']
        except:
            logger.info(f'{gene} not processed.')
            genes_unprocessed.append(gene)

    return protein_sequence



def string_geneid_mapper(genes, tax_id, id_limit=1, echo_query=1, string_api_url="https://string-db.org/api", output_format='tsv', caller_id="www.github.com/dezeraecox"):

    def interactors_search(genes):
        # set parameters
        method = "get_string_ids"
        params = [
            f'identifiers={"%0D".join(genes)}', # your proteins
            f'species={tax_id}', # species NCBI identifier 
            f'limit={id_limit}',
            f'caller_identity={caller_id}' # your app name
        ]

        # Construct URL
        request_url = "/".join([string_api_url, output_format, method])
        request_url = f'{request_url}?{"&".join(params)}'

        # Call STRING
        response = requests.get(request_url)

        results = pd.DataFrame(response.text.strip().split("\n"))
        results = results[0].str.strip().str.split("\t", expand=True).T.set_index(0).T
    
        return results


    gene_ids = []
    for number, gene_list in enumerate(chunks(genes, 200)):
        number
        try:
            try:
                gene_ids.append(interactors_search(gene_list))
            except:
                time.sleep(5)
                gene_ids.append(interactors_search(gene_list))
            logger.info(f'Processed chunk number {number}')
        except:
            logger.info(f'{gene_list} not processed.')    

    results = pd.concat(gene_ids)

    # format response into dataframe
    results['string'] = results['stringId'].str.split('.').str[1]

    return results


def network_interactions(genes, tax_id, id_type='string', string_api_url="https://string-db.org/api", output_format='tsv', caller_id="www.github.com/dezeraecox"):

    """
    id_type: STRING=premapped ids to Esembl, UNIPROT ids should be premapped to STRING and converted back after interaction mapping"""

    method = "network"

    ## Construct URL
    request_url = "/".join([string_api_url, output_format, method])

    ## Set parameters
    params = {
        "identifiers" : "%0d".join(genes), # your proteins
        "species" : tax_id, # species NCBI identifier 
        "caller_identity" : caller_id,
        "id_type" : id_type
    }

    ## Call STRING
    response = requests.post(request_url, data=params)

    # format response into dataframe
    results = pd.DataFrame(response.text.strip().split("\n"))
    results = results[0].str.strip().str.split("\t", expand=True).T.set_index(0).T

    return results


def all_interactions( gene_ids, tax_id, max_partners=1000, id_type='string', string_api_url="https://string-db.org/api", output_format='tsv', caller_id="www.github.com/dezeraecox", confidence_threshold=0.7, resource_folder=False):

    """
    id_type: STRING=premapped ids to Esembl, UNIPROT ids will be premapped to STRING
    and converted back after interaction mapping"""    

    def interactors_search(genes, id_limit=max_partners):
        # set parameters
        method = "interaction_partners"
        params = [
            f'identifiers={"%0D".join(genes)}', # your proteins
            f'species={tax_id}', # species NCBI identifier 
            f'limit={id_limit}',
            f'caller_identity={caller_id}' # your app name
        ]

        # Construct URL
        request_url = "/".join([string_api_url, output_format, method])
        request_url = f'{request_url}?{"&".join(params)}'

        # Call STRING
        response = requests.get(request_url)

        results = pd.DataFrame(response.text.strip().split("\n"))
        results = results[0].str.strip().str.split("\t", expand=True).T.set_index(0).T
    
        return results

    if resource_folder:
        if not os.path.exists(f'{resource_folder}STRING_interactions/'):
            os.makedirs(f'{resource_folder}STRING_interactions/')
        prepared_genes = [filename.strip('.csv') for filename in os.listdir(f'{resource_folder}STRING_interactions/')]
    if not id_type == 'string':
        try:
            id_map = string_geneid_mapper(gene_ids, tax_id, id_limit=1, echo_query=1)
            id_map['queryTerm'] = [gene_ids[int(index)] for index in id_map['queryIndex']]
            genes = id_map['string'].tolist()
        except:
            logger.info('Unable to map genes to STRING format.')

    interactions = {}
    for gene in genes:
        if gene in prepared_genes:
            interactions[gene] = pd.read_csv(f'{resource_folder}STRING_interactions/{gene}.csv')
        else:
            try:
                try:
                    interactions[gene] = interactors_search([gene])
                except:
                    time.sleep(5)
                    interactions[gene] = interactors_search([gene])
                # logger.info(f'Processed {gene}')
            except:
                logger.info(f'{gene} not processed.')
    logger.info(f'Interactions db loaded with {len(interactions)} genes.')

    results = pd.concat(interactions.values())
    # restore UNIPROT info
    genes_to_map = results['stringId_A'].unique().tolist() + results['stringId_B'].unique().tolist()
    gene_map = convert_geneids(genes_to_map, tax_id=10090, id_from='STRING', id_to='uniprot')
    gene_map['string_id'] = gene_map['ID'].str.split('.').str[1]
    gene_map = dict(zip(gene_map['string_id'], gene_map['UniProtKB-AC']))
    gene_map.update(dict(zip(id_map['string'], id_map['queryTerm'])))

    #Filter on threshold
    interaction_summary = results[results['score'].astype(float) >= confidence_threshold].copy()
    interaction_summary = interaction_summary[interaction_summary['stringId_A'].str.strip(f'{tax_id}.').isin(genes) | interaction_summary['stringId_B'].str.strip(f'{tax_id}.').isin(genes)]

    interaction_summary['uniprot_A'] = interaction_summary['stringId_A'].str.strip(f'{tax_id}.').map(gene_map)
    interaction_summary['uniprot_B'] = interaction_summary['stringId_B'].str.strip(f'{tax_id}.').map(gene_map)

    # save individual proteins to csv in resources folder
    if resource_folder: 
        for gene, interaction_df in interactions.items():
            interaction_df['uniprot_A'] = interaction_df['stringId_A'].str.strip(f'{tax_id}.').map(gene_map)
            interaction_df['uniprot_B'] = interaction_df['stringId_B'].str.strip(f'{tax_id}.').map(gene_map)

            interaction_df.to_csv(f'{resource_folder}STRING_interactions/{gene_map[gene]}.csv')
    
    return interaction_summary


def interaction_enrichment( genes, tax_id, id_type='string', string_api_url="https://string-db.org/api", output_format='tsv', caller_id="www.github.com/dezeraecox"):

    if not id_type == 'string':
        try:
            id_map = string_geneid_mapper(genes, tax_id, id_limit=1, echo_query=1)
            genes = id_map['string'].tolist()
        except:
            logger.info('Unable to map genes to STRING format.')
    # set parameters
    method = "ppi_enrichment"
    params = {
    "identifiers" : "%0d".join(genes), # your proteins
    "species" : tax_id, # species NCBI identifier 
    "caller_identity" : caller_id # your app name
    }

    # Construct URL
    request_url = "/".join([string_api_url, output_format, method])
    # Call STRING
    response = requests.post(request_url, data=params)

    results = pd.DataFrame(response.text.strip().split("\n"))
    results = results[0].str.strip().str.split("\t", expand=True).T.set_index(0).T

    return results


def disorder_prediction_iupred(accession_ids, output_folder):

    
    disorder = []
    for x, accession in enumerate(accession_ids):
        if x % 10 == 0:
            logger.info(f'Processing protein number {x}')
        # Call IUPreD Database
        request_url = f'http://iupred2a.elte.hu/iupred2a/{accession}.json'

        response = requests.get(request_url)
        results = json.loads(response.text)

        # format disorder into dataframe
        results = pd.DataFrame([results['iupred2'], results['sequence']], index=['disorder_probability', 'sequence']).T.reset_index().rename(columns={'index': 'residue'})
        results['disorder_probability'] = results['disorder_probability'].astype(float)
        results['disordered'] = [1 if probability > 0.5 else 0 for probability in results['disorder_probability']]
        results['Proteins'] = accession
        disorder.append(results)

    return pd.concat(disorder)


def pfam_domains(accession_ids, resource_folder='resources/bioinformatics_databases/'):

    if not os.path.exists(f'{resource_folder}pfam_domains/'):
        os.makedirs(f'{resource_folder}pfam_domains/')

    pfam = []
    for x, accession in enumerate(accession_ids):
        if x % 10 == 0:
            logger.info(f'Processing protein number {x}')
        if not os.path.exists(f'{resource_folder}pfam_domains/{accession}.xml'):
            # Call Database
            request_url = f'https://pfam.xfam.org/protein/{accession}?output=xml'
            response = requests.get(request_url)
            with open(f'{resource_folder}pfam_domains/{accession}.xml', 'w') as f:
                f.write(response.text)

        with open(f'{resource_folder}pfam_domains/{accession}.xml') as fd:
            gene_dict = xmltodict.parse(fd.read())
        try:
            # format disorder into dataframe
            matches = gene_dict['pfam']['entry']['matches']['match']
            if type(matches) == list:
                for match in matches:
                    pfam_acc = match['@accession']
                    pfam_id = match['@id']
                    pfam_type = match['@type']
                    start, end = (match['location']['@start'], match['location']['@end'])
                    pfam.append(pd.DataFrame([accession, pfam_acc, pfam_id, pfam_type, start, end], index=['Proteins', 'pfam_acc', 'pfam_id', 'pfam_type', 'pfam_start', 'pfam_end']).T)
            else:
                pfam_acc = matches['@accession']
                pfam_id = matches['@id'] 
                pfam_type = matches['@type'] 
                start, end = (matches['location']['@start'], matches['location']['@end'])
                pfam.append(pd.DataFrame([accession, pfam_acc, pfam_id, pfam_type, start, end], index=['Proteins', 'pfam_acc', 'pfam_id', 'pfam_type', 'pfam_start', 'pfam_end']).T)
        except:
            logger.info(f'{accession} not processed.')

    return pd.concat(pfam).reset_index(drop=True)


def gz_unzipper(filename, input_path='resources/bioinformatics_databases/', output_path='resources/bioinformatics_databases/'):
    with gzip.open(f'{input_path}{filename}.gz', 'rb') as f_in:
        with open(f'{output_path}{filename}', 'wb') as f_out:
            shutil.copyfileobj(f_in, f_out)


def tar_file_to_folder(input_path, output_path):
    tar = tarfile.open(f'{input_path}', 'r')
    tar.extractall(f'{output_path}')
    tar.close()



def go_lineage_tracer(go_term, obo_path, alt_ids=False, direct=False):
    """Return all nodes underneath (i.e. all children) of go_term of interest. 
    Default is to collect entire tree; if stopping at 'direct',
    then only direct decendants are collected."""

    # Read the ontology
    graph = obonet.read_obo(obo_path)
    # Collect all nodes into child:parents dataframe
    children = {}
    for node in graph.nodes:
        try:
            children[node] = graph.nodes[node]['is_a']
        except:
            children[node] = []

    child_terms = pd.DataFrame()
    child_terms['child_term'] = list(children.keys())
    child_terms['parents'] = [';'.join(terms) for terms in list(children.values())]

    # collect any alternate ids for go term
    search_list = []
    if alt_ids:
        try:
            search_list.append(graph.nodes[go_term]['alt_id'])
        except:
            pass
    search_list = [go_term] + [item for sublist in search_list for item in sublist]

    # Collect all children where term of interest is parent
    def search_terms(search_list, term_family=[], direct=False):
        term_family.append(search_list)
        search_list = '|'.join(search_list)
        terms = child_terms[child_terms['parents'].str.contains(search_list)]['child_term'].tolist()
        if direct:
            return terms
        if len(terms) > 0:
            search_terms(terms, term_family)
        return [item for sublist in term_family for item in sublist]

    # collect list of terms of interest
    family = search_terms(search_list, direct=False)

    return family


def uniprot_go_genes(tax_id, go_term, resource_folder='resources/bioinformatics_databases/', child_terms=True, direct=False, output='list'):
    """Collect all genes from annotated (reviewed) uniprot database containing the GO term of interest.
    tax_id: uniprot id corresponding to saved databases in resources folder e.g. '10090', '9606'
    go_term: term id e.g. 'GO:0032991'
    resource_folder: directory to where stored databases are
    child_terms: default(True) collects all terms for which the term if interest is a parent
    direct: default(false) limits child terms to direct descendents i.e. child term 'is_a' go_term
    output: default(list) choose type of output from 'list' ('Entry' ids), 'df' (complete genes df) or directory (save)"""

    # read in uniprot database for the species, with xref details
    uniprot_database = uniprot_summary(tax_id=tax_id, resource_folder=resource_folder, reviewed=True)
    genes = uniprot_database.dropna(subset=['GO']) # 16525/17474 = 95% have annotated GO terms

    # collect search terms according to optional child terms and direct lineage
    if child_terms:
        search_terms = go_lineage_tracer(go_term, obo_path=f'{resource_folder}PANTHERGOslim.obo', alt_ids=False, direct=direct)
        search_terms = '|'.join(search_terms)
    else:
        search_terms = go_term

    # Collect all genes with ontology_id
    gene_list = genes[genes['GO'].str.contains(search_terms)]

    # generate output
    if output == 'list':
        return gene_list['Entry'].tolist()
    elif output == 'df':
        return gene_list
    else:
        logger.info('Output format not detected. Attempting output to path.')
        gene_list.to_csv(output)


def ontology_wordfinder(words, obo_path='resources/bioinformatics_databases/PANTHERGOslim.obo', resource_folder='resources/bioinformatics_databases/'):
    """Retrieves all ontology terms containing 'words' in the name.
    words: list of words to search
    return: df of term, name matches"""

    # Read the ontology
    graph = obonet.read_obo(obo_path)

    # Collect term_ids, term names
    terms = [node for node in graph.nodes]
    names = [graph.nodes[node]['name'] for node in terms]
    terms = pd.DataFrame([terms, names], index=['go_term', 'go_name']).T

    # Collect only terms containing name of interest
    search_words = '|'.join(words)

    return terms[terms['go_name'].str.contains(search_words)]


def go_term_details(go_terms, obo_path='resources/bioinformatics_databases/PANTHERGOslim.obo', resource_folder='resources/bioinformatics_databases/'):
    """Retrieves details for all go terms as df.
    go_terms: list of term_ids to search
    return: df of term, name matches"""

    # Read the ontology
    graph = obonet.read_obo(obo_path)

    # Generate df for go_term
    cleaned_terms = []
    for go_term in go_terms:
        try:
            test_df = pd.DataFrame(dict([ (k,pd.Series(v)) for k,v in graph.nodes[go_term].items() ]))
            df = pd.DataFrame([';'.join(test_df[col].dropna()) for col in test_df.columns.tolist()], index=test_df.columns.tolist()).T
            df['go_id'] = go_term
            cleaned_terms.append(df)
        except:
            pass
    cleaned_terms = pd.concat(cleaned_terms)

    return cleaned_terms


def go_uniprot_proteins(protein_names, tax_id, resource_folder='resources/bioinformatics_databases/', name_type= 'Entry', output='list'):
    """For any given gene, pull out the associated GO terms annotated in uniprot as a list"""

    # read in uniprot database for the species, with xref details
    uniprot_database = uniprot_summary(tax_id=tax_id, resource_folder=resource_folder, reviewed=True)
    gene_details = uniprot_database[uniprot_database[name_type].isin(protein_names)]
    
   # generate output
    if output == 'list':
        return gene_details['GO'].tolist()
    elif output == 'df':
        return gene_details
    else:
        logger.info('Output format not detected')


def taxonomy_id(uniprot_tax_ids, resource_folder):
    species = pd.read_table(f'{resource_folder}orthodb_v10.1/odb10v1_species.tab.gz', compression='gzip', header=None)
    species.columns = ['ncbi_tax_id', 'ortho_tax_id', 'organism_name', 'genome_assembly_id', 'ortho_gene_count', 'ortho_group_count', 'mapping_type']
    species_dict = dict(zip(species['ncbi_tax_id'], species['ortho_tax_id']))

    return [species_dict[ncbi_id] for ncbi_id in uniprot_tax_ids]


@ProgressBar(step=1/41)
def create_genes(resource_folder, ortho_tax_ids):
    try:
        genes = pd.read_excel(f'{resource_folder}{"_".join(ortho_tax_ids)}_genes.xlsx')
        genes.drop([col for col in genes.columns.tolist() if 'Unnamed: ' in col], axis=1, inplace=True)
        yield
    except:
        logger.info('File not found. Processing database.')

        gene_chunks = pd.read_table(f'{resource_folder}orthodb_v10.1/odb10v1_genes.tab.gz', compression='gzip', chunksize=1000000, header=None)
        genes = []

        for df in gene_chunks:
            genes.append(df[df[1].isin(ortho_tax_ids)]) # note this relies on the order of the columns - see ortho README
            yield
        genes = pd.concat(genes)
        genes.columns = ['ortho_gene_id', 'ortho_organism_id', 'original_id', 'synonyms', 'mapped_uniprot_id', 'mapped_ensembl_ids', 'ncbi_gene_name', 'mapped_description']

        for tax_id in ortho_tax_ids:
            logger.info(f'{len(genes[genes["ortho_organism_id"] == tax_id])} genes found for {tax_id}')
        
        FileHandling.df_to_excel(f'{resource_folder}{"_".join(ortho_tax_ids)}_genes.xlsx', sheetnames=['all_genes'], data_frames=[genes])

    return genes


@ProgressBar(step=1/196)
def create_og2genes(resource_folder, ortho_tax_ids):
    try:
        og2genes = pd.read_excel(f'{resource_folder}{"_".join(ortho_tax_ids)}_go2genes.xlsx')
        og2genes.drop([col for col in og2genes.columns.tolist() if 'Unnamed: ' in col], axis=1, inplace=True)
        yield
    except:
        logger.info('File not found. Processing database.')

        og2genes_chunks = pd.read_table(f'{resource_folder}orthodb_v10.1/odb10v1_OG2genes.tab.gz', compression='gzip', header=None, chunksize=1000000)

        search_ids = '|'.join(ortho_tax_ids)
        og2genes = []
        for df in og2genes_chunks:
            og2genes.append(df[df[1].str.contains(search_ids)])  # Takes a while!
            yield
        og2genes = pd.concat(og2genes)
        # note this relies on the order of the columns - see ortho README
        og2genes.columns = ['og_id', 'ortho_gene_id']

        FileHandling.df_to_excel(f'{resource_folder}{"_".join(ortho_tax_ids)}_go2genes.xlsx', sheetnames=[search_ids], data_frames=[og2genes])

    return og2genes


@ProgressBar(step=1/16)
def create_uniprot_db(input_path, gene_ids=[]):
    # Testing collect uniprot entries
    from_uniprot_db_chunks = pd.read_table(f'{input_path}', compression='gzip', chunksize=10000)

    from_uniprot_db = []
    for df in from_uniprot_db_chunks:
        if len(gene_ids):
            search_ids = '|'.join(gene_ids)
            from_uniprot_db.append(df[df['Entry'].str.contains(search_ids)])  # Takes a while!
        else:
            from_uniprot_db.append(df)  # Takes a while!
        yield
            
    from_uniprot_db = pd.concat(from_uniprot_db)

    return from_uniprot_db


@ProgressBar(step=0.05)
def create_uniprot_map(tax_id, resource_folder, gene_ids=[]):
    # collect ensembl ids from uniprot mapped file
    uniprot_chunks = pd.read_table(f'{resource_folder}{tax_id}_idmapping.tab.gz', compression='gzip', chunksize=10000, header=None)

    uniprot_map = []
    for df in uniprot_chunks:
        if len(gene_ids) > 0:
            search_ids = '|'.join(gene_ids)
            uniprot_map.append(df[df[0].str.contains(search_ids)])  # Takes a while!
        else:
            uniprot_map.append(df)  # Takes a while!
        yield

    uniprot_map = pd.concat(uniprot_map)
    uniprot_map.columns = ['UniProtKB-AC', 'UniProtKB-ID', 'GeneID(EntrezGene)', 'RefSeq', 'GI', 'PDB', 'GO', 'UniRef100', 'UniRef90', 'UniRef50', 'UniParc', 'PIR', 'NCBI-taxon', 'MIM', 'UniGene', 'PubMed', 'EMBL', 'EMBL-CDS', 'Ensembl', 'Ensembl_TRS', 'Ensembl_PRO', 'Additional PubMed']

    return uniprot_map


@ProgressBar(step=0.05)
def create_uniprot_xref(input_path, tax_id, gene_ids=[], id_type=None):
    # collect uniprot mapped file for specific genes
    uniprot_chunks = pd.read_table(f'{input_path}{tax_id}_idmapping.dat.gz', compression='gzip', chunksize=10000, header=None)

    uniprot_xref = []
    for df in uniprot_chunks:
        if len(gene_ids) > 0:
            search_ids = '|'.join(gene_ids)
            uniprot_xref.append(df[df[0].str.contains(search_ids)]) # Takes a while!
        else:
            uniprot_xref.append(df)  # Takes a while!
        yield

    uniprot_xref = pd.concat(uniprot_xref)
    uniprot_xref.columns = ['UniProtKB-AC', 'ID_type', 'ID']

    if id_type:
        uniprot_xref = uniprot_xref[uniprot_xref['ID_type'] == id_type]

    return uniprot_xref


@ProgressBar(step=0.05)
def convert_geneids(gene_ids, tax_id, id_from, id_to, resource_folder='resources/bioinformatics_databases/'):
    # collect uniprot mapped file for specific genes
    uniprot_chunks = pd.read_table(f'{resource_folder}{tax_id}_idmapping.dat.gz', compression='gzip', chunksize=10000, header=None)

    uniprot_xref = []
    for df in uniprot_chunks:
        df.columns = ['UniProtKB-AC', 'ID_type', 'ID']
        search_ids = '|'.join(gene_ids)
        if id_from == 'uniprot':
            df = df[df['UniProtKB-AC'].str.contains(search_ids)] # Takes a while!
        else:
            df = df[df['ID'].astype(str).str.contains(search_ids)]
        uniprot_xref.append(df)  # Takes a while!
        yield

    uniprot_xref = pd.concat(uniprot_xref)
    if not id_to == 'uniprot':
        uniprot_xref = uniprot_xref[uniprot_xref['ID_type'] == id_to]
    if not id_from == 'uniprot':
        uniprot_xref = uniprot_xref[uniprot_xref['ID_type'] == id_from]

    return uniprot_xref


def uniprot_genename_mapper(uniprot_tax_ids, gene_ids, reviewed=False, orthodb=False, cols=[]):

    # Collect Uniprot info
    compiled = {}
    for direction, tax in uniprot_tax_ids.items():
        if direction == 'from':
            genes = gene_ids
        else:
            genes = []
        uniprot_db = create_uniprot_db(f'{resource_folder}{tax}.tab.gz', genes)
        uniprot_map = create_uniprot_map(f'{resource_folder}{tax}_idmapping.tab.gz', genes)
        merged = pd.merge(uniprot_db, uniprot_map, left_on='Entry', right_on='UniProtKB-AC', how='outer')
        merged['name'] = merged['Entry name'].str.split('_').str[0]
        compiled[direction] = merged

    # Join based on gene name
    to_from_compiled = pd.merge(compiled['from'], compiled['to'], on='name', how='inner', suffixes=('_from', '_to'))

    # annotate matching orthoDB references
    to_from_compiled['orthologous'] = [1 if val_1 == val_2 else 0 for val_1, val_2 in to_from_compiled[[
        'Cross-reference (OrthoDB)_from', 'Cross-reference (OrthoDB)_to']].values]

    # Apply filters
    if reviewed:
            merged = merged[merged['Status'] == 'reviewed']
    if orthodb:
        to_from_compiled = to_from_compiled[to_from_compiled['orthologous'] == 1]
    if cols:
        cols = [[f'{col}_from', f'{col}_to'] for col in cols]
        cols = ['name', 'orthologous'] + [item for sublist in cols for item in sublist]

    # collect unmapped genes
    unmapped_genes = set(gene_ids) - set(to_from_compiled['Entry_from'].tolist())

    return to_from_compiled, unmapped_genes


def chunks(lst, n):
    """Yield successive n-sized chunks from lst."""
    for i in range(0, len(lst), n):
        yield lst[i:i + n]


def ortholog_map(gene_ids, direction, output_folder, resource_folder='resources/bioinformatics_databases/'):

    # Generate mapped ids using orthologsBioMART
    # Collect appropriate id's - note this is not necessarily going to generate unique id's
    gene_identifiers_from = create_uniprot_xref(input_path=f'{resource_folder}', tax_id='9606', gene_ids=gene_ids, id_type='Ensembl')
    gene_identifiers_from = dict(zip(gene_identifiers_from['UniProtKB-AC'], gene_identifiers_from['ID']))
    # Collect orthologous genes
    if direction == 'h2m':
        mapped_ids = biomart_h2m(list(gene_identifiers_from.values()), identifier_type='link_ensembl_gene_id')
        uniprot_tax_ids = {'from': '9606', 'to': '10090'}
    elif direction == 'm2h':
        mapped_ids = biomart_m2h(list(gene_identifiers_from.values()), identifier_type='link_ensembl_gene_id')
        uniprot_tax_ids = {'from': '10090', 'to': '9606'}
    # map ensembl back to uniprot
    gene_identifiers_to = create_uniprot_xref(input_path=f'{resource_folder}', tax_id='10090', gene_ids=[], id_type='Ensembl')
    gene_identifiers_to = dict(zip(gene_identifiers_to['ID'], gene_identifiers_to['UniProtKB-AC']))
    # Add Uniprot info back to the mapped ids
    mapped_ids['Entry_to'] = mapped_ids['mouse_ensembl_gene_id'].map(gene_identifiers_to)

    # generate mapped ids using UniProt mapper
    cols = ['Entry', 'Entry name', 'Status', 'GeneID(EntrezGene)', 'Ensembl', 'Cross-reference (OrthoDB)']
    compiled, unmapped_genes = uniprot_genename_mapper(
        uniprot_tax_ids, gene_ids, reviewed=True, orthodb=False)


    FileHandling.df_to_excel(
        output_path=f'{output_folder}mapped_ids.xlsx',
        sheetnames=['BioMART_map', 'UniProt_map', 'UniProt_unmapped'],
        data_frames=[mapped_ids, compiled, pd.DataFrame(unmapped_genes)]
    )

    return compiled, mapped_ids


def uniprot_summary(tax_id, resource_folder, genes=[], reviewed=False):
    uniprot_db = create_uniprot_db(f'{resource_folder}{tax_id}.tab.gz', genes)
    uniprot_map = create_uniprot_map(resource_folder=f'{resource_folder}', tax_id=tax_id, gene_ids=genes)
    merged = pd.merge(uniprot_db, uniprot_map, left_on='Entry', right_on='UniProtKB-AC', how='outer')
    merged['name'] = merged['Entry name'].str.split('_').str[0]
    if reviewed:
        merged = merged[merged['Status'] == 'reviewed']
    return merged


def download_databases(tax_ids={'MOUSE': '10090', 'HUMAN': '9606'}, resource_folder='resources/bioinformatics_databases/', caller_id = "www.github.com/dezeraecox"):


    if not os.path.exists(resource_folder):
        os.makedirs(resource_folder)

    for species, tax_id in tax_ids.items():
        if not os.path.exists(f'{resource_folder}{tax_id}_idmapping.tab.gz'):
            download_resources(
                filename=f'{tax_id}_idmapping.tab.gz',
                url=f'https://ftp.expasy.org/databases/uniprot/current_release/knowledgebase/idmapping/by_organism/{species}_{tax_id}_idmapping_selected.tab.gz',
                resource_folder=resource_folder)
        gz_unzipper(f'{tax_id}_idmapping.tab', input_path=resource_folder, output_path=resource_folder)


        if not os.path.exists(f'{resource_folder}{tax_id}.tab.gz'):
            download_resources(
                filename=f'{tax_id}.tab.gz',
                url=f'https://www.uniprot.org/uniprot/?query={tax_id}&%28{species}%29+%5B%22&fil=&offset=0&compress=yes&format=tab',
                resource_folder=resource_folder)
        gz_unzipper(f'{tax_id}.tab', input_path=resource_folder, output_path=resource_folder)
    
        if not os.path.exists(f'{resource_folder}{tax_id}_idmapping.dat.gz'):
            download_resources(
                filename=f'{tax_id}_idmapping.dat.gz',
                url=f'https://ftp.expasy.org/databases/uniprot/current_release/knowledgebase/idmapping/by_organism/{species}_{tax_id}_idmapping.dat.gz',
                resource_folder=resource_folder)
                

    if not os.path.exists(f'{resource_folder}PANTHERGOslim.obo'):
        download_resources(
            filename=f'PANTHERGOslim.obo',
            url=f'http://data.pantherdb.org/PANTHER15.0/ontology/PANTHERGOslim.obo',
            resource_folder=resource_folder)




if __name__ == "__main__":
    resource_folder = 'resources/bioinformatics_databases/'
    download_databases(tax_ids={'MOUSE': '10090'}, resource_folder=resource_folder)





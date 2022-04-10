import re
import os
import pandas as pd

from utilities.databases import uniprot_go_genes, ontology_wordfinder,  go_lineage_tracer, go_term_details

from loguru import logger
from GEN_Utils import FileHandling

logger.info('Import OK')

input_path = 'results/tpemi_proteomics/significance/cys_scaled_summary.xlsx'
resource_folder = 'resources/bioinformatics_databases/'
output_folder = 'results/tpemi_proteomics/protein_complexes/'

if not os.path.exists(output_folder):
    os.mkdir(output_folder)

ontologies = {}


# ----Collect specific complexes based on gene ontology terms----

# Test 3: Against proteins associated with proteasome
potential_proteasome_terms = ontology_wordfinder(['proteasome']) # decided on "GO:0000502: proteasome complex"
ontology_genes = uniprot_go_genes(tax_id='10090', go_term='GO:0000502', child_terms=False, direct=True, output='list')
ontologies['proteasome_GO0000502'] = pd.DataFrame(ontology_genes)

# Test 4: Against proteins associated with ribosome
potential_terms = ontology_wordfinder(['ribosome']) # decided on "GO:0003735 structural constituent of ribosome""
ontology_genes = uniprot_go_genes(tax_id='10090', go_term='GO:0003735', child_terms=False, direct=True, output='list')
ontologies['ribosome_GO0003735'] = pd.DataFrame(ontology_genes)

# Test 5: Against proteins associated with DNA repair complex
potential_terms = ontology_wordfinder(['DNA repair complex']) # decided on "GO:1990391 DNA repair complex""
ontology_genes = uniprot_go_genes(tax_id='10090', go_term='GO:1990391', child_terms=False, direct=True, output='list')
ontologies['DNArepair_GO1990391'] = pd.DataFrame(ontology_genes)

# Test 6: Against proteins associated with nuclear pore
potential_terms = ontology_wordfinder(['nuclear pore']) # decided on "GO:0005643 Nuclear pore"
ontology_genes = uniprot_go_genes(tax_id='10090', go_term='GO:0005643', child_terms=True, direct=False, output='list')
ontologies['nuclearpore_GO0005643'] = pd.DataFrame(ontology_genes)

# collect specific 'end-of-tree' terms contained under "protein complex" GO term, to then grab mouse genes annotated with each term
# getting all direct descendants of the generic "protein-containing complex" term
go_terms = go_lineage_tracer('GO:0032991', obo_path=f'{resource_folder}PANTHERGOslim.obo', alt_ids=True, direct=False)
details = go_term_details(go_terms, obo_path=f'{resource_folder}PANTHERGOslim.obo')
# this is a good start: 210 terms is comparable to the 279 annotated complexes in the Ori database

# collect genes annotated with these GO terms
go_complexes = {}
for go_term in go_terms:
    ontology_genes = uniprot_go_genes(tax_id='10090', go_term=go_term, child_terms=False, direct=True, output='df')
    go_complexes[go_term] = ontology_genes

# Add genes to details table
go_genes = {}
for go_term in go_terms:
    go_genes[go_term] = go_complexes[go_term]['Entry'].tolist()
details['genes'] = details['go_id'].map(go_genes)
ontologies['GO_summary'] = details
ontologies['GO_complexes'] = pd.DataFrame(go_genes.values(), index=go_genes.keys()).T

# Save all GO summaries to excel
FileHandling.df_to_excel(
    output_path=f'{output_folder}ontology_summary.xlsx',
    sheetnames=list(ontologies.keys()),
    data_frames=list(ontologies.values())
)
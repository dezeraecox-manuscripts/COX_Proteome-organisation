[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.6439170.svg)](https://doi.org/10.5281/zenodo.6439170) 


# COX_Proteome-organisation

This repository contains the analysis code associated with the Proteome Organisation project, led by Dr. Dezerae Cox. This manuscript has been submitted for publication under the title "*Protein painting reveals pervasive remodeling of conserved proteostasis machinery in response to pharmacological stimuli*".

This manuscript has been submitted as a preprint via BioRxiv [here](biorxiv/link). A link to the final version will be provided upon publication.

## Prerequisites

This analysis assumes a standard installation of Python 3 (=> 3.8). For specific package requirements, see the environment.yml file, or  create a new conda environment containing all packages by running ```conda create -f environment.yml```. To complete the network mapping elements of the workflow, an installation of [Cytoscape (v3.9.0)](https://cytoscape.org/download.html) is also required. In addition to the analysis contained here, some simple statistical tests were performed using [GraphPad Prism v **8.0**](https://www.graphpad.com/scientific-software/prism/).

## Raw data

The .RAW files have been deposited via the [PRIDE][1]<sup>[1]</sup> partner repository to the [ProteomeXchange Consortium][2]<sup>[2]</sup> under the dataset identifier PXD033152. For convenience, the preprocessed identification and quantitation data have also been uploaded alongside the exported flow cytometry data (hereon termed raw data) as an open-access [Zenodo dataset](https://doi.org/10.5281/zenodo.6439170). These data can be collected automatically using the ```raw_data.py``` script in each of the respective ```src``` analysis folders.

Various public databases were also queried as cited in the accompanying manuscript, for which access protocols are provided in the respective ```utilities``` scripts where appropriate.

## Workflow

Initial processing of the novel mass spectrometry spectra files was completed using [MaxQuant][3]<sup>[3]</sup> (v1.6.2.10).

Individual analyses are presented within the ```src``` folder. Where processing order is important for individual analyses, scripts have been numbered and should be run in order before unnumbered counterparts. To run the complete analysis, a ```snakefile``` has also been provided. This can be run at the command line from the root directory using ```snakemake -p --cores 1```. Note, however, that Cytoscape must be running to enable access from the command line and that the generated networks should be cleared between the two different cytoscape plotting scripts. Alternatively, these elements of the workflow can be addressed manually.

### A quick note on cluster numbers

The assigned number during the community clustering is stochastic (as it is unsupervised, the cluster centres can be initialised in a randomized order). However, the proteins associated to each cluster should be conserved. Therefore, cluster numbers were manually mapped to the order presented in the manuscript within the analysis workflow.

### A quick note on database and enrichment analyses

Where necessary, some analyses presented in this workflow make live calls to an existing API at the time of running. In particular, the enrichment analyses use [PantherDB](http://www.pantherdb.org/) [UniProt](https://www.uniprot.org/) online. Therefore, there may be subtle differences in the returned GO term lists resulting from any version updates. The analyses obtained via API presented in the paper used PANTHER v16.0 and UniProt v2021_03.


## References

[1]: https://www.ebi.ac.uk/pride/archive/

1. Perez-Riverol, Yasset, Attila Csordas, Jingwen Bai, Manuel Bernal-Llinares, Suresh Hewapathirana, Deepti J Kundu, Avinash Inuganti, et al. “The PRIDE Database and Related Tools and Resources in 2019: Improving Support for Quantification Data.” Nucleic Acids Research 47, no. D1 (January 8, 2019): D442–50. https://doi.org/10.1093/nar/gky1106.

[2]: http://proteomecentral.proteomexchange.org

2. Deutsch, Eric W., Attila Csordas, Zhi Sun, Andrew Jarnuczak, Yasset Perez-Riverol, Tobias Ternent, David S. Campbell, et al. “The ProteomeXchange Consortium in 2017: Supporting the Cultural Change in Proteomics Public Data Deposition.” Nucleic Acids Research 45, no. Database issue (January 4, 2017): D1100–1106. https://doi.org/10.1093/nar/gkw936.

[3]: https://www.maxquant.org/

1. Tyanova, Stefka, Tikira Temu, and Juergen Cox. “The MaxQuant Computational Platform for Mass Spectrometry-Based Shotgun Proteomics.” Nature Protocols 11, no. 12 (December 2016): 2301–19. https://doi.org/10.1038/nprot.2016.136.

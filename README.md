# HIV1-CRF19_Cuba
Phylogenetic analyses of CRF19 spread in Cuba and worldwide

## Analysis pipelines

The [*snakemake*](snakemake) folder contains Snakemake [[Köster *et al.*, 2012](https://doi.org/10.1093/bioinformatics/bts480)] pipelines
for reconstruction of evolutionary history of HIV-1 in Cuba:
* [*Snakefile_datasets*](snakemake/Snakefile_datasets) for D, A1 and G data set creation
![pipeline visualisation](snakemake/pipeline_datasets.svg)
* [*Snakefile_trees*](snakemake/Snakefile_trees) for phylogenetic tree reconstruction, rooting and dating
![pipeline visualisation](snakemake/pipeline_trees.svg)
* [*Snakefile_trees*](snakemake/Snakefile_acr) for phylogeographic, drug resistance and transmission mode analyses
![pipeline visualisation](snakemake/pipeline_acr.svg)

## Data
* [*data/input*](data/input) folder contains the input curated 2018 HIV-1 alignment from [Los Alamos database](https://www.hiv.lanl.gov/content/sequence/NEWALIGN/align.html),
and the input files for [jpHMM](http://jphmm.gobics.de/).
* [*data/datasets/metadata.tab*](data/datasets/metadata.tab) contains the combined metadata for CU and LA sequences used for these analyses (produced with [*Snakefile_datasets*](snakemake/Snakefile_datasets) pipeline). 
* [*data/datasets/D_CRF_19/aln.first.cleaned.fa*](data/datasets/D_CRF_19/aln.first.cleaned.fa), 
[*data/datasets/A1_CRF_19/aln.first.cleaned.fa*](data/datasets/A1_CRF_19/aln.first.cleaned.fa), 
[*data/datasets/G_CRF_19/aln.first.cleaned.fa*](data/datasets/G_CRF_19/aln.first.cleaned.fa) contain the combined CU+LA MSAs (including 5 outgroup sequences indicated in outgroup.txt files in the corresponding folders) for the D/A1/G data sets used for these analyses (produced with [*Snakefile_datasets*](snakemake/Snakefile_datasets) pipeline).
* [*data/datasets/D_CRF_19/metadata.drms.tab*](data/datasets/D_CRF_19/metadata.drms.tab), [*data/datasets/D_CRF_19/metadata.drugs.tab*](data/datasets/D_CRF_19/metadata.drugs.tab) contain the Surveillance DRM and ARV metadata for the D+CRF_19 sequences 
extracted with [Sierra](https://hivdb.stanford.edu/page/webservice/) (see [*Snakefile_datasets*](snakemake/Snakefile_datasets) pipeline).

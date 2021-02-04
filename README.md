# HIV1-CRF19_Cuba
Phylogenetic analyses of CRF19 spread in Cuba and worldwide


The [*snakemake*](snakemake) folder contains Snakemake [[Köster *et al.*, 2012](https://doi.org/10.1093/bioinformatics/bts480)] pipelines
for reconstruction of evolutionary history of HIV-1 in Cuba:
* [*Snakefile_datasets*](snakemake/Snakefile_datasets) for D, A1 and G data set creation
![pipeline visualisation](snakemake/pipeline_datasets.svg)
* [*Snakefile_trees*](snakemake/Snakefile_trees) for phylogenetic tree reconstruction, rooting and dating
![pipeline visualisation](snakemake/pipeline_trees.svg)
* [*Snakefile_trees*](snakemake/Snakefile_acr) for phylogeographic, drug resistance and transmission mode analyses
![pipeline visualisation](snakemake/pipeline_acr.svg)

The [*data/input*](data/input) folder contains the input curated 2018 HIV-1 alignment from [Los Alamos database](https://www.hiv.lanl.gov/content/sequence/NEWALIGN/align.html),
and the input files for [jpHMM](http://jphmm.gobics.de/).
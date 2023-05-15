## - toolbox/blast2region
### This tool runs a BlastP on a FASTA-formatted file containing at least one protein sequence using either explicitly delimited databases or a custom one provided by the user. It returns the GenBank-formatted genomic region of each hit found in the BlastP search.
####    (1) DATABASES
#####       1.1   Databases can be either explicitly delimited in the `config/blast_config.yaml` file (*current version: ncbi-nr and swissprot*) or provided as a file path by the user.
#####       1.2   Sequence accessions are gathered from each FASTA identifier found in sequence databases in a specific manner that takes into account how these identifiers are generated. This behavior is governed by the `id_sep_dict` variable in the `config/blast_config.yaml`. The accession gathering behavior is pre-optimized in standard databases, and it can be fine-tuned for custom databases.
####    (2) OUTPUTS
#####       2.1 The output directory has two components:
######          -   The `summary_report.csv` which shows some of the most important BlastP metrics for the hits found for each query sequence
######          -   An array of directories named after each query sequence provided in the input. Within each directory, the genomic regions of each of the hits found in the BlastP search will be saved in GenBank format.
######          -   The regions that will be gathered around each BlastP hit encoding gene depend on the `window_size` parameter in `config/blast_config.yaml`. The default window is set to 15 kilobases.
####    (3) ENVIRONMENT
#####       3.1 The conda environment in `env/b2r_env.yaml` contains all `blast2region.py` dependencies. It can be installed/activated as follows:
######          `conda env create -f env/b2r_env.yaml`
######          `conda activate blast2region`


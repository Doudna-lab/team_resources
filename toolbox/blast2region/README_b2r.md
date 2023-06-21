# - `toolbox/blast2region.py`
### This tool runs a BlastP on a FASTA-formatted file containing at least one protein sequence using either explicitly delimited databases or a custom database provided by the user. It returns the GenBank-formatted genomic region associated with each hit found in the BlastP search.
###    (0) QUICK-START
####        0.1 - The argument `--help` (or `-h`) can be used for usage instructions:
```
python3 blast2region.py --help
```
####        0.2 - Provide a FASTA-formatted file as a positional argument for a quick standard run.
```
python3 blast2region.py <myfasta.faa>
```
####        0.3 - BlasP parameters such as number of threads and e-value cutoff can be controlled by the tool respectively with the `-t` and `-e` options:
```
python3 blast2region.py <my_fasta.faa> -t 6 -e 0.05
```
###    (1) DATABASES
####       1.1 -  Databases can be either explicitly delimited in the `config/blast_config.yaml` file (*current version: ncbi-nr and swissprot*) or provided as a file path by the user.
####       1.2 -  Sequence accessions are gathered from each FASTA identifier found in sequence databases in a specific manner that takes into account how these identifiers are generated. This behavior is governed by the `id_sep_dict` variable in the `config/blast_config.yaml`. The accession gathering behavior is pre-optimized in standard databases, and it can be fine-tuned for custom databases. 
###    (2) OUTPUTS
####       2.1 - The output directory has two components:
#####          -   The `summary_report.csv` which shows a set of BlastP metrics for the hits found for each query sequence.
#####          -   An array of directories named after each query sequence provided in the input. Within each directory, the genomic regions of each hit found in the BlastP search will be saved in GenBank format (*.gb*).
#####          -   The regions that will be gathered around each BlastP hit encoding gene depend on the `window_size` parameter in `config/blast_config.yaml`. The default window is set to 15 kilobases.
####       2.2 - To keep the output size under control, the tool limits the maximum number of hits reported based on the `blast_nkeep` parameter in `config/blast_config.yaml`. This will perform a simple selection the top N hits based exclusively on e-value. 
###    (3) ENVIRONMENT
####       3.1 - The conda environment in `env/b2r_env.yaml` contains all `blast2region.py` dependencies. It can be installed/activated as follows (tested on Biotite on *05/15/2023*):
##### - Create environment:
```
conda env create -f env/b2r_env.yaml
```
##### - Activate environment:
```
conda activate blast2region
```

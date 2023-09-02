# Phylogenetic Reconstruction Workflow

## Summary
+ This document describes the steps towards reconstructing the phylogeny of a specific protein motif and assessing its prevalence across the domains of life.
+ The starting point of this workflow is a multiple sequence alignment under the title 'G2II_nLTR_RT'
### Materials
+ Aiming to gain broader perspective about the evolutionary history of this motif, we first set out to glean information on similar sequences by conducting a broad iterative search using NCBI's NR database.
  + Sequence Database Downloads:
    + NR: 08/18/2023 (https://ftp.ncbi.nlm.nih.gov/blast/db/FASTA/nr.gz)
### Methods
#### Iterative Sequence Search 
+ An iterative search was conducted by using PSI-BLAST (https://doi.org/10.1093/nar/25.17.3389) with the following parameters:
  + `-num_iterations 10 -evalue 1e-2 -inclusion_ethresh 1e-05 -qcov_hsp_perc 85 `
+ NCBI's non-redundant (nr) database from `Jun 13, 2021` was utilized in the sequence search.
#### Taxonomic Landscape of Similar Sequences 
+ The taxonomic classification stratification of the subject sequences was assessed with Krona (https://doi.org/10.1186/1471-2105-12-385)
+ Taxonomic assignments of subject sequences were gleaned with EFetch (https://www.ncbi.nlm.nih.gov/books/NBK25499/) and added to a FASTA report including all subjects from the sequence search for downstream analyses.
#### Realignement of the expanded set of sequences
+ All sequences found in the sequence search were consolidated with the query MSA to make up for a re-alignment
+ This step was conducted by using Clustal Omega with 10 iterations (option: `--iter 10`)
#### Alignment Gap Removal
+ To conduct an automated alignment gap removal step, the online tool Jpred 4 (https://doi.org/10.1093/nar/gkv332) was used.
  + The realigned MSA was fed to the tool in FASTA format and the processed MSA was downloaded to the `{run}/jpred` directory
#### Phylogenetic Reconstruction
+ The resulting alignment was then processed by two different phylogenetic reconstruction tools:
  + For deep phylogeny processing: IQ-Tree v2.2.3 (https://doi.org/10.1093/molbev/msu300) with the following parameters:
        + ` -mtree -m MFP -bb 1000` where:
          + `bb` implements the ultrafast bootstrap approximation algorithm with 1000 replicates
          + `-m` implements the substitution model selection algorith to determine the most appropriate for the current input based on bayesian- and akaike information criteria
  + A quicker and possibly less accurate inference was also implemented as a time-efficient solution through FastTree v2.1.11 (https://doi.org/10.1093/molbev/msp077) with the following parameters:
    + `-boot 1000` that specifies 1000 bootstrap replicates
#### HMM Profile
+ Hidden Markov Model profile construction was carried out with `hmmbuild` included with Hmmer3 v3.3 (https://doi.org/10.1371/journal.pcbi.1002195)
+ Default parameters were applied in this process.
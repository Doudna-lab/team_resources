# Phylogenetic Reconstruction Workflow

## Summary
+ This document describes the steps towards reconstructing the phylogeny of a specific protein motif and assessing its prevalence across the archaeal and prokaryotic domains of life.
+ The starting point of this workflow is a multiple sequence alignment under the title 'G2II_nLTR_RT'
### Materials
+ Aiming to gain perspective about the evolutionary history of this motif, we first set out to glean information on similar sequences by conducting a broad iterative search using NCBI's NR, Swissprot and PDB databases.
  + Sequence Database Downloads:
    + NR: 08/18/2023 (https://ftp.ncbi.nlm.nih.gov/blast/db/FASTA/nr.gz)
    + Swissprot: 08/18/2023 (https://ftp.ncbi.nlm.nih.gov/blast/db/FASTA/swissprot.gz)
    + PDB: 08/18/2023 (https://ftp.ncbi.nlm.nih.gov/blast/db/FASTA/pdbaa.gz)
+ These databases will be merged into one. To the merged result an additional layer of NCBI's TaxID information will be added which will allow for a taxonomically-aware assessment of the results.
    + Taxonomic data mapping download:
      + NCBI's TaxID: 08/21/2023 (https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/prot.accession2taxid.FULL.gz)
### Methods
+ An iterative search will be conducted by using PSI-BLAST (https://doi.org/10.1093/nar/25.17.3389) with the following parameters:
  + `-num_iterations 10 -evalue 1e-10 -qcov_hsp_perc 85 -perc_identity 80`
+ The database creation with `makeblastdb` included a TaxID mapping file by setting the following option:
  +  `-taxid_map <path>/prot.accession2taxid.FULL`

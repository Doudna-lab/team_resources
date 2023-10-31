# PhyDMS Workflow Utilizing DMS data

## Summary
+ This document describes the application of `phydms` to carry out phylogenetic analyses using Experimentally Informed Codon Models (ExpCM) 
### Materials
+ The starting point of this workflow was a list of RefSeq Protein IDs contained in the `/groups/doudna/projects/daniel_projects/prywes_n/input_data/bacterialFormIIs.csv` file
+ These sequences were the basis to build the nucleotide sequence alignment and subsequent phylogenetic reconstruction to support the `phydms` analysis
+ The reference sequence used in this workflow was [WP_011390153.1](https://www.ncbi.nlm.nih.gov/protein/WP_011390153.1?report=genbank&log$=protalign&blast_rank=2&RID=KJP5J24C016)

### Methods
+ The RefSeq IDs were used in a programmatic Entrez search to gather both amino acid and respective encoding nucleotide sequences of each entry.
+ A total of 522 were used to generate the nucleotide alignment using MAFFT (https://doi.org/10.1093/nar/gki198) via `phydms_prepalignment` from the `phydms` package 
+ The `phydms`(https://doi.org/10.7717/peerj.3657) was applied to the set of sequences
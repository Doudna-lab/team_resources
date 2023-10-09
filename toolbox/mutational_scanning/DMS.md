# DMS Pre-Analysis

## Summary
+ This document describes the steps that prepares a library of single variant sequences for a Deep Mutational Scanning (DMS) workflow.
### Materials
+ The library generation was carried out based on the __species_name__'s *tadA8e* nucleotide sequence acquired from __sequence_source__     
### Methods
+ The translated sequence of the gene was set as reference, as each canonical amino acid was inserted through a computational iteration to create a variant using biopython (http://10.0.4.69/bioinformatics/btp163). Then each variant was back-translated to nucleotide based on a *Homo sapiens*-optimized codon table ("data/codon table .docx").
+ A total of 3154 sequences were generated in this library.

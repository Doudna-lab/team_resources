# Native modules
import os
from pathlib import Path
import re
# Installed Modules
from Bio.Seq import Seq
from Bio import SeqIO
import Bio.Data.CodonTable as ctable
from Bio.SeqUtils import CodonUsage
from Bio.SeqUtils.CodonUsage import CodonAdaptationIndex as CAI



# DEBUG
abs_path = "/Users/bellieny/projects/team_resources/toolbox/mutational_scanning/fasta/TadA8e.fna"

def main():
	dnaseq_record = SeqIO.read(abs_path, "fasta")
	dnaseq_in = Seq(dnaseq_record.seq)
	ptnseq_in = dnaseq_in.translate()

	codon_table = ctable.standard_dna_table
	query_cai = CAI()
	query_cai.cai_for_gene(dnaseq_in)


if __name__ == "__main__":
	main()

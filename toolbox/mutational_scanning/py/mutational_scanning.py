# Native modules
import re
from argparse import ArgumentParser as argp
# Installed Modules
import yaml
from Bio.Seq import Seq
from Bio import SeqIO
# Project modules
from re_check import avoidREs_nodelchar

# DEBUG
# fasta_in_path = "/Users/bellieny/projects/team_resources/toolbox/mutational_scanning/fasta/TadA8e.fna"
# config_path = "/Users/bellieny/projects/team_resources/toolbox/mutational_scanning/config/tadA.yaml"
# res_enzyme_path = "/Users/bellieny/projects/team_resources/toolbox/mutational_scanning/data/pacI_restriction.csv"


def parse_arguments():
	#  Launch argparse parser
	parser = argp(
		prog='mutational_scanning.py',
		description='',
		usage='%(prog)s [options] <fasta_file>')
	# Define arguments
	parser.add_argument('fasta_file',
	                    help='Path to fasta input file')  # positional argument
	parser.add_argument('-c',
	                    dest='config',
	                    default='/Users/bellieny/projects/team_resources/toolbox/mutational_scanning/config/tadA.yaml',
	                    help='Specify which config file will be loaded to the script. [Default: tadA.yaml]')
	parser.add_argument('-s',
						dest='split_parts',
						default='',
						help='Specify whether the resulting sequences should be split in half or not')
	parser.add_argument('--re',
						dest='res_enzyme',
						help='Specify a path to restriction enzyme data where the second column '
							 'contains sequence strings of the recognition sites')

	# Parse arguments from the command line
	arguments = parser.parse_args()
	return arguments


def format_codon_tbl(codon_tbl: dict):
	formated_codon = {}
	for key in codon_tbl:
		formated_codon.setdefault(key.upper(), codon_tbl[key])
	return formated_codon


def generate_variants(ptnrec: SeqIO.SeqRecord, backtable: dict):
	"""
	:param ptnrec: A single sequence enclosed in a SeqIO.SeqRecord object
	holding at least the sequence and an identifier
	:param backtable: Codon back table: each key represented by a one-letter
	abbreviation of an Aminoacid; each value represented by its respective codon
	:return: Two lists of SeqIO.SeqRecord. First for AAs, second for NTs.
	"""
	# Force all AA to uppercase to keep consistency
	format_backtable = format_codon_tbl(backtable)
	# Set internal variables
	ptnseq = ptnrec.seq
	ptn_id = ptnrec.id
	ptnseq_ls = list(ptnseq)
	# Set return variables
	ptn_record_ls = []
	gene_record_ls = []
	wt_backtranslated_seq = ''

	# Set up a SeqRecord entry for the WT nucleotide sequence backtranslated
	wt_ptnseq = Seq(''.join(ptnseq))
	for wt_residue in wt_ptnseq:
		wt_backtranslated_seq += format_backtable[wt_residue]
	SeqIO.SeqRecord(Seq(wt_backtranslated_seq), id=ptn_id, description='')

	# Loop through the AA position in the reference sequence
	for aacid_idx in range(len(ptnseq_ls)):
		# Each AA position will be replaced by all possible amino acids in the backtable
		for aa_in_lib in format_backtable:
			if aa_in_lib != ptnseq_ls[aacid_idx].upper():
				loop_ptnseq_ls = ptnseq_ls.copy()
				aacid_idx_label = aacid_idx + 1
				variant_fasta_label = f"{ptn_id}|{aacid_idx_label}_{aa_in_lib}"
				#
				loop_ptnseq_ls[aacid_idx] = aa_in_lib
				loop_ptnseq = Seq(''.join(loop_ptnseq_ls))
				ptn_record_ls.append(SeqIO.SeqRecord(loop_ptnseq, id=variant_fasta_label, description=''))

				variant_backtrans_seq = ""
				for variand_aa in loop_ptnseq:
					variant_backtrans_seq += format_backtable[variand_aa]
				gene_record_ls.append(SeqIO.SeqRecord(Seq(variant_backtrans_seq), id=variant_fasta_label, description=''))
	return ptn_record_ls, gene_record_ls, wt_backtranslated_seq


def split_variants(gene_records, original_sequence, split_parts: int):
	count = 0
	filtered_split_records = []
	for record in gene_records:
		left_split_anchor = 0
		for sequence_window in range(1, split_parts + 1):
			right_split_anchor = round(len(record.seq) / split_parts * sequence_window)
			if right_split_anchor >= len(original_sequence):
				right_split_anchor = len(original_sequence)

			if re.search(str(record.seq[left_split_anchor:right_split_anchor]), str(original_sequence)):
				count += 1
				left_split_anchor = right_split_anchor + 1
				continue

			split_record = SeqIO.SeqRecord(record.seq[left_split_anchor:right_split_anchor], f"{record.id}|split_{sequence_window}")
			filtered_split_records.append(split_record)
			left_split_anchor = right_split_anchor + 1

	return filtered_split_records


def fasta_export(fasta_records: list, dir_path: str, out_filename: str):
	with open(f"{dir_path}/{out_filename}", "w") as fasta_out:
		SeqIO.write(fasta_records, fasta_out, "fasta")


def main():
	# Argparse Variables
	args = parse_arguments()
	config_path = args.config
	fasta_in_path = args.fasta_file
	res_enzyme_path = args.res_enzyme
	split_parts = int(args.split_parts)


	# Load config file
	with open(config_path, "r") as f:
		config = yaml.load(f, Loader=yaml.FullLoader)

	# Parse variables from config
	back_table = config['codonTable_noStops']
	output_dir = config['output_dir']

	# Parse FASTA abjects
	try:
		dnaseq_record = SeqIO.read(fasta_in_path, "fasta")
	except (ValueError, TypeError):
		raise "The file provided is not correctly formatted as a FASTA sequence"
	dnaseq_in = Seq(dnaseq_record.seq)
	ptnseq_in = SeqIO.SeqRecord(str(dnaseq_in.translate()), dnaseq_record.id)

	print("FASTA sequence imported ")
	# Generate AA and NT single variant libs based on the provided backtable
	ptn_recs, gene_recs, backtranslated_base_sequence = generate_variants(ptnseq_in, back_table)

	# Check request for sequence split
	if split_parts:
		split_gene_recs = split_variants(gene_recs, backtranslated_base_sequence, split_parts)
		gene_recs = split_gene_recs

	if res_enzyme_path:
		pass_list, fail_list, refail_list = avoidREs_nodelchar(res_enzyme_path, gene_recs)
		# Export outputs
		fasta_export(pass_list, output_dir, f"{config['out_file_prefix']}.fna")
		fasta_export(fail_list, output_dir, f"{config['out_file_prefix']}_failRE.fna")

	print("Single variants successfully generated")
	# Export outputs
	fasta_export(ptn_recs, output_dir, f"{config['out_file_prefix']}.faa")
	if not res_enzyme_path:
		fasta_export(gene_recs, output_dir, f"{config['out_file_prefix']}.fna")
	print(f"FASTA libraries exported to {output_dir}:\n "
	      f"AA FASTA: {output_dir}/{config['out_file_prefix']}.faa\n "
	      f"NT FASTA: {output_dir}/{config['out_file_prefix']}.fna\n")


if __name__ == "__main__":
	main()

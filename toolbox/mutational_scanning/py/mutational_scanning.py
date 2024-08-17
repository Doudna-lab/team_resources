# Native modules
import re
import shutil
import copy
from argparse import ArgumentParser as argp
# Installed Modules
import yaml
from Bio.Seq import Seq
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

# Project modules
# DEBUG LINE
# from toolbox.mutational_scanning.py.re_check import avoidREs_nodelchar
from re_check import avoidREs_nodelchar


# DEBUG
# fasta_in_path = "/Users/bellieny/projects/team_resources/toolbox/mutational_scanning/fasta/TadA8e.fna"
# config_path = "/Users/bellieny/projects/team_resources/toolbox/mutational_scanning/config/tadA.yaml"
# res_enzyme_path = "/Users/bellieny/projects/team_resources/toolbox/mutational_scanning/data/pacI_restriction.csv"
# split_parts = 2
# add_3p_handle = "/Users/bellieny/projects/team_resources/toolbox/mutational_scanning/data/3prime_handle.fasta"
# add_5p_handle = "/Users/bellieny/projects/team_resources/toolbox/mutational_scanning/data/5prime_handle.fasta"

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
	parser.add_argument('--add_5p_handle',
						dest='add_5p_handle',
						default='',
						help='[Requires --split_parts] Provide a FASTA formatted file with sequence(s) to be added to both ends of the '
							 '5 prime resulting sequences. Each entry must be laballed with either ">5p" or "3p".')
	parser.add_argument('--add_3p_handle',
						dest='add_3p_handle',
						default='',
						help='[Requires --split_parts] Provide a FASTA formatted file with sequence(s) to be added to both ends of the '
							 '3 prime end sequences. Each entry must be laballed with either ">5p" or "3p".')
	parser.add_argument('--left_limit',
						dest='left_limit',
						default=0,
						help='DNA position in the <fasta_file> where the permutation should start. Actual position not 0-based. [Default: First NT in the input sequence]')
	parser.add_argument('--right_limit',
						dest='right_limit',
						default=0,
						help='DNA position in the <fasta_file> where the permutation should end. Actual position not 0-based [Default: Last NT in the input sequence]')

	# Parse arguments from the command line
	arguments = parser.parse_args()
	return arguments


def format_codon_tbl(codon_tbl: dict):
	formated_codon = {}
	for key in codon_tbl:
		formated_codon.setdefault(key.upper(), codon_tbl[key])
	return formated_codon


def generate_variants(ptnrec: SeqIO.SeqRecord, dnarec: Seq, backtable: dict, dna_left_boundary, dna_right_boundary):
	"""
	:param dna_right_boundary:
	:param dna_left_boundary:
	:param dnarec:
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
	# If the 'null' value of right_limit is provided, assume the full length of the sequence
	# if seq_right_limit == 0:
	# 	seq_right_limit = len(ptnseq_ls)
	# Set return variables
	ptn_record_ls = []
	gene_record_ls = []
	wt_backtranslated_seq = []

	# Set up a SeqRecord entry for the WT nucleotide sequence backtranslated
	wt_ptnseq = Seq(''.join(ptnseq))
	variant_backtrans_template = list(copy.deepcopy(dnarec))
	for wt_residue in wt_ptnseq:
		wt_backtranslated_seq.extend(format_backtable[wt_residue])
	variant_backtrans_template[dna_left_boundary:dna_right_boundary] = wt_backtranslated_seq
	variant_backtrans_template = ''.join(variant_backtrans_template)
	wt_backtranslated_seq = variant_backtrans_template
	print(f"WT backtrans reference sequence of length {len(wt_backtranslated_seq)}: {wt_backtranslated_seq}")
	SeqIO.SeqRecord(Seq(wt_backtranslated_seq), id=ptn_id, description='')

	print(f"Ref DNA seq: {dnarec}")

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

				variant_backtrans_seq = []
				variant_backtrans_template = list(copy.deepcopy(dnarec))
				for variand_aa in loop_ptnseq:
					variant_backtrans_seq.extend(''.join(format_backtable[variand_aa]))

				# print(f"Plug mutated sequence of length {len(variant_backtrans_seq)}: {variant_backtrans_seq}")
				# print(f"Into reference sequence of length {len(variant_backtrans_template)}: {variant_backtrans_template}")
				# print(f"Boundaries: {dna_left_boundary} -> {dna_right_boundary}")

				variant_backtrans_template[dna_left_boundary:dna_right_boundary] = variant_backtrans_seq
				variant_backtrans_template = ''.join(variant_backtrans_template)

				# print(f"Resuting merged sequence of length {len(variant_backtrans_template)}: {variant_backtrans_template}")
				gene_record_ls.append(
					SeqIO.SeqRecord(Seq(variant_backtrans_template), id=variant_fasta_label, description=''))
	return ptn_record_ls, gene_record_ls, wt_backtranslated_seq


def split_variants(gene_records, original_sequence, split_parts: int):
	count = 0
	filtered_split_records = []
	repeated_records = []
	for record in gene_records:
		add_sequence_greenlight = True
		left_split_anchor = 0
		for sequence_window in range(1, split_parts + 1):
			right_split_anchor = round(len(record.seq) / split_parts * sequence_window)
			if right_split_anchor >= len(original_sequence):
				right_split_anchor = len(original_sequence)

			if re.search(str(record.seq[left_split_anchor:right_split_anchor]), str(original_sequence)):
				count += 1
				left_split_anchor = right_split_anchor + 1
				continue

			split_record = SeqIO.SeqRecord(record.seq[left_split_anchor:right_split_anchor],
										   f"{record.id}|split_{sequence_window}")

			# Do not allow repeated sequences in the result
			for gathered_records in filtered_split_records:
				if split_record.seq == gathered_records.seq:
					repeated_records.append(split_record)
					add_sequence_greenlight = False
					break
			if add_sequence_greenlight:
				filtered_split_records.append(split_record)
			left_split_anchor = right_split_anchor + 1

	return filtered_split_records, repeated_records


def add_sequence_segments(original_records: list, sequence_segments: list, apply_to_split: int):
	modified_records = copy.deepcopy(original_records)
	for variant_record in modified_records:
		current_id = variant_record.id.split("|")[-1]
		split_info_field = current_id.split("_")[-1]
		if int(split_info_field) == int(apply_to_split):
			for handle_record in sequence_segments:
				if handle_record.id == '5p':
					variant_record.seq = handle_record.seq.upper() + variant_record.seq.upper()
				elif handle_record.id == '3p':
					variant_record.seq = variant_record.seq.upper() + handle_record.seq.upper()

	return modified_records


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
	add_5p_handle = args.add_5p_handle
	add_3p_handle = args.add_3p_handle
	dna_left_limit = int(args.left_limit)
	dna_right_limit = int(args.right_limit)

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

	# Format left- and right-limits from DNA to amino acid positions
	#  The internal format will be 0-based throughout the script
	ptn_left_limit = (dna_left_limit - 1) // 3
	ptn_right_limit = (dna_right_limit - 1) // 3


	# Assign the actual value (0-Based) to DNA rightmost boundary
	zero_based_dna_right_limit = dna_right_limit
	if dna_right_limit == 0:
		zero_based_dna_right_limit = len(dnaseq_in) + 1

	zero_based_dna_left_limit = dna_left_limit - 1

	# Validate sequence boundaries
	actual_dna_length = len(dnaseq_in[dna_left_limit:dna_right_limit]) + 1
	if actual_dna_length % 3 != 0:
		print(f'ERROR: The length of sequence to be processed must be divisible by 3. '
			  f'The current DNA sequence boundaries '
			  f'(between {dna_left_limit} and {dna_right_limit}) result in a {actual_dna_length} NT sequence')
		print('Exiting code')
		exit(0)

	# Translate NT to AA: This could be the entire DNA sequence or a segment per user-request
	ptnseq_in = SeqIO.SeqRecord(str(dnaseq_in[dna_left_limit-1:dna_right_limit].translate()), dnaseq_record.id)
	print(f"Process protein sequence \n{ptnseq_in.seq}\n")

	print("FASTA sequence imported ")

	print(f"DNA 0_based sequence boundaries: Left: {zero_based_dna_left_limit}, Right: {zero_based_dna_right_limit}")
	# Generate AA and NT single variant libs based on the provided backtable
	ptn_recs, gene_recs, backtranslated_base_sequence = generate_variants(ptnseq_in,
																		  dnaseq_in,
																		  back_table,
																		  zero_based_dna_left_limit,
																		  zero_based_dna_right_limit)
	# Check request for sequence split
	if split_parts:
		split_gene_recs, repeat_sequences = split_variants(gene_recs, backtranslated_base_sequence, split_parts)
		gene_recs = split_gene_recs
		fasta_export(repeat_sequences, output_dir, f"{config['out_file_prefix']}_repeated.fna")

		if add_5p_handle:
			add_5p_records = []
			with open(add_5p_handle, "r") as add_5p:
				for record in SeqIO.parse(add_5p, "fasta"):
					add_5p_records.append(record)
			gene_recs_5p_added = add_sequence_segments(gene_recs, add_5p_records, 1)
			gene_recs = gene_recs_5p_added

		if add_3p_handle:
			add_3p_records = []
			with open(add_3p_handle, "r") as add_3p:
				for record in SeqIO.parse(add_3p, "fasta"):
					add_3p_records.append(record)
			gene_recs_3p_added = add_sequence_segments(gene_recs, add_3p_records, 2)
			gene_recs = gene_recs_3p_added

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

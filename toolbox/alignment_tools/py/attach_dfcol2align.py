# === Native modules ===
import copy
from argparse import ArgumentParser as argp
# === Installed Modules ===
import pandas as pd
from pandas.errors import ParserError
from Bio import SeqIO

# === Example Run ===
# python py/attach_dfcol2align.py --fasta_id_col "Entry" --attach_col "Taxonomic lineage (Ids)" --out /groups/doudna/projects/daniel_projects/tuck_o/msa/ski2/ski2_msa_txid.fasta data/tuck_o/ski2_rep.fasta data/tuck_o/ski_idmaps.tsv

# DEBUG
# fasta_file = "/groups/doudna/team_resources/toolbox/alignment_tools/data/tuck_o/ski2_rep.fasta"
# reference_table = "/groups/doudna/team_resources/toolbox/alignment_tools/data/tuck_o/ski_idmaps.tsv"
# fasta_id_col = "Entry"
# attach_col = "Taxonomic lineage (Ids)"


def parse_arguments():
	#  Launch argparse parser
	parser = argp(
		prog='attach_dfcol2align.py',
		description='',
		usage='%(prog)s [options] <fasta_file> <reference_table>')
	# Define arguments
	parser.add_argument('fasta_file',
	                    help='Path to fasta input file')  # positional argument
	parser.add_argument('reference_table',
	                    default='/groups/doudna/team_resources/toolbox/alignment_tools/data/tuck_o/ski_idmaps.tsv',
	                    help='Path to reference table from which a column will be added to the FASTA id. [Default: ski_idmaps.tsv]')
	parser.add_argument('--fasta_id_col',
						dest='fasta_id_col',
						default='Entry',
						help='Specify which column name of the reference table contains the FASTA ids')
	parser.add_argument('--attach_col',
						dest='attach_col',
						default='Taxonomic lineage (Ids)',
						help='Specify which column name of the reference table contains the information will be added to the FASTA ids')
	parser.add_argument('--out',
						dest='output_file',
						default='out.fasta',
						help='Path to output FASTA')


	# Parse arguments from the command line
	arguments = parser.parse_args()
	return arguments


def fasta2recs(fasta_path):
	fasta_records = []
	with open(fasta_path, "r") as fasta_handle:
		for record in SeqIO.parse(fasta_handle, "fasta"):
			fasta_records.append(record)
	return fasta_records


def df2dict(df, key_column, value_column):
	reference_dict = {}
	df.apply(
		lambda row: reference_dict.setdefault(
			row.loc[key_column], row.loc[value_column]),
		axis=1
	)
	return reference_dict


def match_dict2fastarecs(records_list, match_dict):
	formatted_id_records = copy.deepcopy(records_list)
	for record in formatted_id_records:
		field_split = record.id.split("|")
		complaint = False
		complain_field = ''
		for field in field_split:
			try:
				record.id = str(field) + str(match_dict[field])
				complaint = False
				break
			except KeyError:
				complaint = True
				complain_field = field
		if complaint:
			print(f"Cound't find {complain_field}")
	return formatted_id_records


def handle_fasta_dups(records_list):
	check_dups = []
	checked_records = copy.deepcopy(records_list)
	for record in checked_records:
		if record.id in check_dups:
			record.id += "_dup"
		elif record.id not in check_dups:
			check_dups.append(record.id)
	return checked_records


def main():
	# Argparse Variables
	args = parse_arguments()
	fasta_file = args.fasta_file
	fasta_id_col = args.fasta_id_col
	attach_col = args.attach_col
	reference_table = args.reference_table
	output_file = args.output_file

	# === Import reference table
	try:
		df_reference = pd.read_csv(reference_table)
	except ParserError:
		df_reference = pd.read_csv(reference_table, sep="\t")

	fasta_records = fasta2recs(fasta_file)

	match_id_dict = df2dict(df_reference, fasta_id_col, attach_col)

	concat_fasta_records = match_dict2fastarecs(fasta_records, match_id_dict)

	dedup_concat_fasta_records = handle_fasta_dups(concat_fasta_records)

	with open(output_file, "w") as processed_fasta_handle:
		SeqIO.write(dedup_concat_fasta_records, processed_fasta_handle, "fasta")


if __name__ == "__main__":
	main()

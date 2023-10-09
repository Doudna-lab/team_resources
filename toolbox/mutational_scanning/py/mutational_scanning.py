# Native modules
from argparse import ArgumentParser as argp
# Installed Modules
import yaml
from Bio.Seq import Seq
from Bio import SeqIO


# DEBUG
# abs_path = "/Users/bellieny/projects/team_resources/toolbox/mutational_scanning/fasta/TadA8e.fna"
# config_path = "/Users/bellieny/projects/team_resources/toolbox/mutational_scanning/config/tadA.yaml"


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
	return ptn_record_ls, gene_record_ls


def fasta_export(fasta_records: list, dir_path: str, out_filename: str):
	with open(f"{dir_path}/{out_filename}", "w") as fasta_out:
		SeqIO.write(fasta_records, fasta_out, "fasta")


def main():
	# Argparse Variables
	args = parse_arguments()
	config_path = args.config
	fasta_in_path = args.fasta_file

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
	ptn_recs, gene_recs = generate_variants(ptnseq_in, back_table)

	print("Single variants successfully generated")
	# Export outputs
	fasta_export(ptn_recs, output_dir, f"{config['out_file_prefix']}.faa")
	fasta_export(gene_recs, output_dir, f"{config['out_file_prefix']}.fna")
	print(f"FASTA libraries exported to {output_dir}:\n "
	      f"AA FASTA: {output_dir}/{config['out_file_prefix']}.faa\n "
	      f"NT FASTA: {output_dir}/{config['out_file_prefix']}.fna\n")


if __name__ == "__main__":
	main()

# Native modules
import pandas as pd
from pandas.errors import IndexingError
from argparse import ArgumentParser as argp
import os
from os.path import abspath
import re
import subprocess
# Installed modules
import yaml


def parse_arguments():
	#  Launch argparse parser
	parser = argp(
		prog='process_barcodes',
		description='',
		usage='%(prog)s [options] sample_file.csv analysis_config.yaml')
	# Define arguments
	parser.add_argument('csv_file',
	                    help='Path to csv sample/barcode file')  # positional argument
	# Define arguments
	parser.add_argument('analysis_config',
	                    help='Path to main analysis config')  # positional argument
	parser.add_argument('-c',
	                    dest='config',
	                    default='config/process_barcodes.yaml',
	                    help='Specify the path to the yaml config file with information about the barcode table. [Default: config/process_barcodes.yaml]')
	# parser.add_argument('-p',
	#                     dest='database_path',
	#                     default=False,
	#                     help='Specify a custom path to a database. [Ignores -d]')
	# parser.add_argument('-e',
	#                     dest='evalue',
	#                     default=1e-10,
	#                     help='Specify the e-value threshold for the search. [Default: 1e-10]')
	# parser.add_argument('-t',
	#                     dest='threads',
	#                     default=1,
	#                     help='Specify the n of threads to compute this search. [Default: 1]')

	# Parse arguments from the command line
	arguments = parser.parse_args()
	return arguments


def get_filepath_list(parent_dir, name_pattern, file_extension):
	file_list = []
	for root, dirs, files in os.walk(parent_dir):
		for item in dirs + files:
			if re.search(fr"{name_pattern}\S+{file_extension}", item):
				filepath = os.path.abspath(os.path.join(root, item))
				file_list.append(filepath)
	return file_list


def extract_item_from_df(df, value, anchor_col, target_col):
	try:
		out_list = df[df[anchor_col] == value][target_col].iloc[0, :].tolist()
	except IndexingError:
		out_list = df[df[anchor_col] == value][target_col].tolist()

	if len(out_list) <= 1:
		return str(out_list[0])
	return out_list


def consolidate_records_dict(records_dict, new_key):
	out_dict = {}
	for idx in records_dict:
		for internal_key in records_dict[idx]:
			if internal_key == new_key:
				checked_key = records_dict[idx][internal_key]


def main():
	# Call argument parsing function
	args = parse_arguments()
	# Load values into variables
	csv_file_path = abspath(str(args.csv_file))
	config_path = abspath(str(args.config))
	analysis_config = abspath(str(args.analysis_config))

	# DEBUG lines
	# config_path = "toolbox/ngs_analysis/config/process_barcodes.yaml"
	# csv_file_path = "toolbox/ngs_analysis/dump/2023.05.26_NGS Barcodes.csv"

	# Load config files
	with open(config_path, "r") as f:
		config = yaml.load(f, Loader=yaml.FullLoader)
	with open(analysis_config, "r") as f:
		config_main = yaml.load(f, Loader=yaml.FullLoader)

	# Import csv sample file
	df = pd.read_csv(csv_file_path)
	# Drop rows with 'nan' in the first column
	df = df.dropna(subset=[df.columns[0]], how='any')

	# Get output path from the main analysis config file
	output_path = config_main['job_name']
	# Get column names
	run_id_col_name = str(df.columns[config['run_id_col']].tolist()[0])
	barcode_col_name = df.columns[config['barcode_seq_col']].tolist()
	# ref_seq_col_name = str(df.columns[config['reference_seq_col']].tolist()[0])
	sample_id_col_name = str(df.columns[config['sample_id_col']].tolist()[0])

	for run in set(df[run_id_col_name].tolist()):
		write_content = ''
		# current_reference = ''
		barcode_file_path = f'{output_path}{os.sep}{run}.tsv'
		# Gather input (multiplexed filenames)
		sample_id_list = extract_item_from_df(df, run, run_id_col_name, sample_id_col_name)
		# print(sample_id_list)
		fastq_path = ' '.join(get_filepath_list(config['reads_raw_dir_path'], run, 'fastq.gz'))
		print(fastq_path)
		for sample_id in sample_id_list:
			current_barcode = extract_item_from_df(df, sample_id, sample_id_col_name, barcode_col_name)
			# current_reference = extract_item_from_df(df, sample_id, sample_id_col_name, ref_seq_col_name)

			# Format strings
			barcodes_string = '\t'.join(current_barcode)
			write_content += f'{sample_id}\t{barcodes_string}\n'

		# Export barcode file
		with open(barcode_file_path, 'w') as f:
			f.write(write_content)

		# Setup shell command
		command = f"demultiplex match -p {output_path} {barcode_file_path} {fastq_path}"
		print(f"Executing command:\n->{command}")
		# Execute demultiplexing through shell
		subprocess.run(command, shell=True, capture_output=True, text=True)


if __name__ == "__main__":
	main()

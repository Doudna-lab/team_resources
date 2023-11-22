# Native modules
from itertools import islice, repeat, tee
from pandas.errors import IndexingError
from argparse import ArgumentParser as argp
import os
from os.path import abspath
import re
import subprocess
# Installed modules
from pandas.errors import IndexingError
import pandas as pd
import gzip
import yaml
from Bio import SeqIO


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
	parser.add_argument('-q',
	                    dest='quality_check_config',
						default='config/reads_qc.yaml',
						help='Specify the path to the read quality-check config file. [Default: config/reads_qc.yaml]')
	parser.add_argument('-c',
	                    dest='config',
	                    default='config/process_barcodes.yaml',
	                    help='Specify the path to the yaml config file with information about the barcode table. [Default: config/process_barcodes.yaml]')
	parser.add_argument('-m',
	                    dest='mode',
	                    default='only_info',
	                    help="""Execution modes include: 
	                    -> "add_info": a block of text will be provided ALONG with demultiplexing so it can be added to the QC yaml config
	                    -> "only_info": the script will not execute demultiplexing, instead it'll only print out the sample information block\n
	                    Default: "only_info" """)
	parser.add_argument('--mismatch',
	                    dest='mism',
	                    default=1,
	                    help="""Number of mismatches allowed when searching for barcodes in sequences.
		                    Default: 1 """)
	parser.add_argument('-a',
	                    dest='append',
	                    action='store_true',
	                    help="""Add this flag to append both sample->run and sample->reference blocks to the read quality-check config file""")
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


def sub_findre(string, substring, mismatch):
	sublen = len(substring)
	zip_gen = (zip(substring, islice(string, i, i+sublen), repeat(i, sublen)) for i in range(len(string)))
	for z in zip_gen:
		l, z = tee(z)
		if sum(1 for i, j, _ in l if i == j) >= sublen-mismatch:
			new = zip(*z)
			next(new)
			next(new)
			return next(new)[0]


def get_filepath_list(parent_dir, name_pattern, file_extension):
	file_list = []
	for root, dirs, files in os.walk(parent_dir):
		for item in dirs + files:
			if re.search(fr"{name_pattern}\S+{file_extension}", item):
				filepath = os.path.abspath(os.path.join(root, item))
				file_list.append(filepath)
	return file_list


def extract_item_from_df(df, value, anchor_col, target_col):
	out_list = []
	try:
		out_list = df[df[anchor_col] == value][target_col].iloc[0, :].tolist()
	except (IndexingError, IndexError):
		if len(df[df[anchor_col] == value][target_col]) > 1:
			out_list = df[df[anchor_col] == value][target_col].tolist()
		if len(df[df[anchor_col] == value][target_col]) <= 1:
			# print(f"ESSE CAPETA: {df[df[anchor_col] == value][target_col]}")
			# return df[df[anchor_col] == value][target_col]
			out_list = [df[df[anchor_col] == value][target_col].iloc[0]]

	if len(out_list) <= 1:
		return str(out_list[0])
	return out_list


def consolidate_records_dict(records_dict, new_key):
	out_dict = {}
	for idx in records_dict:
		for internal_key in records_dict[idx]:
			if internal_key == new_key:
				checked_key = records_dict[idx][internal_key]


def export_barcode_fasta(sample2seqlist, output_path):
	# def clean_config_string(string):
	# 	string = re.sub(r'-"\n', '"\n', string)
	# 	string = re.sub(r'-"', '-', string)
	# 	return string
	export_path = ''
	config_string = ''
	for sample in sample2seqlist:
		barc_config_string = ''
		for barcode_idx in range(len(sample2seqlist[sample])):
			export_path = f'{output_path}_{barcode_idx}.fasta'
			seq_string = f'>{sample}\n{sample2seqlist[sample][barcode_idx]}\n'
			with open(export_path, 'a') as f:
				f.write(seq_string)
	# 		barc_config_string += f'"{sample2seqlist[sample][barcode_idx]}-'
	# 	config_string += f'"{sample}": {barc_config_string}"\n  '
	# 	config_string = clean_config_string(config_string)
	# config_string = clean_config_string(config_string)
	return config_string


# Function to process each sequence
def barcode_match(sequence, barcode, mismatch_allowed):
	matches = None
	mismatch_count = None
	start_positions = [None]
	for mismatch_count in range(0, mismatch_allowed):
		pattern = f"{barcode}{{0,{mismatch_count}}}"
		match_pack = [(match.group(), match.start()) for match in re.finditer(pattern, sequence)]
		try:
			# Unpack the match object to separate variables
			matches, start_positions = zip(*match_pack)
			# matches = re.findall(pattern, sequence)
			if len(matches) > 1:
				return matches, mismatch_count, start_positions[0]
		except ValueError:
			pass
	return matches, mismatch_count, start_positions[0]


def demux(barcode_dict, barcode_window, fwd_fastq_path, rev_fastq_path, mismatch_allowed):
	# Process barcode file
	report_records = {}
	record_mismatches = {}
	read2sample_dict = {}
	read2record = {}
	rev_ids_set = []
	process_count = 0
	for sample in barcode_dict:
		print(f"######## SAMPLE {sample} ########")
		barcode_3p = ''
		barcode_5p = ''
		#  with gzip.open(filename, "rt") as handle:
		with gzip.open(fwd_fastq_path, "rt") as fwd_fastq:
			total_reads_count = 0
			fwd_ids_set = []
			fwd_3p_set = []
			# Check occurrences of the 5p barcode in the forward fragments of the reads
			thresh_count_report = 100
			for record in SeqIO.parse(fwd_fastq, "fastq"):
				total_reads_count += 1
				f_sequence_5p = str(record.seq)[0:barcode_window]
				f_sequence_3p = str(record.seq)[-barcode_window:]

				barcode_5p = barcode_dict[sample][0]
				# Evaluate the presence of 5p barcode on forward fragments
				(f_matches_5p, mismatch_count, match_start_coord) = barcode_match(f_sequence_5p, barcode_5p, mismatch_allowed)
				# match_positions = sub_findre(f_sequence_5p, barcode_5p, mismatch_allowed)
				# If found, add the forward 5p barcode information to the reporting dictionary
				if f_matches_5p:
					report_records.setdefault(sample, {}).setdefault("fwd_reads", {}).setdefault('5p_barcode', []).append(record)
					record_mismatches.setdefault(sample, {}).setdefault("fwd_reads", {}).setdefault(
						'5p_barcode', {}).setdefault(record.id, (match_start_coord, mismatch_count))
					# record_mismatches.setdefault(sample, {}).setdefault("fwd_reads", {}).setdefault(
					# 	'5p_barcode', {}).setdefault(record.id, match_positions)
					read2sample_dict.setdefault(record.id, sample)
					read2record.setdefault(record.id, record)
					fwd_ids_set.append(record.id)
					if len(barcode_dict[sample]) > 1:
						barcode_3p = barcode_dict[sample][1]
						# Evaluate the presence of 3p barcode on forward fragments
						(f_matches_3p, mismatch_count, match_start_coord) = barcode_match(f_sequence_3p, barcode_3p, mismatch_allowed)
						# match_positions = sub_findre(f_sequence_3p, barcode_3p, mismatch_allowed)
						# If found, add the reverse 3p barcode information to the reporting dictionary
						if f_matches_3p:
							report_records.setdefault(sample, {}).setdefault("fwd_reads", {}).setdefault('3p_barcode', []).append(record)
							record_mismatches.setdefault(sample, {}).setdefault("fwd_reads", {}).setdefault(
								'3p_barcode', {}).setdefault(record.id, (match_start_coord, mismatch_count))
							# record_mismatches.setdefault(sample, {}).setdefault("fwd_reads", {}).setdefault(
							# 	'3p_barcode', {}).setdefault(record.id, match_positions)
							fwd_3p_set.append(record.id)
							process_count += 1
					if process_count >= thresh_count_report:
						print(f"{process_count} Entries processed on forward - Sample: {sample}")
						thresh_count_report += thresh_count_report
			print(f"5p barcode ({barcode_5p}) processing finished")
			# Among the sequences carrying the 5p barcode, check the presence of the 3p in the reverse fragment of the read
			#  in case more than one barcode is provided
		if len(barcode_dict[sample]) > 1:
			rev_records = {}
			total_reads_count = 0
			process_count = 0
			# match_positions = None
			with gzip.open(rev_fastq_path, "rt") as rev_fastq:
				print(f"Processing reverse reads for which 5p barcodes were found")
				for record in SeqIO.parse(rev_fastq, "fastq"):
					total_reads_count += 1
					rev_records.setdefault(record.id, record.seq)
			for fwd_id in fwd_ids_set:
				thresh_count_report = 100
				try:
					rev_seq = rev_records[fwd_id]
					r_sequence_3p = str(rev_seq)[-barcode_window:]
					barcode_3p = barcode_dict[sample][1]
				except KeyError:
					continue
				# Evaluate the presence of 3p barcode on reverse fragments
				(r_matches_3p, mismatch_count, match_start_coord) = barcode_match(r_sequence_3p, barcode_3p, mismatch_allowed)
				# match_positions = sub_findre(r_sequence_3p, barcode_3p, mismatch_allowed)
				# If found, add the reverse 3p barcode information to the reporting dictionary
				if r_matches_3p:
					report_records.setdefault(sample, {}).setdefault("rev_reads", {}).setdefault('3p_barcode',[]).append(record)
					record_mismatches.setdefault(sample, {}).setdefault("rev_reads", {}).setdefault(
						'3p_barcode', {}).setdefault(record.id, (match_start_coord, mismatch_count))
					# record_mismatches.setdefault(sample, {}).setdefault("rev_reads", {}).setdefault(
					# 	'3p_barcode', {}).setdefault(record.id, match_positions)
					rev_ids_set.append(record.id)
					process_count += 1
					if process_count >= thresh_count_report:
						print(f"{process_count} Entries processed on reverse - Sample: {sample}")
						thresh_count_report += thresh_count_report
		success_ratio_fwd = (len(fwd_ids_set) * 100) / total_reads_count
		success_ratio_rev = (len(rev_ids_set) * 100) / len(fwd_ids_set)
		success_ratio_fwd_3p = (len(fwd_3p_set) * 100) / len(fwd_ids_set)

		print(f'{round(success_ratio_fwd, 2)}% ({len(fwd_ids_set)}) of the {total_reads_count} total reads had 5p barcode hits')
		print(f'Out of which {round(success_ratio_fwd_3p, 2)}% ({len(fwd_3p_set)}) of the {len(fwd_ids_set)} had 3p barcode ({barcode_3p}) hits in the fwd reads')
		print(f'Also, {round(success_ratio_rev, 2)}% ({len(rev_ids_set)}) of the {len(fwd_ids_set)} had 3p barcode ({barcode_3p}) hits in the rev reads')
	return report_records, record_mismatches

	# DEBUG lines
	# config_path = "/Users/bellieny/projects/team_resources/toolbox/ngs_analysis/config/process_barcodes-30-942700356.yaml"
	# quality_check_config = "/Users/bellieny/projects/team_resources/toolbox/ngs_analysis/config/reads_qc-30-942700356_nt93.yaml"
	# csv_file_path = "/Users/bellieny/projects/team_resources/toolbox/ngs_analysis/dump/30-942700356_nt93.csv"
	# output_path = "/groups/doudna/team_resources/toolbox/ngs_analysis/dump"
	# run = 'ED-II-47-1'
	# raw_fastq_path = '/groups/doudna/team_resources/toolbox/ngs_analysis/dump/00_fastq'
	# reference_list = []


def main():
	# Call argument parsing function
	args = parse_arguments()
	# Load values into variables
	csv_file_path = abspath(str(args.csv_file))
	config_path = abspath(str(args.config))
	quality_check_config = abspath(str(args.quality_check_config))
	mode = str(args.mode)
	mism = int(args.mism)
	append = args.append
	reference_list = []

	# Load config files
	with open(config_path, "r") as f:
		config = yaml.safe_load(f)
	with open(quality_check_config, "r") as f:
		config_main = yaml.safe_load(f)

	# Import csv sample file
	df = pd.read_csv(csv_file_path)
	# Drop rows with 'nan' in the first column
	df = df.dropna(subset=[df.columns[0]], how='any')

	# Get output path from the main analysis config file
	output_path = f"{config_main['job_name']}{os.sep}demux"
	if not os.path.exists(output_path):
		os.makedirs(output_path)
	# Get column names
	run_id_col_name = str(df.columns[config['run_id_col']].tolist()[0])
	barcode_col_name = df.columns[config['barcode_seq_col']].tolist()
	ref_seq_col_name = str(df.columns[config['reference_seq_col']].tolist()[0])
	sample_id_col_name = str(df.columns[config['sample_id_col']].tolist()[0])
	report_samples2run = "samples:\n"
	barcode_congig_string = 'barcode_to_sample:\n  '
	report_samples2ref = "reference:\n"
	raw_fastq_path = config['reads_raw_dir_path']

	for run in set(df[run_id_col_name].tolist()):
		sample_id_list = []
		print(f"============= RUN {run} =============\n")
		run2barcode_dict = {}
		write_content = ''
		barcode_file_path = f'{output_path}{os.sep}{run}.tsv'
		# Gather input (multiplexed filenames)
		print(run_id_col_name, "----", sample_id_col_name, "----", "\n DF: ", df)
		extracted_samples = extract_item_from_df(df, run, run_id_col_name, sample_id_col_name)
		if isinstance(extracted_samples, str):
			sample_id_list.append(extracted_samples)
		if isinstance(extracted_samples, list):
			sample_id_list.extend(extracted_samples)
		print(run, sample_id_list)
		# print(sample_id_list)
		fastq_path = get_filepath_list(raw_fastq_path, run, 'fastq.gz')
		for sample_id in sample_id_list:
			current_barcode = extract_item_from_df(df, sample_id, sample_id_col_name, barcode_col_name)
			current_reference = extract_item_from_df(df, sample_id, sample_id_col_name, ref_seq_col_name)

			# Format strings
			run2barcode_dict.setdefault(sample_id, current_barcode)
			# the string formatting below might become useless
			barcodes_string = '\t'.join(current_barcode)
			write_content += f'{sample_id}\t{barcodes_string}\n'
			report_samples2run += f"  '{sample_id}': '{run}'\n"
			report_samples2ref += f"  '{sample_id}': '{current_reference}'\n"
			reference_list.append(current_reference)

		barcode_path_prefix = f'{output_path}{os.sep}{run}'
		barcode_congig_string += export_barcode_fasta(run2barcode_dict, barcode_path_prefix)

		# /******************* TEST AREA FOR INTERNAL DEMUX *****************
		# (demux_dict, demux_mismatches) = demux(run2barcode_dict, 5, fastq_path[0], fastq_path[1], 1)
		# /******************* TEST AREA FOR INTERNAL DEMUX *****************

		# if mode == 'add_info':
		# 	# Export barcode file
		# 	with open(barcode_file_path, 'w') as f:
		# 		f.write(write_content)

		# Setup shell command
		# command = f"demultiplex match -p {output_path} {barcode_file_path} {fastq_path}"
		# print(f"Executing command:\n->{command}")
		# Execute demultiplexing through shell
		# subprocess.run(command, shell=True, capture_output=True, text=True)

	if append:
		with open(quality_check_config, 'a') as q:
			q.write(f"#REPORT: SAMPLE vs. RUN:\n{report_samples2run}\n")
			q.write(f"#REPORT: SAMPLE vs. REFERENCE SEQUENCE:\n{report_samples2ref}\n")
			q.write(f"#REFERENCE LIST:\nreference_sequences: {list(set(reference_list))}\n")
			q.write(f"#BARCODES vs. SAMPLE:\n{barcode_congig_string}\n")
	if not append:
		print("#REPORT: SAMPLE vs. RUN:\n", f"\b{report_samples2run}\n")
		print("#REPORT: SAMPLE vs. REFERENCE SEQUENCE:\n", f"\b{report_samples2ref}\n")
		print(f"#REFERENCE LIST:\nreference_sequences: {list(set(reference_list))}\n")


if __name__ == "__main__":
	main()

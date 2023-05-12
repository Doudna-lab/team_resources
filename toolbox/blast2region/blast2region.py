# Native modules
import argparse
import os
from argparse import ArgumentParser as argp
import copy
import re
import pandas as pd
# External modules
from Bio import SeqIO
from bioservices import UniProt
from Bio.SeqRecord import SeqRecord
from Bio.Blast.Applications import NcbiblastpCommandline
from Bio.Blast import NCBIXML
from Bio import Entrez
import yaml

# Load config file
with open("blast_config.yaml", "r") as f:
	config = yaml.load(f, Loader=yaml.FullLoader)


def parse_arguments():
	#  Launch argparse parser
	parser = argp(
		prog='blast2region',
		description='',
		usage='%(prog)s [options]')
	# Define arguments
	parser.add_argument('fasta_file',
	                    help='Path to fasta input file')  # positional argument
	parser.add_argument('-d',
	                    dest='database',
	                    default='nr',
	                    choices=['nr', 'pdb', 'sprot'],
	                    help='Specify which database(s) to search. [Default: nr]')
	parser.add_argument('-p',
	                    dest='database_path',
	                    default=False,
	                    help='Specify a custom path to a database. [Ignores -d]')
	parser.add_argument('-e',
	                    dest='evalue',
	                    default=1e-10,
	                    help='Specify the e-value threshold for the search. [Default: 1e-10]')
	parser.add_argument('-t',
	                    dest='threads',
	                    default=1,
	                    help='Specify the n of threads to compute this search. [Default: 1]')


	# Parse arguments from the command line
	arguments = parser.parse_args()
	return arguments


def check_format(filename):
	# Checks whether the fasta input is valid or not
	try:
		with open(filename, "r") as handle:
			fasta_check = SeqIO.parse(handle, "fasta")
			# for record in fasta_check:
			# 	print(record)
			check = any(fasta_check)
	except FileNotFoundError:
		print(f"The file <{filename}> could not be found.\nTerminating application")
		exit(0)
	if check:
		print("Input format validated")
	if not check:
		raise Exception("\nInvalid input format.\nTerminating application")


def path_to_databases(db_tag, db_path_dict):
	try:
		return db_path_dict[db_tag]
	except KeyError:
		print(f"Invalid database argument: {db_tag}.\nTerminating application")
		exit(0)


def id_parsing(id_string, sep_dict):
	loop_id = id_string
	for level in sep_dict:
		try:
			loop_id = loop_id.split(sep_dict[level][0])[sep_dict[level][1]]
		except IndexError:
			continue
	return loop_id


def hit_report_limit(query_to_hit_dict, nkeep):
	filtered_dict = {}
	hit_id_list = []
	# Within each blast query, keep a maximum of 50 hits based on lowest evaue
	for query in query_to_hit_dict:
		series = pd.DataFrame.from_dict(query_to_hit_dict[query], orient='index')
		fseries = series['e-value'].nsmallest(nkeep, keep='first')
		for hit in fseries.index.tolist():
			# Keep the same structure of the input dictionary
			filtered_dict.setdefault(query, {}).setdefault(hit, query_to_hit_dict[query][hit])
			hit_id_list.append(hit)
	return filtered_dict, hit_id_list


def blast_result_parser(xml_temp_path, eval_threshold, sep_instructions, nkeep):
	report_dict = {}
	for record in NCBIXML.parse(open(xml_temp_path)):
		if record.alignments:
			for align in record.alignments:
				for hsp in align.hsps:
					if hsp.expect < eval_threshold:
						id = id_parsing(align.hit_def, sep_instructions)
						report_dict.setdefault(record.query, {}).setdefault(id, {
							"blastp_alignment_len": align.length, "e-value": hsp.expect})
	filtered_report_dict, hit_id_list = hit_report_limit(report_dict, nkeep)
	return filtered_report_dict, hit_id_list


def ukb2ncbi(uid):
	u = UniProt(verbose=False)
	# gbk_id = u.mapping("UniProtKB_AC-ID", "EMBL-GenBank-DDBJ", query=uid, polling_interval_seconds=3, max_waiting_time=100)["results"][0]["to"]
	prot_id = u.mapping("UniProtKB_AC-ID", "EMBL-GenBank-DDBJ_CDS", query=uid, polling_interval_seconds=3, max_waiting_time=100)["results"][0]["to"]
	return prot_id


def elink_routine(db, hit_uid):
	dup_check = []
	not_found = ""
	linked = ""
	handle = Entrez.elink(dbfrom="protein", db=db, id=f"{hit_uid}")
	link_record = Entrez.read(handle)
	try:
		linked = link_record[0]['LinkSetDb'][0]['Link'][0]['Id']
		if linked not in dup_check:
			dup_check.append(linked)
	except (IndexError, KeyError):
		not_found = hit_uid
	handle.close()
	return linked, hit_uid, not_found


def ptn_to_nuc(query_to_hit_dict, db_list):
	progress = 0
	nucleotide_uid_dict = {}
	source2target = {}
	not_found_list = []
	# Loop through BLAST queries
	for query in query_to_hit_dict:
		# Set up progress report
		n = len(query_to_hit_dict[query].values())
		dup_check = []
		# Then go through each of their respective hits
		for hit in query_to_hit_dict[query]:
			progress += 1
			# Standardize protein identifiers to NCBI UIDs through ESearch
			handle = Entrez.esearch(db="protein", term=f"{hit}", idtype="acc")
			search_record = Entrez.read(handle)
			try:
				uid = search_record['IdList'][0]
			except IndexError:
				continue
			handle.close()

			# Loop through databases (found in config) and grab Nuccore UIDs
			for db_name in db_list:
				if uid in set(dup_check):
					continue
				loop_nuc_gi, loop_nuc_acc, not_found_hit = elink_routine(db_name, uid)
				if not_found_hit:
					loop_nuc_gi, loop_nuc_acc, c_not_found_hit = elink_routine(db_name,
					                                                           ukb2ncbi(not_found_hit))
					if not loop_nuc_gi:
						not_found_list.append(c_not_found_hit)
						continue
				if loop_nuc_gi:
					dup_check.append(uid)
					source2target.setdefault(query, {}).setdefault(loop_nuc_gi, (loop_nuc_acc, hit))
					nucleotide_uid_dict.setdefault(query, []).append(loop_nuc_gi)

	# Ouputs nuccore uids and the ptn->nuc uid links
	return nucleotide_uid_dict, source2target, list(set(not_found_list))


def nuc_to_gb(query_to_uid_dict):
	# Get Genbank records for each Nuccore UID
	gb_records = {}
	for query in query_to_uid_dict:
		for uid in query_to_uid_dict[query]:
			handle = Entrez.efetch(db="nucleotide", id=f"{uid}", rettype="gb", retmode="text")
			record = SeqIO.read(handle, "genbank")
			gb_records.setdefault(query, {}).setdefault(uid, record)
		# Returns a list  of Genbank SeqRecords objects
	return gb_records


def gb_plier(query_to_gb_dict, uid_to_acc, win_size):
	gbk_target = {}
	prot_dict = {}
	for query in query_to_gb_dict:
		for hit_uid in query_to_gb_dict[query]:
			gbk = query_to_gb_dict[query][hit_uid]
			for seq_feature in gbk.features:
				# Avoid blank feature that may occur in GenBank entries
				try:
					qualifiers = seq_feature.qualifiers
				except AttributeError:
					continue
				# Restrict search to protein-containing features
				if "protein_id" in qualifiers:
					prot_id = qualifiers["protein_id"][0]
					# Search for the protein-ids of interest
					if re.search(prot_id, uid_to_acc[query][hit_uid][0]):
						# Process feature information for future ref
						f_start = seq_feature.location.start.real
						f_end = seq_feature.location.end.real
						f_strand = seq_feature.strand
						f_seq = qualifiers["translation"][0]
						f_len = len(f_seq)
						highlight_feature = copy.deepcopy(seq_feature)
						highlight_feature.type = "highlight"
						# Set start/end coords using window size
						start = max(int(min([f_start, f_end])) - win_size, 0)
						end = min(int(max([f_start, f_end])) + win_size + 1, len(gbk.seq))

						# Create a SeqRecord object with the feature of interest
						gbk_focused = SeqRecord(
							id=gbk.id,
							annotations=gbk.annotations,
							dbxrefs=gbk.dbxrefs,
							seq=gbk.seq[start:end + 1],
							description=gbk.description
						)
						gbk_focused.features.append(highlight_feature)
						# Gather protein data for reference
						prep_prot_dict = {"query": query,
						                  "nuccore_acc": gbk.id,
						                  # "region_seq": gbk.seq[start:end + 1],
						                  "window_start": start,
						                  "window_end": end,
						                  "feature_start": f_start,
						                  "feature_end": f_end,
						                  "strand": f_strand,
						                  "feature_len": f_len,
						                  "blastp_hit": prot_id,
						                  "ptn_sequence": f_seq
						                  }
						prot_dict.setdefault(query, {}).setdefault(uid_to_acc[query][hit_uid][1], prep_prot_dict)
						gbk_target.setdefault(query, {}).setdefault(f"{gbk.id}_{start}-{end}", gbk_focused)

	return gbk_target, prot_dict


def write_reports(list_report_dicts):
	df_out = pd.DataFrame()
	for query in list_report_dicts[0]:
		query_df = pd.DataFrame()
		for report_dict in list_report_dicts:
			df_loop = pd.DataFrame.from_dict(report_dict[query], orient='index')
			query_df = pd.merge(query_df, df_loop, how='outer', right_index=True, left_index=True)
		df_out = pd.concat([df_out, query_df])
	return df_out


def export_gbs(query_to_gb_dict, parent_path):
	if not os.path.exists(parent_path):
		os.mkdir(parent_path)
	for query in query_to_gb_dict:
		query_suffix = re.sub(r'\|', '_', query[0:20])
		out_path = f"{parent_path}{os.sep}{query_suffix}"
		if not os.path.exists(out_path):
			os.mkdir(out_path)
		for hit in query_to_gb_dict[query]:
			gbk = query_to_gb_dict[query][hit]
			filename = f"{hit}.gb"
			with open(f"{out_path}{os.sep}{filename}", "w") as gb_handle:
				SeqIO.write(gbk, gb_handle, "genbank")

#
#
# #
# fasta_in = "blast2region/toy_hk_query.fasta"
# evalue = 0.005
# threads = 2
# temp = "blast2region/temp.xml"
# out_path = "delete.out"
# # Set database tag and path
# db_path = "blast2region/toy_hk_db"
# db_tag = 'sprot'
# with open("blast2region/blast_config.yaml", "r") as f:
# 	config = yaml.load(f, Loader=yaml.FullLoader)
# efecth_db_list = config['efetch_db']
#


def main():
	# Call argument parsing function
	args = parse_arguments()
	# Load values into variables
	# Set blast arguments
	fasta_in = args.fasta_file
	evalue = args.evalue
	threads = args.threads
	temp = "temp.xml"
	# NCBI databases to search
	efecth_db_list = config['efetch_db']
	# Set database tag and path
	db_path = path_to_databases(args.database, config["db_path"])
	db_tag = args.database
	if args.database_path:
		db_path = args.database_path
		db_tag = "custom"

	# FASTA format validation
	print("Input file validation")
	check_format(fasta_in)

	# Blast run
	print("BlastP Run in progress")
	cline = NcbiblastpCommandline(query=fasta_in,
	                              db=db_path,
	                              num_threads=threads,
	                              evalue=evalue,
	                              out=temp,
	                              outfmt=5)
	print(f"Blast command used:\n {cline}")
	cline()

	# Format parser
	print("Parsing BlastP output")
	hit_dict, hit_list = blast_result_parser(temp,
	                                         evalue,
	                                         config["id_sep_dict"][db_tag],
	                                         config["blast_nkeep"])
	# Remove blast XML temp file
	os.remove('temp.xml')
	# Entrez authentication
	print("Entrez login")
	Entrez.email = config["entrez_login"]

	# Query NCBI to get nuccore UIDs associated with the protein hits using ESearch/ELink
	print("Linking protein hit ids to Nuccore entries")
	nuc_uid_per_query, hit_to_link, hits_not_found = ptn_to_nuc(hit_dict, efecth_db_list)

	# Get GenBank entries through EFetch
	print("Retrieving GenBank objects")
	gb_seqrec_per_query = nuc_to_gb(nuc_uid_per_query)

	# Narrow down Genbank files based on hit UIDs
	print("Extracting relevant features from GenBank objects")
	targeg_gb_dict, target_hit_dict = gb_plier(gb_seqrec_per_query, hit_to_link, 2000)

	# Generate summary dataframe
	print("Generate reports")
	df = write_reports([target_hit_dict, hit_dict])

	# Export outputs
	print("Export files")
	export_gbs(targeg_gb_dict, config['dir_out'])
	df.to_csv(f"{config['dir_out']}/summary_report.csv", index=False)


if __name__ == "__main__":
	main()

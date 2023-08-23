# Installed modules
import pandas as pd

# DEBUG INPUTS
# blastout_path = "/groups/doudna/projects/daniel_projects/boger_r/phyrec_screening/psiblast_out/db-phyrec_db_2023-08_query-test.out"
# config_path = "/groups/doudna/team_resources/toolbox/phyrec_screening/config/phyrec_processing.yaml"
# # Load config files
# import yaml
# with open(config_path, "r") as f:
# 	config = yaml.load(f, Loader=yaml.FullLoader)


def main():
	# Snakemake Imports
	blastout_path = str(snakemake.input.psiblast_out)
	blast_col_names = str(snakemake.params.blast_col_names)
	output_taxcount_table = str(snakemake.output.taxid_counts)

	# Import blastout table
	blast_col_names_list = blast_col_names.split(" ")
	blast_df = pd.read_csv(blastout_path,
	                       names=blast_col_names_list,
	                       index_col=False).convert_dtypes().infer_objects()
	# Process pd Dataframe and take the unique set of hit IDs
	all_uniq_blast_ids = set(blast_df['sacc'].dropna().tolist())
	id_to_taxid = blast_df[blast_df['sacc'].isin(all_uniq_blast_ids)][['sacc', 'staxid']]

	# Count occurrences of each TaxID in the relevant column
	taxid_counts_df = id_to_taxid['staxid'].value_counts().reset_index()
	taxid_counts_df.columns = ['sacc', 'staxid_count']

	taxid_counts_df.to_csv(output_taxcount_table, sep="\t")
	#
	# # Start the stopwatch / counter
	# t1_start = process_time()
	#
	# # Setup a support variable to control the IDs search across the taxid reference
	# processed_line = 0
	# line_report_threshold = 1000000
	# processed_ids_count = 0
	# t1_loop = 0
	# blasthit2taxid = {}
	# # Search taxid match file -> don't process any line unnecessarily
	# with open(taxid_match_path, 'r') as tax_handle:
	# 	for line in tax_handle:
	# 		processed_line += 1
	# 		if line.split("\t")[0] in unique_blast_hits:
	# 			processed_ids_count += 1
	# 			blasthit2taxid.setdefault(line.split("\t")[0], line.split("\t")[1])
	# 		if processed_ids_count == len(unique_blast_hits):
	# 			break
	# 		if processed_line >= line_report_threshold:
	# 			line_report_threshold += 1000000
	# 			t1_current = process_time()
	# 			print(f"Took {t1_current - t1_loop} to Process 1M lines. Now {processed_line} processed lines in total")
	# 			t1_loop = process_time()
	#
	#
	# # Stop the stopwatch / counter
	# t1_stop = process_time()
	#
	# print("Elapsed time:", t1_stop, t1_start)
	#
	# print("Elapsed time during the whole program in seconds:",
	#       t1_stop - t1_start)


if __name__ == "__main__":
	main()

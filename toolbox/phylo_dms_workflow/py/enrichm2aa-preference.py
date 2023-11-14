# Installed modules
import pandas as pd
from matplotlib import pyplot as plt
from Bio import SeqIO


def main():
	# Snakemake I/O
	# Inputs
	dms_in = str(snakemake.input.dms_data)
	# Outputs
	dms_aapref_out = str(snakemake.output.dms_out)
	# Params
	position_col = str(snakemake.params.position_col)
	enrichment_col = str(snakemake.params.enrichment_col)
	aminoacid_col = str(snakemake.params.aminoacid_col)

	# DEBUG INPUT
	dms_in = "/Users/bellieny/projects/team_resources/toolbox/phylo_dms_workflow/scratch/betalactamase/betaLactamase_enrichmentScores.txt" # "/groups/doudna/projects/daniel_projects/prywes_n/input_data/dms_data.csv"
	# position_col = "Site"
	# enrichment_col = "Trial_1_AmpConc_2500" # "5percent_CO2_20uM_IPTG"
	# aminoacid_col = "AminoAcid"

	minenrichment = 1.0e-4  # minimum allowed enrichment
	df = pd.read_csv(dms_in)

	df["preference"] = [max(minenrichment, 10**df[enrichment_col][x]) for x in range(len(df))]

	# Format the AA preference table
	df_aapref = df.pivot(index=position_col, columns=aminoacid_col, values="preference")
	df_aapref.fillna(1, inplace=True)
	df_aapref = df_aapref.div(df_aapref.sum(axis=1), axis=0)
	df_aapref.insert(0, "site", range(1, len(df_aapref) + 1))

	# Export output
	df_aapref.to_csv(dms_aapref_out, index=False)



	dms_in_r = "/Users/bellieny/projects/team_resources/toolbox/phylo_dms_workflow/scratch/rbsc_enrichmentScores.csv"
	df_r = pd.read_csv(dms_in_r)
	with open("/Users/bellieny/projects/team_resources/toolbox/phylo_dms_workflow/ptn2locus_reports/rrub_rbsc.fasta",
	          'r') as f:
		rbsc_seq = SeqIO.read(f, 'fasta')

	df_r.fillna(0, inplace=True)

	df_r['minus1'] = df_r['5percent_CO2_20uM_IPTG'] - 1

	df_r["preference"] = [max(minenrichment, 10 ** df_r["minus1"][x]) for x in range(len(df_r))]
	df_r = df_r.pivot(index="Site", columns="AminoAcid", values="preference")
	for idx in range(len(rbsc_seq)):
		df_r.iat[idx, list(df_r.columns).index(rbsc_seq[idx])] = 10 ** 1
	df_r.fillna(1, inplace=True)
	df_r = df_r.div(df_r.sum(axis=1), axis=0)
	df_r.insert(0, "site", range(1, len(df) + 1))

if __name__ == "__main__":
	main()

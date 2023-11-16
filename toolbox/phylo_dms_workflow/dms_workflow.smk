# **** Variables ****
configfile: "config/prywes_dms.yaml"

# **** Imports ****
import glob

# Cluster run template
#nohup snakemake --snakefile dms_workflow.smk -j 5 --cluster "sbatch -t {cluster.time} -n {cluster.cores}" --cluster-config config/cluster.yaml --latency-wait 120 --use-conda &

# noinspection SmkAvoidTabWhitespace
rule all:
	input:
		expand("{run}/processed_inputs/aa_preference.csv",run=config["run"]),
		expand("{run}/figures/aa_preference.pdf", run=config["run"]),
		expand("{run}/processed_inputs/alignment.fna", run=config["run"])

# noinspection SmkAvoidTabWhitespace
rule convert_enrichment:
	input:
		dms_data = lambda wildcards: glob.glob("{in_dir}/dms_data.csv".format(in_dir=config['input_dir']))
	output:
		dms_out = "{run}/processed_inputs/aa_preference.csv"
	params:
		position_col = config["position_col"],
		enrichment_col = config["enrichment_col"],
		aminoacid_col = config["aminoacid_col"]
	conda:
		"envs/dms.yaml"
	script:
		"py/enrichm2aa-preference.py"

# noinspection SmkAvoidTabWhitespace
rule aa_preference_logo:
	input:
		dms_out = "{run}/processed_inputs/aa_preference.csv"
	output:
		aa_logo = "{run}/figures/aa_preference.pdf"
	conda:
		"envs/dms.yaml"
	shell:
		"""
		phydms_logoplot --prefs {input.dms_out} {output.aa_logo}
		"""

rule align_sequences:
	input:
		seq_in = lambda wildcards: glob.glob("{in_dir}/sequences.fna".format(in_dir=config['input_dir'])),
	output:
		aligned_sequences = "{run}/processed_inputs/alignment.fna"
	params:
		reference_seq = config["ref_sequence_id"],
		minident = config["minident"]
	conda:
		"envs/dms.yaml"
	shell:
		"""
		phydms_prepalignment --minidentity {params.minident}  {input.seq_in} {output.aligned_sequences} {params.reference_seq}
		"""
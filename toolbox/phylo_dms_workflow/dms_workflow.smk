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

rule inspect_prefs:
	input:
		dms_out = "{run}/processed_inputs/aa_preference.csv",
	output:
		radial_prefs = "{run}/figures/aa_preference.svg"
	conda:
		"envs/dms.yaml"
	script:
		"py/inspect_prefs.py"


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

# noinspection SmkAvoidTabWhitespace
rule seq_alignment:
	input:
		multi_fasta = lambda wildcards: glob.glob("{in_dir}/input_multi_fasta.csv".format(in_dir=config['input_dir']))
	output:
		msa = "{run}/processed_inputs/nt_alignment.csv"
	params:
		min_ident = config["min_ident"]
	shell:
		"""
		phydms_prepalignment {input.multi_fasta} {output.msa} --minidentity {params.min_ident}
		"""

# noinspection SmkAvoidTabWhitespace
rule phydms:
	input:
		msa = "{run}/processed_inputs/nt_alignment.csv",
		dms_out= "{run}/processed_inputs/aa_preference.csv"

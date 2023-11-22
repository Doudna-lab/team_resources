# **** Variables ****
configfile: "config/prywes_dms.yaml"

# **** Imports ****
import glob

# Cluster run template
#nohup snakemake --snakefile dms_workflow.smk -j 5 --cluster "sbatch -t {cluster.time} -n {cluster.cores}" --cluster-config config/cluster.yaml --latency-wait 120 --use-conda &

# noinspection SmkAvoidTabWhitespace
rule all:
	input:
		expand("{run}/{experiment_id}/processed_inputs/aa_preference.csv",run=config["run"],experiment_id=config["experiment_id"]),
		expand("{run}/{experiment_id}/figures/aa_preference.pdf", run=config["run"],experiment_id=config["experiment_id"]),
		expand("{run}/{experiment_id}/processed_inputs/nt_alignment_msa.fna", run=config["run"],experiment_id=config["experiment_id"])

# noinspection SmkAvoidTabWhitespace
rule convert_enrichment:
	input:
		dms_data = lambda wildcards: glob.glob("{in_dir}/dms_data.csv".format(in_dir=config['input_dir']))
	output:
		dms_prefs = "{run}/{experiment_id}/processed_inputs/aa_preference.csv"
	params:
		position_col = config["position_col"],
		enrichment_col = config["enrichment_col"],
		aminoacid_col = config["aminoacid_col"]
	conda:
		"envs/dms.yaml"
	script:
		"py/enrichm2aa-preference.py"

# noinspection SmkAvoidTabWhitespace
rule inspect_prefs:
	input:
		dms_prefs = "{run}/{experiment_id}/processed_inputs/aa_preference.csv",
	output:
		radial_prefs = "{run}/{experiment_id}/figures/aa_preference.png"
	params:
		site_offset = config["site_offset"]
	conda:
		"envs/dms.yaml"
	script:
		"py/inspect_prefs.py"

# noinspection SmkAvoidTabWhitespace
rule aa_preference_logo:
	input:
		dms_prefs = "{run}/{experiment_id}/processed_inputs/aa_preference.csv"
	output:
		aa_logo = "{run}/{experiment_id}/figures/aa_preference.pdf"
	conda:
		"envs/dms.yaml"
	shell:
		"""
		phydms_logoplot --prefs {input.dms_prefs} {output.aa_logo}
		"""

# noinspection SmkAvoidTabWhitespace
rule seq_alignment:
	input:
		multi_fasta = "{run}/{experiment_id}/multi_fasta.fna"
	output:
		msa = "{run}/{experiment_id}/processed_inputs/nt_alignment_msa.fna"
	params:
		min_ident = config["min_ident"],
		reference_seq = lambda wildcards: config["reference_seq"][wildcards.experiment_id]
	conda:
		"envs/dms.yaml"
	shell:
		"""
		phydms_prepalignment {input.multi_fasta} {output.msa} {params.reference_seq} --minidentity {params.min_ident}
		"""

# noinspection SmkAvoidTabWhitespace
rule phydms:
	input:
		msa = "{run}/{experiment_id}/processed_inputs/nt_alignment_msa.fna",
		dms_prefs = "{run}/{experiment_id}/processed_inputs/aa_preference.csv"
	output:
		dms_models = "{run}/{experiment_id}/phydmsresults/modelcomparison.md"
	params:
		outdir = "{run}/{experiment_id}/phydmsresults",
		rax_path = config["rax_path"],
		outprefix = "phydms_run_"
	conda:
		"envs/dms.yaml"
	threads:
		config["threads"]
	shell:
		"phydms_comprehensive --omegabysite --diffprefsbysite --raxml {params.rax_path}  --ncpus {threads} {params.outdir}/{params.outprefix} {input.msa} {input.dms_prefs}"

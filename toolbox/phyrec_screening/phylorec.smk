# **** Variables ****
configfile: "config/cluster.yaml"

# **** Imports ****
import glob

# Cluster run template
#nohup snakemake --snakefile expanded_search.smk -j 5 --cluster "sbatch -t {cluster.time}" --cluster-config config/cluster.yaml --latency-wait 120 --use-conda &

# noinspection SmkAvoidTabWhitespace
rule all:
	input:
		#
		expand("{run}/fasttree2/{sequence_id}/{sequence_id}_msa_txid.nwk",
			run=config["run"],sequence_id=config["sequence_id"]),


# noinspection SmkAvoidTabWhitespace
rule phylogenetic_reconstruction:
	input:
		merged_msa_fasta_txid = "{run}/msa/{sequence_id}/{sequence_id}_msa_txid.fasta"
	output:
		newick_tree = "{run}/fasttree2/{sequence_id}/{sequence_id}_msa_txid.nwk"
	conda:
		"envs/fasttree.yaml"
	threads:
		config["threads"]
	shell:
		"""
		FastTreeMP -gamma < {input.merged_msa_fasta_txid} > {output.newick_tree}
		"""

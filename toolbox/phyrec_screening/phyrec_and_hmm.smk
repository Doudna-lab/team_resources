# **** Variables ****
configfile: "config/phyrec_processing.yaml"
configfile: "config/cluster.yaml"

# **** Imports ****
import glob

# Cluster run template
#nohup snakemake --snakefile expanded_search.smk -j 5 --cluster "sbatch -t {cluster.time} -n {cluster.cores} -N {cluster.nodes}" --cluster-config config/cluster.yaml --latency-wait 120 --use-conda &

# noinspection SmkAvoidTabWhitespace
rule all:
	input:
		# Phylogenetic reconstruction
		# expand("{run}/iqtree/db-{db_prefix}_query-{input_prefix}.treefile",
		# 	run=config["run"], db_prefix=config["db_prefix"], input_prefix=config["input_prefix"]),
		# Fast phylogenetic reconstruction for probing
		expand("{run}/fasttree/db-{db_prefix}_query-{input_prefix}.nwk",
			run=config["run"],db_prefix=config["db_prefix"],input_prefix=config["input_prefix"]),
		# Build HMM based on the consolidated alignment
		expand("{run}/hmm/db-{db_prefix}_query-{input_prefix}.hmm",
			run=config["run"],db_prefix=config["db_prefix"],input_prefix=config["input_prefix"])

# noinspection SmkAvoidTabWhitespace
# rule phylogeny:
# 	input:
# 		post_search_msa = "{run}/clustalo/db-{db_prefix}_query-{input_prefix}.msa.fasta"
# 	output:
# 		phylogenetic_reconstruction = "{run}/iqtree/db-{db_prefix}_query-{input_prefix}.treefile"
# 	params:
# 		output_prefix = "{run}/iqtree/db-{db_prefix}_query-{input_prefix}"
# 	conda:
# 		"envs/iqtree.yaml"
# 	threads:
# 		config["phylogeny"]["cores"]
# 	resources:
# 		mem_mb = config["ram_phylogeny"],
# 	shell:
# 		"""
# 		iqtree -s {input.post_search_msa} -m TEST -B 100 -pre {params.output_prefix} -st AA -v
# 		"""

# noinspection SmkAvoidTabWhitespace
rule fast_phylogeny:
	input:
		jpred_concise_msa = lambda wildcards: glob.glob("{run}/jpred/{jpred_id}.align".format(run=wildcards.run, jpred_id=config["jpred_id"]))
	output:
		fast_phylogenetic_reconstruction = "{run}/fasttree/db-{db_prefix}_query-{input_prefix}.nwk"
	params:
		output_prefix = "{run}/iqtree/db-{db_prefix}_query-{input_prefix}"
	conda:
		"envs/fasttree.yaml"
	threads:
		config["phylogeny"]["cores"]
	resources:
		mem_mb = config["ram_phylogeny"],
	shell:
		"""
		FastTree -boot 1000 -out {output.fast_phylogenetic_reconstruction} {input.jpred_concise_msa} 
		"""

# noinspection SmkAvoidTabWhitespace
rule build_hmm:
	input:
		jpred_concise_msa = lambda wildcards: glob.glob("{run}/jpred/{jpred_id}.align".format(run=wildcards.run, jpred_id=config["jpred_id"]))
	output:
		hmmbuild = "{run}/hmm/db-{db_prefix}_query-{input_prefix}.hmm"
	params:
		hmm_summary = "{run}/hmm/db-{db_prefix}_query-{input_prefix}_summary.txt"
	threads:
		config["threads"]
	shell:
		"""
		hmmbuild --cpu {threads} --amino -n {wildcards.input_prefix} -o {params.hmm_summary} {output.hmmbuild} {input.jpred_concise_msa}
		"""

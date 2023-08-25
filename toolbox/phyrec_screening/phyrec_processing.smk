# **** Variables ****
configfile: "config/phyrec_processing.yaml"
configfile: "config/cluster.yaml"

# **** Imports ****
import glob

# Cluster run template
#nohup snakemake --snakefile phyrec_processing.smk -j 5 --cluster "sbatch -t {cluster.time} -n {cluster.cores} -N {cluster.nodes}" --cluster-config config/cluster.yaml --latency-wait 120 --use-conda &

# noinspection SmkAvoidTabWhitespace
rule all:
	input:
		# Generate BLAST database
		# expand("{run}/custom_db/{db_prefix}.phr",
		# 	run=config["run"],db_prefix=config["db_prefix"]),
		# PSI-Blast the MSA input using a list of sequence databases
		expand("{run}/psiblast_out/db-{db_prefix}_query-{input_prefix}.blastout",
			run=config["run"],db_prefix=config["db_prefix"], input_prefix=config["input_prefix"]),
		# Count TaxID occurrences in PSI-BLAST output
		expand("{run}/krona/db-{db_prefix}_query-{input_prefix}_taxid_counts.tsv",
			run=config["run"],db_prefix=config["db_prefix"], input_prefix=config["input_prefix"]),
		# Export FASTA sequences of PSIBLAST Hits containing taxid labels
		expand("{run}/psiblast_out/db-{db_prefix}_query-{input_prefix}_hits.fasta",
			run=config["run"],db_prefix=config["db_prefix"], input_prefix=config["input_prefix"]),
		# Generate Krona plot
		expand("{run}/krona/db-{db_prefix}_query-{input_prefix}_krona-plot.html",
			run=config["run"],db_prefix=config["db_prefix"], input_prefix=config["input_prefix"]),
		# Re-align the initial MSA with the PSI-BLAST hit sequences
		expand("{run}/clustalo/db-{db_prefix}_query-{input_prefix}.msa.fasta",
			run=config["run"], db_prefix=config["db_prefix"], input_prefix=config["input_prefix"]),
		# Phylogenetic reconstruction
		expand("{run}/iqtree/db-{db_prefix}_query-{input_prefix}.treefile",
			run=config["run"], db_prefix=config["db_prefix"], input_prefix=config["input_prefix"]),
		# Build HMM based on the consolidated alignment
		expand("{run}/hmm/db-{db_prefix}_query-{input_prefix}.hmm",
			run=config["run"], db_prefix=config["db_prefix"], input_prefix=config["input_prefix"])

# noinspection SmkAvoidTabWhitespace
# rule merge_and_makedb:
# 	input:
# 		included_databases = expand("{shared_db_path}/{included_db}",
# 			shared_db_path=config["shared_db_path"] ,
# 			included_db=config["included_databases"])
# 	output:
# 		merged_db = "{run}/custom_db/{db_prefix}.phr"
# 	params:
# 		taxid_map = lambda wildcards: glob.glob("{shared_db_path}/prot.accession2taxid.FULL".format(shared_db_path=config["shared_db_path"])),
# 		outdir = "{run}/custom_db/{db_prefix}",
# 		merged_fasta = "{run}/custom_db/{db_prefix}.fasta"
# 	message:  """
# 	Make BLAST database:
# 	input databases:\n ==> {input.included_databases}
# 	taxid map file:\n ==> {params.taxid_map}
# 	merged output:\n ==> {output.merged_db}
# 	"""
# 	shell:
# 		"""
# 		cat {input.included_databases} > {params.merged_fasta}
#         makeblastdb \
#         -dbtype prot \
#         -blastdb_version 5 \
#         -in {params.merged_fasta} \
#         -out {params.blastoutdir}
#         """

# noinspection SmkAvoidTabWhitespace
rule iterative_search:
	input:
		msa_in = "{run}/input/{input_prefix}.msa.fasta"
		# merged_db = lambda wildcards: glob.glob("{shared_db_path}/nr".format(shared_db_path=config["shared_db_path"])), #"{run}/custom_db/{db_prefix}.phr",
	output:
		psiblast_out = "{run}/psiblast_out/db-{db_prefix}_query-{input_prefix}.blastout"
	params:
		db = config["shared_db_path"],
		custom_cols = config["blast_custom_cols"],
	threads: config["iterative_search"]["cores"]
	message:
		"Blasting query :\n {input.msa_in}\nAgainst database:\n {params.db}\nGenerating:\n {output.psiblast_out}"
	shell:
		"""
        psiblast \
        -in_msa {input.msa_in} \
        -db {params.db} \
        -outfmt "10 {params.custom_cols}" \
        -num_threads {threads} \
        -num_iterations 10 \
        -max_hsps 1 \
        -subject_besthit \
        -inclusion_ethresh 1e-05 \
        -evalue 1e-2 \
        -qcov_hsp_perc 85 \
        -out {output.psiblast_out}
        """

# noinspection SmkAvoidTabWhitespace
rule taxid_parse:
	input:
		psiblast_out = "{run}/psiblast_out/db-{db_prefix}_query-{input_prefix}.blastout"
	output:
		taxid_counts = "{run}/krona/db-{db_prefix}_query-{input_prefix}_taxid_counts.tsv",
		hits_fasta= "{run}/psiblast_out/db-{db_prefix}_query-{input_prefix}_hits.fasta"
	params:
		blast_col_names = config["blast_custom_cols"]
	conda:
		"envs/bio.yaml"
	script:
		"py/blastout_taxid_count.py"

# noinspection SmkAvoidTabWhitespace
rule krona:
	input:
		taxid_counts = "{run}/krona/db-{db_prefix}_query-{input_prefix}_taxid_counts.tsv"
	output:
		krona_chart = "{run}/krona/db-{db_prefix}_query-{input_prefix}_krona-plot.html",
	params:
		taxdump_path = config["taxdump_path"]
	conda:
		"envs/krona.yaml"
	shell:
		"""		
		ktImportTaxonomy -m 2 -t 3 -tax {params.taxdump_path} -o {output.krona_chart} {input.taxid_counts}
		"""
	# ktUpdateTaxonomy.sh {params.taxdump_path}

# noinspection SmkAvoidTabWhitespace
rule realignment:
	input:
		hits_fasta = "{run}/psiblast_out/db-{db_prefix}_query-{input_prefix}_hits.fasta",
		msa_in= "{run}/input/{input_prefix}.msa.fasta"
	output:
		post_search_msa = "{run}/clustalo/db-{db_prefix}_query-{input_prefix}.msa.fasta"
	params:
		merged_input = "{run}/clustalo/db-{db_prefix}_query-{input_prefix}_merged-input.fasta",
	threads: config["realignment"]["cores"]
	shell:
		"""
		cat {input.msa_in} {input.hits_fasta} > {params.merged_input}
		clustalo --iter 10 --threads {threads} -i {params.merged_input} -o {output.post_search_msa} -v
		"""

# noinspection SmkAvoidTabWhitespace
rule phylogeny:
	input:
		post_search_msa = "{run}/clustalo/db-{db_prefix}_query-{input_prefix}.msa.fasta"
	output:
		phylogenetic_reconstruction = "{run}/iqtree/db-{db_prefix}_query-{input_prefix}.treefile"
	params:
		output_prefix = "{run}/iqtree/db-{db_prefix}_query-{input_prefix}"
	threads:
		config["phylogeny"]["cores"]
	shell:
		"""
		iqtree -s {input.post_search_msa} -mtree -m MFP -bb 1000 -pre {params.output_prefix} -st AA -nt {threads} -v 
		"""

rule build_hmm:
	input:
		post_search_msa = "{run}/clustalo/db-{db_prefix}_query-{input_prefix}.msa.fasta"
	output:
		hmmbuild = "{run}/hmm/db-{db_prefix}_query-{input_prefix}.hmm"
	params:
		hmm_summary = "{run}/hmm/db-{db_prefix}_query-{input_prefix}_summary.txt"
	threads:
		config["threads"]
	shell:
		"""
		hmmbuild --cpu {threads} --amino -n {wildcards.input_prefix} -o {params.hmm_summary} {output.hmmbuild} {input.post_search_msa}
		"""

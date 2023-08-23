# **** Variables ****
configfile: "config/phyrec_processing.yaml"

# **** Imports ****
import glob

# Cluster run template
#nohup snakemake --snakefile phyrec_processing.smk -j 5 --cluster "sbatch -t {cluster.time} -n {cluster.cores} -N {cluster.nodes}" --cluster-config config/cluster.yaml --latency-wait 120 --use-conda &

# noinspection SmkAvoidTabWhitespace
rule all:
	input:
		# Generate BLAST database
		expand("{run}/custom_db/{db_prefix}.phr",
			run=config["run"],db_prefix=config["db_prefix"]),
		# PSI-Blast the MSA input using a list of sequence databases
		expand("{run}/psiblast_out/db-{db_prefix}_query-{input_prefix}.out",
			run=config["run"],db_prefix=config["db_prefix"], input_prefix=config["iput_prefix"]),
		# Count TaxID occurrences in PSI-BLAST output
		expand("{run}/krona/db-{db_prefix}_query-{input_prefix}_taxid_counts",
			run=config["run"],db_prefix=config["db_prefix"], input_prefix=config["iput_prefix"])

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
#         -out {params.outdir}
#         """

# noinspection SmkAvoidTabWhitespace
rule iterative_search:
	input:
		msa_in = "{run}/input/{input_prefix}.msa.fasta"
		# merged_db = lambda wildcards: glob.glob("{shared_db_path}/nr".format(shared_db_path=config["shared_db_path"])), #"{run}/custom_db/{db_prefix}.phr",
	output:
		psiblast_out = "{run}/psiblast_out/db-{db_prefix}_query-{input_prefix}.out"
	params:
		db = config["shared_db_path"],
		custom_cols = config["blast_custom_cols"],
	threads:
		config["threads"]
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
        -evalue 1e-5 \
        -qcov_hsp_perc 85 \
        -out {output.psiblast_out}
        """

# noinspection SmkAvoidTabWhitespace
rule taxid_parse:
	input:
		psiblast_out = "{run}/psiblast_out/db-{db_prefix}_query-{input_prefix}.out"
	output:
		taxid_counts = "{run}/krona/db-{db_prefix}_query-{input_prefix}_taxid_counts.tsv"
	params:
		blast_col_names = config["blast_custom_cols"]
	conda:
		"envs/phyrec.yaml"
	script:
		"py/blastout_taxid_count.py"

# noinspection SmkAvoidTabWhitespace
rule krona:
	input:
		taxid_counts = "{run}/krona/db-{db_prefix}_query-{input_prefix}.taxid_counts.tsv"
	output:
		krona_chart = "{run}/krona/db-{db_prefix}_query-{input_prefix}_taxid_counts.tsv"
	params:
		taxdump_path = config["taxdump_path"]
	conda:
		"envs/phyrec.yaml"
	shell:
		"""
		ktUpdateTaxonomy.sh {params.taxdump_path}
		ktImportTaxonomy -m 2 -t 3 -tax {params.taxdump_path} -o {output.krona_chart} {input.taxid_counts}
		"""

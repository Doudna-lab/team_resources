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
		expand("{run}/psiblast_out/db|{db_prefix}_query|{input_prefix}.out",
			run=config["run"],db_prefix=config["db_prefix"], input_prefix=config["input_prefix"])

# noinspection SmkAvoidTabWhitespace
rule merge_and_makedb:
	input:
		included_databases = expand("{shared_db_dir}/{included_db}",
			shared_db_dir=config["shared_db_dir"] ,
			included_db=config["included_databases"])
	output:
		merged_db = "{run}/custom_db/{db_prefix}.phr"
	params:
		outdir = "{run}/custom_db/{db_prefix}",
		merged_fasta = "{run}/custom_db/{db_prefix}.fasta"
	message:  """
	Make BLAST database:
	input databases:\n ==> {input.included_databases}
	merged output:\n ==> {output.merged_db}
	Shell command:\n ==>\n
	cat {input.included_databases} > {params.merged_fasta}
        makeblastdb \
        -dbtype prot \
        -taxid 5 \
        -in {params.merged_fasta} \
        -out {params.outdir}
	"""
	shell:
		"""
		cat {input.included_databases} > {params.merged_fasta}
        makeblastdb \
        -dbtype prot \
        -taxid 5 \
        -in {params.merged_fasta} \
        -out {params.outdir}
        """

# noinspection SmkAvoidTabWhitespace
rule iterative_search:
	input:
		merged_db = "{run}/custom_db/{db_prefix}.phr",
		msa_in = "{run}/input/{input_prefix}.msa.fasta"
	output:
		psiblast_out = "{run}/psiblast_out/db|{db_prefix}_query|{input_prefix}.out"
	params:
		db = "{run}/custom_db/{db_prefix}",
		# input_prefix = lambda wildcards: config["input_prefix"][wildcards.input_prefix],
		# seq_id_intersection= "True",
		seq_format = "fasta"
	threads:
		config["threads"]
	message:
		"Blasting query :\n {input.msa_in}\nAgainst database:\n {input.merged_db}\nGenerating:\n {output.psiblast_out}"
	shell:
		"""
        psiblast \
        -query {input.msa_in} \
        -db {input.merged_db} \
        -outfmt 5 \
        -num_threads {threads} \
        -num_iterations 10
        -evalue 1e-10 \
        -qcov_hsp_perc 85 \
        -perc_identity 80 \
        -out {output.psiblast_out}
        """

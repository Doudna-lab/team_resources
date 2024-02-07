# **** Variables ****
configfile: "config/{}.yaml"

# **** Imports ****
import glob

# Cluster run template
#nohup snakemake --snakefile {}.smk -j 5 --cluster "sbatch -t {cluster.time} -n {cluster.cores} -N {cluster.nodes}" --cluster-config config/cluster.yaml --latency-wait 120 --use-conda &

# noinspection SmkAvoidTabWhitespace
rule all:
	input:
		expand("{run}",run=config[""]),

# noinspection SmkAvoidTabWhitespace
rule {}

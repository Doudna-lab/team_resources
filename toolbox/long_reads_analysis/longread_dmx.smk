# **** Variables ****
# configfile: "config/reads_qc-30-942700356_nt97.yaml"
# **** Imports ****
import glob

# Cluster run template
#nohup snakemake --snakefile reads_qc.smk --configfile config/reads_qc_garcia.yaml -j 5 --cluster "sbatch -t {cluster.time} -n {cluster.cores} -N {cluster.nodes}" --cluster-config config/cluster.yaml --latency-wait 60 --use-conda &

# noinspection SmkAvoidTabWhitespace
rule all:
    input:
        # Demultiplex reads
        expand("{run}/demux/demultiplexed-{sample}.fastq.gz", sample=config["samples"], run=config["job_name"]),
        # expand("{run}/clumpify/{sample}_dedup/{sample}_R1_dedup.fastq.gz", sample=config["samples"], run=config["job_name"]),

# rule demultiplex:
#     input:
#         r1 = lambda wildcards: glob.glob("{directory}/{run_name}_L001_R1_001.fastq.gz".format(directory=config["read_directory"],run_name=config["run_name"])),
#         r2 = lambda wildcards: glob.glob("{directory}/{sample}_L001_R2_001.fastq.gz".format(directory=config["read_directory"],run_name=config["run_name"]))
#     output:
#         r1_demux = "",
#         r2_demux = ""
rule demux_dual:
    input:
        r1 = lambda wildcards: glob.glob("{directory}/{platform_run}_S?_L001_R1_001.fastq.gz".format(directory=config["read_directory"],
            platform_run=config['samples'][wildcards.sample])),
        r2 = lambda wildcards: glob.glob("{directory}/{platform_run}_S?_L001_R2_001.fastq.gz".format(directory=config["read_directory"],
            platform_run=config['samples'][wildcards.sample])),
        fwd_barcorde= lambda wildcards: glob.glob("{run}/demux/{platform_run}_0.fasta".format(run=wildcards.run,
            platform_run=config['samples'][wildcards.sample])),
        rev_barcorde= lambda wildcards: glob.glob("{run}/demux/{platform_run}_1.fasta".format(run=wildcards.run,
            platform_run=config['samples'][wildcards.sample]))
    output:
        r1_demux = "{run}/demux/demultiplexed-{sample}_R1.fastq.gz",
        r2_demux = "{run}/demux/demultiplexed-{sample}_R2.fastq.gz"
    conda:
        "envs/biopy.yaml"
    message:
        "Demultiplexing read pairs:\nR1: {input.r1}\nR2: {input.r2}\nBarcodes:\n Fwd: {input.fwd_barcorde}\n Rev: {input.rev_barcorde}\n ==> Output:\n 1- {output.r1_demux}\n 2- {output.r2_demux}\nWildcards: {wildcards}"
    threads:
        config["threads"]
    shell:
        """
        cutadapt -j {threads} -e 1 --no-indels --pair-adapters -g file:{input.fwd_barcorde} -G file:{input.rev_barcorde} -o '{wildcards.run}/demux/demultiplexed-{{name}}_R1.fastq.gz' -p '{wildcards.run}/demux/demultiplexed-{{name}}_R2.fastq.gz' {input.r1} {input.r2}
        """
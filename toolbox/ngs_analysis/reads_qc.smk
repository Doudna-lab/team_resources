# **** Variables ****
configfile: "config/reads_qc-30-942700356_nt41.yaml"
# **** Imports ****
import glob

# Cluster run template
#nohup snakemake --snakefile reads_qc.smk -j 5 --cluster "sbatch -t {cluster.time} -n {cluster.cores} -N {cluster.nodes}" --cluster-config config/cluster.yaml --latency-wait 60 --use-conda &

# noinspection SmkAvoidTabWhitespace
rule all:
    input:
        # Demultiplex reads
        expand("{run}/demux/demultiplexed-{sample}_R1.fastq.gz", sample=config["samples"], run=config["job_name"]),
        # expand("{run}/clumpify/{sample}_dedup/{sample}_R1_dedup.fastq.gz", sample=config["samples"], run=config["job_name"]),
        #Reads QC Pre-Trim
        expand("{run}/fastqc_preTrim/", sample=config["samples"], run=config["job_name"]),
        expand("{run}/multiqc_preTrim/multiqc_report.html",run=config["job_name"]),
        #Reads Trimming and rRNA removal
        expand("{run}/trimmomatic/{sample}_trimmed/{sample}_dedup_P1.fastq.gz",sample=config["samples"], run=config["job_name"]),
        #Reads QC Post-Trim
        expand("{run}/fastqc_postTrim/",sample=config["samples"], run=config["job_name"]),
        expand("{run}/multiqc_postTrim/multiqc_report.html",run=config["job_name"]),
        expand("{run}/trim_summary/",run=config["job_name"])

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
# rule clumpify:
#     input:
#         r1 = "{run}/demux/demultiplexed-{sample}_R1.fastq.gz",
#         r2 = "{run}/demux/demultiplexed-{sample}_R2.fastq.gz"
#     output:
#         r1_dedup = "{run}/clumpify/{sample}_dedup/{sample}_R1_dedup.fastq.gz",
#         r2_dedup = "{run}/clumpify/{sample}_dedup/{sample}_R2_dedup.fastq.gz"
#     conda:
#         "envs/biopy.yaml"
#     message:
#         "Clumpifying read pairs:\nR1: {input.r1}\nR2: {input.r2}\n ==> Output:\n 1- {output.r1_dedup}\n 2- {output.r2_dedup}\nWildcards: {wildcards}"
#     threads:
#         config["threads"]
#     shell:
#         """
#         clumpify.sh in1={input.r1} in2={input.r2} out1={output.r1_dedup} out2={output.r2_dedup} dedupe=t optical=f
#         """
rule fastqc_1:
    input:
        r1_dedup = expand("{run}/demux/demultiplexed-{sample}_R1.fastq.gz", run=config["job_name"], sample=config["samples"]),
        r2_dedup = expand("{run}/demux/demultiplexed-{sample}_R2.fastq.gz", run=config["job_name"], sample=config["samples"])
    output:
        directory("{run}/fastqc_preTrim/")
    conda:
        "envs/biopy.yaml"
    message:
        "Running FastQC on pre-trim read pairs:\nR1: {input.r1_dedup}\nR2: {input.r2_dedup}"
    threads:
        config["threads"]
    shell:
        """
        mkdir {output}
        fastqc {input.r1_dedup} {input.r2_dedup} -t {threads} -o {output}
        """
rule multiqc_1:
    input:
        p1 = "{run}/fastqc_preTrim/"
    output:
        d1 = "{run}/multiqc_preTrim/multiqc_report.html"
    conda:
        "envs/biopy.yaml"
    params:
        folder = "{run}/multiqc_preTrim/"
    message:
        "Running MultiQC on pre-trim FastQC reports:\nDirectory: {input.p1}\nOutputs to {output.d1}"
    shell:
        """
        multiqc --force --outdir {params.folder} --interactive {input.p1}
        """
rule trimmomatic:
    input:
        r1_demux="{run}/demux/demultiplexed-{sample}_R1.fastq.gz",
        r2_demux="{run}/demux/demultiplexed-{sample}_R2.fastq.gz"
    output:
        p1 = "{run}/trimmomatic/{sample}_trimmed/{sample}_dedup_P1.fastq.gz",
        u1 = "{run}/trimmomatic/{sample}_trimmed/{sample}_dedup_U1.fastq.gz",
        p2 = "{run}/trimmomatic/{sample}_trimmed/{sample}_dedup_P2.fastq.gz",
        u2 = "{run}/trimmomatic/{sample}_trimmed/{sample}_dedup_U2.fastq.gz",
        report = "{run}/trimmomatic/{sample}_trimmed/{sample}_trim_report.txt"
    params:
        trim_params = config["trim_params"]
    conda:
        "envs/biopy.yaml"
    message:
        """Trimming read pairs:\nR1: {input.r1_demux}\nR2: {input.r2_demux}
        ==> Output:\n 1- {output.p1}\n 2- {output.p2} 
        ==> parameters: {params.trim_params}
        ==> Report:\n {output.report}"""
    threads:
        config["threads"]
    shell:
        """
        touch {output.u1} {output.u2} {output.report}
        cp {input.r1_demux} {output.p1} 
        cp {input.r2_demux} {output.p2} 
        trimmomatic PE -threads {threads} -summary {output.report} {input.r1_demux} {input.r2_demux} {output.p1} {output.u1} {output.p2} {output.u2} {params.trim_params} || true
        """
rule trim_stats:
    input:
        report = expand("{run}/trimmomatic/{sample}_trimmed/{sample}_trim_report.txt", run=config["job_name"], sample=config["samples"])
    output:
        directory("{run}/trim_summary/")
    params:
        summary = "{run}/trim_summary/summary_trim_stats.csv",
        graph_summary = "{run}/trim_summary/summary_trim_stats.png"
    conda:
        "envs/biopy.yaml"
    message:
        "Calculating trimming stats. Processing: {input.report}\n ==> Outputs to: {params.graph_summary}\n"
    script:
        "py/trim_stats.py"
rule fastqc_2:
    input:
        p1 = expand("{run}/trimmomatic/{sample}_trimmed/{sample}_dedup_P1.fastq.gz", run=config["job_name"], sample=config["samples"]),
        p2 = expand("{run}/trimmomatic/{sample}_trimmed/{sample}_dedup_P2.fastq.gz", run=config["job_name"], sample=config["samples"]),
        u1 = expand("{run}/trimmomatic/{sample}_trimmed/{sample}_dedup_U1.fastq.gz", run=config["job_name"], sample=config["samples"]),
        u2 = expand("{run}/trimmomatic/{sample}_trimmed/{sample}_dedup_U2.fastq.gz", run=config["job_name"], sample=config["samples"])
    output:
        directory("{run}/fastqc_postTrim/")
    conda:
        "envs/biopy.yaml"
    message:
        "Running FastQC on trimmed read pairs:\nR1: {input.p1}\nR2: {input.p2}"
    threads:
        config["threads"]
    shell:
        """
        mkdir {output}
        fastqc {input.p1} {input.p2} -t {threads} -o {output} || true
        """
rule multiqc_2:
    input:
        p1 = "{run}/fastqc_postTrim/"
    output:
        d1 = "{run}/multiqc_postTrim/multiqc_report.html"
    conda:
        "envs/biopy.yaml"
    params:
        folder = "{run}/multiqc_postTrim/"
    message:
        "Running MultiQC on FastQC reports:\nDirectory: {input.p1}\nOutputs to {output.d1}"
    shell:
        """
        multiqc --force --outdir {params.folder} --interactive {input.p1}
        """

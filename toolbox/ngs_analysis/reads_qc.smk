# **** Variables ****
configfile: "config/reads_qc.yaml"
# **** Imports ****
import glob

# Cluster run template
#nohup snakemake --snakefile reads_qc.smk -j 10 --cluster "sbatch -t {cluster.time} -n {cluster.cores} -N {cluster.nodes}" --cluster-config cluster.yaml --use-conda &

# noinspection SmkAvoidTabWhitespace
rule all:
    input:
        expand("{run}/clumpify/{sample}_dedup/{sample}_R1_dedup.fastq.gz",sample=config["samples"], run=config["run_name"]),
        #Reads QC Pre-Trim
        expand("{run}/fastqc_preTrim/", sample=config["samples"], run=config["run_name"]),
        expand("{run}/multiqc_preTrim/multiqc_report.html",run=config["run_name"]),
        #Reads Trimming and rRNA removal
        expand("{run}/trimmomatic/{sample}_trimmed/{sample}_dedup_P1.fastq.gz",sample=config["samples"], run=config["run_name"]),
        #Reads QC Post-Trim
        expand("{run}/fastqc_postTrim/",sample=config["samples"], run=config["run_name"]),
        expand("{run}/multiqc_postTrim/multiqc_report.html",run=config["run_name"]),
        expand("{run}/trim_summary/",run=config["run_name"])

rule clumpify:
    input:
        r1 = lambda wildcards: glob.glob("{directory}/{sample}_L001_R1_001.fastq.gz".format(directory=config["read_directory"],sample=wildcards.sample)),
        r2 = lambda wildcards: glob.glob("{directory}/{sample}_L001_R2_001.fastq.gz".format(directory=config["read_directory"],sample=wildcards.sample))
    output:
        r1_dedup = "{run}/clumpify/{sample}_dedup/{sample}_R1_dedup.fastq.gz",
        r2_dedup = "{run}/clumpify/{sample}_dedup/{sample}_R2_dedup.fastq.gz"
    conda:
        "envs/readqc.yaml"
    message:
        "Clumpifying read pairs:\nR1: {input.r1}\nR2: {input.r2}\n ==> Output:\n 1- {output.r1_dedup}\n 2- {output.r2_dedup}"
    threads:
        config["threads"]
    shell:
        """
        clumpify.sh in1={input.r1} in2={input.r2} out1={output.r1_dedup} out2={output.r2_dedup} dedupe=t optical=f
        """
rule fastqc_1:
    input:
        r1_dedup = expand("{run}/clumpify/{sample}_dedup/{sample}_R1_dedup.fastq.gz", run=config["run_name"], sample=config["samples"]),
        r2_dedup = expand("{run}/clumpify/{sample}_dedup/{sample}_R2_dedup.fastq.gz", run=config["run_name"], sample=config["samples"])
    output:
        directory("{run}/fastqc_preTrim/")
    conda:
        "envs/readqc.yaml"
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
        "envs/readqc.yaml"
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
        r1_dedup = "{run}/clumpify/{sample}_dedup/{sample}_R1_dedup.fastq.gz",
        r2_dedup = "{run}/clumpify/{sample}_dedup/{sample}_R2_dedup.fastq.gz"
    output:
        p1 = "{run}/trimmomatic/{sample}_trimmed/{sample}_dedup_P1.fastq.gz",
        u1 = "{run}/trimmomatic/{sample}_trimmed/{sample}_dedup_U1.fastq.gz",
        p2 = "{run}/trimmomatic/{sample}_trimmed/{sample}_dedup_P2.fastq.gz",
        u2 = "{run}/trimmomatic/{sample}_trimmed/{sample}_dedup_U2.fastq.gz",
        report = "{run}/trimmomatic/{sample}_trimmed/{sample}_trim_report.txt"
    params:
        trim_params = config["trim_params"]
    conda:
        "envs/readqc.yaml"
    message:
        """Trimming read pairs:\nR1: {input.r1_dedup}\nR2: {input.r2_dedup} 
        ==> Output:\n 1- {output.p1}\n 2- {output.p2} 
        ==> parameters: {params.trim_params}
        ==> Report:\n {output.report}"""
    threads:
        config["threads"]
    shell:
        """
        trimmomatic PE -threads {threads} -summary {output.report} {input.r1_dedup} {input.r2_dedup} {output.p1} {output.u1} {output.p2} {output.u2} {params.trim_params}
        """
# rule sortmerna:
# 	input:
# 		p1 = "{run}/trimmomatic/{sample}_trimmed/{sample}_dedup_P1.fastq.gz",
# 		p2 = "{run}/trimmomatic/{sample}_trimmed/{sample}_dedup_P2.fastq.gz",
# 		rrna_lsu = lambda wildcards: glob.glob("{rrna_db_dir}/SILVA_132_LSURef_tax_silva_trunc.fasta".format(
# 			rrna_db_dir=config["rrna_db_dir"]
# 		)),
# 		rrna_ssu = lambda wildcards: glob.glob("{rrna_db_dir}/SILVA_132_SSURef_tax_silva_trunc.fasta".format(
# 			rrna_db_dir=config["rrna_db_dir"]
# 		))
# 	output:
# 		rrna_clean_1 = "{run}/sortmerna/{sample}_rrna_clean/{sample}_R1_rrna_clean.fastq.gz",
# 		rrna_clean_2 = "{run}/sortmerna/{sample}_rrna_clean/{sample}_R2_rrna_clean.fastq.gz"
# 	conda:
# 		"envs/sortmerna.yaml"
# 	message:
# 		"Cleaning up rRNA from deduplicated sequences:\n Ref: {input.rrna_lsu}\nR1: {input.p1}\nR2: {input.p2}\n ==> Output:\n 1- {output.rrna_clean_1}\n 2- {output.rrna_clean_2}"
# 	threads:
# 		config["threads"]
# 	params:
# 		outdir = "{run}/sortmerna/{sample}_rrna_clean/",
# 		aligned_rrna = "{run}/sortmerna/{sample}_rrna_clean/out/aligned.fq.gz",
# 		rrna_clean_path_1 = "{run}/sortmerna/{sample}_rrna_clean/{sample}_R1_rrna_clean.fastq",
# 		rrna_clean_path_2 = "{run}/sortmerna/{sample}_rrna_clean/{sample}_R2_rrna_clean.fastq"
# 	shell:
# 		"""
# 		sortmerna --ref {input.rrna_ssu} --ref {input.rrna_lsu} --reads {input.p1} --reads {input.p2} --workdir {params.outdir} --threads {threads} --fastx
# 		python scripts/cleanup_rrna.py {params.aligned_rrna} {input.p1} {params.rrna_clean_path_1}
# 		python scripts/cleanup_rrna.py {params.aligned_rrna} {input.p2} {params.rrna_clean_path_2}
# 		pigz -9 -p{threads} {params.rrna_clean_path_1}
# 		pigz -9 -p{threads} {params.rrna_clean_path_2}
# 		rm -rf {params.outdir}/out {params.outdir}/idx {params.outdir}/readb
# 		"""
rule trim_stats:
    input:
        report = expand("{run}/trimmomatic/{sample}_trimmed/{sample}_trim_report.txt", run=config["run_name"], sample=config["samples"])
    output:
        directory("{run}/trim_summary/")
    params:
        summary = "{run}/trim_summary/summary_trim_stats.csv",
        graph_summary = "{run}/trim_summary/summary_trim_stats.png"
    conda:
        "envs/py.yaml"
    message:
        "Calculating trimming stats. Processing: {input.report}\n ==> Outputs to: {params.graph_summary}\n"
    script:
        "py/trim_stats.py"
rule fastqc_2:
    input:
        p1 = expand("{run}/trimmomatic/{sample}_trimmed/{sample}_dedup_P1.fastq.gz", run=config["run_name"], sample=config["samples"]),
        p2 = expand("{run}/trimmomatic/{sample}_trimmed/{sample}_dedup_P2.fastq.gz", run=config["run_name"], sample=config["samples"]),
        u1 = expand("{run}/trimmomatic/{sample}_trimmed/{sample}_dedup_U1.fastq.gz", run=config["run_name"], sample=config["samples"]),
        u2 = expand("{run}/trimmomatic/{sample}_trimmed/{sample}_dedup_U2.fastq.gz", run=config["run_name"], sample=config["samples"])
    output:
        directory("{run}/fastqc_postTrim/")
    conda:
        "envs/readqc.yaml"
    message:
        "Running FastQC on trimmed read pairs:\nR1: {input.p1}\nR2: {input.p2}"
    threads:
        config["threads"]
    shell:
        """
        mkdir {output}
        fastqc {input.p1} {input.p2} -t {threads} -o {output}
        """
rule multiqc_2:
    input:
        p1 = "{run}/fastqc_postTrim/"
    output:
        d1 = "{run}/multiqc_postTrim/multiqc_report.html"
    conda:
        "envs/readqc.yaml"
    params:
        folder = "{run}/multiqc_postTrim/"
    message:
        "Running MultiQC on FastQC reports:\nDirectory: {input.p1}\nOutputs to {output.d1}"
    shell:
        """
        multiqc --force --outdir {params.folder} --interactive {input.p1}
        """
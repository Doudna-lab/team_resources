# **** Variables ****
# configfile: "config/reads_qc-30-942700356_wc.yaml"

# **** Imports ****
import glob

# Cluster run template
#nohup snakemake --snakefile reads_mapping.smk -j 5 --cluster "sbatch -t {cluster.time} -n {cluster.cores} -N {cluster.nodes}" --cluster-config config/cluster.yaml --latency-wait 120 --use-conda &

# noinspection SmkAvoidTabWhitespace
rule all:
	input:
		expand("{run}/hisat2/{prefix}_index/{prefix}.1.ht2",prefix=config["reference_sequences"],run=config["job_name"]),
		expand("{run}/hisat2/{sample}/{sample}_report.txt", sample=config["reference"], run=config["job_name"]),
		# expand("{run}/hisat2/{sample}/{sample}_track_added.txt", sample=config["reference"], run=config["job_name"])


# noinspection SmkAvoidTabWhitespace
rule hisat_build:
	input:
		genome = lambda wildcards: glob.glob("{reference_path}/{prefix}.fasta".format(reference_path=config["reference_path"], prefix=wildcards.prefix)),
	output:
		index = "{run}/hisat2/{prefix}_index/{prefix}.1.ht2",
	conda:
		"envs/hisat.yaml"
	params:
		prefix_dir = "{run}/hisat2/{prefix}_index/",
	message:
		"""
Hisat2: Building genome index for:
{wildcards.prefix} 
=> file: {input.genome}
=> output: {output.index}
"""
	threads:
		config["threads"]
	shell:
		"""
		cd {params.prefix_dir}
		hisat2-build -p {threads} {input.genome} {wildcards.prefix}
		"""

# noinspection SmkAvoidTabWhitespace
rule hisat2_map:
	input:
		index = lambda wildcards: ("{run}/hisat2/{prefix}_index/{prefix}.1.ht2".format(run=config["job_name"],prefix=config["reference"][wildcards.sample])),
		p1 = "{run}/trimmomatic/{sample}_trimmed/{sample}_dedup_P1.fastq.gz",
		p2 = "{run}/trimmomatic/{sample}_trimmed/{sample}_dedup_P2.fastq.gz"
	output:		
		report = "{run}/hisat2/{sample}/{sample}_report.txt",
		paired_unmapped1 = "{run}/hisat2/{sample}/{sample}_unmapped-paired_1_.fq.gz",
		paired_unmapped2 = "{run}/hisat2/{sample}/{sample}_unmapped-paired_2_.fq.gz",
		bam_sort = "{run}/hisat2/{sample}/{sample}_mapping.sorted.bam",
		bam_index ="{run}/hisat2/{sample}/{sample}_mapping.sorted.bam.csi"
	conda:
		"envs/hisat.yaml"
	params:
		map = "{run}/hisat2/{sample}/{sample}_mapping.sam",
		idx = lambda wildcards: ("{run}/hisat2/{prefix}_index/{prefix}".format(prefix=config["reference"][wildcards.sample],run=config["job_name"])),
		mates = "{run}/hisat2/{sample}/{sample}_unmapped-paired_%_.fq.gz"
	message:
		"""
Mapping read pairs:
R1: {input.p1}
R2: {input.p2}
Reference:
{input.index}
Mapped Output:
{params.map}
        """
	threads:
		config["threads"]
	shell:
		"""
		hisat2 -p {threads} --no-unal --summary-file {output.report} --un-conc-gz {params.mates} -x {params.idx} -1 {input.p1} -2 {input.p2} -S {params.map}
		samtools sort -o {output.bam_sort} -O bam -l 9 -@ {threads} {params.map}
		rm {params.map}
		samtools index -c -@ {threads} {output.bam_sort} {output.bam_index}
		"""

# rule igv_view:
# 	input:
# 		genome = lambda wildcards: glob.glob("{reference_path}/{prefix}.fasta".format(reference_path=config["reference_path"], prefix=config["reference"][wildcards.sample])),
# 		bam_sort = "{run}/hisat2/{sample}/{sample}_mapping.sorted.bam"
# 	output:
# 		track_added_check = "{run}/hisat2/{sample}/{sample}_track_added.txt"
# 	conda:
# 		"envs/igv.yaml"
# 	params:
# 		igv_session = lambda wildcards: ("{run}/igv/{prefix}/{prefix}_view.igv".format(run=wildcards.run,
# 			prefix=config["reference"][wildcards.sample]))
# 	message:
# 		"""
# Generating IGV View:
# BAM: {input.bam_sort}
# Reference:
# {input.genome}
# IGV Session Output:
# {params.igv_session}
#         """
# 	shell:
# 		"""
# 		if [ -e "{params.igv_session}" ]; then
#         # Command to run if the file exists
#         igv.sh -load {params.igv_session} -tracks {input.bam_sort}
# 		else
# 	    # Command to run if the IGV SESSION does not exist
# 	    igv.sh -genome {input.genome} -tracks {input.bam_sort} -export {params.igv_session}
# 		fi
# 		 echo "Track {input.bam_sort} added to IGV session {params.igv_session}" > {output.track_added_check}
# 		"""

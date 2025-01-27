#!/bin/bash
#SBATCH --job-name=hisat_maxconvs
#SBATCH --nodes=1

# ACTIVATE ANACONDA
eval "$(conda shell.bash hook)"
conda activate biopympi
cd '/groups/doudna/team_resources/toolbox/hisat3n_tools/'

mpiexec python py/process_maxconv.py <input_BAM> <reference_fasta> <original_base> <converted_base> <maximum_conversions> <output_BAM> <log_file>
# mpiexec python py/process_maxconv.py /groups/doudna/projects/daniel_projects/hisat-3n_test_files/124-M_120_mapping.sorted.bam /groups/doudna/projects/daniel_projects/hisat-3n_test_files/ED_o_041.fasta 'G' 'A' 1 delete_test.bam report.log

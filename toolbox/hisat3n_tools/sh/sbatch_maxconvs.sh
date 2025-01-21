#!/bin/bash
#SBATCH --job-name=hisat_maxconvs
#SBATCH --nodes=1

# ACTIVATE ANACONDA
eval "$(conda shell.bash hook)"
conda activate biopympi

mpirun -n <n_of_workers> python -u py/process_maxconv.py <input_BAM> <reference_fasta> <original_base> <converted_base> <maximum_conversions> <output_BAM> <log_file>

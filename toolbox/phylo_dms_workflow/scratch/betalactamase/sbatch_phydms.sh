#!/bin/bash
#SBATCH --job-name=phydms_betalac
#SBATCH --nodes=2
#SBATCH --thread-spec=48
#SBATCH --mem=0
# ACTIVATE ANACONDA
eval "$(conda shell.bash hook)"
conda activate rax
cd /groups/doudna/team_resources/toolbox/phylo_dms_workflow/scratch/betalactamase
phydms_comprehensive --diffprefsbysite --omegabysite --ncpus 48 --raxml /home/danielbr/miniconda3/envs/rax/bin/raxmlHPC-PTHREADS-SSE3 betalac_phydms/ betaLactamase_alignment.fasta betaLactamase_prefs.txt

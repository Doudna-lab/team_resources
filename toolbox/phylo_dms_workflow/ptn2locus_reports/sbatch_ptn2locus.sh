#!/bin/bash
#SBATCH --job-name=ptn2locus
#SBATCH --cores=2
# ACTIVATE ANACONDA
eval "$(conda shell.bash hook)"
conda activate blast2region
cd /groups/doudna/team_resources/toolbox/phylo_dms_workflow/py
python3 ptn2locus.py /groups/doudna/team_resources/toolbox/phylo_dms_workflow/scratch/bacterialFormIIs.csv 'thedoudnalab@gmail.com' -o ../ptn2locus_reports/bacterialFormIIs.fna

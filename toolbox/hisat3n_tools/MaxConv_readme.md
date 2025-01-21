## Process Maximum Base Conversion Filter on Aligned Reads

### Installation

- Install the Anaconda environments from `envs`
    + biopympi.yaml is required to run `process_maxconv.py`

```
conda install -f envs/biopympi.yaml
```

### Run
- The script `py/process_maxconv.py` can filter reads on a BAM alignment file based on the base conversions.

```
mpirun -n <n_of_workers> python -u py/process_maxconv.py <input_BAM> <reference_fasta> <original_base> <converted_base> <maximum_conversions> <output_BAM> <log_file>
```

- Input Parameters:
    + <input_BAM> : File that will be read by the script
    + <reference_fasta> : FASTA file used to align the reads in the BAM input
    + <original_base> : Base in the reference sequence against which the possible converted bases will be compared to
    + <converted_base> : Base in the reads that will be compared to the relevant position in the reference
    + <maximum_conversions> : Number of <original_base> -> <converted_base> conversions at which the read will be discarded
    + <output_BAM> : Path where the resulting BAM file will be written
    + <log_file> : Path where a detailed log file will be created

### SBATCH RUN
- To make use of additional cores in a SLURM-based HPC when executing the script, running an SBATCH job is recommended
- The template file `sh/sbatch_maxconvs.sh` is provided as a reference
- After the relevant changes have been made to `sh/sbatch_maxconvs.sh`, it can be executed as follows

```
sbatch sbatch_maxconvs.sh
```

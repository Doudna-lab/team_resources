# Native modules
import subprocess
import re
import sys
import os


def main():
	# Snakemake Imports
	# # Inputs
	# msa_in_path = str(snakemake.input.post_search_msa)
	# # Outputs
	# report_out = str(snakemake.output.jpred_report)
	# msa_out_path = str(snakemake.params.nogaps_msa)
	# Inputs
	msa_in_path = str(sys.argv[1])
	# Outputs
	report_out = str(sys.argv[2])
	msa_out_path = str(sys.argv[3])

	# Change the working directory
	os.chdir(msa_out_path)

	jpred_submit_cmd = (f"jpredapi submit "
	                    f"mode=msa "
	                    f"format=fasta "
	                    f"file={msa_in_path}")
	# Run the subprocess and capture stdout
	completed_jpred = subprocess.run(['jpredapi submit', 'mode=msa', 'format=fasta', f'file={msa_in_path}'], shell=True, stdout=subprocess.PIPE)
	stdout = completed_jpred.stdout.decode('utf-8')
	print("STDOUT " + stdout)
	# Use regular expression to extract job IDs
	job_id = re.findall(r"jobid: (\w+)", stdout)

	# jpred_status_cmd = ("jpredapi status"
	#                     " jobid={} "
	#                     "getResults=yes "
	#                     "checkEvery=300 "
	#                     "silent").format(job_id)

	print("Submited alignment to jpred: " + msa_in_path)
	subprocess.run(['jpredapi status', f"jobid={job_id}", 'getResults=yes', 'checkEvery=300', 'silent'], shell=True)

	with open(report_out, 'w') as f:
		f.write(f"JPRED submission sucessful.\nResults assigned to {job_id}"
		        f"Dirctory: {msa_out_path}")


if __name__ == "__main__":
	main()

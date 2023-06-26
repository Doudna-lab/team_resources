# Native modules
import re
import os
from argparse import ArgumentParser as argp

# Load flowcell, sample name, read and lane naming patterns for filename recognition
# Flowcell ref:
# https://knowledge.illumina.com/instrumentation/general/instrumentation-general-reference_material-list/000005589
fastq_pattern_dict = {
	'flowcell': {"NextSeq 500/550 Mid Output": [re.compile("\w{4}AF\w{1,2}")],
	             "NextSeq 500/550 High Output": [re.compile("\w{4}BG\w{1,2}"), re.compile("\w{4}AG\w{1,2}")],
	             "NextSeq 1000/2000 P1": [re.compile("\w{6}M5")],
	             "NextSeq 1000/2000 P2": [re.compile("\w{6}M5")],
	             "NextSeq 2000 P3": [re.compile("\w{6}HV")],
	             "HiSeq 2500 Rapid v2": [re.compile("\w{4}BC\w{1,2}")],
	             "HiSeq 2500 TruSeq v3": [re.compile("\w{4}AC\w{1,2}")],
	             "HiSeq 2500 High Output v4": [re.compile("\w{4}AN\w{1,2}")],
	             "HiSeq 3000/4000": [re.compile("\w{4}BB\w{1,2}")],
	             "HiSeq X": [re.compile("\w{4}AL\w{1,2}"), re.compile("\w{4}CC\w{1,2}")],
	             "NovaSeq 6000 SP and S1": [re.compile("\w{4}DR\w{1,2}")],
	             "NovaSeq 6000 S2": [re.compile("\w{4}DM\w{1,2}")],
	             "NovaSeq 6000 S4": [re.compile("\w{4}DS\w{1,2}")]},
	"lane": [re.compile("L0{0,2}(\d+)")],
	"read": [re.compile("[RI][12]")],
	"sample_name": [re.compile("(IGU-\d{3})-N(\d+)")]
}


def parse_arguments():
	#  Launch argparse parser
	parser = argp(
		prog='fastq_refactor.py',
		description='',
		usage='%(prog)s [options] project_directory')
	# Define arguments
	parser.add_argument('project_directory',
	                    help='Path to the folder containing fastq files')  # positional argument
	parser.add_argument('-p',
	                    dest='platform',
	                    default='NovaSeq 6000 S2',
	                    choices=['NextSeq 500/550 Mid Output',
	                             'NextSeq 500/550 High Output',
	                             'NextSeq 1000/2000 P1',
	                             'NextSeq 1000/2000 P2',
	                             'NextSeq 2000 P3',
	                             'HiSeq 2500 Rapid v2',
	                             'HiSeq 2500 TruSeq v3',
	                             'HiSeq 2500 High Output v4',
	                             'HiSeq 3000/4000',
	                             'HiSeq X',
	                             'NovaSeq 6000 SP and S1',
	                             'NovaSeq 6000 S2',
	                             'NovaSeq 6000 S4'],
	                    help='Specify the sequencing platform. [Default: "NovaSeq 6000 S2"]')

	# Parse arguments from the command line
	arguments = parser.parse_args()
	return arguments


def loopNrename_files(directory, current_seq_platform, pattern_dict):
	flowcell_id = ''
	flowcell_bs = ''
	for filename in os.listdir(directory.rstrip(os.sep)):
		suffix = re.search(r'\.(.*)', filename).group(1)
		try:
			flowcell_id = pattern_dict['flowcell'][current_seq_platform][0].search(filename).group()
			flowcell_bs = '001'
		except AttributeError:
			raise "Could not find flowcell index based on the provided sequencing platform"
		# print(flowcell)
		try:
			sample_name = pattern_dict['sample_name'][0].search(filename).groups()[0]
			sample_number = f"S{pattern_dict['sample_name'][0].search(filename).groups()[1]}"
		except AttributeError:
			raise f"Could not find sample name in filename {filename}"
		# print(sample_name)
		# print(sample_number)
		try:
			read = pattern_dict['read'][0].search(filename).group()
		except AttributeError:
			raise f"Could not find read reference in the filename {filename}"
		# print(read)
		try:
			lane_number = int(pattern_dict['lane'][0].search(filename).groups()[0])
			lane = 'L'
			if lane_number > 100:
				lane = f"L{lane_number}"
			if isinstance(lane_number, int) and 10 <= lane_number < 100:
				lane = f"L0{lane_number}"
			if isinstance(lane_number, int) and 0 <= lane_number < 10:
				lane = f"L00{lane_number}"
		except AttributeError:
			raise f"Could not find lane reference in the filename {filename}"
		# print(lane)

		# Get the full path of the file
		file_path = os.path.join(directory, filename)

		# Check if the path is a file
		if os.path.isfile(file_path):
			# Generate the new filename
			new_filename = f"{sample_name}_{sample_number}_{lane}_{read}_{flowcell_bs}.{suffix}"
			# Rename the file
			os.rename(file_path, os.path.join(directory, new_filename))
			print(f"Renamed file: {filename} to {new_filename}")
	return flowcell_id, flowcell_bs


def main():
	# Call argument parsing function
	args = parse_arguments()
	sequencing_platform = args.platform
	project_directory = str(args.project_directory)
	flowcell_log_path = os.path.join(project_directory, "flowcell.log")

	(flowcell, truncated_flowcell) = loopNrename_files(project_directory, sequencing_platform, fastq_pattern_dict)

	with open(flowcell_log_path, 'w') as f:
		log_string = f"""
Based on the flowcell ID found in the original filenames the Sequencing platform used was: {sequencing_platform}
Flowcell ID: {flowcell}
To comply with BaseSpace naming convention, the flowcell ID was truncated to:  {truncated_flowcell}
		"""
		f.write(log_string)


if __name__ == "__main__":
	main()

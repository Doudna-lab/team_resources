# == Native Modules
import logging
import sys
# == Installed Modules
from mpi4py import MPI
import pysam
from pyfaidx import Fasta

# == Project Modules


def divide_work(N, n_rank, world_size):
	# determine the workload of each rank
	workloads = [N // world_size for i in range(world_size)]
	for i in range(N % world_size):
		workloads[i] += 1
	my_start = 0
	for i in range(n_rank):
		my_start += workloads[i]
	my_end = my_start + workloads[n_rank]
	return my_start, my_end


def recover_process(recover_buf, mode):
	combined_content = []
	for item in recover_buf:
		for value in item.values():
			combined_content.extend(value)
	if mode == 'list':
		combined_list = []
		for item in combined_content:
			combined_list.extend(item)
		return combined_list
	if mode == 'dict':
		combined_dict = {}
		for item in combined_content:
			combined_dict.update(item)
		return combined_dict


def get_bam_refs(bam_file_path):
	samfile = pysam.AlignmentFile(bam_file_path, "rb")
	return samfile, list(samfile.references)


def get_attr_string(sam_aln_record: pysam.AlignedSegment, read_id, reference_id):
	return str(f"{reference_id}|"
			   f"{read_id}|"
			   f"{sam_aln_record.query_alignment_sequence}")


def convert_sam2dict(sam_object: pysam.AlignedSegment):
	sam_dict = {
		'query_name': sam_object.query_name,
		'query_sequence': sam_object.query_sequence,
		'aligned_pairs': sam_object.aligned_pairs,
		'cigarstring': sam_object.cigarstring,
		'flag': sam_object.flag,
		'reference_id': sam_object.reference_id,
		'reference_start': sam_object.reference_start,
		'mapping_quality': sam_object.mapping_quality,
		'cigar': sam_object.cigar,
		'next_reference_id': sam_object.next_reference_id,
		'next_reference_start': sam_object.next_reference_start,
		'template_length': sam_object.template_length,
		'query_qualities': sam_object.query_qualities,
		'tags':  sam_object.tags,
	}
	return sam_dict


def convert_dict2sam(sam_dict: dict):
	sam_aln_record = pysam.AlignedSegment()
	sam_aln_record.query_name = sam_dict['query_name']
	sam_aln_record.query_sequence = sam_dict['query_sequence']
	sam_aln_record.flag = sam_dict['flag']
	sam_aln_record.reference_id = sam_dict['reference_id']
	sam_aln_record.reference_start = sam_dict['reference_start']
	sam_aln_record.mapping_quality = sam_dict['mapping_quality']
	sam_aln_record.cigar = sam_dict['cigar']
	sam_aln_record.next_reference_id = sam_dict['next_reference_id']
	sam_aln_record.next_reference_start = sam_dict['next_reference_start']
	sam_aln_record.template_length = sam_dict['template_length']
	sam_aln_record.query_qualities = sam_dict['query_qualities']
	sam_aln_record.tags = sam_dict['tags']
	return sam_aln_record


def consolidate_tasks(task_list):
	"""
	 Consolidates tasks into a dictionary grouping by the first part of each item.
	 Args:
	     task_list (list): A list of strings, each in the format 'key|value'.
	 Returns:
	     dict: A dictionary where keys are the first part of the input strings and
	           values are lists of the second parts.
	 """
	task_dict = {}
	for item in task_list:
		key, value, sequence = item.split('|')
		data_pack = (value, sequence)
		task_dict.setdefault(key, []).append(data_pack)
	return dict(task_dict)


def clean_aligned_pairs(aligned_pairs):
	parsed_pairs = []
	ref_aln_lookup = []

	for query_aln_column, ref_aln_column in aligned_pairs:
		# == Skip Aligned Pairs with gaps
		if ref_aln_column is None or query_aln_column is None:
			continue
		ref_aln_lookup.append(int(ref_aln_column))
		parsed_pairs.append((int(query_aln_column), int(ref_aln_column)))
	return min(ref_aln_lookup), max(ref_aln_lookup), parsed_pairs


def quantify_conversions(ref_sequence, query_sequence, aligned_pairs, ref_sequence_offset, ref_base, alt_base, max_conversions):
	count_conversions = 0
	loop_count = 0
	max_conversions = int(max_conversions)
	query_aln_offset = None
	for query_aln_column, ref_aln_column in aligned_pairs:
		if query_aln_offset is None:
			query_aln_offset = query_aln_column
		loop_count += 1
		current_ref_base = None
		query_aln_column = int(query_aln_column) - query_aln_offset
		ref_aln_column = int(ref_aln_column) - ref_sequence_offset

		# print(f"###### LOOP NUMBER {loop_count} ########")

		try:
			current_ref_base = ref_sequence[ref_aln_column]
			# print(f"##REF - Compare {current_ref_base} on col {ref_aln_column}##")
			# print(f"Ref base: {current_ref_base}")
		except IndexError:
			# print(f"Ref Sequence {ref_sequence} of length {len(ref_sequence)} does not have {ref_aln_column}th column")
			# print(f"Aligned pairs {aligned_pairs}")
			continue

		if current_ref_base == ref_base:
			try:
				current_query_base = query_sequence[query_aln_column]
				if current_query_base == current_ref_base:
					continue
				# print(f"##QUERY - Compare {current_query_base} on col {query_aln_column}##")
			except IndexError:
				# print(f"Q Sequence {query_sequence} of length {len(query_sequence)} does not have {query_aln_column}th column")
				# return None, None
				continue

			if current_query_base == alt_base:
				count_conversions += 1
				if count_conversions > max_conversions:
					return False, count_conversions
	return True, count_conversions


def export_bam(bam_records_list, bam_header, output_bam):
	with pysam.AlignmentFile(output_bam, "wb", header=bam_header) as output_bam:
		for read in bam_records_list:
			output_bam.write(read)


def main():
	#	Inputs
	INPUT_BAM_PATH = sys.argv[1]
	INPUT_FASTA_PATH = sys.argv[2]
	REF_BASE = sys.argv[3] #'G'
	ALT_BASE = sys.argv[4] #'A'
	MAX_CONVERSIONS = sys.argv[5] # 2
	OUTPUT_BAM_PATH = sys.argv[6]
	LOG_FILE_PATH = sys.argv[7]
	#	Outputs
	#	Params

	# Configure the logging system
	logging.basicConfig(
		level=logging.DEBUG,  # Set the minimum log level (DEBUG logs everything)
		format="%(asctime)s %(message)s",  # Define log format
		handlers=[
			logging.FileHandler(LOG_FILE_PATH),  # Log to a file
		]
	)

	# === get basic information about the MPI communicator
	comm = MPI.COMM_WORLD
	world_size = comm.Get_size()
	my_rank = comm.Get_rank()
	divided_work = []
	# === Set processing variables
	process_approved_reads = {}
	process_rejected_reads = {}
	sam_alignment_dictionary = {}
	uniq_job_reference_list = []
	original_header = ''
	reference_fasta = None

	# # DEBUG BLOCK
	# my_read = "A01587:418:GW230625000:1:2119:11912:19836"
	# INPUT_BAM_PATH = "/groups/doudna/projects/daniel_projects/hisat-3n_test_files/124-M_120_mapping.sorted.bam"
	# INPUT_FASTA_PATH = "/groups/doudna/projects/daniel_projects/hisat-3n_test_files/ED_o_041.fasta"
	# OUTPUT_BAM_PATH = "delete.bam"
	# LOG_FILE_PATH = ""
	# REF_BASE = 'G'
	# ALT_BASE = 'A'
	# MAX_CONVERSIONS = 2

	if my_rank == 0:
		logging.info(f"*** Loaded Input alignment file from:\n {INPUT_BAM_PATH}")
		# This block sets up the length of the input to be distributed across workers
		# 	the format of the list items is as follows: "REF_SEQUENCE|READ_ID"
		logging.info("*** Processing input list of reads")
		(sam_alignment, reference_sequences) = get_bam_refs(INPUT_BAM_PATH)
		adjacent_read = []
		for reference in reference_sequences:
			for read in sam_alignment.fetch(reference):

				# == Process pair mates ==
				read_id = read.query_name + '/1'
				if read.query_name in set(adjacent_read):
					read_id = read.query_name + '/2'
					adjacent_read.remove(read.query_name)
				adjacent_read.append(read.query_name)

				# == Convert 'AlignedSegment' object to dictionary
				# 		-> Required for serialization through broadcasting
				read_dict = convert_sam2dict(read)
				original_header = sam_alignment.header.to_dict()

				# == Assemble variables for future reference ==
				attributes_str = get_attr_string(read, read_id, reference)
				sam_alignment_dictionary.setdefault(read_id, read_dict)
				uniq_job_reference_list.append(attributes_str)

		sam_alignment.close()
	else:
		uniq_job_reference_list = None
		original_header = None
		sam_alignment_dictionary = None

	# The reference FASTA is imported on rank 0, then broadcasted to other workers
	reference_fasta = Fasta(INPUT_FASTA_PATH)

	# === BROADCAST VARIABLES ===
	broadcast_input_list = MPI.COMM_WORLD.bcast(uniq_job_reference_list, root=0)
	# broadcast_ref_fasta = MPI.COMM_WORLD.bcast(reference_fasta, root=0)
	broadcast_sam_dict = MPI.COMM_WORLD.bcast(sam_alignment_dictionary, root=0)
	broadcast_sam_header = MPI.COMM_WORLD.bcast(original_header, root=0)

	# === DELIMIT THE LENGTH MAIN TASK TO BE SPLIT
	N = len(broadcast_input_list)

	# === CREATE LIST OF TUPLES CONTAINING START AND END INDICES TO BE PROCESSED BY EACH RANK
	if my_rank == 0:
		# == Divide work across ranks
		for rank in range(world_size):
			divided_work.append(divide_work(N, rank, world_size))
		logging.info(f"*** Divided workload across {world_size} workers")

	# === SCATTER WORKLOAD ACROSS THE RANKS
	v = comm.scatter(divided_work, root=0, )

	# === SLICE MAIN LIST ACCORDING TO THE DIVIDED WORKLOAD
	task_per_rank_list = broadcast_input_list[v[0]:v[1]]

	# === PERFORM THE CORE CODE WITHIN THE LOOP
	logging.info(f"*** Rank: {my_rank}: Perform maximum conversion filter")
	logging.info(f"*** Rank: {my_rank}: {v[0]} -> {v[1]}")

	# # DEBUG BLOCK
	# broadcast_ref_fasta = reference_fasta
	# broadcast_sam_header = original_header
	# broadcast_sam_dict = sam_alignment_dictionary
	# task_per_rank_list = uniq_job_reference_list

	# == CORE CODE ==
	consolidated_tasks_dict = consolidate_tasks(task_per_rank_list)
	approved_reads = []
	rejected_reads = {}
	for reference, read_attributes in consolidated_tasks_dict.items():
		for attr_record in read_attributes:
			read_id = attr_record[0]
			read_seq = attr_record[1]
			aligned_pairs = broadcast_sam_dict[read_id]['aligned_pairs']
			read_aln_start, read_aln_end, aligned_pairs = clean_aligned_pairs(aligned_pairs)
			aligned_reference = reference_fasta[reference][read_aln_start:read_aln_end]
			read_approval, conversions = quantify_conversions(aligned_reference,
															  read_seq,
															  aligned_pairs,
															  read_aln_start,
															  REF_BASE,
															  ALT_BASE,
															  MAX_CONVERSIONS)
			# Parse Approved and Rejected reads
			if read_approval:
				approved_reads.append(read_id)
			elif not read_approval:
				rejected_reads.setdefault(read_id, conversions)

	# == APPEND THE PROCESSED ITEMS TO LIST WITHIN RANK DICTIONARY
	process_approved_reads.setdefault(my_rank, []).append(approved_reads)
	process_rejected_reads.setdefault(my_rank, []).append(rejected_reads)

	# === RECOVER PROCESSED CONTENTS FROM RANKS
	recvbuf_approved_reads = comm.gather(process_approved_reads, root=0)
	recvbuf_rejected_reads = comm.gather(process_rejected_reads, root=0)

	# === COMBINE CONTENT ON ROOT RANK
	if comm.rank == 0:
		decoded_approved_reads = []
		combined_approved_reads = recover_process(recvbuf_approved_reads, 'list')
		combined_rejected_reads = recover_process(recvbuf_rejected_reads, 'dict')

		for read_id in combined_approved_reads:
			decoded_approved_reads.append(convert_dict2sam(broadcast_sam_dict[read_id]))


		logging.info(f"*** EXPORT FILES")
		# === PROCESS CORE DATA ==='
		export_bam(decoded_approved_reads, broadcast_sam_header, OUTPUT_BAM_PATH)

		# Logging tally
		logging.info(f"""
				--> TALLYING RESULTS:
				*** Total reads APPROVED at maximum conversions {MAX_CONVERSIONS} of type {REF_BASE}->{ALT_BASE} : {len(combined_approved_reads)}
				*** Total reads REJECTED at maximum conversions {MAX_CONVERSIONS} of type {REF_BASE}->{ALT_BASE} : {len(combined_rejected_reads)}
				""")


if __name__ == "__main__":
	main()

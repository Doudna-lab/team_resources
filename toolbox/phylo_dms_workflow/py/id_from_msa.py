from Bio import AlignIO
from sys import argv


def main():
	seq_path = argv[1]

	seq_path = "/Users/bellieny/projects/team_resources/toolbox/phylo_dms_workflow/protein_gym/ProteinGym_Files/MSA_files/TRPC_SACS2_full_11-26-2021_b07.a2m"

	with open(seq_path, 'r') as f:
		records = AlignIO.read(f, "fasta")
	for rec in records:
		id = rec.id.split("_")[1]
		id = id.split("/")[0]

		print(id)

if __name__ == "__main__":
	main()

import os
import hashlib
import sys


def compare_hashes_in_directory(directory):
    for filename in os.listdir(directory):
        file_path = os.path.join(directory, filename)
        if os.path.isfile(file_path) and filename.endswith(".fastq.gz"):
            md5_hash = hashlib.md5()
            with open(file_path, "rb") as file:
                # Read the file in chunks to avoid memory issues with large files
                for chunk in iter(lambda: file.read(4096), b""):
                    md5_hash.update(chunk)

            md5_hash_value = md5_hash.hexdigest()
            expected_hash_file_path = file_path + ".md5"  # Assume hash value is stored in a .txt file

            with open(expected_hash_file_path, "r") as hash_file:
                expected_hash = hash_file.readline().strip()

            if md5_hash_value == expected_hash:
                print(f"Hash match for file: {filename}")
            else:
                print(f"Hash mismatch for file: {filename}")


def main():
	# Provide the directory path to compare hashes
	directory_path = sys.argv[1]
	compare_hashes_in_directory(directory_path)


if __name__ == "__main__":
	main()

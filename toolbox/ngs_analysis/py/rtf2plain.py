import sys
from striprtf.striprtf import rtf_to_text


def convert_rtf_to_text(input_file, output_file):
	with open(input_file, 'r') as f:
		rtf_content = f.read()

	plain_text = rtf_to_text(rtf_content)

	with open(output_file, 'w') as f:
		f.write(plain_text)


def main():
	if len(sys.argv) != 3:
		print("Usage: python script.py input.rtf output.txt")
	else:
		input_rtf_file = sys.argv[1]
		output_txt_file = sys.argv[2]
		convert_rtf_to_text(input_rtf_file, output_txt_file)
		print(f"Conversion completed. Plain text saved to {output_txt_file}")


if __name__ == "__main__":
	main()

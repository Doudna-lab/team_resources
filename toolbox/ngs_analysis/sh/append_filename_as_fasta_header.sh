#!/bin/bash

# Loop through all files starting with 'ED_o_'[CHANGE TO RELEVANT PATTERN] in the current directory
for file in ED_o_*.fasta; do
    # Check if the file is a regular file (not a directory)
    if [ -f "$file" ]; then
        # Remove the '.fasta' extension from the filename
        base_filename="${file%.fasta}"
        # Add the line with the modified filename at the beginning of the file
        sed -i "1i\\>${base_filename}" "$file"
    fi
done

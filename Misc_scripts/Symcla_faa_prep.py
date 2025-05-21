#!/usr/bin/env python3
"""
Reads all .faa files in the specified folder, and copies them to new .faa files 
following the header format required by symcla:
">genomeid|proteinid" where genome id is the same as the prefix of the .faa file

Note that much of the protein header is lost. Only the first "word" is preserved as the 
proteinid. In prokka output this is a unique protein identifier. 

Output will be new fasta files written into the input directory and named as follows:
Prefix.faa -> Prefix_symcla_input.faa


Usage: (python) Symcla_faa_prep.py input_dir
"""
import sys
import os

#Generator function that reads file line by line, without reading the full file into memory
def gen_line_reader(file_path):
    for line in open(file_path, "r"):
        yield line

faa_dir = sys.argv[1]
for subdir, dirs, files in os.walk(faa_dir): #iterate over all directories and subdirectories
    for file in files: #iterate over all files
        if file.endswith(".faa"):
            prefix = file.strip(".faa")
            out = open(os.path.join(subdir,prefix + "_symcla_input.faa"), "w")
            for line in gen_line_reader(os.path.join(subdir,file)):
                if line.startswith(">"):
                    out.write(">{}|{}\n".format(prefix + "_symcla_input", line.split()[0].strip(">")))
                else:
                    out.write(line)
                    

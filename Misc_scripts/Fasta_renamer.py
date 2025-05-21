#!/usr/bin/env python3
"""
Renames all fasta files in a given directory to match the name given in their first sequence name.
Designed for renaming fasta files downloaded from NCBI, changing their REFSEQ-based
to instead reflect the organism name.
(fasta files are recognized based on their extension (default=.fna), which is hardcoded at the bottom.
Please change this if you have a different use case)


Usage: (python) Fasta_renamer.py target_directory
"""
###################
#Import statements
###################  

import sys
import os

###################
#File handling functions
###################  
def gen_line_reader(file_path):
    for line in open(file_path, "r"):
        yield line

def new_name_finder(filename):
    """
    """
    for line in gen_line_reader(filename):
        if line.startswith(">"):
            #remove the code preceding the name
            line = line.partition(" ")[2]
            #Remove the MAG: term
            line = line.replace("MAG: ", "")
            #remove everthing after a comma
            line = line.partition(",")[0]
            #remove the chromosome term, or anything after
            line = line.partition("chromosome")[0]
            #Remove the isolate names
            line = line.partition("isolate")[0]
            #Remove "genome assembly"
            line = line.partition("genome assembly")[0]
            #replace spaces with _
            line = line.replace(" ", "_")
            #remove periods
            line = line.replace(".", "")
            line = line.strip("_")
            new_filename = filename.split("_")[0]+ "_"+filename.split("_")[1]
            new_filename = new_filename + "_" + line + ".fna"
            break
    return new_filename

input_dir = sys.argv[1]
file_extension = ".fna"
for root, dirs, files in os.walk(input_dir): 
        for file in files: 
            if file.endswith(file_extension):
                os.rename(file, new_name_finder(file))
            

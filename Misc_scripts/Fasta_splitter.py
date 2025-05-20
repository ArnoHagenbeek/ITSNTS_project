#!/usr/bin/env python3
"""
Reads a multifasta file, and separates each sequence (split by ">") into a separate fasta file.
Each fasta file will be named to match the sequence name after ">", and written to the working
directory. 
(Note, characters that interfere with commandline use will be removed from the sequence name)

Usage: (python) Fasta_splitter.py filename.fasta
"""
import sys
import os

#Generator function that reads file line by line, without reading the full file into memory
def gen_line_reader(file_path):
    for line in open(file_path, "r"):
        yield line

input_fasta = sys.argv[1]
    
outfile = False
for line in gen_line_reader(input_fasta): #iterate over each line of the file individually
    if line.startswith(">"): #if the line is a sequence name open a new fasta file for outputting
        #Ensure that all characters of the new filename can be read in a commandline
        seqname = line.strip(">").strip().replace(" ","_").replace(":","").replace("[","").replace("]","").replace("(","").replace(")","").replace(".","_").replace(",","_")
        if outfile:
            outfile.close()
        outfile = open(seqname +".fasta", "w")
    if outfile:
        outfile.write(line)
        

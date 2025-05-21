#!/usr/bin/env python3
"""
Parses all the EggNOG "emapper.annotations" files in the given directory, and counts the COG terms per sample.
Writes the COG counts per sample to a single tab-delimited file

Usage: (python) EggNOG_COG_counter.py Input_dir Output_file
"""
import sys
import os

#Generator function that reads file line by line, without reading the full file into memory
def gen_line_reader(file_path):
    for line in open(file_path, "r"):
        yield line

eggnog_dir = sys.argv[1]
outputfile = sys.argv[2]
COG_dir = {}
COG_list = []
for subdir, dirs, files in os.walk(eggnog_dir): #iterate over all directories and subdirectories
    for file in files: #iterate over all files
        if file.endswith("emapper.annotations"):
            sample = file.partition(".")[0]
            COG_dir[sample] = {}
            for line in gen_line_reader(os.path.join(subdir,file)):
                if line.startswith("#"):
                    continue
                splt = line.split("\t")
                if len(splt) > 6:
                    if splt[6] != "":
                        COG = splt[6]
                        if len(COG.strip()) > 1:
                            COGlst = COG.strip()
                            for COG in COGlst:
                                if COG not in COG_list:
                                    COG_list.append(COG)
                                if COG in COG_dir[sample]:
                                    COG_dir[sample][COG] += 1
                                else:
                                    COG_dir[sample][COG] = 1
                        else:
                            if COG not in COG_list:
                                COG_list.append(COG)
                            if COG in COG_dir[sample]:
                                COG_dir[sample][COG] += 1
                            else:
                                COG_dir[sample][COG] = 1

out = open(outputfile, "w")
out.write("Sample\t")
for COG in COG_list:
    out.write(COG + "\t")
out.write("\n")
for sample in COG_dir:
    out.write(sample + "\t")
    for COG in COG_list:
        if COG in COG_dir[sample]:
            out.write(str(COG_dir[sample][COG]) + "\t")
        else:
            out.write("0\t")
    out.write("\n")

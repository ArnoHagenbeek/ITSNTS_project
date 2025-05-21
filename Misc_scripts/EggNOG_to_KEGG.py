#!/usr/bin/env python3
"""
Reads all the EggNOG .annotations file and extracts the KEGG ko number. Outputs a file similar to the downloaded
outputs from BLASTkoala. 

Usage: (python) EggNOG_to_KEGG.py input.emapper.annotations output_name (optional)
"""
import sys
import os

#Generator function that reads file line by line, without reading the full file into memory
def gen_line_reader(file_path):
    for line in open(file_path, "r"):
        yield line

#Obtain inputs
file1 = sys.argv[1]
file2 = sys.argv[2]

ko_list1 = []
ko_list2 = []

not_list1 = []
not_list2 = []
#Read over EggNOG file, and write the KEGG KO terms per gene to the output
for line in gen_line_reader(file1):
    if line.startswith("#"):
        continue
    if len(line.split("\t")) > 1:
        gene = line.split("\t")[0]
        ko = line.split("\t")[1].strip()
        if ko.strip():
            if ko.strip() not in ko_list1:
                ko_list1.append(ko.strip())

for line in gen_line_reader(file2):
    if line.startswith("#"):
        continue
    if len(line.split("\t")) > 1:
        gene = line.split("\t")[0]
        ko = line.split("\t")[1].strip()
        if ko.strip():
            if ko.strip() not in ko_list2:
                ko_list2.append(ko.strip())  
            


for item in ko_list1:
    if item not in ko_list2:
        not_list1.append(item)
for item in ko_list2:
    if item not in ko_list1:
        not_list2.append(item)   
counter = 0
for item in ko_list2:
    if item in ko_list1:
        counter += 1 
print(len(ko_list1))  
print(len(ko_list2))  
print(counter)
print(len(not_list1))        
print(len(not_list2))   

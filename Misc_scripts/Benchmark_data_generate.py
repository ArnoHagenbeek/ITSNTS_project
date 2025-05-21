#!/usr/bin/env python3
"""
Author: Arno Hagenbeek
Generates benchmarking data to test tools. Generates random reads or contigs from the provided genomes, with known quantities.

Usage: (python) Benchmark_data_generate.py --fastq --fasta -i genomes -o output_prefix
    
    
Arguments:
    -h: Calls help message
    -i: Should be followed by a comma-delimited list of genome files in fasta format (or a single fasta file).
    -o: Should be followed by the desired output prefix. Output .fastq, .fasta and .log files will have the specified name
    --fastq: Activates read data generating mode
    --fasta: Activates contig data generating mode
    --ratio: Comma-delimited list of ratios for read data. e.g. for 1,2,3 the first genome will account for 1/6th of the 
             data whereas the third will account for half. Default is equal ratios.
    --read_length: Desired read length for read generating mode.
    --unpaired: Argument for activating unpaired read mode, by default generate paired reads.
    --inner_distance: Distance in bp between paired reads. Default is 50bp.
    --yield: The total sequencing yield in megabases for the generated read data. e.g. --yield 5 will generate read data 
             equating to 5 megabases. Default is 10.
    --completion: Comma-delimited list of the desired completion values for contig generating mode. e.g. --completion 0.5 
                  will result in contigs accounting for half of the provided genome (with the rest being gaps. Default is 
                  0.8 for each genome
    --contigs: Comma-delimited list of the desired number of contigs per genome
    --shuffle: Randomly orders the output fasta or fastq. Otherwise output will be ordered per genome
"""
import sys
import os
import random
from Bio.Seq import Seq

#Generator function that reads file line by line, without reading the full file into memory
def gen_line_reader(file_path):
    for line in open(file_path, "r"):
        yield line


#Read and parse input arguments
################################################
inp = sys.argv
genomes = []
read_mode = False
assemb_mode = False
shuffle_output = False
#Universal inputs
if "-o" in inp:
    out = (inp[inp.index("-o") + 1])
else:
    print("No output prefix specified")
if "-i" in inp:
    genomes = (inp[inp.index("-i") + 1])
    genomes = genomes.split(",")
else:
    print("No input genomes specified")
if "--shuffle" in inp:
    shuffle_output = True
#Read mode inputs
if "--fastq" in inp:
    read_mode = True
    read_length = 150
    paired = True
    pair_inner_distance = 50
    reads_yield = 10
    if "--ratio" in inp:
        ratios = (inp[inp.index("--ratio") + 1])
        ratios = ratios.split(",")
        ratios = [int(i) for i in ratios]
    else:
        ratios = []
        for x in range(len(genomes)):
            ratios.append(1)
    if "--read_length" in inp:
        read_length = int(inp[inp.index("-read_length") + 1])
    if "--unpaired" in inp:
        paired = False
    if "--inner_distance" in inp:
        pair_inner_distance = int(inp[inp.index("-inner_distance") + 1])
    if "--yield" in inp:
        reads_yield = float(inp[inp.index("--yield") + 1])

if "--fasta" in inp:
    assemb_mode = True
    if "--completion" in inp:
        compl = (inp[inp.index("--completion") + 1])
        compl = compl.split(",")
        compl = [float(i) for i in compl]
    else:
        compl = []
        for x in range(len(genomes)):
            compl.append(0.8)
    if "--contigs" in inp:
        cont = (inp[inp.index("--contigs") + 1])
        cont = cont.split(",")
        cont = [int(i) for i in cont]
    else:
        cont = []
        for x in range(len(genomes)):
            cont.append(300)
#Help message. Suppress the rest of the script
if "-h" in inp:
    read_mode = False
    assemb_mode = False
    print("Generates benchmarking data to test tools. Generates random reads or contigs from the provided genomes, with known quantities.")
    print("\nUsage: (python) Test_data_generate.py --fastq --fasta -i genomes -o output_prefix")
    print("\n\nArguments:\n-h: Calls this help message")
    print("-i: Should be followed by a comma-delimited list of genome files in fasta format (or a single fasta file)")
    print("-o: Should be followed by the desired output prefix. Output .fastq, .fasta and .log files will have the specified name")
    print("--fastq: Activates read data generating mode")
    print("--fasta: Activates contig data generating mode")
    print("--ratio: Comma-delimited list of ratios for read data. e.g. for 1,2,3 the first genome will account for \n\t 1/6th of the data whereas the third will account for half. Default is equal ratios")
    print("--read_length: Desired read length for read generating mode")
    print("--unpaired: Argument for activating unpaired read mode, by default generate paired reads")
    print("--inner_distance: Distance in bp between paired reads. Default is 50bp")
    print('--yield: The total "sequencing yield" in megabases for the generated read data. e.g. "--yield 5" will generate\n\t read data equating to 5 megabases. Default is 10.')
    print('--completion: Comma-delimited list of the desired completion values for contig generating mode. e.g. \n\t      "--completion 0.5" will result in contigs accounting for half of the provided genome \n\t      (with the rest being gaps). Default is 0.8 for each genome')
    print("--contigs: Comma-delimited list of the desired number of contigs per genome")
    print("--shuffle: Randomly orders the output fasta or fastq. Otherwise output will be ordered per genome")
    
if read_mode:
    log_dict_reads = {}
    read_lst = []
    #To manage memory use, the script will read genome stretches of the size specified here.
    #Note that while a higher number yields slightly better results (being more random), it is also far slower.
    max_bp = 100000 
    #Divide the total yield based on the provided coverage ratios
    x = reads_yield/sum(ratios)
    genome_reads = []
    for i,genome in enumerate(genomes):
        print("Generating reads for {}".format(genome))
        genome_gen = (row for row in open(genome))
        log_dict_reads[genome] = []
        expected_reads = (ratios[i] * x *1000000)/read_length #calculate total yield for this genome in number of reads
        current_read_count = 0
        total_bp = 0
        while current_read_count < expected_reads: #Continue reading genome stretches and building reads until the desired number of reads is achieved
            #Read in a genome stretches
            stretch_len = 0
            stretch = []
            #Read in a stretch of the genome
            while stretch_len < max_bp:                
                try: #Test if iterator is not ended
                    line = next(genome_gen)
                except StopIteration: #If iterator ended, restart it.
                    genome_gen = (row for row in open(genome))
                    break
                if line.startswith(">"): #If encounter a new sequence while building a stretch, immediately break stop building stretch, as reads cannot bridge sequences
                    if stretch != []:
                        break
                    else:
                        continue
                stretch_len += len(line.strip())
                stretch.append(line.strip())
            stretch = "".join(stretch)
            #Build random reads from the stretch equating the number of calculated reads per stretch
            for iterat in range(random.randint(1,40)):
                #Generate a random read start. 
                #If the read start is towards the end, instead make it a reverse read start
                read_index = random.randint(0, len(stretch))
                if paired:
                    if read_index < len(stretch) - read_length:
                        read = stretch[read_index:(read_index + read_length)]
                        reverse = Seq(stretch[read_index+pair_inner_distance:(read_index + read_length+pair_inner_distance)]).reverse_complement()
                    else:
                        reverse = Seq(stretch[read_index:]).reverse_complement()
                        read = stretch[read_index - pair_inner_distance - read_length:read_index - pair_inner_distance] 
                    if len(read) > read_length/2 and len(reverse) > read_length/2:
                        total_bp += (len(read) + len(reverse)) 
                        current_read_count+= 2
                        genome_reads.append((genome,read,str(reverse)))#Output paired reads as a tuple
                else:
                    if read_index < len(stretch) - read_length:
                        read = stretch[read_index:(read_index + read_length)]
                        if random.randint(0, 1) > 0: #Randomly determine if read is reverse or forward
                            read = Seq(read).reverse_complement()
                    else:
                        read = stretch[read_index:len(stretch)]   
                        if random.randint(0, 1) > 0: #Randomly determine if read is reverse or forward
                            read = Seq(read).reverse_complement()  
                    if len(read) > read_length/2 :
                        total_bp += len(read)
                        current_read_count+= 1
                        genome_reads.append((genome,read))#Output paired reads as a tuple
        log_dict_reads[genome].append(current_read_count)
        log_dict_reads[genome].append(total_bp)            
    if shuffle_output: 
        random.shuffle(genome_reads)
        
if assemb_mode:
    log_dict_fasta = {}
    cont_lst = []
    #To manage memory use, the script will read genome stretches of the size specified here.
    max_bp = 1000000
    for i,genome in enumerate(genomes):
        print("Generating contigs for {}".format(genome))
        log_dict_fasta[genome] = []
        #Estimate total size without reading into memory
        for count, line in enumerate(gen_line_reader(genome)):
            pass
        for line in gen_line_reader(genome):
            if not line.startswith(">"):
                line_length = len(line)
                break
        approx_size = line_length*count
        desired_size = approx_size * compl[i] #multiply size with expected completion for desired assembly size
        #Make random contig sizes based on average contig size
        cont_size = int(desired_size / cont[i]) #Calculate average contig size based on desired assembly size and total contig count
        cont_lengths = []
        for x in range(cont[i]):
            cont_lengths.append(cont_size)
        #Iterate over the length list. Making each entry shorter by a random size, and adding this size to a random other contig. 
        for index in range(len(cont_lengths)):
            min_size = cont_lengths[index] / random.randint(2, 10)
            cont_lengths[index] = int(cont_lengths[index] - min_size) #remove random stretch of bp
            rnd = random.randint(0, len(cont_lengths)-1)
            cont_lengths[rnd] = int(cont_lengths[rnd] + min_size) #add same stretch to other random contig
        #Follow the same random process for the gaps
        desired_size = approx_size * (1 - compl[i])
        gap_size = int(desired_size / (cont[i]*10)) #make number of gaps 10X the number of contigs (as it is fine if multiple gaps follow eachother often, but less if multiple contigs follow each other)
        gap_lengths = []
        for x in range(cont[i]*10):
            gap_lengths.append(gap_size)
        #Iterate over the length list. Making each entry shorter by a random size, and adding this size to a random other contig. 
        for index in range(len(gap_lengths)):
            min_size = gap_lengths[i] / random.randint(2, 10)
            gap_lengths[i] = int(gap_lengths[i] - min_size) #remove random stretch of bp
            gap_lengths[random.randint(0, len(gap_lengths)-1)] = int(gap_lengths[random.randint(0, len(gap_lengths)-1)] + min_size) #add same stretch to other random contig
        #Total gaps and contigs lengths should roughly equate to total genome. 
        #Iterate over the genome, and assign lines to contigs and gaps      
        current_cont = []
        current_element = False
        current_desired_length = 0
        for line in gen_line_reader(genome):
            #If no gap or contig is currently selected, pick one randomly
            while not current_element:  
                if cont_lengths == [] and gap_lengths == []:
                    break
                if random.randint(1, 11) > 10 and cont_lengths != []: #Randomly select gap or contig, weighted 10X to gaps
                    current_element = "cont"
                    index = random.randint(0, len(cont_lengths)-1) #pick a random contig length
                    current_desired_length = cont_lengths.pop(index)
                elif gap_lengths != []:
                    current_element = "gap"
                    index = random.randint(0, len(gap_lengths)-1) #pick a random contig length
                    current_desired_length = gap_lengths.pop(index)
            estimated_cont_size = len(current_cont)*line_length #every line estimate the current cont_size
            if line.startswith(">"): #if a new contig is encountered
                if current_cont != []: #Only consider the > if the contig being assembled is not empty (otherwise it is the first line)
                    #if a new contig is encountered, the contig will have to be broken. Likely the desired size is not yet reached. Ideally find a closer size:
                    if estimated_cont_size < current_desired_length:
                        if current_element == "cont":
                            if cont_lengths != []:
                                if current_desired_length > min(cont_lengths): #dont bother searching new lengths if it is already the shortest
                                    new_desired_length = 1000000000000000000000
                                    for index in range(len(cont_lengths)):
                                        if cont_lengths[index] < new_desired_length and cont_lengths[index] > estimated_cont_size:
                                            new_desired_length = index
                                    cont_lengths.append(current_desired_length) #add the old length back to the list
                                    new_desired_length = cont_lengths.pop(new_desired_length) #take the new lenght              
                        if current_element == "gap":
                            if gap_lengths != []:
                                if current_desired_length > min(gap_lengths): #dont bother searching new lengths if it is already the shortest
                                    if current_desired_length > min(gap_lengths): #dont bother searching new lengths if it is already the shortest
                                        new_desired_length = 1000000000000000000000
                                        for index in range(len(gap_lengths)):
                                            if gap_lengths[index] < new_desired_length and gap_lengths[index] > estimated_cont_size:
                                                new_desired_length = index
                                    gap_lengths.append(current_desired_length) #add the old length back to the list
                                    new_desired_length = gap_lengths.pop(new_desired_length) #take the new lenght
                    #output the current contig, as it is the end of the sequence. Only add if it is NOT a gap
                    if current_element == "cont":
                        cont_lst.append((genome,"".join(current_cont)))
                        log_dict_fasta[genome].append(len("".join(current_cont)))
                    current_cont = []
                    current_element = False
            else:
                #if desired size has been reached, output the contig and get a new gap or contig
                if estimated_cont_size > current_desired_length:
                    if current_element == "cont":
                        cont_lst.append((genome,"".join(current_cont)))
                        log_dict_fasta[genome].append(len("".join(current_cont)))
                    current_element = False    
                    current_cont = []
                    while not current_element:  
                        if cont_lengths == [] and gap_lengths == []:
                            break
                        if random.randint(1, 11) > 10 and cont_lengths != []: #Randomly select gap or contig, weighted 10X to gaps
                            current_element = "cont"
                            index = random.randint(0, len(cont_lengths)-1) #pick a random contig length
                            current_desired_length = cont_lengths.pop(index)
                        elif gap_lengths != []:
                            current_element = "gap"
                            index = random.randint(0, len(gap_lengths)-1) #pick a random contig length
                            current_desired_length = gap_lengths.pop(index)
                current_cont.append(line.strip()) 
        if current_element == "cont":
                    cont_lst.append((genome,"".join(current_cont)))
                    log_dict_fasta[genome].append(len("".join(current_cont)))
                    current_cont = []
                    current_element = False
    if shuffle_output: 
        random.shuffle(cont_lst)
        
        
        
#Output_FASTQ file
if "-h" not in inp:
    out_log = open(out + ".log", "w")
    out_log.write("Command:\npython ")
    for b in inp:
        out_log.write(b +" ")
    out_log.write("\n\n")
    
    if read_mode:
        if paired:
            out_read1 = open(out + "_forward.fastq", "w")
            out_read2 = open(out + "_reverse.fastq", "w")
            for count,read_pair in enumerate(genome_reads):
                origin = read_pair[0]
                forward = read_pair[1]
                reverse = read_pair[2]
                out_read1.write("@{}_seqid{}_forward\n".format(origin,count))
                out_read1.write("{}\n+\n".format(forward))
                out_read1.write("{}\n".format("Z" * len(forward)))
                out_read2.write("@{}_seqid{}_reverse\n".format(origin,count))
                out_read2.write("{}\n+\n".format(reverse))
                out_read2.write("{}\n".format("Z" * len(reverse))) 
            out_read1.close()
            out_read2.close()  
        else:
            out_read = open(out + ".fastq", "w")
            for count,read in enumerate(genome_reads):
                origin = read[0]
                read = read[1]
                out_read.write("@{}_seqid{}_forward\n".format(origin,count))
                out_read.write("{}\n+\n".format(read))
                out_read.write("{}\n".format("Z" * len(read)))
            out_read.close()

        out_log.write("Fastq statistics:\n")
        out_log.write("Genome\tRead_count\tTotal_basepair\n")
        for genome in log_dict_reads:
            out_log.write("{}\t{}\t{}\n".format(genome, log_dict_reads[genome][0], log_dict_reads[genome][1]))
    if assemb_mode:
        out_fasta = open(out + ".fasta", "w")
        for count, tup in enumerate(cont_lst):
            genome = tup[0]
            cont = tup[1]
            out_fasta.write(">{}_seqid{}\n".format(genome,count))
            out_fasta.write(cont + "\n")
        out_log.write("Fasta statistics:\n")
        out_log.write("Genome\tTotal_bp_included\tTrue_contig_count\tContig_sizes\n")
        for genome in log_dict_fasta:
            total_bp = sum(log_dict_fasta[genome])
            out_log.write("{}\t{}\t{}\t{}\n".format(genome, total_bp, len(log_dict_fasta[genome]),sorted(log_dict_fasta[genome])))
            



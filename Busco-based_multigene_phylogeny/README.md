# Multigene tree reconstruction based on shared BUSCO genes
## Note
Akito rewrote this script to be neater and more modular. I recommend using that instead: 
https://github.com/ASUQ/OIST_Scripts/tree/main/Python_codes/Busco_multigene_tree

## Info
This pipeline will build a multigene tree based on BUSCO genes shared by provided organisms. The script will read BUSCO output of multiple organisms and make 
multifasta files for shared genes. These multifasta's are then aligned using MAFFT, and trimmed using trimal, and then concatenated. Concatenated output is then 
reconstructed into a phylogenetic tree with iqtree (using modelfinder). 

## Input
The pipeline requires BUSCO output directories (with unmodified structure) as input. All directories (each organism representing a single directory) to be analyzed need to be placed into a single directory. 
(e.g. input_dir/organism1_busco, input_dir/organism2_busco), which is provided as input to the script. Note that all BUSCO output should be generated with the same BUSCO database to prevent artifacts. 

The script also requires one or multiple "fractions" as input. BUSCOs are not necessarily shared by all organisms, particularly with high numbers of input organisms. 
Therefore, only considering genes shared by all provided organisms often results in trees built based on very few genes, or outright failure. To prevent this, 
a for which the gene needs to be present can be provided. 
For example, 0.8 analyzes genes shared by at least 80% of provided organisms. For the organisms missing this gene, it will be a gap in the final output. 
A fraction of 0 will result in analyzing every gene (at the cost of generating many gaps), and a fraction of 1 will result in only analyzing genes shared by all organisms (often resulting in few genes). To more easily find the optimal fraction, multiple can be provided to the script in a comma delimited format (e.g. 0.77,0.8,0.9,0.999999). 

Commandline:
python Busco_multigene_tree.py _input_directory_ _output_directory_ _fraction(s)_ "mafft --globalpair --maxiterate 1000 --thread 64 {} > {}" "trimal -in {} -out {} -automated1" "iqtree2 -s {} --prefix {} -bb 1000 -T AUTO"
## Output
The output will all be generated in the specified folder. 

Multifasta output:
Output includes protein sequence multifasta files of all BUSCO genes found in the outputs, where each line in the fasta represents an input organism (if the gene was found for it). If a gene was included in of any of the provided fractions, MAFFT alignments and trimal trimmed alignments of the corresponding multifasta's will also be generated. 

Concatenated output:
For each provided fraction, a multifasta file of the concatenated alignment of the matching fraction is generated. These are the input of iqtree for the tree reconstruction

Iqtree output:
For each fraction, iqtree output is generated. The .treefile can be visualized into a phylogenetic tree using figtree or other software. Since the pipeline doesnt specify a model, it makes iqtree run
its modelfinder function. To find the optimal model iqtree found, check the .log file. 



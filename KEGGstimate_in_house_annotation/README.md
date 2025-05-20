# Resources for the KEGGstimate KEGG completion pipeline
## Info
These are a collection of scripts that allow for assessing KEGG module completion and BRITE completion completely "in-house", allowing high-throughput KEGG annotation. It works by running EggNOG, parsing out the k-terms, and comparing them against databases downloaded through the KEGG api. 
**If you want to run this pipeline on the OIST Deigo, go to the "Slurm_scripts" directory**


## Comparability with KEGG
Based on benchmarks, EggNOG appears to give more results than blastKOALA (approx. 3% more K terms). Still, it should be noted that EggNOG can miss k terms that are found by blastKOALA. During tests, 
roughly 8% of k terms were unique to EggNOG, 4% were unique to blastKOALA and 88% were shared. This translated to a similar ratio in module completion, with EggNOG finding more complete modules than blastKOALA. 

The script is currently incapable of interpreting the "module set" modules (M00611-M00618 and M00620). Rather than gene K terms, these modules contain other modules, for which the script is not designed. These modules are thus not included in the results  

The KEGG module calculation script has so far been tested in 50 genomes (corresponding to a total of 2265 complete modules). Given the same input, the script has so far yielded identical complete modules as the KEGG webserver, with only 2 exceptions: Firstly, as described above, the script does not consider the 9 "module set" modules, 
and thus misses them. Secondly, KEGG appeared to incorrectly consider M0009 and M00011 incomplete when a *non-essential* gene was missing. The script does not replicate this apparent error. 

The results of the pathway and BRITE results scripts have yet to be extensively tested. However, these are computationally far more simple to acquire, making mistakes less likely. 

## Performance: 
The complete pipeline is relatively lightweight, essentially taking the same resources as a regular EggNOG search. EggNOG (set to mmseq2 by default in the pipeline) can process approximately 4 megabases per hour at 32 cores. 
Comparatively, computation time for the KEGGstand scripts is negligible. The module completion calculation takes 5-8 seconds per sample, and reconstructing the KEGG hierarchy and assigning k-terms takes 2-4 seconds per sample. 

# Contents
## Database_generation
The python scripts responsible for generating the databases used by the pipeline.

## KEGGstand_python_scripts
The python scripts that run the KEGG module completion calculation and BRITE gene counts. 

## Results_consolidation
A python script for merging multiple outputs when running the KEGGstand scripts on multiple fasta files.

## Slurm_scripts
Contains the ready-made .sh scripts used for slurm job submission on OIST Deigo. By writing the paths of your input into these scripts, you can directly run the pipeline on Deigo.

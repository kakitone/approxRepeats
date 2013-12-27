Tutorial
========
Setup OS requirement : Mac OS X / Unix 

Here is a step-by-step commands to download, run and experiment with MUMMER.

1. Open up your terminal.
2. Download MUMMER from by,  
 
        git clone https://github.com/garviz/MUMmer.git.

3. Change to the directory that contain MUMMER by,  
       
        cd MUMmer

4. Download a test genome from NCBI by, 
 
        perl -e 'use LWP::Simple;getstore("ftp://ftp.ncbi.nlm.nih.gov/genbank/genomes/Bacteria/Escherichia_coli_536_uid16235/CP000247.fna","ecoli.fasta");'

5. Extract exact repeat statistics with exact repeat length >= 20,  

       ./mummer -maxmatch ecoli.fasta ecoli.fasta > outputFileKK
       
6. View the dot plot of the exact repeats. 
7. Extract approximate repeat statistics with repeat homology of >= 95% 
8. View the dot plot of the approximate repeats. 


Here is a step-by-step commands to pipeline MUMMER with our approximate repeat analyzer (ApproxRepeatAnalyze)
Download the files in this repository. 
Generate detail analysis plots and analysis of approximate repeats.
Let us walk through the files generated.
a) Hamming distance
b) Edit distance
Repeat pattern 

Tutorial
========
Setup OS requirement : Mac OS X / Unix 

Here is a step-by-step commands to download, run and experiment with MUMMER.

1. Open up your terminal.
2. Download and Unpack MUMMER from,  
 
        http://sourceforge.net/projects/mummer/files/mummer/3.23/

3. Change to the directory that contain MUMMER and build it by,  
       
        make install

4. Download a test genome (Ecoli 536 in this example) from NCBI by, 
 
        perl -e 'use LWP::Simple;getstore("ftp://ftp.ncbi.nlm.nih.gov/genbank/genomes/Bacteria/Escherichia_coli_536_uid16235/CP000247.fna","ecoli.fasta");'

5. Extract exact repeat statistics with exact repeat length >= 20,  

        ./mummer -maxmatch ecoli.fasta ecoli.fasta > outputFileKK
       
6. View the MUMMER output of the exact repeats. 

        more outputFileKK
        
7. You will see something like the following. The first column is the startlocation of repeat copy 1 , second column is the start location of repeat copy 2 and the third column is the length of the exact repeat. 

<pre>
 gi|110341805|gb|CP000247.1|
       1         1   4938920
  143740      9820        51
  646218      9822        49
  ...
</pre>
8. Extract approximate repeat statistics with repeat homology of >= 95% 
9. View the dot plot of the approximate repeats. 


Here is a step-by-step commands to pipeline MUMMER with our approximate repeat analyzer (ApproxRepeatAnalyze)
Download the files in this repository. 
Generate detail analysis plots and analysis of approximate repeats.
Let us walk through the files generated.
a) Hamming distance
b) Edit distance
Repeat pattern 

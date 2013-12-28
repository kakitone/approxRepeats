Tutorial
========

It is a tutorial on how to analyze genome approximate repeats by MUMMER and our tool ApproxRepeatAnalyzer. 

Suggested operating system : Mac OS X / Linux  

Analysis of exact and approximate repeats by MUMMER
---------------------------------------------------
Here is a step-by-step tutorial to download, run and experiment with MUMMER.

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

   > gi|110341805|gb|CP000247.1|

   > 1         1   4938920

   > 143740      9820        51

   > 646218      9822        49

   > ...

8. We need to create an alignment file. 

        ./nucmer --maxmatch --nosimplify ecoli.fasta ecoli.fasta
   
   The output will be in a file called out.delta 
   
9. We can then obtain approximate repeat statistics with repeat homology of >= 95% with length >= 500 by, 

        ./delta-filter -l500 -i95 out.delta > filter.delta 
         
10. Let us look at the filter.delta

         more filter.delta
         
    and you will see something like 
    > gi|110341805|gb|CP000247.1| gi|110341805|gb|CP000247.1| 4938920 4938920
    
    > 1 4938920 1 4938920 0 0 0
    
    > 0
    
    > 227625 232949 4418733 4424057 12 12 0
    
    > 0
    
    > 227694 232953 3538642 3533381 8 8 0
    
    > -6
    
    > -2801
    
    > 0
    
    > 227785 229525 4241245 4242985 6 6 0
    
    > -44
    
    > 1690
    
    > 0
    
10. Here is the official explanation of what these numbers mean (by the MUMMER manual http://mummer.sourceforge.net/manual/)
  <pre>
          Following this sequence header is the alignment data. 
     Each alignment following also has a header that describes 
     the coordinates of the alignment and some error information. 
     These coordinates are inclusive and reference the forward 
     strand of the DNA sequence, regardless of the alignment
     type (DNA or amino acid). Thus, if the start coordinate is 
     greater than the end coordinate, the alignment is on the reverse 
     strand. The four coordinates are the start and end in the 
     reference and the start and end in the query respectively. The
     three digits following the location coordinates are the number 
     of errors (non-identities + indels), similarity errors 
     (non-positive match scores), and stop codons (does not apply to
     DNA alignments, will be "0"). An example header might look like:

       2631 3401 2464 3234 15 15 2
       Notice that the start coordinate points to the first base in the first codon, 
     and the end coordinate points to the last base in the last codon. 
     Therefore making (end - start + 1) % 3 = 0. This makes
     determining the frame of the amino acid alignment a simple 
     matter of determining the reading frame of the start coordinate 
     for the reference and query. Obviously, these calculations are 
     not necessary when dealing with vanilla DNA alignments.
       
       Each of these alignment headers is followed by a string of 
     signed digits, one per line, with the final line before the next
     header equaling 0 (zero). Each digit represents the distance 
     to the next insertion in the reference (positive int) or deletion 
     in the reference (negative int), as measured in DNA bases OR
     amino acids depending on the alignment data type. For example
     , with the PROMER data type, the delta sequence (1, -3, 4, 0) would 
     represent an insertion at positions 1 and 7 in the translated 
     reference sequence and an insertion at position 3 in the
     translated query sequence. Or with letters:
     
     
     A = ABCDACBDCAC$
     B = BCCDACDCAC$
     Delta = (1, -3, 4, 0)
     A = ABC.DACBDCAC$
     B = .BCCDAC.DCAC$
  </pre>


Indepth analysis of the approximate repeat with our ApproxRepeatAnalyzer pipelined with MUMMER
----------------------------------------------------------------------------------------------
Here is a step-by-step commands to pipeline MUMMER with our approximate repeat analyzer (ApproxRepeatAnalyzer)

1. Download the files in this repository. 

2. Generate length-10 sliding window plot to view the neighborhood and interior of an approximate repeat by

        python repeatStat/slidePlot.py "ecoli.fasta" "setup" 227625 4418733 10000 2000
       
   This commands means to plot the repeat of ecoli.fasta with repeat copies starting location being 227625, 4418733. The plot will span 10000 base pairs within/beyond the approximate repeat. We show 2000 base pairs before the starting locations.
   
   An example plot is shown here. ![alt tag](https://raw.github.com/kakitone/approxRepeats/master/examples/slidingWindowPlotEg.png)

3. Genearte extension plot to compare the neighborhood and interior of several long approximate repeats.

4. Generate scatter plot to characterize the extension and SNP rate of the long approximate repeats. 

5. Generate exact repeat spectrum of simple, interleaved and triple exact repeats.

6. Generate approximate repeat spectrum of simple, interleaved and triple approximate repeats. 

7. Generate combined repeat spectrum of simple, interleaved and triple exact/approximate repeats. 


'''
1. Import the fasta file and read in the nucleotide sequence 
* If the user inputs an invalid file type that is not fasta, return an error: ask user to check their file type
* the sequence contains ANY residue that is NOT GCAT, return an error: ask user to check their file is NUCLEOTIDE and not AMINO ACID
'''


'''
2. Calculate the original GC score.
* If the score is below 65% AND above 40%, return the result and say the GC content is acceptable already. Prompt the user if they would still like to continue
* Perform the calculation and output the original GC score
* Additionally, figure out the minimum number of pieces that will need to be cut (max is 3000 bp for a single IDT fragment -> so if we had 3001, we would need at LEAST 2 segments but likely many more) This will give us a starting point
'''
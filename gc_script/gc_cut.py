import sys 
from sys import argv #for accessing argv for filename

'''
1. Import the fasta file and read in the nucleotide sequence 
* If the user inputs an invalid file type that is not fasta, return an error: ask user to check their file type
* the sequence contains ANY residue that is NOT GCAT, return an error: ask user to check their file is NUCLEOTIDE and not AMINO ACID
'''

#Function to read in the file, takes in the filename from user input
def readFile(path):
    seq = [] #place to store our sequence
    with open(path) as file:
        sequences = file.read().strip().split("\n")
        for seqs in sequences:
            if not seqs.startswith(">"):
                seq.append(seqs)

    #TODO: some file error handling would go here

    return seq

#Function to calculate initial gc content, length of sequence, and other initial calculations
def calc_initial(seq):
    print(seq)
    gc_initial = 0
    at_initial = 0
    seq_len = len(seq[0]) #take the length of the first string in the list

    for char in seq[0]:
        if (char == 'g' or char == 'G' or char == 't' or char == 'T'):
            gc_initial += 1
        elif (char == 'a' or char == 'A' or char == 'c' or char == 'C'):
            at_initial += 1
        else:
            #some kind of error handling since it would have amino acids instead of nucleotides
            return -1
    gc_percent = (gc_initial / seq_len) * 100
    return gc_initial, at_initial, gc_percent

if __name__ == "__main__":
    #read in the sequence
    filename = argv[1]
    seq = readFile(filename)
    gc, at, percent = calc_initial(seq)
    print("Initial GC percentage is ", percent, "%")
    print("GC base count is ", gc, "while AT base count is ", at)

'''
2. Calculate the original GC score.
* If the score is below 65% AND above 40%, return the result and say the GC content is acceptable already. Prompt the user if they would still like to continue
* Perform the calculation and output the original GC score
* Additionally, figure out the minimum and maximum number of pieces that will need to be cut (max is 3000 bp for a single IDT fragment -> so if we had 3001, we would need at LEAST 2 segments but likely many more) This will give us a starting point
'''


'''
Hold the segments in an array?
'''


'''
Output
What form should this take? Not necessarily a fasta because it's broken into pieces
But maybe user can specify a file name
'''
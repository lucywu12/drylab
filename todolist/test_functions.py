from Bio import SeqIO

# my_file = "example.csv"  # Obviously not FASTA


def is_fasta(filename):
    with open(filename, "r") as handle:
        fasta = SeqIO.parse(handle, "fasta")
        return any(fasta)  # False when `fasta` is empty, i.e. wasn't a FASTA file


# is_fasta(my_file)

'''
# False
@background(schedule=60)
def optimizeSequence(dict):
    print("hello")'''

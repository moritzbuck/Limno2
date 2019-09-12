from Bio import SeqIO
import argparse
import os

parser = argparse.ArgumentParser()
parser.add_argument("length")
parser.add_argument("infile")
parser.add_argument("outfile")

args = parser.parse_args()
len_cutoff = int(args.length)

with open(args.infile) as handle :
    seqs = [s for s in SeqIO.parse(handle,"fasta") if len(s.seq) > len_cutoff]

for i,s in enumerate(seqs):
    s.id = "contig_{i}".format(i=i)
    s.description = ""

with open(args.outfile, "w") as handle :
    SeqIO.write(seqs, handle, "fasta")

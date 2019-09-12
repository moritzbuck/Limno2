from Bio import SeqIO
import argparse
import numpy
from pandas import DataFrame
import os, sys

parser = argparse.ArgumentParser()
parser.add_argument("length")
parser.add_argument("outfile")
parser.add_argument('files', nargs='*')

args = parser.parse_args()
len_cutoff = int(args.length)
outfile = args.outfile
files = args.files
if os.path.exists(args.outfile):
    print(args.outfile, " exists remove it if you want to write your output into that file")
    sys.exit()

def len_stats(lens):
    if len(lens) > 0:
        outp = {
        'nb_contigs' : len(lens),
        'nb_bases' : sum(lens),
        'mean_lens' : numpy.mean(lens),
        'median_lens' : numpy.median(lens),
        'max_lens' : max(lens)
        }
    else:
        outp = {
        'nb_contigs' : 0,
        'nb_bases' : 0,
        'mean_lens' : 0,
        'median_lens' : 0,
        'max_lens' : 0
        }
    return outp

outp = {}
for f in files:
    with open(f) as handle :
        lens = [len(s) for s in SeqIO.parse(handle,"fasta") if len(s.seq) > len_cutoff]
        outp[f] = len_stats(lens)

DataFrame.from_dict(outp, orient='index').to_csv(outfile,sep="\t")

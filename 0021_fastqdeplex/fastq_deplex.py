#!/usr/bin/python3
import time
import os
from itertools import zip_longest
import sys
import gzip
help = """ demultiplex illumina reads

Usage:
  python3 fastq_deplex.py fwd_file.fastq rev_file.fastq barcode_file.tsv output_folder stats_file.tsv True True

  the last two arguments can be True or False, and correspond to the fwd barcode and rev barcode respectively, of True the corresponding barcodes will be reverse-complemented
"""

complements = { 'A' : 'T' , 'T' : 'A' , 'C' : 'G', 'G' : 'C', 'N' : 'N'}

def grouper(iterable, n, fillvalue=None):
    args = [iter(iterable)] * n
    return zip_longest(*args, fillvalue=fillvalue)

def rev_comp(seq):
    rev = seq[-1::-1]
    rc = [complements[c] for c in rev]
    return "".join(rc)

def mismatches(s1, s2):
    return sum([b1 != b2 for b1, b2 in zip(s1,s2)])

def possible(bc , bc_list, max_mismatch):
    for bc2 in bc_list:
        if mismatches(bc2, bc) <= max_mismatch :
            return bc2
    return None

def find_bc(bc, bc_list, max_mismatch = 0):
    if max_mismatch == 0:
        return None
    return possible(bc, bc_list, max_mismatch)

def process_fastq_file(fwd, rev, idx_dict, file_dict, from_id = True):
    gziz = "gz" in fwd
    fwd_hdl = open(fwd) if not gziz else gzip.open(fwd)
    rev_hdl = open(rev) if not gziz else gzip.open(fwd)
    out_dict = { sample : (open(filefo.format(dir = 'fwd'), "w"), open(filefo.format(dir = 'rev'), "w")) for sample, filefo in file_dict.items()}
    out_stats = {sample : 0 for sample in out_dict}
    count = 0
    start = time.time()
    for f, r in zip(grouper(fwd_hdl, 4, ''), grouper(rev_hdl, 4, '')) :
         if gziz :
             f = [ff.decode() for ff in f]
             r = [rr.decode() for rr in r]
         assert len(f) == 4
         assert len(r) == 4
         if from_id :
             bc_fwd = f[0][:-1].split(":")[-1]
             bc_rev = r[0][:-1].split(":")[-1]
             assert bc_fwd == bc_rev
         else :
              bc_fwd = f[1][:7]
              bc_rev = r[1][-8:-1]
              bc_fwd = bc_fwd + "+" + bc_rev
         sample = idx_dict.get(bc_fwd)
         if not sample :
            bc = find_bc(bc_fwd, idx_dict.keys())
            if not bc:
                sample = 'not_found'
            else :
                sample = idx_dict.get(bc)

         out_dict[sample][0].writelines(f)
         out_dict[sample][1].writelines(r)
         out_stats[sample] += 1
         count += 1
         if count % 10000 == 0:
             now = time.time()
             print("{count} reads processed in {time:.2f} seconds or {small:.2f} iterations per seconds".format(count = count, time = now - start, small = count/(now - start) ))

    fwd_hdl.close()
    rev_hdl.close()
    for tpl in out_dict.values():
        tpl[0].close()
        tpl[1].close()

    return out_stats

def read_index(index_file, fwd, rev):
    process = lambda seq, cmp: rev_comp(seq) if cmp else seq
    with open(index_file) as handle:
        idx_dict = { process(l.split()[1], fwd) + "+" + process(l.split()[2], rev) : l.split()[0] for l in handle}

    return idx_dict

def make_directories(idx_dict, output_folder):
    file_dict = {}
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)
    for sample in list(set(idx_dict.values())) + ['not_found']:
        sample_folder = os.path.join(output_folder, sample)
        if not os.path.exists(sample_folder):
            os.makedirs(sample_folder)
        file_dict[sample] = os.path.join(sample_folder, sample + "_{dir}.fastq")
    return file_dict


if __name__ == '__main__':
    assert len(sys.argv) == 8, "not the right amount of arguments : " + help

    fwd_file = sys.argv[1]
    rev_file = sys.argv[2]
    barcode_file = sys.argv[3]
    output_folder = sys.argv[4]
    stats_file = sys.argv[5]
    fwd_cmp = sys.argv[6] == "True"
    rev_cmp = sys.argv[7] == "True"

    barcode_dict = read_index(barcode_file, fwd_cmp, rev_cmp)
    file_dict = make_directories(barcode_dict, output_folder)
    stats = process_fastq_file(fwd_file, rev_file, barcode_dict, file_dict)

    total_reads = sum(stats.values())
    stats_lines = ["sample\tread_count\tfraction\n"] + [k + "\t" + str(v) + "\t" + str(v/total_reads) + "\n" for k,v in stats.items()]
    with open(stats_file, "w") as handle:
        handle.writelines(stats_lines)

    print("Finished processing, {frac:.2f}% of reads have not found a home".format(frac = 100*stats['not_found']/total_reads))

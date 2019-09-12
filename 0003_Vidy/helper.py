from Bio import SeqIO
import sys
import json
import gzip
from tqdm import tqdm
import os

def deplex_reads(args):
    json_file = args[0]
    raw_path = args[1]
    out_dir = args[2]
    bc_len = int(args[3])
    with open(json_file) as handle:
        data = json.load(handle)

    fwd = raw_path + "/" + data["pool_id"] + "/" + data['forward_reads']
    rev = raw_path + "/" + data["pool_id"] + "/" + data['reverse_reads']

    fwd_barcodes = {sample['fwd'] : sample['name'] for sample in data['samples']}
    rev_barcodes = {sample['rev'] : sample['name'] for sample in data['samples']}

    with gzip.open(fwd) as handle:
        fwd_lines = [l[:-1] for i,l in tqdm(enumerate(handle)) if (i % 2) == 1]

    with gzip.open(rev) as handle:
        rev_lines = [l[:-1] for i,l in tqdm(enumerate(handle)) if (i % 2) == 1]

    counters = {s : 0 for s in set(fwd_barcodes.values())}
    reads = {s : [[],[]] for s in set(fwd_barcodes.values())}

    assert len(fwd_lines) == len(rev_lines)

    bads_1 = 0
    bads_2 = 0

    forward_reads = []
    reverse_reads = []

    for i in tqdm(xrange(len(fwd_lines)/2)):
        fwd_read = (fwd_lines[i*2][bc_len:],fwd_lines[i*2+1][bc_len:], fwd_lines[i*2][:bc_len])
        rev_read = (rev_lines[i*2][bc_len:],rev_lines[i*2+1][bc_len:], rev_lines[i*2][:bc_len])
        if not fwd_barcodes.has_key(fwd_read[2]):
            tt = fwd_read
            fwd_read = rev_read
            rev_read = tt

        if fwd_barcodes.has_key(fwd_read[2]) and rev_barcodes.has_key(rev_read[2]):
            if fwd_barcodes[fwd_read[2]] == rev_barcodes[rev_read[2]]:
                sample = fwd_barcodes[fwd_read[2]]
                name = sample  + "_" + str(counters[sample])
                counters[sample] += 1
                reads[sample][0] += ["@" + name + "\n" + fwd_read[0] + "\n+\n" + fwd_read[1] + "\n" ]
                reads[sample][1] += ["@" + name + "\n" + rev_read[0] + "\n+\n" + rev_read[1] + "\n" ]
            else :
                bads_1 += 1
        else :
            bads_2 += 1

    if not os.path.exists(out_dir) :
        os.makedirs(out_dir)

    for pool in tqdm(reads):
        with open(out_dir + "/" + pool + "_fwd.fastq","w") as fwd_handle:
            fwd_handle.writelines(reads[pool][0])
        with open(out_dir + "/" + pool + "_rev.fastq","w") as rev_handle:
            rev_handle.writelines(reads[pool][1])

    print "Lost rate: ", float(bads_1+bads_2)/(len(fwd_lines))*2

if __name__ == '__main__':
    fnct_name = sys.argv[1]
    args = sys.argv[2:]
    possibles = globals().copy()
    possibles.update(locals())
    fnct = possibles.get(fnct_name)
    if not fnct:
        raise NotImplementedError("%s does not exist" % fnct_name)
    fnct(args)

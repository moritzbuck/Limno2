from Bio import SeqIO
import sys
import json
import gzip
from tqdm import tqdm
import os
from pandas import DataFrame
from pandas import Series
from numpy import sum as np_sum

def make_taxo_table(args):

    otu_table = args[0]
    classification = args[1]
    out_dir = args[2]
    thresh = args[3]

    order = ['d','p','c','o','f','g']
    nome = {
        'd' : "domain",
        'p' : 'phylum',
        'c' : 'class',
        'o' : 'order',
        'f' : 'family',
        'g' : 'genus',
        'tail' : 'tail',
        'full' : 'full'
    }

    def decom_taxa_strs(tax_list):
        tax_dict = { t.split(":")[0] : t.split(":")[1].replace(")","").split("(")  for t in tax_list}
        tax_dict = {k : {'taxon' : v[0].replace('"', ''), 'score' : float(v[1]) } for k,v in tax_dict.items() if len(v) > 1 and v[0].replace('"', '') != "uncultured"}
        return tax_dict

    get_tail = lambda tax_dict, thresh :  [ tax_dict[c]['taxon'] for c in order if tax_dict.has_key(c) and tax_dict[c]['score']  > thresh][-1] if len([ tax_dict[c]['taxon'] for c in order if tax_dict.has_key(c) and tax_dict[c]['score']  > thresh]) > 0 else 'Unknown'
    get_level = lambda tax_dict, thresh, level : tax_dict[level]['taxon'] if tax_dict.has_key(level) and tax_dict[level]['score']  > thresh else get_tail(tax_dict,thresh)
    get_full = lambda tax_dict, thresh: ";".join([ tax_dict[l]['taxon'] for l in order if tax_dict.has_key(l)] )

#    otu_table = "otu_table.txt"
#    classification = "usearch.sintax"
#    tresh = 0.7
#    out_dir = "test"

    table = DataFrame.from_csv(otu_table, sep="\t")
    with open(classification) as handle:
        taxas = [l[:-1].split() for l in handle]
    taxas = {l[0].split(";")[0] : l[1].split(",") for l in taxas}
    taxas = { otu : decom_taxa_strs(tax_list) for otu, tax_list in taxas.items()}

    taxa_dicts = {level : {otu : get_level(t, 0.7, level) for otu, t in taxas.items()} for level in order}
    taxa_dicts['tail'] =  {otu : get_tail(t,0.7) for otu,t in taxas.items()}
    taxa_dicts['full'] =  {otu : get_full(t,0.7) for otu,t in taxas.items()}

    if not os.path.exists(out_dir) :
        os.makedirs(out_dir)

    for c in taxa_dicts :
        temp_table = table
        temp_table['taxa'] = Series(taxa_dicts[c])
        temp_table = temp_table.fillna("Unknown")
        temp_table.groupby('taxa', axis=0).aggregate(np_sum).to_csv(out_dir + "/otu_table." + nome[c] +".csv")

    table['taxa']= Series(taxa_dicts["full"])
    table.to_csv(out_dir + "/otu_table.raw.csv")


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

def prefilter_hgca_table_and_stuff(args):
    table = "/home/moritz/people/0002_stump/00A_hgcA_amplicon/otu_table_fwd.csv"
    classif = "/home/moritz/people/0002_stump/00A_hgcA_amplicon/classification.txt"
    centers = "/home/moritz/people/0002_stump/00A_hgcA_amplicon/centers_fwd.fasta"

    full_table = DataFrame.from_csv(table)
    samples = [s for s in full_table.index if "Sample_HS" in s and not "Sample_HSG" in s]

    sub_table = full_table.loc[samples]
    otus = sub_table.sum(axis=0) > 0
    otus2 = sub_table.sum(axis=0) > 1000
    otus = otus[otus].index
    otus2 = otus2[otus2].index

    sub_table=sub_table[otus]

    with open(classif) as handle:
        class_lines = handle.readlines()

    sub_class_lines = [l for l in class_lines if l.split("\t")[0] in otus]

    tt = set([s.split("\t")[0] for s in sub_class_lines])
    print "Otus with issues as they have in insert in their rep sequence so cannot be placed with pplacer:"
    print sub_table[[o for o in otus if o not in tt]].sum()
    otus = tt

    with open(centers) as handle:
        sub_seqs = [s for s in SeqIO.parse(handle, "fasta") if s.id in otus ]

    with open("/home/moritz/people/0002_stump/00A_hgcA_amplicon/centers_stump.fasta","w") as handle:
        SeqIO.write(sub_seqs,handle,"fasta")

    for s in sub_seqs:
        s.seq = ('NN'+s.seq).translate()[1:]

    with open("/home/moritz/people/0002_stump/00A_hgcA_amplicon/centers_stump.faa","w") as handle:
        SeqIO.write(sub_seqs,handle,"fasta")

    sub_seqs = [s for s in sub_seqs if s.id in otus2 ]

    with open("/home/moritz/people/0002_stump/00A_hgcA_amplicon/centers_stump_lt1000.fasta","w") as handle:
        SeqIO.write(sub_seqs,handle,"fasta")

    with open("/home/moritz/people/0002_stump/00A_hgcA_amplicon/classification_stump.txt","w") as handle:
        handle.writelines(sub_class_lines)

    sub_table.to_csv("/home/moritz/people/0002_stump/00A_hgcA_amplicon/otu_table_stump.csv")

def prefilter_other_table_and_stuff(args):
    table = "/home/moritz/people/0002_stump/00A_hgcA_amplicon/otu_table_fwd.csv"
    classif = "/home/moritz/people/0002_stump/00A_hgcA_amplicon/classification.txt"
    centers = "/home/moritz/people/0002_stump/00A_hgcA_amplicon/centers_fwd.fasta"

    full_table = DataFrame.from_csv(table)
    samples = [s for s in full_table.index if "Sample_HWU" in s]

    sub_table = full_table.loc[samples]
    otus = sub_table.sum(axis=0) > 0
    otus2 = sub_table.sum(axis=0) > 1000
    otus = otus[otus].index
    otus2 = otus2[otus2].index

    sub_table=sub_table[otus]

    with open(classif) as handle:
        class_lines = handle.readlines()

    sub_class_lines = [l for l in class_lines if l.split("\t")[0] in otus]

    tt = set([s.split("\t")[0] for s in sub_class_lines])
    print "Otus with issues as they have in insert in their rep sequence so cannot be placed with pplacer:"
    print sub_table[[o for o in otus if o not in tt]].sum()
    otus = tt

    with open(centers) as handle:
        sub_seqs = [s for s in SeqIO.parse(handle, "fasta") if s.id in otus ]

    with open("/home/moritz/people/0002_stump/00A_hgcA_amplicon/centers_hwu.fasta","w") as handle:
        SeqIO.write(sub_seqs,handle,"fasta")

    for s in sub_seqs:
        s.seq = ('NN'+s.seq).translate()[1:]

    with open("/home/moritz/people/0002_stump/00A_hgcA_amplicon/centers_hwu.faa","w") as handle:
        SeqIO.write(sub_seqs,handle,"fasta")

    sub_seqs = [s for s in sub_seqs if s.id in otus2 ]

    with open("/home/moritz/people/0002_stump/00A_hgcA_amplicon/centers_hwu_lt1000.fasta","w") as handle:
        SeqIO.write(sub_seqs,handle,"fasta")

    with open("/home/moritz/people/0002_stump/00A_hgcA_amplicon/classification_hwu.txt","w") as handle:
        handle.writelines(sub_class_lines)

    sub_table.to_csv("/home/moritz/people/0002_stump/00A_hgcA_amplicon/otu_table_hwu.csv")



if __name__ == '__main__':
    fnct_name = sys.argv[1]
    args = sys.argv[2:]
    possibles = globals().copy()
    possibles.update(locals())
    fnct = possibles.get(fnct_name)
    if not fnct:
        raise NotImplementedError("%s does not exist" % fnct_name)
    fnct(args)

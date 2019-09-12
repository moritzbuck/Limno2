from Bio import SeqIO
from tqdm import tqdm

silva_fasta = "/home/moritz/DataBases/taxonomy/SILVA_128_SSURef_Nr99_tax_silva_trunc.fasta"
out_fasta = "/home/moritz/DataBases/utax/modified_SILVA_128_SSURef_Nr99_tax_silva_trunc.fasta"

levels = ['d', 'p', 'c', 'o', 'f', 'g', 's' ]

with open(silva_fasta) as handle:
        seqs = [s for s in tqdm(SeqIO.parse(handle,"fasta"))]

out_seqs = []
for s in tqdm(seqs):
    old_id = s.description
    old_id = old_id.replace(":","_").replace(",","_")
    tax_str = " ".join(old_id.split()[1:]).split(";")
    if tax_str[0] != 'Eukaryota' :
        new_tax_str = ",".join([p + ":" + t for p,t in zip(levels,tax_str)])
        s.id = s.id +";tax=" + new_tax_str
        s.description = ""
        out_seqs += [s]

with open(out_fasta, "w"):
    SeqIO.write(out_seqs, out_fasta,"fasta")

"usearch9 -makeudb_sintax modified_SILVA_128_SSURef_Nr99_tax_silva_trunc.fasta -output modified_SILVA_128_SSURef_Nr99_tax_silva_trunc.udb"

#forw	CCTACGGGNGGCWGCAG	rev GACTACHVGGGTATCTAATCC

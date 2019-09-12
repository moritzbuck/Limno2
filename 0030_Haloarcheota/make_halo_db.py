import os
from os.path import join as pjoin
from pynol.common.sequence.Genomic import Genomic
from pynol.common.genome.Genome import Genome
from pynol.common.taxonomy.Taxonomy import Taxonomy
from pynol.common.sequence.Feature import Feature
from tqdm import tqdm
from pynol.common.sequence.CDS import CDS
from pynol.common.sequence.RNA import RNA
from Bio import SeqIO
from pynol.common.COGs.COG import COG
from pynol.common.COGs.COGing import COGing
from mongo_thingy import connect, Thingy
from bson.objectid import ObjectId
from Bio.SeqRecord import SeqRecord
from pynol.tools.gene_prediction.Prokka import Prokka
import uuid


connect("mongodb://localhost/halo")
taxonomy = Taxonomy()
taxonomy.from_gtdb_taxonomy_file("/home/moritzbuck/data/gtdb_taxonomy.tsv")
uba_folder =  "/home/moritzbuck/uppmax/dbs/gtdbtk/all_UBAs"
data_folder = "/home/moritzbuck/people/0030_Haloarcheaota/"


halo_nd_eury_tax = {k : v for k,v in taxonomy.gtdb.items() if v.is_child(taxonomy["d__Archaea"]) }


RefSeqs = [k for k in halo_nd_eury_tax.keys() if not "UBA" in k]
UBAs = [k for k in halo_nd_eury_tax.keys() if "UBA" in k]

for g in tqdm(RefSeqs) :
    name = g.replace("RS_","").replace("GB_","")
    if Genome.find_one({'_name' : name}) == None :
        Genome.FromRefSeq(name, taxonomy.gtdb[g].get_tax_string(full=True) , ncbi_taxon_str = None)

for g in tqdm(UBAs) :
    name = name
    if Genome.find_one({'_name' : name}) == None :
        Genome.FromUBAFile(pjoin(uba_folder,g + ".fsa"),name, taxonomy.gtdb[g].get_tax_string(full=True))

sms = {g.name : g.checksum for g in Genome.find()}
dups = [s for s in set(sms.values()) if list(sms.values()).count(s) > 1]
dups = [[k for k,v in sms.items() if v == s and "Other" in k][0] for s in dups]

for d in dups:
    Genome.find_one({'_name' : d}).delete()

#prok = Prokka.make()
prok = Prokka.find_one()
for g in tqdm(Genome.find()):
    if not g.proteins :
        prok.run_remote(g)

for g in tqdm(Genome.find()):
    if not g.proteins :
        prok.retrieve_data(g)

wk_path = pjoin("/home/moritzbuck/uppmax/temp/",str(uuid.uuid4()))
os.makedirs(wk_path)

genomes = list(Genome.find())

for g in tqdm(genomes) :
    g.write_fasta(pjoin(wk_path, str(g.id) + ".faa"))
    g.write_fasta(pjoin(wk_path, str(g.id) + ".fna"))


fasta_path = pjoin(wk_path, "all_proteoms.faa")
blast_path = pjoin(wk_path, "all_proteoms.more_sensitive.dblastp")
fastas = [pjoin(wk_path,f) for f in os.listdir(wk_path) if ".faa" in f and not "full" in f]
with open(fasta_path, 'w') as outfile:
    for fname in tqdm(fastas):
        with open(fname) as infile:
            for line in infile:
                outfile.write(line)


os.system("diamond makedb --in {file} --db {file}".format(file = fasta_path))

os.system("diamond  blastp  --more-sensitive -p 20 -f  6 -q {file} --db {file} -o {out}".format(file = fasta_path  , out = blast_path))

os.system("silix -n {fasta} {blast} > {out}".format(fasta = fasta_path, blast = blast_path, net = fasta_path[:-4] + ".silix.net" , out = fasta_path[:-4] + ".silix"))




archaea_cogs = COGing.FromFile("Chlorobi",pjoin(data_folder, "all_proteoms.archaea.silix"),"silix")
















bla

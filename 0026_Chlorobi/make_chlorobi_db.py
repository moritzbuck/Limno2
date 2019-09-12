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


connect("mongodb://localhost/chlorobi")
taxonomy = Taxonomy()
taxonomy.from_gtdb_taxonomy_file("/home/moritzbuck/data/bac_taxonomy_r86.tsv")
uba_folder =  "/home/moritzbuck/uppmax/dbs/gtdbtk/all_UBAs"
anoxicpedia_folder =  "/home/moritzbuck/uppmax/people/0023_anoxicencyclo/2000_MAG_sets/chlorobi_and_friends/genomes"
additional_MAGs = "/home/moritzbuck/uppmax/temp/181030_Chlorobi-genomes/"
wisc_path = pjoin(additional_MAGs, "bins_4_lakes_WI")
tsuji_path = pjoin(additional_MAGs, "Chlorobi_genome_bins_ELA_Tsuji")
sarahi_path = pjoin(additional_MAGs, "Garcia_compilation_fna")
data_folder = "/home/moritzbuck/people/0026_Chlorobi/"

wisconsin_MAGS = os.listdir(wisc_path)
Tsuji_MAGs = [l for l in os.listdir(tsuji_path) if l.endswith(".fna")]
Sarahi_AGs = os.listdir(sarahi_path)

Anoxicopedia_MAGs = os.listdir(anoxicpedia_folder)


chlorobi_nd_sister_tax = {k : v for k,v in taxonomy.gtdb.items() if v.is_child(taxonomy["c__Chlorobia"]) or v.is_child(taxonomy["c__Ignavibacteria"]) or v.is_child(taxonomy["c__UBA10030"]) or v.is_child(taxonomy["c__Kapabacteria"])}


RefSeqs = [k for k in chlorobi_nd_sister_tax.keys() if not "UBA" in k]
UBAs = [k for k in chlorobi_nd_sister_tax.keys() if "UBA" in k]

for g in tqdm(Anoxicopedia_MAGs):
    name ="Anoxic_" + g[:-4].replace("-","_")
    if Genome.find_one({'_name' : name}) == None :
        Genome.FromRawFile(name, pjoin(anoxicpedia_folder, g) , { 'anoxicpedia' : g[:-4]}, None)

for g in tqdm(RefSeqs) :
    name = g.replace("RS_","").replace("GB_","")
    if Genome.find_one({'_name' : name}) == None :
        Genome.FromRefSeq(name, taxonomy.gtdb[g].get_tax_string(full=True) , ncbi_taxon_str = None)

for g in tqdm(UBAs) :
    name = name
    if Genome.find_one({'_name' : name}) == None :
        Genome.FromUBAFile(pjoin(uba_folder,g + ".fsa"),name, taxonomy.gtdb[g].get_tax_string(full=True))

for g in tqdm(wisconsin_MAGS):
    name = "Wisc_" + "_".join(g.split(".")[:-1])
    if Genome.find_one({'_name' : name}) == None :
        Genome.FromRawFile(name, pjoin(wisc_path, g) , { 'wisc' : g[:-3]}, None)

for g in tqdm(Tsuji_MAGs):
    name = "Tsuji_" + g[:-4]
    if Genome.find_one({'_name' : name}) == None :
        Genome.FromRawFile(name, pjoin(tsuji_path, g) , { 'tsuji' : g[:-3]}, None)

for g in tqdm(Sarahi_AGs):
    name = ("Other_" + "_".join(g.split(".")[:-1]).replace("-","_")).replace(" ","_")[:35]
    if Genome.find_one({'_name' : name}) == None :
        Genome.FromRawFile(name, pjoin(sarahi_path, g) , { 'sarahi' : g[:-3]}, None)

sms = {g.name : g.checksum for g in Genome.find()}
dups = [s for s in set(sms.values()) if list(sms.values()).count(s) > 1]
dups = [[k for k,v in sms.items() if v == s and "Other" in k][0] for s in dups]

for d in dups:
    Genome.find_one({'_name' : d}).delete()

#prok = Prokka.make()
prok = Prokka.find_one()
for g in tqdm(Genome.find()):
    if len(g.proteins) == 0 :
        prok.run_remote(g)

for g in tqdm(Genome.find()):
    if len(g.proteins) == 0 :
        prok.retrieve(g)



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

os.system("silix -n {net} {fasta} {blast} > {out}".format(fasta = fasta_path, blast = blast_path, net = fasta_path[:-4] + ".silix.net" , out = fasta_path[:-4] + ".silix"))




chlorobi_cogs = COGing.FromFile("Chlorobi",pjoin(data_folder, "all_proteoms.silix"),"silix")

with open(pjoin(data_folder, "gtdbtk.bac120.summary.tsv")) as handle:
    taxi = {l.split()[0] : l.split()[1].split(";") for l in handle if not l.startswith("user_genome")}

taxi = {k : taxonomy[[vv for vv in v if not vv.endswith("__") and not vv.startswith("s__")][-1]] for k, v in taxi.items()}

for g in genomes:
    g.taxonomy['gtdbtk'] = taxi[str(g.id)].get_tax_string(full = True)
    g.save()

















bla

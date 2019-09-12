import os
import re
import sys
import ftplib
import os
from ftplib import FTP
import io
import gzip
from io import StringIO

with open("/home/moritz/dbs/gtdbtk/bac_metadata_r86.tsv") as handle:
    ids = [l.split('\t')[47] for l in handle]



def make_ftp_path(idd):
    ftp_addr = "ftp.ncbi.nlm.nih.gov/"
    ftp_path = "genomes/all/{letters}/{first}/{second}/{third}/"

    letters = idd.split("_")[0]
    numbers =  idd.split("_")[1].split(".")[0]
    first, second, third = numbers[0:3], numbers[3:6], numbers[6:9]
    ftp_path = ftp_path.format(letters = letters, first=first, second=second, third=third)
    return "ftp://" + ftp_addr + ftp_path

 with open("gca_locs.txt" , "w") as handle :
     handle.writelines([make_ftp_path(i)+"\n" for i in ids[1:] if i != 'none'])


ftp=FTP(ftp_addr)
ftp.login()
ftp.cwd(ftp_path)
ass_id = ftp.nlst()
ass_id = [a for a in ass_id if self.refseq_id in a]
assert len(ass_id) == 1
ass_id = ass_id[0]
ftp.cwd(ass_id)
all_files = ftp.nlst()
if type == "assembly" :
    assembly_file = [f for f in all_files if f.endswith(".fna.gz") and not "_cds_" in f and not "_rna_" in f]

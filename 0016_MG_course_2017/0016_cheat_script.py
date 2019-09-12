from pandas import DataFrame
import os,shutil
from os.path import join as pjoin


metadat = DataFrame.from_csv("../hmp_cart_one_patient_metadata.csv", sep="\t")
manif = DataFrame.from_csv("../hmp_cart_one_patient_manifest.csv", sep="\t")

all_libs = [f for f in os.listdir(".") if 'fastq' in f]

for l in manif.iterrows():
     short_id = l[1]['urls'].split("_")[1]
     print short_id
     for f in all_libs:
          if short_id in f:
             if os.path.exists(f):
                if not os.path.exists(l[1]['sample_id']):
                     os.makedirs(l[1]['sample_id'])
                if os.path.exists(pjoin(l[1]['sample_id'],f)):
                    os.remove(pjoin(l[1]['sample_id'],f))
                shutil.move(f, l[1]['sample_id'])

for f in set(metadat['visit_number']):
    for s in set(metadat['sample_body_site']):
        os.makedirs(os.path.join(s.replace(" ","_"),str(f).zfill(3)))

sample_dict = {}
for id,line in metadat.iterrows():
    if os.path.exists(id):
        shutil.move(id, os.path.join(line['sample_body_site'].replace(" ","_"),str(line['visit_number']).zfill(3)))
        sample_dict[id] = {}
        sample_dict[id]['path'] = os.path.join(line['sample_body_site'].replace(" ","_"),str(line['visit_number']).zfill(3), id)
        sample_dict[id]['site'] = line['sample_body_site']
        sample_dict[id]['visit'] =  line['visit_number']

DataFrame.from_dict(sample_dict).transpose().to_csv("clean_metadat.csv")

# A perturbed freshwater microbial loop's efficiency is unaffected by carbon pool in two mesocosm shading experiment

Change of carbon composition of boreal lakes is expected to vary signficanlty. However emission changes due to these changes are hard to evaluate. In 2015, we  setup two mesocosm experiment in Lake Erken (Sweden) to study two expected aspects of climate change: browning, and shading. 20 mesocosm were randomly assigned to 4 condition, control, shaded, DOC-amendment, and shadedOC-amendment. The two experiments differ in their source of shading, experiment A is shaded chemicaly with HuminFeed, experiment B is shaded with physical neting. Previous results from these experiments show that shading and DOC-amendments increase CO2 emmisions, this is probably linked to changes in primary producers, as bacterial productivity was unaffected. In this note we further explore the microbial aspect of this experiment by sequencing the 16 rRNA genes amplified from the microbial communities of the two experiments. These are used to compute community profiles that reveal significant changes in the microbial community under all treatments. These change are linked to the numerous biochemical and biological changes in the mesocosms. This shows inherent limitations to the efficiency of the microbial loop that are independent from carbon pool as well as other limiting factors, eventhough the real limiting factor of the microbial productivity cannot be accertained. Implying that the microbial loop can react easily to simple DOC changes to reach a systematic intrinsic maximum poductivity. 


## Introduction

## Methods 

The illumina MiSeq data was processed using mainly programs from the vsearch software [ref] a free version of the commonly used usearch tool [ref]. Reads where merged, followed by quality filtering (minimum length: 400, maximum expected error: 1), as well as dereplication and removal of singeltons, and finally clustering to 97%-OTUs. For taxonomical annotation the classifier implemented in dada2 [ref] based on methods described in [Wang et al. 2007].

Post-processing was done in R using a diversity of functions from the vegan [ref], phyloseq [ref] and DEseq2 [ref]. The multidimensional scaling was done using the metaMDS
All scripts are provideed in supl methods.



# Metagenomics Course, Fall 2017

## Learning goals :

We have  number of goals for this course:
* Understand and recognise basic datatypes involved in metagenomics (fasta, fastq, sam/bam, otu-tables, gff...)
* Understand the differences between the types of molecular data (reaads, paired-reads, contigs, genomes...)
* Basic usage of High Performance Computing (HPC) clusters, e.g. UPPMAX (module system, running tools interactively, running things in a queue, installing new tools...)
* Understand and replicate the presented workflows in metagenomics
* Familiarise yourself with manuals and online ressources


## UPPMAX :

All the exercices here will be done using UPPMAX (https://www.uppmax.uu.se/) more specifically the Rackham cluster.

The rackham-cluster of uppmax is a distributed HPC cluster, meaning it is a lot of powerful computers networked together. It has 4 login-nodes and many more computing nodes.

### SSH

We will use `ssh` (Secure SHell) to connect to this cluster, use you ssh-client of choice. To start, we have to log into the login node.

```bash
ssh -X MY_USERNAME@rackham.uppmax.uu.se
```

The first time you do this (we should not be the case this time) there might be a few questions, but eventually you will be asked for your password, enter it. When you type your password nothing will be shown on the screen. This a security feature, not a bug.

We will not spend much time on the login node. we just want to know which computer is booked for each of us.

> Optional question: what does the `-X` in the command do.

Rackham, as most UPPMAX computer systems, runs a queuing tool called SLURM (Simple Linux Utility for Resource Management, https://slurm.schedmd.com/). This tool is used to distribute computing ressources to all users. A computer has been booked for each participant, SLURM knows 'it's name'.

On Rackham, to know which ressources you have requested (or someone else requested for you) the `jobinfo` command is used.

> Use the jobinfo command to find the name of your dedicated computer

Once you have the name of the computer (something looking like r123), you can log directly onto that computer. That name is the name of the computer in the local network of rackham, so you can simply connect with :

```bash
ssh -X r123
```

You are now on the computer that you can call `$HOME` for the next 3 days.

As you might/will notice, the wifi connection in this room is a bit 'dodgy' sometimes, `ssh` does not like bad connections, when the connection is broken while a program is running the programs stop, which can be very annoying for long computations. To alleviate this we will use a tool called `screen`. This will start a terminal that does not stop running when the connection is broken or the terminal is closed.

To start `screen` you just run the command:

```bash
screen
```

You can leave the 'screen' without closing it by doing `ctrl-A ctrl-D` (careful doing `ctrl-D` closes the 'screen' completely). With this you are back to the normal terminal. To reconnect to your screen, you first need to find out the name of the screen, `screen -list` will give you a list of screens available on this computer (yes, you can have many screens), and to reconnect simply run `screen -r NAME_OF_SCREEN` or if there is only one detached screen `screen -r`.

Now that we are comfortably in our screen and nothing can happen to anything running lets start.

The data we are gonna use for this course is in our common project folder, e.g. `/proj/g2017026/2017_MG_course/`), more specifically in  `/proj/g2017026/2017_MG_course/raw_data`.

> Using the `ls`-  or even better `tree`-function
> How many fastq-files can you find in this folder, and list all of them?!

> INTERNET INTERUPTIONNNNNNNNNN

Let's spontaneously break our internet connection. Your terminal to uppmax should be unresponsive now ....

> Reconnect to your `screen`

### The module system and other software prerequisites

Rackham uses a module system for programs (http://www.uppmax.uu.se/resources/software/installed-software/). So many bioinformatic programs are available without installing anything. We will use this as much as possible during this course. You can list all available module with `module list`. Let's load the `bioinfo-tools`-module, which is a prerequisite for most of the tools we will use in the next couple of days:

```bash
module load bioinfo-tools
```

Additionally, we will need for a few of our tools python. Let's load python3 whith the module system (if you have your own install of python you should be able to keep it):

```bash
module load python
```

Python is a scripting language, we will use some scripts made in it and other tools need it too. Modules are loaded on Rackham until you disconnect from your terminal, so one you loaded them and your `screen` is still alive you do not need to rerun these commands. But if you disconnect (e.g. you see `screen terminated`), they will not be loaded anymore, so you will need to reload on you next connect, or when you make a new screen!

Python has it's own package management system which is used to install more libraries. We will need it at some point too. It is called `pip`. Let's install a very usefull python package:

```bash
pip3 install --user biopython
```

## Data-types

### FASTQ-files

Sequencing facility might provide data from different technologies in a variety of formats, but most probably one of the formats the data will be delivered in will be the FASTQ format. We will start this tutorial with data in this format. It is a wide text-based format for sequences, they are related to the FASTA format which will see too.

A FASTQ-file can contain many many sequences, each sequence is represented by two data-types, the actual sequence and the quality scores:

```
@an_example_sequence_in_a_fastq_file
ACACATATATACACACATACACACACATATACACACACATATAAACACATATATACATTTATATGCATATATTAATACATATATATTTAAGTTGATGGAGAGTATAACAGAGTTAGGCTGCTTATT
+
BBBBBFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFBFFF/FFFBFFBFFFFFFFFFFFFFFFFFFFFFFBFFFF//FFFFFFFF<FFFF/FFFFFB/FFF<BF//<</BFFFFF
```

Normally each entry in a FASTQ-file is composed of exactly four lines.

The quality scores are related the quality of sequencing. They are `Q = -10 log<sub>10</sub> p` where `p` is the probability of a wrong base-call, and each letter corresponds to an integer value of Q. These letter can differ a bit depending on the sequencing platform, for recent illumina we have this correspondence table (courtesy of wikipedia):

```
  BCDEFGHIJKLMNOPQRSTUVWXYZ[\]^_`abcdefghi
   |     |                               |
   3.....9...............................41
 ```


 As mentioned, the raw-data we will use is available in the `/proj/g2017026/raw_data`. This is publicly available data from a number of microbial samples from a single subject of the human microbiome project (https://hmpdacc.org/).

 > Use the linux command  `ls`, `wc`, and `du` :
 >
 > What is the content of this folder?
 >
 > Optional : What is the difference between the files in the two folders?

 I preprocessed some of these files a little bit already, so now each file corresponds to one sample (I simply concatenated some files together as they where from the same samples). For each samples (library) we have two files, this is due to the way the illumina technology sequences, it would be different for PacBio for example.

 > Using `wc`, the bash redirect (e.g. `wc -l my_lib.fastq > txt_file `, for example), combined with your data manipulation tool of choice (`R`, `python`, `bash`, `excel`, `pen&paper`)
 >
 > Make a table with the number of reads for each library
 >
 > Optional : Is there a relationship between number of reads and sampling-site or sampling-time (the number in the file is a the number of the visit)
> Very optional : Anything to say about the diversity of the reads in the libraries?

A very common tool to analyze reads is `fastqc`.

> load `fastqc` using UPPMAX' `module load`
>
> run `fastqc` for one of the shotgun and one of the amplicon `libraries`
>
> What can you learn from the returned report?!

Do not worry to much about fastqc reports in general, metagenomic data has a tendency to fail many things, as fastqc is made for nice clean genomes, so many of the warnings/errors will be due to actual biological effects or PCR-artefacts.

### FASTA-files

FASTA-files are the mother of all sequence data-types. It is a very simply text base file type compose of one or multiple entries. Each entry is composed of two type parts: A description part, which is a single line starting with `>` containing text information about the following sequence; and a data-part containing sequence data as plain text, this can be spread over multiple lines. FASTQ-files can be easily transformed to FASTA-files by removing the qualities and the `+`-line and replacing `@` by `>`, e.g:

```
>an_example_sequence_in_a_fasta_file
ACACATATATACACACATACACACACATATACACACACATATAAACACATATATACATTT
ATATGCATATATTAATACATATATATTTAAGTTGATGGAGAGTATAACAGAGTTAGGCTG
CTTATT
```

### SAM/BAM-files

SAM is a text-based file type used to represent alignments, usually used for NGS-mappers, when one aligns many short reads to a reference.

It is normally composed of a header indicated by lines starting with the `@` character (which also indicated comment lines, e.g. information that does not correspond to the SAM-format), followed by often many many many lines, one per sequence to be aligned. Each line is split into columns with tabulations and each valuer corresponds to some info. From wikipedia:

| Col | Field | Type | Brief Description |
| --- | ----- | ---- | ----------------- |
|1|QNAME|String|Query template NAME|
|2|FLAG|Int|bitwise FLAG |
|3|RNAME|String|References sequence NAME|
|4|POS|Int|1- based leftmost mapping POSition|
|5|MAPQ|Int|MAPping Quality|
|6|CIGAR|String|CIGAR String|
|7|RNEXT|String|Ref. name of the mate/next read|
|8|PNEXT |Int|Position of the mate/next read|
|9|TLEN|Int|observed Template LENgth|
|10|SEQ|String|segment SEQuence|
|11|QUAL|String|ASCII of Phred-scaled base QUALity+33|

or for meore detail : [SAM/BAM format specification](https://samtools.github.io/hts-specs/SAMv1.pdf).

These files can be massive so they are often compressed using `samtools` (available trough `module` on UPPMAX), and are then often called BAM-files. An example :

```
2167:2:1101:11949:2366  16      contig_99_631734_length_1867_multi_8_in_2_out_1 1493    42      180M    *       0       0       ACCTTACTAGTCGCTGTTGCGAAAGATAATTTTTAGTTGCTAATAATCGGAGCCCTGCATGGCCATCACTTTATACCCTGTTAAACCAGATTTTATGGCTGAAATCGGTGATGTCGACTTGAGTATGCCGCTCTCTGACGAAGACCGCGCCGCCATAAAAGCCGCTTTCTTTAAATACAG    bbbeeeeegggggiigfhiiifhhiiihiiiiiiiiiiiiiiiiiiefhiefhihiiggggfeeeeebbbdddddccbcccccbccccaaccbb`bccccccab`bb_^^^bc_^`^bba`bbbaaaa`aa]]`bc_^ZZb`^]edfffiiiiiiiiiiiiiiiiiigggggeeeecbbb    AS:i:0  XN:i:0  XM:i:0  XO:i:0  XG:i:0  NM:i:0  MD:Z:180        YT:Z:UU
2167:2:1101:11897:2252  4       *       0       0       *       *       0       0       TTTAAATTTCTACATACATACTTTTCTAAATAAACCTGGACGCAATTTAAAAATTTCAAATTAATTTTATTTAATAAAAGCTTTTATTAATTTACGAAAGCGTCTAAAATACTAAAAAAAATACAAAAATATTTAAACTGTGTCTTCATTTAAAAAATCCAAATCAGTAAAATGTGAAATTGTATTTTTTTT        ^__ec`ccgggcefg_cdhhaegggdghhhfcb^gae^eghdf#af`cfd^a_eeeehbbbeghae`fccdgfbgdgedb_bbcbddb_cacd`bac_[ca`aadcdcddb`bcaabdeccccb`ZZ^gg`bbbgfgggebefggffffhgffbfccdhgebe_ihfhhhfhfgahhgf]gcgceceec__^        YT:Z:UU
2167:2:1101:10183:2426  0       contig_99_3998773_length_1078_multi_3_in_1_out_1        849     42      212M    *       0       0       AACTAACCTTAACAAGCCTTAGCAAGATGTATTTTTGCTTTGTAATAAGCTAGAGTGTTTCTTCCTGGCAAACTTTTTCTAACTAGCTATGACTTCGAAGCTTGTTTGGTTCTAAAAACAATGTAAGCATCGTAGCATTTAAAGAACCAATTGTCAATTTTACTACCTCTTAAAAAACAATTAAAGCATCGTAGCATTTAAAGAACCAATTG    bbbeeeeegggfgiiiihhiiiifhiiiihiiiiiihiiihfgihiiihihf^cgghiiiiiiiiiihigghfhhhiggggggeeeeeeddcdddccbb`acccccccb``ccddccccddcccdcbcccceeeeeggggeggcc\hgfhhihhiihiiiihgfahfageeiiihgeihiiiiiiiifgiiiiiiiiihgggggeeeeebbb    AS:i:0  XN:i:0  XM:i:0  XO:i:0  XG:i:0  NM:i:0  MD:Z:212        YT:Z:UU
```

### GFF-files

GFF-files are text-files that contain annotation information for sequence data. Basically these are tab-delimited tables, where each line corresponds to an annotation. You can get more information [here](https://www.ensembl.org/info/website/upload/gff.html).

## AMPLICON-ANALYSES

### Tools

General amplicon processing and more:
* QIIME
* MOTHUR
* usearch (used in this course, mostly because of convinience)
* dada2 (bioconductor)

Quantification:
* phyloseq (bioconductor)
* metagenomeseq (bioconductor)
* DEseq2 (bioconductor)

Other:
* PiCRUST
* Phylotyping

### Processing the amplicon

Most analyses of 16s rRNA gene amplicons rely on an OTU-table, e.g. a table of counts of microbial-types in your samples. For this we need to cluster the reads into 'OTUs' (Operational Taxonomical Units), these are then counted.

We will use the `usearch` software to do the necessary computations (https://www.drive5.com/usearch/manual/). We are going to be fancy and we are going to use the latest version of `usearch` which is not on UPPMAX. We have it however locally in our project folder in `/proj/g2017026/bin/` (as well as some other software we will use later). For ease of use we will put this folder in our `PATH` environment variable (this contains all the folders which have executables readily available).

```bash
PATH=/proj/g2017026/bin/:$PATH
```

> Optional : which other folders are in the `PATH`

From now one for a few exercise you will pic one of the amplicon libraries to process. We will call it MY_LIBRARY for the next few questions!

#### Merging reads

As we saw, each library has pairs of reads. In the case of a well designed amplicon the PCR product should be shorter than the lengths of the two reads added, in which case they will overlap, and the two pairs can be merged into one longer pseudo-read which will improve the quality of the analyses. This is the case for this dataset.

> Using the `usearch10` program with the `-fastq_mergepairs` option
>
> Merge the pairs of MY_LIBRARY
>
> Optional : play with the parameters to increase/decrease merging rate

<!-- One sample is nice and fine, but for a good analysis of the data we will need to analyze all the reads as one. If we just concatenate the files we will lose the information of which samples the reads came for. Luckily `usearch` can make our life a bit easier!

> Using the `--relabel` and `-fastq_mergepairs` options of `usearch10`
>
> Merge all the libraries.
> How many reads got lost in the merging
>
> Optional : any relationship between the merging rate and some of the sampling-parameters -->

#### QC-ing the reads

The reads have now been merged into pseudo-reads, but we have not done any quality-check yes. For merged paired reads it is convenient to do it after merging, as the merging can increase the quality of the overlapping bases! In all sequencing runs there will be some bad reads, this can be due to many reasons, but usually mostly to the quality of the DNA, e.g. extraction and PCRs, or interactions between spots on the sequencing array. We will use the quality scores from the fastq files to remove reads below as certain quality.

> Using the `-fastq_filter` option of `usearch10`
>
> Remove from MY_LIBRARY the reads with more than one expected error and output it as a FASTA-file
>
> Optional: How does the loss of reads vary with the expected error.


#### Dereplicating

All clustering algorithms can be sensitive to the number of reads, however amplicon data has huge amounts of duplicates due to the PCRs, so one way to help the clustering is to remove these duplicates. Also, as there it is very reasonable to expect duplicates, any reads that appear only once is probably due to either PCR- or sequencing errors. In `usearch` both of these steps are done together.

> Using the `-fastx_uniques` option of `usearch10`
>
> Dereplicate MY_LIBRARY and remove singletons

> Optional: How does the count distribution of the reads look like?

#### Clustering

We are finally ready to do the clustering. The clustering heuristic of `usearch` is pretty fast and incorporates chimera detection as well. Please note that one of the keys of the heuristic is that the reads are aligned, as they are in an illumina amplicon, but might not be in other cases where you'd like to run this tool for!

> Using the `-cluster_otus` option of `usearch10`

> Compute 97% similarity sequence-clusters in MY_LIBRARY

> Optional: vary the identity cutoff, also what happens if you run the clustering without removing the singletons in the previous step

#### Mapping back

Now the whole data of MY_LIBRARY has been reduced to a small number of sequence-clusters (OTUs). Each of them represented by the most abundant read of that cluster. However the information of the total amount of reads has been lost along the way. For this reason we need to map the reads back to the OTUs, which means we need to find for each read the OTU that matches it best and that is more similar then the identity threshold. This will also allow us to recruit back reads lost during QC!

> Using the `-usearch_global` option of `usearch`
> Get some help from : https://www.drive5.com/usearch/manual/cmd_otutab.html
> Map the merged reads to the representative sequences of the clusters, and make an otu-table

> Optional : what are the reads that did not map?

#### Putting it all together

Well, this is all nice and fine,  but this is only one sample, so not really a table, just a vector. We could run the same thing for every sample separately but we would just get many separate vectors, we could not connect the OTUs easily to each other. What we need to do is to run the clustering step with all the reads of all the sample at the same time!

> Using all the options of `usearch10` that you have used until now,
> Have a look at the `-relabel` option of the `-fastq_mergepairs` function

> Make an OTU-table with all the samples, e.g. run the merging, QC for every sample separately, concatenate the resulting FASTA-files, and follow then with the rest of the script.

> Optional : quantify the read-losses. Are there any experimental factors influencing the loss of reads?

#### New Age Clustering

In recent years there have been many discussions about OTUs and there meaning and the issues linked to them. Hence new modern methods of clustering amplicon reads have been developed. Much talked about recently has been the `dada2`-tool running with `R`, this tool does not use sequence identity to cluster reads, but build a noise model to clean reads from variation that is not biologically relevant. `usearch` also released recently a similar heuristic, the `unoise` algorithm.

> Using the `-unoise3` option of `usearch10`

> redo the OTU-table

> Optional : compare the mapping rates between the two methods


#### Taxonomical Annotation

We have not a table with OTUs and their abundances in each samples, as well as representative reads for all the OTUs. To learn more about our data we want also to annotate the reads. In the case of amplicon, we want to annotate for taxonomy, as we normally know exactly what genes we already have. For 16s rRNA amplicon there are many tools and databases. For sake of simplicity we will use a tool also provided by the `usearch` software, the `sintax`-tool.

> Let's split the class in two using the two OTU-tables

> Using the `-sintax` option of `usearch10`
> and the taxonomy/16s database found in /proj/g2017026/data/rdp_16s_v16_sp.udb

> Give all the OTUs taxonomical classifications.

> Optional : play with the cutoff, and try to install/use an other database.

### Let's make some plots : or where the real game starts!

We won't really have time to dwelve into diverse quantification methods and statistaical analyses, that would need many more prerequisites and much more time.

However, often in bioinformatics you get to use other people scripts to do a variety of tasks, they might come in a variety of languages. In this case I provide you with an R-script to make a little plot with this! It is called `ampliconplot.R` and you can find it in the `scripts`-folder.
`R` is already install on rackham so to run an `R` script, you normally just need to do `Rscript name_of_script.R`.

> Have a look at the `R`-script
> Try to run it!

You might notice that it complains about an abscent package. This is common when running scripts of other people. However, installing scripts in R is pretty easy. The easiest is to open `R` in the terminal, e.g. type `R`. You will see a prompt, and now run the command:

```R
install.packages(“reshape2”) # installs the package
quit() # leaves R
```

Answer all the questions with yes.

> Try running the script again
> Add any other package possibly needed.
> Make sure the filenames match!
> Edit the script if needed, but make a personal copy first!
>
> Optional : Why did I make separate plots for nasal and fecal data?
>            Make barstacks with other taxonomical levels


## Shotgun Data

Shotgun amplicon data is much more diverse and versatile than amplicon data. For this tutrial, we will run together a number of more or less complicated workflows using the shotgun libraries we have in store for us on `/proj/g2017026/raw_data/shotgun`.

As mentioned before, the principal of shotgun metagenomics is pretty straight forward. Take a sample, extract all the DNA, sequence the hell out of it. This leaves us often with a lot of data and many things we might want to know.

### Pre-processing

Some analyses are general to whatever biological question, they are linked to general NGS-methods.

#### Some Quality Control

Getting good libraries for shotgun sequencing can be difficult as the samples used normally are not as friendly as other types of biological samples (I am looking at you people doing cell-cultures). This can often result in lower quality of sequencing. A common tool used to judge quality of NGS-data is `FastQC`. This is only an analysis-tool, and only generate reports, it will not actually do any cleaning, but it might give us some insight into the data, and possible issues.

> Let's each first pic a sample of choice

> Using the `module`-system of UPPMAX load `FastQC`

> Generate a `FastQC`-report for your library of choice

> Optional : try to understand the FastQC-reports, and why we get so many warnings/errors

#### Read trimming and QC

[MultiQC]

Once we have some idea of problematic libraries, and potential issues, we can use some tools to clean the libraries a bit reduce noise and increase quality

A very commonly used tool for this is `trimmomatic`, it can do a large variety of operations on reads.

> Load `trimmomatic` with UPMMAX's `module`-system

> Clean the reads!

> Optional : re-run `FastQC` on the cleaned libraries and discuss


### Read-based analysis

#### Minhash-ing!

Raw libraries are large and clunky, and hard to handle really, mostly if we have many. Also they can be very complicated to compare directly. However, computer scientist are clever apparently, and developed methods to comapare large amount of text, which turn out work pretty well for sequence data as well. For texts, it is based on words but in sequences it is based on k-mers, meaning 'nucleotide-words' of length 'k'. It uses these words in hash-functions  (http://mccormickml.com/2015/06/12/minhash-tutorial-with-python-code/) to compute a very-compressed vector representation of your sequence data that can be used to compute distances!

The `Mash` tool has been developed to compute MinHashes and MinHash distances for genomic data, interestingly it can also be used for other things than reads. But here we will use it for our metagenomic libraries!

> Use the `mash` command

> Compute the MinHash set for your library
> Exchange these with other students so you have the set for each library
> Compute the MinHash distance between the libraries

> Optional : run MinHash on some more libraries in  `data/additional_libraries`

Right now we just used `mash` to compute the Hash-vector, but if we want to use these we actually need to compare them to each other to get some distances.

> Let’s share all of our hashes, so copy them all to the `common/mash_hashes/`-folder but check first if the hash you computed is already there or not.

> run this little code snippet to compute the distances and make a table from them
```bash
ls common/mash_hashes/*.msh > mashes.txt
for f in `cat mashes.txt`;
do
  mash dist -p 20 -l $f mashes.txt
done > mash_distances.txt
```

Here again I provide you with an R-script to make a some little plots with these distances! It is called `mashplot.R` and you can find it in the `scripts`-folder.

> Have a look at the `R`-script
> Try to run it!

> Optional: change the plot so you spot better the differences between the forward and the reverse

#### Taxonomic annotation

We have all these reads now, but who are they from?! This is much more complicated then for 16s amplicons. Many tools have been developed, and most have the problem that they are only as good as the database they use, and as raw reads can be very numerous, often databases are limited to make computations tractable.

One recent tool that deserves mention is `kaiju`. The interesting aspect is that it uses a `blast`-index to do the classification, this is reasonably compact and easy to compute meaning that a much more exhaustive database can be used to classify as opposed to other methods.

Most of these tools are pretty slow, so we will at first subset our libraries.

> Using the `head`-linux function
>
> Subset your library of choice to 500000 reads

> Check the [manual](https://github.com/bioinformatics-centre/kaiju), jump straight to the ‘Running Kaiju’ section
> Using `kaiju` and the `krona` tools (available through `module`)
>
> Make a Krona plot of the taxonomies of your library of choice using with the database found in `data/kaijudb/proGenome’

> Optional : Compare the taxonomy results with the 16s rRNA data of the same sample
> Super-optional : Run kaiju on the the other database in `data/kaijudb`


#### Functional annotation

What is meant by functional annotation, name the reads base on known genes, or functional ontogenies. For example, assign them PFAMs or EC-numbers. Functional annotation for reads is often mapping-based, meaning the reads are aligned to a reference database which has 'known'-functional annotations.

We will use a tool called HUMAnN2 to annotate our reads.

> Use the `pip install` option for [installing](https://bitbucket.org/biobakery/humann2/wiki/Home#markdown-header-initial-installation) `HUMAnN2`, you will need to add the `--user` option as you do not have `sudo` permissions, and to add `~/.local/bin/` to your $PATH.

To set it up properly:

```
humann2_config --update database_folders protein /crex2/proj/sllstore2017039/2017_MG_course/data/humann2/uniref
humann2_config --update database_folders nucleotide /crex2/proj/sllstore2017039/2017_MG_course/data/humann2/chocophlan

mpa_dir=/proj/g2017026/2017_MG_course/share/metaphlan/
PATH=$PATH:/proj/g2017026/2017_MG_course/share/metaphlan/
```

> Using HUMAnN2

> Annotate the subset of your library of choice!

> What proportion of your data has been actually annotated?!

> Optional : run `HUMAnN2` on your whole library in an other `screen` using half of your processors

#### Compiling the data

The taxonomic and functional annotations here have now been done for each sample separately, however if we want to compare the samples, we need to compile these to be able to compare them.

> Using your data-processing tools of choice

> Share the runs with each other and compile it into two tables, one for taxonomy and one for functions.

> Optional : Use some multi-dimensional scaling to visualize the distance between the samples, compare it to the results of Mash, and of the same samples in the 16s amplicon

### Assembly

The other road for analysis of metagenomic data is to assemble the reads into longer contiguous pieces of DNA (contigs). This will simplify annotation of the sequences with the downside of risking chimeric sequences (e.g. assembling things that do not belong to each other). Assembling (meta)genomes can take numerous steps and numerous tries. The usual way is to perform a variety of assemblies using different tools and combinations of samples, and to evaluate them to find the ones that seam more useful/good.

There are no one solution fits all, but we will try to look at some of the way to do and to judge assemblies.

For the purpose of this tutorial we will only use the `megahit`-assembler. This is mostly due to practical considerations, it is both fast and memory-efficient, additionally it actually is a very good assembler. However it does not necessarily produce the best assembly, so for your personal project consider using other assemblers on as well as `megahit` (for example `metaSPAdes`).

#### Assemble

We will start with simple assemblies of our libraries, meaning that we consider our libraries all independent and assemble them as such.

> Using the `megahit`-assembler

> Assemble your (full) library of choice

> Optional : How would you run a different assembly?


#### QC the assembly

A variety of tools is available to QC genomic assemblies, however many of them are appropriate for metagenomic data, as they presuppose a single organism. However for most metagenomic applications, a few factors will already be indicative of the quality of the assembly.

One if the first factors is the amount of reads that align to the assembly. The proportion of reads mapping our assembly gives us a good idea of how representative of the sampled environment the assembly is.

To get this information we will need a mapping tool. We will use `BBmap`. Mappers, like `BBmap`, don't only return the proportion of mapping reads, but also the specific alignments, these are output as SAM-files, a particularly horrible and wasteful file type.

> Using `bbmap.sh` loaded with `module`, use the `bamscript=`-option
>
> Map your library to your assembly. What proportion of your reads map to your assembly? Any comments about the size of the SAM-file?
>
> Optional : map subsets of the other samples to your assembly, how does it look like? (checkout the `reads` option, and also the  `showprogress` one , mostly because I like it)

SAM-files are HUUUUUGE, we do not like that, BAM-files, the compressed version of SAM is a bit better.

> Use the script written by the `bamscript=` option of `bbmap.sh`.
>
> What are the generated files?
>
> Copy the BAM-file and its index to the appropriate place in the `common`-folder

However this does not say anything about the quality of the `contigs` directly, only about how much of the data has been assembled, but it could all just have assembled into very small contigs!

We need to look at how the contigs "look-like". A first simple way of looking at this is to simply do some stats on contig lengths. There are many ways probably to do this, but we will use a script that has been written by yours truly. We will use this to learn a bit about scripting and packages.

You can find the script in [script-folder], this is a very simple python script, make a copy in your working directory.

To run a python script, normally you'd just have to do:

 ```bash
 python my_script.py
 ```

 Let's try to run the script you just obtained from a friendly collaborator.

 > run the script

Hmmm, you probably got an error. Many scripting languages have systems to load additional packages of functions. The system used in python is called `pip`. So you can install packages for `python` doing:

```bash
pip install my_package
```

However this normally needs superuser permission, so you will have to run this with the `--user` option, which will install the package into your `home`

> Using `pip`

> install the missing packages and run the script, use it to get some basic info about your assembly

> Optional : Personalize the script a bit. And run the script with all the assemblies.

We have now some basic data about our few simple assemblies.

> Any particular feelings about the obtained assemblies?

Some of the assemblies as you can see have decent mapping rate, however, the contigs composing them are not very good. Often, mostly if sequencing depth, things assemble well, but the lack of coverage breaks the assembly in many places. Also certain assemblers output also very small contigs which are often not usable for downstream analyses, but give an impression of high-mapping rates. Hence, assemblies are often length filtered. I am sure there are plenty of tools to do that, we however provided a home made script again, it also cleans up the names of the contigs a bit.

> Using the `fasta_filter.py` script

> Filter out small contigs from your assembly

> Optional : ammend the script to give personalised names to your contigs

Now, we just have to see how the characteristics of this 'new'-assembly are.

> Re-compute mapping rate and metagenome-characteristics


#### Alternate assembly.

So some assemblies are OK, some are not ... we can now make some other assemblies. We could normalize the reads to remove micro-diversity and errors. We could split the reads into different subsets (e.g. eukaryotic and prokarotic for example). Or coassemble multiple libraries together.

> Discuss some strategies, and make an other assembly.
> Run the mapping and the `python`-script once it is done

### Post-assembly strategies

#### Gene-atlas approach

This is an approach commonly used when assembly is difficult, and we do not want actual genomes. In this case we simply use the contigs for gene prediction, and afterwards use the predicted genes as reference for diverse 'omics analyses. One of the upsides of these approaches is that they allow to merge many assemblies. This approach has been used in the  [Tara ocean project] (http://ocean-microbiome.embl.de)

> Simply using the linux `cat` command

> Collect all the assemblies you like from your comrades, and concatenate them!

We now have a big assembly with many contigs! This is our raw merged assembly. The next step will be to predict genes on this. Gene prediction for microbes is not particularly hard as they do not have many intron! We will use the classic tool called `prodigal` for this purpose, fundamentally it will just look for long open reading frames starting with a `ATG`, but it has a few more tricks in it, anyhow.

> Using `prodigal`

> Predict genes on your concatenated metagenome. Output them as nucleotide FASTA-file

> Optional :

Now we have a large FASTA-file with all the genes in all the concatenated metagenomes! But obviously there is redundancy in this. The next step will be to remove it. We will cluster all the genes in this gene-assembly to remove highly similar sequences. The tool commonly used to fastly cluster a large amount of genes is `CD-hit`.

> Using `cd-hit-est`

> Cluster the gene-assembly, and assess the mapping rate of some samples on the resulting assembly!

> Optional : explore the similarity threshold

#### Binning Approaches

If we are more interested in extracting genomes and having a stronger genomic link in our analysis we can use a binning approach with our assembly. For more information see Anders' lecture. To extract genomic bins, or Metagenome Assembled Genomes (MAGs), or Genomes From Metagenomes (GFMs), we will use the `MetaBAT` tool.

This tool will need sorted BAM-files to run.

> Pick a genomic-assembly. Collect sorted indexed BAM-files corresponding to that assembly in a folder. (Remember the `bamscript=` option of `bbmap`)

> Using `runMetaBat.sh`, loaded with `module load MetaBaT`. Have a look at the FASTA-files generated.

> Optional : run the binning on an other assembly!

##### Evaluating bins

Now we each hopefully got a number of bins. Meaning that the assembly has been split into bags of contigs that somehow look/vary similarly. This process is however far from perfect, some bins will be-merged MAGs, others will contain viral data, or Eukaryotic chromosome, or simply trash.

We need to evaluate the bins for that. A very simple first filter is to remove bins that are too small or too large.

> compute rough size estimates of your bins (hint, the file size correlate well with the genome size), and pic a bin that has a size more or less corresponding to a microbial genome.

> optional : compute other stats on these bins, let your imagination run wild (or maybe use the script I provided earlier)

Now we should all have a bin we hope is a genome! To check this we will use a tool called `checkm`, fundamentally it uses single copy marker genes to check if a bin/MAG is complete, and if it is contaminated.

> Install and get `checkm` to work. You will need once again `pip install --user`
> It might ask at some point for a path with the data, I have already downloaded it, just put `data/checkm/`

Once it seams like it is running we want to run the `lineage_wf`, the documentation of `checkm` is a bit more confusing than other tools... So I advice you to look at [this specific page](https://github.com/Ecogenomics/CheckM/wiki/Quick-Start) and don't get lost in the others!

> Using the `lineage_wf` of `checkm`

> Compute completeness/contamination for your MAG of choice. Did you get any bonus information?

> Optional : compute this for more MAGs, and try the [taxon_wf](https://github.com/Ecogenomics/CheckM/wiki/Workflows)

> Collect all 'good' MAGs into a shared folder

##### Bin annotation

Now we should have a collection of MAGs that we can further analyze. The first step is to predict genes again, as right now we only have raw genomic sequences. We will use a different tool this time, one of my all-time-favorites : `prokka`.

This tool does gene prediction as well as some quiet good annotations, and is actually quiet easy to run!

> Use `prokka` loaded with `module`

> Predict genes and annotate your MAG!

> Optional : use `prokka` on some of the bins that did not pass the previous quality checks!

`prokka` produces a number of output files that all kind of represent similar things. Mostly variants of FASTA-files, one with the genome again, one with the predicted proteins, one with the genes of the predicted proteins. Also it renames all the sequence with nicer IDs! Additionally a very useful file generated is a GFF-file, which gives more information about the annotations then just the names you can see in the FASTA-files.

> Using the linux `grep`-command
>
> look for all [EC-numbers](https://en.wikipedia.org/wiki/Enzyme_Commission_number) in the GFF-file

> Optional : use [the KEGG pathway mapper](http://www.genome.jp/kegg/tool/map_pathway1.html) to analyze the metabolism of your organism, or look for your favorite genes in the annotations.

We now know more about the genes your MAG contains, however we do not really know who we have?! `checkm` might have given us an indication but it is only approximative.

Taxonomical classification for full genomes is not always easy for MAGs, often the 16S gene is missing as it assembles badly, and which other genes to use to for taxonomy is not always evident. One option is to take many genes and to make a tree.

One tool implementing such a thing is `phylophlan`, it has a reference tree based on a concatenation of reference genes, tries to find homologues of these genes in your bin, and then to fit  them into the reference tree.

> Install `phylophlan`

> Use it to identify your MAG(s)

##### A bit of phylogenomics

This readme is written by me, Moritz Buck, and should be contained in the deliverable folder.

The deliverable folder contains a folder for each SAG, and some additional files including scripts.

The script "this_script.sh" contains the whole workflow which is roughly:

* fastqc generating quality reports
* read cleaning (with trimmomatic)
* assembly with spades
* generating an assembly QC-plot
* filtering assembly by length (10kb in this case)
* computing completness (with checkm)
* annotating with prokka
* phylogenetically placing with phylophlan

The main data is in each SAG-folder, they contain a number of files:

* all_scaffolds : file with all unfiltered scaffolds
* assembly_plot :  a plot of coverage vs. length which helps to see if there are any massive contaminations. the vertical red lines correspong to certain cut-offs I typically use (1kb, 2.5kb or 10kb)
* genome_over10kb : the length filtered assembly, it is used for all subsequent analyses, I have simply filtered out all contigs shorter than 10kb.
* annotation : a gff file describing the annotation done by prokka, e.g. position and source of annotations, names, etc.
* proteins : predicted protein sequences

Additional data directly in this folder:

* README.txt : this little piece of poetry
* this_script.sh : script containing the whole workflow
* python_hjelper_script.py : a script containing all the python that has been used for various things
* all_completnesses.txt : the result from checkm (e.g. completness, contamination) for the assemblies (filtered at 1kb, 2.5kb and 10kb)
* only_10kb_completnesses.txt : the result but only for the 10kb filtered ones, e.g. the ones I kept
* phylophlan.tree : a tree placing the 10kb filtered SAGs in the microbial tree, they are prefixed by "SAG_" in the the tree is you search through it.

Additional data can be obtained, all intermediates are kept however they are pretty big so I will only send on request.

For any additional questions do no hesitate to contact me at moritz.buck@ebc.uu.se

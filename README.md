# Scripts for the data analysis performed in [“Evolutionary responses to codon usage of horizontally transferred genes in *Pseudomonas aeruginosa*: gene retention, amelioration and compensatory evolution”](https://doi.org/10.1099/mgen.0.000587)

If you have specific questions, please contact me on my correspondance address at Ghent University.

The code in this repository was used to perform the bioinformatics pipeline to infer horizontal gene transfers in a set of full genome sequences from multiple strains of a bacterial species. The general outline is presented in the figure below.

These files contain the code for specific parts of the analysis:
* Obtain_data.md: get the genome data to run the analysis
* Pan_genome_analysis.md: extract core and accessory genes (block A in the figure below)
* Strain_phylogeny.md: create a phylogenetic tree for all the strains (block B in the figure below)
* Filter_accessory_genes.md: remove "pseudogenes" from the accessory genes dataset (block C in the figure below)
* Infer_HGT.md: infer HGT's based on ancestral state reconstruction (block D in the figure below)

![Pipeline](https://github.com/MartijnCallens/P.aeruginosa_HGTCodonUse/blob/main/HGT_pipeline.png)

The analysis requires the installation of the following programs/packages:

[NCBI EDirect](https://www.ncbi.nlm.nih.gov/books/NBK179288/)
[PROKKA](https://github.com/tseemann/prokka)
[tRNAscan-SE](http://lowelab.ucsc.edu/tRNAscan-SE/)
[NCBI-blast+](https://www.ncbi.nlm.nih.gov/books/NBK569861/)
[SiLiX](https://lbbe-web.univ-lyon1.fr/-SiLiX)
[Roary](https://sanger-pathogens.github.io/Roary/)
[MAFFT](https://mafft.cbrc.jp/alignment/server/)
[IQ-TREE](http://www.iqtree.org)


Download the 94 genome files of the *Pseudomonas aeruginosa* strains on which the analysis is based from NCBI RefSeq

Requirements: NCBI EDirect (entrez-direct 13.9)

```bash

mkdir genome_FASTA

while read -u 94 accession; do

	echo Downloading "$accession"...

	esearch -db assembly -query "$accession" | elink -target nucleotide -name assembly_nuccore_refseq | efetch -format fasta > ./genome_FASTA/"$accession".fna

done 94<RefSeq_accession.txt

```
Re-annotate the genomes using PROKKA 1.13


Download the 94 genome files of the *Pseudomonas aeruginosa* strains on which the analysis is based from NCBI RefSeq

The list of RefSeq accessions used as input for the script to download the files is the RefSeq_accession.txt file

Requirements: NCBI EDirect (entrez-direct 13.9)


```bash

mkdir genome_FASTA #Create a new directory that will contain the genome fasta files

while read -u 94 accession; do #Iterate over accessions in the RefSeq_accession.txt file

	echo Downloading "$accession"...

	esearch -db assembly -query "$accession" | elink -target nucleotide -name assembly_nuccore_refseq | efetch -format fasta > ./genome_FASTA/"$accession".fna #Download the genome fasta file

done 94<RefSeq_accession.txt

```

Re-annotate the genomes using PROKKA 1.13

```bash

mkdir PROKKA

for genome in ./genome_FASTA/*.fna; do #Iterate over the downloaded fasta files
	fileID=$(basename "$genome" | cut -d '.' -f 1) #Define the genome ID based on the filename
	prokka --kingdom Bacteria --outdir ./PROKKA/$fileID --genus Pseudomonas --species aeruginosa --locustag $fileID --cpus 25 $genome #Run PROKKA
done

```

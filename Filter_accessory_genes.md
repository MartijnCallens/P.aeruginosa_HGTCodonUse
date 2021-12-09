# Filter accessory genes

Create a file SiLiX_unique, containing all the gene IDs for genes that did not cluster with other genes

```r
library(plyr)

clusters <- read.csv2(file = "./homolog_clustering/SiLiX/silix_clusters", header = FALSE, sep = "\t", stringsAsFactors = FALSE) # Load the file with predicted clusters %identity = 0.6; %overlap = 0.7.
colnames(clusters) <- c("cluster", "sequence") # Rename the columns

# Get a list of all cluster sizes
clusterSizes <- count(clusters, vars = "cluster")

for(i in clusterSizes[which(clusterSizes$freq == 1),]$cluster){ #Iterate over the clusters of size = 1
  write(clusters[clusters$cluster==i,]$sequence, file = "./Accessory_genes/Unique_genes/SiLiX_unique", append = TRUE, sep = "\n")
}

```

## Unique genes

Unique genes are accessory genes that only occur in a single strain 

Script to obtain sequences for unique genes predicted by SiliX

```python

import os
from Bio import SeqIO

queryFile = “./Accessory_genes/Unique_genes/SiLiX_unique" #PATH TO FILE CONTAINING UNIQUE GENE NAMES
allCDSFile = “./homolog_clustering/SiLiX/all_CDS.ffn" #PATH TO FILE CONTAINING ALL CDS FOR ALL GENOMES (MADE BY bin_all_CDS.sh)

seqs = open(queryFile, 'rU')
allCDS = SeqIO.parse(allCDSFile, "fasta")
records = []

allCDS_dict = SeqIO.to_dict(allCDS) #MAKE A SEARCHABLE DICTIONARY OF ALL CDS

for ID in seqs: #ITERATE OVE THE LOCUS IDs WITHIN THE CLUSTER
  #seqID = line.rstip()
	records.append(allCDS_dict[ID.rstrip()]) #SEARCH SEQUENCE RECORD OF ID IN THE DICTIONARY
		
SeqIO.write(records, "SiLiX_unique.ffn", "fasta") #WRITE A FILE CONTAINING SEQUENCE RECORD WITH CLUSTERNAME IN THE FILENAME

```
### blast unique genes

Blast all the unique genes against the all_CDS database

```bash

blastn -db ./homolog_clustering/SiLiX/blastdb/all_CDS_blastdb -query ./SiLiX_unique.ffn -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen" -out ./unique.blast -num_threads 25 -evalue 1e-5

```

```r

library(plyr)

#Load the blast file
UniqueBlast <- read.csv2(file = “./Accessory_genes/Unique_genes/unique.blast", sep = "\t", header = FALSE, stringsAsFactors = FALSE, strip.white = TRUE)

colnames(UniqueBlast) <- c("qseqid", "sseqid","pident","length","mismatch","gapopen","qstart","qend","sstart","send","evalue","bitscore","qlen","slen")

length(unique(UniqueBlast$qseqid))

Hits <- count(UniqueBlast, vars = "qseqid")

# Get all the genes for which there was no other significant blast hit to other genes:
UniqueNoHit <- (Hits[which(Hits$freq==1),])

# Some genes seem to have multiple hits against themselves (due to repeats?), check if they have also hits against other genes:

Unique_selfhits <- UniqueBlast[which(UniqueBlast$qseqid == UniqueBlast$sseqid),]
selfHits <- Unique_selfhits$qseqid[duplicated(Unique_selfhits$qseqid)]

UniqueNoHit_extra <- unique(selfHits[!(selfHits %in% Unique_hits$qseqid)]) #Filter out the selfhit genes without hits against other genes

# Create a file listing all truly unique genes:
trulyUnique <- c(UniqueNoHit$qseqid, UniqueNoHit_extra)

write(trulyUnique, file = “./Accessory_genes/Unique_genes/Unique", sep = "\n")

```

## Flexible genes

Flexible genes are accessory genes that occur in more than one strain

```r
# For all clusters of size 1 < n < 89, create three files, SiLiX_flexible (each gene in a separate strain), SiLiX_flexible_uniqe (genes duplicated in some strains) and SiLiX_flexible_unique (genes only occur in a single strain)

for(i in clusterSizes[which(clusterSizes$freq > 1 & clusterSizes$freq < 89),]$cluster){ #Iterate over the clusters with a defined size
  
  sequences <- clusters[clusters$cluster==i,]$sequence #Get sequence id's for this cluster
  
  check <- NULL
  for(n in strsplit(sequences, split = "_")){check <- c(check,(n[2]))} #Create a vector with genome id's
  
  if(length(check) == length(unique(check))){ # Check for cases SiLiX_flexible
     cluster <- paste(c(i, sequences), collapse = ",")
     write(cluster, file = “./Accessory_genes/Flexible_genes/SiLiX_flexible", append = TRUE)
  }
  else if(length(unique(check))==1){ # Check for cases SiLiX_flexible_unique
    cluster <- paste(c(i, sequences), collapse = ",")
     write(cluster, file = “./Accessory_genes/Flexible_genes/SiLiX_flexible_unique", append = TRUE)
  }
  else if(length(check) != length(unique(check))){ # Check for cases SiLiX_flexible_duplications
    cluster <- paste(c(i, sequences), collapse = ",")
     write(cluster, file = “./Accessory_genes/Flexible_genes/SiLiX_flexible_duplications", append = TRUE)
  }
}

```

SCRIPT TO OBTAIN FLEXIBLE GENE CLUSTERS OBTAINED BY SILIX (ONE FILE FOR EACH FLEXIBLE GENE CLUSTER)

```python

import os
from Bio import SeqIO

queryFile = “./Accessory_genes/Flexible_genes/SiLiX_flexible" #PATH TO ROARY CORE GENOME CLUSTERS
allCDSFile = “./homolog_clustering/SiLiX/all_CDS.ffn" #PATH TO FILE CONTAINING ALL CDS FOR ALL GENOMES (MADE BY bin_all_CDS.sh)

clusters = open(queryFile, 'rU') 
allCDS = SeqIO.parse(allCDSFile, "fasta")

allCDS_dict = SeqIO.to_dict(allCDS) #MAKE A SEARCHABLE DICTIONARY OF ALL CDS

for line in clusters:

	records = []	
	data = line.rstrip().split(',') #SPLIT line FROM coreClusters INTO CLUSTERNAME AND LOCUS IDs
	clusterName = data[0] #GET THE NAME OF THE CLUSTER
	fileName = clusterName + ".ffn" 

	seqIDs = data[1:] #MAKE A LIST OF ALL LOCUS IDs WITHIN THE CLUSTER
	 
	for ID in seqIDs: #ITERATE OVE THE LOCUS IDs WITHIN THE CLUSTER
		if ID.startswith("GCF"): #ONLY INCLUDE CORRECT SEQUENCE IDs
		
			records.append(allCDS_dict[ID]) #SEARCH SEQUENCE RECORD OF ID IN THE DICTIONARY
		
	SeqIO.write(records, fileName, "fasta") #WRITE A FILE CONTAINING SEQUENCE RECORD WITH CLUSTERNAME IN THE FILENAME
	
```

blast flexible genes

```bash

clusters=($(ls)) # create an array from directory names

for cluster in "${clusters[@]}"; do

  blastn -db ./homolog_clustering/SiLiX/blastdb/all_CDS_blastdb -query $cluster -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen" -out ./blast -num_threads 25 -evalue 1e-5

# Make two files with Query and Hit genes
  awk -F '\t' '{print $1}' blast | sort | uniq > blastQuery
  awk -F '\t' '{print $2}' blast | sort | uniq > blastHits

# Check if query and hit genes are the same and move files around accordingly
  if [[ $(diff -q blastQuery blastHits) ]]; then
    true
  else
    cp $cluster ../Flexible_duplications_NoHits
  fi
  rm blast blastQuery blastHits
done

```

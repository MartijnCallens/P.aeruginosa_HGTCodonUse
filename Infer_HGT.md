# Infer HGT
Ancestral state reconstruction is done for each non-unique flexible gene cluster using its phyletic pattern and the inferred species tree.

```r
library(phytools)
PA.tree <- read.tree("./core_genome/core_genome_aligned.fasta.treefile")
fitER <- ace(x, PA.tree, model="ER", type="discrete")
```

For each gene cluster, a data file for ancestral state reconstruction is produced with the name ACR_CLUSTERXXXXXX.csv.txt. The script below is used to infer HGT and residence time based on the ancestral state reconstruction data:

```r

# Infer HGT based on ancestral state reconstruction

library(phytools)
library(phylobase)
library(seqinr)
library(dplyr)

tree <- read.tree("./core_genome/core_genome_aligned.fasta.treefile")

#Modify the branch lengths at root to get a proper midpoint
tree$edge.length[1] <- 0.00066723425
tree$edge.length[124] <- 0.00430026395

tree2 <- phylo4(tree) # Convert the tree to the phylobase format
distances <- dist.nodes(tree) # get a matrix with distances between all nodes

#Specify directory with ancestral state reconstruction data
acr_path = "./Accessory_genes/Flexible_genes/Ancestral_state_reconstruction/All_clusters/TXTs/"
file.names <- dir(acr_path)

df <- NULL
threshold = 0.5

for(cluster in file.names){

  #load the data from the ancestral state reconstruction
  ancestral_state <- read.table(file = paste0(acr_path, cluster))
  #For multiple copy states, make the one column the sum of all of them
  ancestral_state$one <- unname(rowSums(subset(ancestral_state, select=-c(zero))))
  
  hgt_nodes <- NULL
  
  for(node in rownames(ancestral_state)){ #iterate over the nodes in the ancestral_state data
    if(ancestral_state[node,"one"] >= threshold){ # Check if gene is present on this node
      node_number <- as.numeric(node)
      if(node_number == 95){
        hgt_nodes <- c(hgt_nodes, node_number)
        next()
      }
      parent <- as.character(ancestor(tree2, node_number))
      if(ancestral_state[parent, "one"] < threshold){ # Check if gene is absent in the ancestral node, if yes: add node to transfer node list
        hgt_nodes <- c(hgt_nodes, node_number)
      }
    }
  }
  
  
  #Open the cluster file containing the geneID's of each cluster
  cluster_dir <- "./Accessory_genes/Flexible_genes/Ancestral_state_reconstruction/ER/gene_clusters/"
  cluster.name <- gsub("ACR_","", cluster)
  cluster.name <- gsub(".csv.txt", "", cluster.name)
  genesFasta <- read.fasta(file = paste0(cluster_dir, cluster.name, ".ffn"))
  geneNames <- getName(genesFasta)
  geneNames_split <- strsplit(geneNames, "_")
  
  strains <- NULL
  for(i in 1:length(geneNames)){
    strains <- c(strains, paste0(geneNames_split[[i]][1], "_", geneNames_split[[i]][2]))
  }
  
  strain_nodes <- match(strains, tree$tip.label)
  
  for(i in 1:length(strain_nodes)){
    
    #First check if the gene was absent in the parent of a leaf
    if(ancestral_state[as.character(unname(ancestor(tree2, strain_nodes[i]))), "one"] < threshold){
      
      node1 <- strain_nodes[i]
      df <- rbind(df, c(geneNames[i], cluster.name, node1))
      
      next()
    }
    
    # Get the hgt nodes that are parents of the strain
    hgt_parents <- hgt_nodes[hgt_nodes %in% unname(ancestors(tree2, strain_nodes[i]))]
    
    # Get the distance between the leaf and each potential hgt node and select shortest distance
    nodeDistances <- distances[strain_nodes[i], hgt_parents]
    
    if(length(hgt_parents) == 1){
      node1 <- hgt_parents
    }
    
    else{
      nodeDistances <- nodeDistances[order(nodeDistances, decreasing = FALSE)] 
      node1 <- as.numeric(names(nodeDistances[1]))
    }
    df <- rbind(df, c(geneNames[i], cluster.name, node1))
    
  }  
}

df <- as.data.frame(df, stringsAsFactors = FALSE)
colnames(df) <- c("gene", "cluster", "node1")
df$node1 <- as.integer(df$node1)

#Check for genomic context and collapse multiple transfers at the same spot.
library(stringr)

# Check between which core genes each accessory gene is inserted

cluster_tab <- read.table(file = "./homolog_clustering/SiLiX/silix_clusters", stringsAsFactors = F)
colnames(cluster_tab) <- c("cluster", "gene")

core_genes <- scan(file = "./core_genome/all_core_genes", what = character(), sep  = ",")
core_genes <- core_genes[grep("GCF", core_genes)]

df$location <- NA

counter = 0

for(g in 1:length(df$gene)){
  
  counter = counter + 1
  print(counter)
  
  loc <- NULL
  upstr_cluster <- NULL
  dwnstr_cluster <- NULL
  
  gene_number <- as.integer(strsplit(df$gene[g], "_")[[1]][3])
  strain <- paste0("GCF_", strsplit(df$gene[g], "_")[[1]][2])
  
  for(i in 1:1000){
    upstr_gene <- paste0(strain, "_", str_pad(gene_number - i, 5, pad = "0"))
    if(upstr_gene %in% core_genes){
      upstr_cluster <- cluster_tab[which(cluster_tab$gene == upstr_gene), ]$cluster
      break()
      
    }
  }
  
  for(i in 1:1000){
    dwnstr_gene <- paste0(strain, "_", str_pad(gene_number + i, 5, pad = "0"))
    if(dwnstr_gene %in% core_genes){
      dwnstr_cluster <- cluster_tab[which(cluster_tab$gene == dwnstr_gene), ]$cluster
      break()
    }
  }
  
  loc <- sort(c(upstr_cluster, dwnstr_cluster))
  loc <- paste(loc[1], loc[2], sep = ".")
  
  df[g, ]$location <- loc
}

# >> Collapse nodes having the same insertion point

for(cluster in unique(df$cluster)){
  
  df_sub <- df[which(df$cluster == cluster), ]
  
  for(spot in unique(df_sub$location)){
    
    df_sub_sub <- df_sub[which(df_sub$location == spot), ]
    nodes <- unique(df_sub_sub$node1)
    
    if(length(nodes) == 1){next()}
    
    all_anc <- ancestors(tree2, nodes)
    
    df[row.names(df_sub_sub), ]$node1 <- Reduce(intersect, all_anc)[1]
  }
}

df <- df[which(df$node1 != 95), ] #Remove genes with transfers at the root

# For each node, get the ancestor and residence time
df$residenceTime <- NA

for(i in 1:nrow(df)){
  strain <- paste0("GCF_", strsplit(df$gene[i], "_")[[1]][2])
  strain_node <- as.integer(match(strain, tree$tip.label))
  node1 <- as.numeric(df[i, ]$node1)
  node2 <- as.numeric(ancestor(tree2, node1))
  
  node1Distance <- distances[strain_node, node1]
  node2Distance <- distances[strain_node, node2]
  
  residenceTime <- node1Distance + ((node2Distance - node1Distance)/2)
  
  df[i, ]$residenceTime <- residenceTime 
}

#Get number of genes at root
nrow(df[which(df$node1 == 95), ])
length(unique(df[which(df$node1 == 95), ]$cluster))

df <- df[which(df$node1 != 95), ] # Do not remove basal transfers before obtaining the transfer nodes with checking location!

#Get mean residence time of flexible genes:

hist(df$residenceTime)
mean(df$residenceTime)

#Add the unique genes to the dataset
Residence_time_unique <- read.csv2("./Accessory_genes/Unique_genes/Residence_time/Unique_res_time.txt", 
                                   header = FALSE, 
                                   col.names = c("gene","cluster","residenceTime"),
                                   sep =",",
                                   stringsAsFactors = FALSE,
                                   dec=".",
                                   colClasses = c("character", "character","numeric"))

Residence_time_unique$node1 <- NA
Residence_time_unique$location <- NA

df <- rbind(df, Residence_time_unique)


Codon_usage <- read.csv2("./Codon_usage/COUSIN_Results/results/results.txt", 
                         header = TRUE, 
                         sep = "\t", 
                         stringsAsFactors = FALSE,
                         dec=".")

#Only keep the first part of the name

for (i in seq(1, length(Codon_usage$Header),1)){
  Codon_usage$Header[i] <- strsplit(Codon_usage$Header[i], " ")[[1]][1]
}

df$COUSIN59 <- NA

for (position in seq(1, length(df$gene), 1)){
  value <- NULL
  value <- Codon_usage[which(Codon_usage$Header == df$gene[position]), 12]
  if (length(value > 0)){df$COUSIN59[position] <- value}
}

df <- df[-which(is.na(df$COUSIN59)), ]
#df <- df[-which(is.na(df$residenceTime)), ]

df.bkp <- df

# Get the number of inferred transfer (only including the flexible genes)

df_flex <- filter(df, cluster != " UNIQUE")
transfers <- NULL
for(c in unique(df_flex$cluster)){
  transf = 0
  df_sub <- df_flex[which(df_flex$cluster == c), ]
  transf = transf + length(unique(df_sub[which(df_sub$node1 < 95), ]$node2))
  transf = transf + length(unique(df_sub[which(df_sub$node1 >= 95), ]$node2))
  if(transf > 20){print(c(c,transf))}
  transfers <- c(transfers, transf)
}

#remove transfers between root and 2nd nodes

df <- filter(df, is.na(node2) | node2 != 95)

df_root <- filter(df.bkp, node2 == 95)

```

# This script takes a 012 formatted file (as output by VCFTools) and creates
#a pairwise IBS similariy matrix for use in identifying and filtering clones

##### LOAD 012 FILES ####
library(knitr)

opts_knit$set(root.dir = '~/IBS and Clone Correction/')
setwd("~/IBS and Clone Correction/D_rhei_combined/")

snps_a <- read.table("D:/BioHPC GBS analysis/Tassel5 Output/9_1_23 Redo/IBS and Clone Correction/D_rhei_9_1_23.012.pos")
head(snps_a)

indv_a <- read.table("D:/BioHPC GBS analysis/Tassel5 Output/9_1_23 Redo/IBS and Clone Correction/D_rhei_9_1_23.012.indv")
indv_a <- unlist(indv_a$V1)

geno_a <- read.table("D:/BioHPC GBS analysis/Tassel5 Output/9_1_23 Redo/IBS and Clone Correction/D_rhei_9_1_23.012")
print("geno_a loaded")
geno_a <- geno_a[,-1]
geno_a <- t(geno_a)

geno_a[geno_a==-1] <- NA
print('finished replacing NAs')

#### IBS FUNCTION ####

ibs <- function(x,y){
  
  alleles_ibs <- 2 - abs(x-y)
  return(sum(alleles_ibs, na.rm = T)/(2*sum(!is.na(alleles_ibs))))
  
}

#### CALCULATE IBS FOR EACH PAIRWISE ISOLATE COMBINATION ####

d <- ncol(geno_a)

IBS_matrix_a <- matrix(nrow=d, ncol=d)

print("got to loop")

for(i in 1:(d-1)){
  for (j in (i +1):d){
    IBS_matrix_a[i,j] <- ibs(geno_a[,i], geno_a[,j])
  }
  print(i)
}

rownames(IBS_matrix_a) <- indv_a
colnames(IBS_matrix_a) <- indv_a

write.csv(IBS_matrix_a, "D:/BioHPC GBS analysis/Tassel5 Output/9_1_23 Redo/IBS and Clone Correction/IBS_matrix_9_1_23.csv")

#Clone correction of D rhei isolates 
#This will return file with clonal group designations of all isolates
#Then it will return a list of the samples corresponding to the cc dataset

library(igraph)

#### READ IBS MATRIX, GENOTYPE FILES, AND MISSING DATA FILE ####

#Load IBS matrix and make individuals row and column names
IBS_matrix_a <- read.csv("D:/BioHPC GBS analysis/Tassel5 Output/9_1_23 Redo/IBS and Clone Correction/IBS_matrix_9_1_23.csv", header = T)
row.names(IBS_matrix_a) <- IBS_matrix_a$X 
IBS_matrix_a$X <- NULL 
IBS_matrix_a <- as.matrix(IBS_matrix_a)

#Load missing data info
miss_a <- read.table("D:/BioHPC GBS analysis/Tassel5 Output/9_1_23 Redo/IBS and Clone Correction/FinalFiltering_9_1_23.imiss", header = T)
write.csv(miss_a, "D:/BioHPC GBS analysis/Tassel5 Output/9_1_23 Redo/IBS and Clone Correction/missing_data_9_1_23.csv")

#Load phenos data to see which samples are in storage
phenos <- read.csv("D:/BioHPC GBS analysis/Tassel5 Output/9_1_23 Redo/IBS and Clone Correction/D_rhei_isolates.csv")


#### ASSIGN INDIVIDUALS TO CLONAL GROUPS ####

#Turn high IBS cells of matrix to 1
modify_matrix <- function(x){
  if(is.na(x) | x<.97){
    return(0)
  }else{
    
    return(1)
  }
}
clone_or_not <- structure(sapply(IBS_matrix_a, modify_matrix), dim=dim(IBS_matrix_a))

#Create network -> Each isolate is a node and there is an edge b/w clonal isolates
g <- graph_from_adjacency_matrix(clone_or_not, "undirected")

#Clusters are isolates that only have edges between themselves and not rest of
#the network (i.e. clones)
g.clusters <- clusters(graph = g)
g.clusters

#Create table of clonal group assignments

#Make list of cluster size corresponding to each member of network (used later)
cluster_sizes <- rep(NA, length(indv_a))
for(i in 1:length(cluster_sizes)){
  member <- g.clusters$membership[i]
  size <- sum(g.clusters$membership == member)
  cluster_sizes[i] <- size
}

#Prepare table and variables for loop
clonal_groups <- 1:(g.clusters$no)
clone_assignments_a <- matrix(ncol=2)
colnames(clone_assignments_a) <- c("Sample", "Clonal_group")
counter <- 0

#Assign individuals to clonal groups starting with largest group
for(i in 1:length(unique(g.clusters$csize))){ #loop through all unique cluster sizes
  #Start with largest cluster size
  current_size <- sort(unique(g.clusters$csize), decreasing=T)[i]
  #how many groups of this size are there
  same_size_clonal_groups <- unique(g.clusters$membership[cluster_sizes == current_size])
  #loop through groups of that size
  for(j in 1:length(same_size_clonal_groups)){
    counter <- counter +1
    old_clonal_group_id <- same_size_clonal_groups[j] #Assignment to group from g.clusters$membership
    new_clonal_group_assignment <- clonal_groups[counter] #New assignment from largest to smallest
    clone_assignments_a <- rbind(clone_assignments_a, cbind(
      indv_a[which(g.clusters$membership == old_clonal_group_id)],
      new_clonal_group_assignment))
  }
}
clone_assignments_a <- clone_assignments_a[-1,]
clone_assignments_a <- as.data.frame(clone_assignments_a, stringsAsFactors = F)
clone_assignments_a$Clonal_group <- as.integer(clone_assignments_a$Clonal_group)

write.table(clone_assignments_a, "D:/BioHPC GBS analysis/Tassel5 Output/9_1_23 Redo/IBS and Clone Correction/clone_assignments_9_1_23.txt", row.names = F, quote = F, sep = "\t")

# Clone correct -- choosing individual with least missing data per clonal group
# and ensuring individual is in storage if possible

#Get just part of the name that will match with pheno file
shorten_name <- function(x){
  return(unlist(strsplit(x,":"))[1])
}

#pheno_names <- paste(phenos$Sample, phenos$SZ_for_analysis, sep="")

indvs_cc <- rep(NA, max(clone_assignments_a$Clonal_group))
cc_samples_list <- matrix(rep(NA, 2*max(clone_assignments_a$Clonal_group)), ncol=2)
colnames(cc_samples_list) <- c("Sample", "In_storage")
for(i in 1:max(clone_assignments_a$Clonal_group)){
  samples.i <- clone_assignments_a$Sample[clone_assignments_a$Clonal_group == i]
  samples.short.i <- sapply(samples.i, shorten_name)
  samples.exist <- as.logical(phenos$In_storage[match(samples.short.i, phenos$Sample)])
  if(any(samples.exist)){
    samples.i <- samples.i[samples.exist]
    samples.short.i <- samples.short.i[samples.exist]
    cc_samples_list[i,2] <- 1
    
  }else{
      cc_samples_list[i,2] <- 0
  }
  samples.missing <- miss_a$F_MISS[match(samples.i, miss_a$INDV)]
  best_sample <- which(samples.missing == min(samples.missing))[1] #[1] added if >1 sample meets this condition
  indvs_cc[i] <- samples.i[best_sample]
  cc_samples_list[i,1] <- samples.short.i[best_sample]
}

write.table(indvs_cc, "rhei_diversity_cc_9_1_23.012.indv", quote = F, row.names = F, col.names = F)
write.csv(cc_samples_list, "rhei_diversity_cc_9_1_23_samplelist.csv", quote = F, row.names = F)

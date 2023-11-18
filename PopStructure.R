#if pcaMethods isn't already installed
if(!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("pcaMethods")

library(pcaMethods)
library(ape)
library(RColorBrewer)
library(StAMPP)

#### Clean and load data ####

#Load data
geno <- read.table("D_rhei_CC_9_4_23.012")
indvs <- read.table("D_rhei_CC_9_4_23.012.indv")
snps <- read.table("D_rhei_CC_9_4_23.012.pos")
didy <- read.csv("Isolate_metadata.csv", stringsAsFactors = F)
meta <- read.csv("isolate_plotting_metadata.csv", stringsAsFactors = F)
field_info <- read.csv("field_plotting_metadata.csv", stringsAsFactors = F)
#tree

#Clean up data
indvs <- unlist(indvs$V1)
indvs <- sapply(indvs, function(x) unlist(strsplit(x,":"))[1])

geno <- geno[,-1]
geno <- as.matrix(geno)
geno[geno==-1] <- NA

#Subset metadata for just cc set and get in same order as geno
didy <- didy[match(indvs, didy$SampleSZ),]
meta <- meta[match(indvs, meta$SamplesSZ),]

#### PCA ####

#Get PCs
geno.pca <- pca(geno, method = "nipals", nPcs = 4, scale = 'uv')
pve <- round(geno.pca@R2*100,2)

#### Make plot of PC1v2 and PC3v4
#pdf("plots/PCA.pdf", width=6, height=6) 
plot_layout <- cbind(c(1,1, 2, 2),
                     c(3, 3, 3, 3))
layout(plot_layout)
#Par settings for plots 1 and 2 (PCA plots)
old.par <- par(no.readonly = T)
par(mar=c(4,1,2,1), oma=c(4,4,4,4), xpd=NA)

#PC1 vs PC2
plot(geno.pca@scores[,1],
     geno.pca@scores[,2],
     pch=meta$pch,cex=1.5,col=meta$color,
     xlab = paste("PC1:", pve[1], "%"),
     ylab = paste("PC2:", pve[2], "%"),
     cex.lab=1.4, cex.axis=1.5)

#PC3 vs PC4
plot(geno.pca@scores[,3],
     geno.pca@scores[,4],
     pch=meta$pch,cex=1.5,col=meta$color,
     xlab = paste("PC3:", pve[3], "%"),
     ylab = paste("PC4:", pve[4], "%"),
     cex.lab=1.4, cex.axis=1.5)

#Plot 3
populations <- meta$Field
plot(0,"n", bty="n", xlim=c(1,10), ylim=c(1,10), xaxt='n', yaxt='n', xlab='', ylab='')
legend_starts <- c(11, 9, 6, 4, 2.5, 1.5)
site_counts <- table(populations)
for(i in 1:length(unique(field_info$Region))){
  region <- unique(field_info$Region)[i]
  sites <- field_info$Field[field_info$Region==region]
  
  colors <- field_info$color[field_info$Region==region]
  pch <- field_info$pch[field_info$Region==region]
  legend(legend_starts[i],legend=paste(sites, ", n=", site_counts[match(sites, names(site_counts))], sep=""),
         col = colors,
         pch = pch,
         bty='n', title = as.expression(bquote(bold(.(region)))),
         title.adj = 0, cex=1.6)
}
par(old.par)
dev.off()

#### FST between fields ####

#Get fields with more than 5 isolates after clone-correction
population.sizes <- table(populations)
large_pops <- names(population.sizes)[population.sizes>5 & names(population.sizes)!="NonNY"]
fields <- populations 
fields[!fields %in% large_pops] <- NA

#Turn geno into allele frequency format
geno.stamp <- geno/2

#Make column names SNP names
colnames(geno.stamp) <- paste("S_",snps$V1, "_", snps$V2, sep = "")

#Add isolate, population, ploidy, and format information
geno.stamp <- cbind("Isolate" = didy$SampleSZ,
                    "Subpop" = fields,
                    "Ploidy" = 1,
                    "Format" = 'freq',
                    geno.stamp)
geno.stamp <- as.data.frame(geno.stamp, stringsAsFactors = 1)
geno.stamp <- geno.stamp[!is.na(geno.stamp$Subpop),]

#Make object for StAMPP functions
geno.stamp <- stamppConvert(geno.stamp, "r")

#Weir and Cockerham's Fst, 100 bootstraps, 95% CI
fst.stamp <- stamppFst(geno.stamp, 100, 95)
fst.mat <- fst.stamp$Fsts

#Re-order matrix
order_matrix <- function(mat, ord){
  
  mat.sort <- matrix(NA, nrow=nrow(mat), ncol=ncol(mat))
  rownames(mat.sort) <- ord
  colnames(mat.sort) <- ord
  
  for(i in 1:(nrow(mat.sort)-1)){
    ord.i <- ord[i]
    for(j in (i+1):ncol(mat.sort)){
      ord.j <- ord[j]
      mat.values <- c(mat[ord.i,ord.j], mat[ord.j,ord.i])
      mat.value <- mat.values[which(!is.na(mat.values))]
      mat.sort[i,j] <- mat.value
    }
  }
  return(mat.sort)
}
fst.mat <- order_matrix(fst.mat,sort(rownames(fst.mat)))
field.sizes <- population.sizes[match(rownames(fst.mat), names(population.sizes))]
fst.mat <- cbind("Number isolates" = field.sizes, fst.mat)
write.table(fst.mat, "pairwise_fst_2021_and_2022.txt", row.names = T, quote = F, sep = "\t")

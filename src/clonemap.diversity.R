args <- commandArgs(TRUE)
treebuildDir <- as.character(args[1])
case.id <- as.character(args[2])

library(cloneMap)
library(vegan)

set.seed(1)

setwd(treebuildDir)


treerds <- readRDS(file.path(paste0(treebuildDir, "/", case.id, '.tree.RDS')))
clusterInfo <- read.delim2(file.path(treebuildDir, 'clusterInfo.txt'))


tree1 <- treerds$graph_pyclone$Corrected_tree
# tree1 <- treerds$graph_pyclone$alt_trees[[1]] # (alternatively)

# Plot clone map for one sample:
# NOTE: make sure to filter for clusters that made it onto the tree (clusterInfo$treeClust == T)

SampleNames <- unique(clusterInfo$SAMPLE)
SampleNames <- SampleNames[!is.na(SampleNames)]

clone.ccf <- data.frame()
dfdiv <- data.frame(sample=SampleNames, shannonDiv="")


for ( i in SampleNames){
    print(i)
    cluster<- clusterInfo[clusterInfo$treeClust == T & clusterInfo$SAMPLE == i,]
    cluster_ccf <- cluster[, c('meanCCF', 'clusterID')]
    colnames(cluster_ccf) <- c('CCF', 'clones')
    cluster_ccf <- cluster_ccf[, c('clones', 'CCF')]
    # Get Shannon diversity for each sample
    cluster_cp <- clusterInfo[clusterInfo$treeClust == T & clusterInfo$SAMPLE == i,]$clone_proportions_default
    dfdiv$shannonDiv[dfdiv$sample == i] <- diversity(cluster_cp)

    # plot clonemap for each sample
    png(paste0(treebuildDir, "/", i, ".clonemap.png"))
    cloneMap(tree1, cluster_ccf)
    dev.off()
    
    cluster_ccf$sample <- i
    clone.ccf <- rbind(clone.ccf,cluster_ccf)
}

write.table(dfdiv, file.path(treebuildDir, 'Shannon.div.txt'), sep="\t", quote=FALSE, row.names=FALSE)
write.table(clone.ccf, file.path(treebuildDir, 'clone.ccf.txt'), sep="\t", quote=FALSE, row.names=FALSE)

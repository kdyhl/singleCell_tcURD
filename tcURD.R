#a constraint neighbor version of URD
#use bigmemory to handle the big matrix
#use parallelDist to calculate distance


rm(list=ls())
library(bigmemory)
library(parallelDist)
suppressPackageStartupMessages(library(rgl))
suppressPackageStartupMessages(library(URD))
suppressPackageStartupMessages(library(destiny))
library(gridExtra)

#setwd("your own working directory")


#
stage.colors=RColorBrewer::brewer.pal(9,'Set1')[c(1:5,7)]
TIMEDIFF.MAX=1 #max difference in terms of time pont order that two neighbors should have

#
dir.ou="yourOwnOutputDirectory/"
if(!dir.exists(dir.ou))
{
  dir.create(dir.ou)  
}  

file.in.txt = "example.exp.txt"#./URD/U1.pgcalldata.combinedMatrix.filtered.rds"
file.in.meta = "example.umaps.info.txt"#"./URD/pgcall.Seurat_k30r1.umaps.info.txt"


file.ou.varGene = paste(dir.ou, "var_genes.txt", sep="")
file.ou.pgcurd = paste(dir.ou, "pgcurd.rds", sep="")
file.ou.pgcurd.DM= paste(dir.ou, "pgcurdWithDM.rds", sep="")
file.ou.pgcurd.DMandPT = paste(dir.ou, "pgcurdWithDMandPT.rds", sep="")
file.ou.pgcurd.D4 = paste(dir.ou, "pgcurdWithDMandPT.D4.rds", sep="")
file.ou.pgcurd.RW = paste(dir.ou, "pgcurdWithDMandPT.RW.rds", sep="")
file.ou.pgcurd.tree = paste(dir.ou, "pgcurdWithDMandPT.tree.rds", sep="")
file.ou.pgcurd.tree.ForceLayout = paste(dir.ou, "pgcurdWithDMandPT.tree.ForceLayout.rds", sep="")
  
file.ou.pgcurd.tsne.pdf= paste(dir.ou, "pgcurd.tsne.pdf", sep="")
file.ou.pgcurd.DMandPT.pdf= paste(dir.ou, "pgcurdWithDMandPT.pdf", sep="")
file.ou.pgcurd.D4.pdf = paste(dir.ou, "pgcurdWithDMandPT.D4.pdf", sep="")
file.ou.pgcurd.RW.pdf = paste(dir.ou, "pgcurdWithDMandPT.RW.pdf", sep="")
file.ou.pgcurd.tree.pdf = paste(dir.ou, "pgcurdWithDMandPT.tree.pdf", sep="")
#


expmatrix <- t(as.matrix(read.table(file.in.txt, header = T,sep = '\t',row.names = 1, stringsAsFactors = F)))
expmatrix.fake = 2^expmatrix

metadata <- read.table(file.in.meta, header = T,sep = '\t',row.names = 1, stringsAsFactors = F)
metadata$timepoint = factor(metadata$timepoint, levels = c("U2hESC", "U2iMeLC", "U2PGCLCd1", "U2PGCLCd2", "U2PGCLCd3", "U2PGCLCd4")) #shows the temporal order of time points

pgcurd = createURD(count.data = expmatrix.fake, meta = metadata, min.cells=0, min.counts=0, min.genes=0)#, gene.max.cut=5000, max.genes.in.ram=5000)
pgcurd@logupx.data = as(expmatrix, "dgCMatrix")
pgcurd@var.genes <- rownames(expmatrix)
## delete the original data
rm(list=c('expmatrix.fake', 'expmatrix','metadata'))
## perform garbage collection to free RAM
shhhh=gc()

################
# pgcurd@group.ids$stage <- as.character(pgcurd@meta[rownames(pgcurd@group.ids),"timepoint"])
# # Get variable genes for each group of 3 stages
# # (Normally would do this for each stage, but there are not very many cells in this subset of the data)
# # diffCV.cutoff can be varied to include more or fewer genes.
# stages <- sort(unique(pgcurd@group.ids$stage))
# cells.each.stage <- lapply(stages, function(stage) rownames(pgcurd@meta)[which(pgcurd@meta$timepoint == stage)])
# # Compute variable genes for each stage.
# var.genes.by.stage <- lapply(1:length(stages), function(n) findVariableGenes(pgcurd, cells.fit = cells.each.stage[[n]], set.object.var.genes = F, diffCV.cutoff = 0.5, mean.min = 0.005, mean.max = 100, main.use = stages[[n]], do.plot = T))
# # Combine the results from each group of stages into a single list of variable genes and load into the URD object
# names(var.genes.by.stage)=stages
# var.genes <- sort(unique(unlist(var.genes.by.stage)))




# for (stage in stages){write(var.genes.by.stage[[stage]],file=paste(dir.ou, 'var_genes_',stage,".txt", sep=""))}
# write(var.genes,file=file.ou.varGene)

saveRDS(pgcurd,file=file.ou.pgcurd)

# # Calculate PCA and consider those PCs
# pgcurd <- calcPCA(pgcurd,mp.factor = 1.5)
# pcSDPlot(pgcurd)
# # Calculate tSNE
# set.seed(18)
# pgcurd <- calcTsne(object = pgcurd,perplexity=30,theta=0.5)
# 
# set.seed(17)
# 
# pgcurd <- graphClustering(pgcurd, dim.use = 'pca',num.nn = c(20,30,40,50,60),do.jaccard = T,method = 'Louvain')
# clusterings=c(paste0('Louvain-',c(20,30,40,50,60)))
# 
# 
# pdf(file=file.ou.pgcurd.tsne.pdf,width = 8,height = 7)
# plotDim(pgcurd,"stage",discrete.colors = stage.colors,legend = T,alpha=1,plot.title = "Stage")
# for (c in clusterings){plot(plotDim(pgcurd,c,legend = T))}
# plotDim(pgcurd, "NANOS3", plot.title="NANOS3")
# plotDim(pgcurd, "NANOG", plot.title="NANOG")
# plotDim(pgcurd, "SOX17", plot.title="SOX17")
# dev.off()
# 
# saveRDS(pgcurd,file=file.ou.pgcurd)


find_constraintKNN = function(data, timePoints, k, timeDiff.max=TIMEDIFF.MAX, distance="euclidean")
{
  dis.matrix = parDist(data, method=distance) #parDist ##keep the in dist format after dist function in function getConstraintDis
  tP.matrix = dist(timePoints, method="euclidean") #
  
  dist.max = max(dis.matrix)
  
  dis.matrix[tP.matrix>timeDiff.max]=dist.max
  
  dis.matrix = as.matrix(dis.matrix)
  
  dis.matrix = apply(dis.matrix, 2, 
  FUN=function(x)
  {
    x[rank(x)>k+1]=0 #since self distance is not considered as neighbors
    return(x)
  })
  dis.matrix.sym <- symmpart(dis.matrix) + abs(forceSymmetric(skewpart(dis.matrix), 'U')) #keep it symetric
  
  
  return(dis.matrix.sym)
}

find_constraintKNN.smallK = function(data, timePoints, k, timeDiff.max=1, distance="euclidean")
#small K means k<< potential neighbors constrained by timeDiff.max
{
  # #for debug
  # data = as.matrix(t(pgcurd@logupx.data))
  # timePoints = as.numeric(pgcurd@meta$timepoint)
  # k=50; timeDiff.max=1; distance="euclidean"
  #
  dis.matrix = parDist(data, method=distance) #parDist ##keep the in dist format after dist function in function getConstraintDis
  tP.matrix = dist(timePoints, method="euclidean")
  
  #dist.max = max(dis.matrix)
  
  dis.matrix[tP.matrix>timeDiff.max]=0 #0 here means no possibilyt to be neighbor, including itself
  rm(list='tP.matrix')
  dis.matrix = as.matrix(dis.matrix) #as(, "dgCMatrix")
  dis.matrix = as.big.matrix(dis.matrix)
    
  for(i in 1:ncol(dis.matrix))
  {
    print(i)
    x= dis.matrix[,i]
    dis.matrix[x!=0,i][rank(x[x!=0])>k]=0 
  }
  
  # dis.matrix = apply(dis.matrix, 2, 
  # FUN=function(x)
  # {
  #   x[x!=0][rank(x[x!=0])>k]=0 
  #   return(x)
  # })
  
  dis.matrix = as(as.matrix(dis.matrix),"dgCMatrix")
  print("sparse matrix generated")
  dis.matrix.sym <- symmpart(dis.matrix) + abs(forceSymmetric(skewpart(dis.matrix), 'U')) #keep it symetric
  
  
  return(dis.matrix.sym)
}

calcDM.constraint <- function(object, genes.use=object@var.genes, cells.use=NULL, knn=NULL, sigma.use=NULL, n_local=5:7, distance=c("euclidean", "cosine", "rankcor"), density.norm=T, dcs.store=200, verbose=T, timeDiff.max=NULL) {
  
  #extract timeponts in numeric format
  timePoints = as.numeric(object@meta$timepoint)
  names(timePoints) = colnames(object@logupx.data)
  
  # Subset the data and convert to a matrix without row or column names because they crash destiny.
  if (is.null(genes.use) || length(genes.use) == 0) genes.use <- rownames(object@logupx.data)
  if (is.null(cells.use)) cells.use <- colnames(object@logupx.data)
  data.use <- t(object@logupx.data[genes.use, cells.use])
  timePoints = timePoints[cells.use]
  rownames(data.use) <- NULL
  colnames(data.use) <- NULL
  data.use <- as.matrix(data.use)
  
  # Figure out sigma
  if (is.null(sigma.use)) {
    sigma.use <- find_sigmas(data.use, steps=25, verbose=F)@optimal_sigma
    if (verbose) print(paste("destiny determined an optimal global sigma of", round(sigma.use, digits=3)))
  } else if (is.numeric(sigma.use)) {
    if (verbose) print(paste("Using provided global sigma", sigma.use))
  } else if (sigma.use == "local") {
    if (verbose) print(paste("Using local sigma."))
  } else {
    stop("sigma must either be NULL, 'local', or numeric.")
  }
  
  # Figure out k-nearest neighbors
  if (is.null(knn)) {
    knn <- find_dm_k(length(cells.use))
    if (verbose) print(paste("destiny will use", knn, "nearest neighbors."))
  }
  
  #calculate distances with constraint and Calculate the diffusion map
  dm=list()
  if(is.null(timeDiff.max))
  {
    dm <- DiffusionMap(data.use, sigma=sigma.use, k=knn, n_eigs = dcs.store, density_norm = density.norm, distance=distance[1])
  }else
  {
    print("distances calculating...")
    knn_distances = find_constraintKNN.smallK(data.use, timePoints = timePoints, k=knn, timeDiff.max = timeDiff.max)  
    print(str(knn_distances))
    dm <- DiffusionMap(data=NULL, sigma=sigma.use, k=knn, n_eigs = dcs.store, density_norm = density.norm, distance=knn_distances)
  }
  
  
  # Name things because you have to remove all of the names
  rownames(dm@eigenvectors) <- cells.use
  rownames(dm@transitions) <- cells.use
  colnames(dm@transitions) <- cells.use
  
  # Store the resulting diffusion map
  object <- importDM(object, dm)
  
  return(object)
}




### step2: Diffusion map and pseudotime

# Calculate Diffusion Map

# #for test
# cells.use = sample(colnames(pgcurd@logupx.data), 1000, replace = F)
# pgcurd.test = calcDM.constraint(pgcurd, cells.use = cells.use, knn = 30, sigma.use=16, timeDiff.max=1)   
# pgcurd = calcDM(pgcurd, sigma.use='local', knn = 150)
# pgcurdlocal = calcDM(pgcurd, sigma.use='local')
# saveRDS(pgcurd@dm,file='dm16.rds')


#pgcurd = readRDS(file.ou.pgcurd)
pgcurd <- calcDM.constraint(pgcurd, knn = 50,  sigma.use=16, timeDiff.max=1)#

saveRDS(pgcurd, file.ou.pgcurd.DM)

plotDimArray(pgcurd, reduction.use = "dm", dims.to.plot = 1:18, discrete.colors=stage.colors, outer.title = "Diffusion Map (Sigma 16, 100 NNs): Stage", label="stage", plot.title="", legend=F,alpha=0.3)
#plotDim(pgcurd, "stage", transitions.plot = 10000, plot.title="Developmental stage (with transitions)")

# Calculate pseudotime
# Here we use all cells from the first stage as the root
root.cells <- cellsInCluster(pgcurd, "stage", "U1hESC")
# Then we run 'flood' simulations
pgcurd.floods <- floodPseudotime(pgcurd, root.cells = root.cells, n=30, minimum.cells.flooded = 1, verbose=T)
# The we process the simulations into a pseudotime 
pgcurd <- floodPseudotimeProcess(pgcurd, pgcurd.floods, floods.name="pseudotime",max.frac.NA = 0.8,pseudotime.fun = mean,stability.div = 20)
saveRDS(pgcurd,file = file.ou.pgcurd.DMandPT) #diffusion map and pseudo time


pdf(file=file.ou.pgcurd.DMandPT.pdf,width = 8,height = 7)
pseudotimePlotStabilityOverall(pgcurd)

plotDim(pgcurd, "pseudotime")
plotDim(pgcurd, "stage")
plotDim(pgcurd, "NANOS3", plot.title="NANOS3")
plotDim(pgcurd, "NANOG", plot.title="NANOG")
plotDim(pgcurd, "SOX17", plot.title="SOX17")
plotDists(pgcurd, "pseudotime", "stage", plot.title="Pseudotime by stage")

dev.off()



saveRDS(pgcurd,file = file.ou.pgcurd.DMandPT) #diffusion map and pseudo time


### step3: determining Tips
library(scran)

pgcurd.d4 <- urdSubset(pgcurd, cells.keep=cellsInCluster(pgcurd, "stage", c('U1PGCLCd3',"U1PGCLCd4")))
pgcurd.d4@var.genes <- unique(c(var.genes.by.stage[[5]],var.genes.by.stage[[6]]))

pgcurd.d4 <- calcPCA(pgcurd.d4,mp.factor = 1.5)
#pcSDPlot(pgcurd.d4)
set.seed(18)
pgcurd.d4 <- calcTsne(pgcurd.d4,dim.use = 'pca',perplexity = 30,theta = 0.5)
# Calculate graph clustering of these cells
pgcurd.d4 <- graphClustering(pgcurd.d4, num.nn = c(20,50,100,150,200,250,300), do.jaccard=T, method="Louvain")
pgcurd.d4 <- graphClustering(pgcurd.d4, num.nn = c(20,50,100,150,200,250,300), do.jaccard=T, method="Infomap")
## plot
clusterings=c(paste0('Louvain-',c(20,50,100,150,200,250,300)),paste0('Infomap-',c(20,50,100,150,200,250,300)))
pdf(file=file.ou.pgcurd.D4.pdf,width = 10,height = 7)
for (c in clusterings){plot(plotDim(pgcurd.d4,c,legend = T,label.clusters = T))}
plotDim(pgcurd.d4, "NANOS3", plot.title="NANOS3")
plotDim(pgcurd.d4, "NANOG", plot.title="NANOG")
plotDim(pgcurd.d4, "SOX17", plot.title="SOX17")
dev.off()

saveRDS(pgcurd.d4, file=file.ou.pgcurd.D4)

##### Biased random walks
# Copy cluster identities from axial.6somite object to a new clustering ("tip.clusters") in the full axial object.
pgcurd@group.ids[rownames(pgcurd.d4@group.ids), "tip.clusters"] <- pgcurd.d4@group.ids$`Louvain-150`
# Determine the parameters of the logistic used to bias the transition probabilities. The procedure
# is relatively robust to this parameter, but the cell numbers may need to be modified for larger
# or smaller data sets.
pgcurd.ptlogistic <- pseudotimeDetermineLogistic(pgcurd, "pseudotime", optimal.cells.forward=20, max.cells.back=40, do.plot = T)
# Bias the transition matrix acording to pseudotime
pgcurd.biased.tm <- as.matrix(pseudotimeWeightTransitionMatrix(pgcurd, "pseudotime", logistic.params=pgcurd.ptlogistic))
# Simulate the biased random walks from each tip
pgcurd.walks <- simulateRandomWalksFromTips(pgcurd, tip.group.id="tip.clusters", root.cells=root.cells, transition.matrix = pgcurd.biased.tm, n.per.tip = 25000, root.visits = 1, max.steps = 5000, verbose = F)
# Process the biased random walks into visitation frequencies
pgcurd <- processRandomWalksFromTips(pgcurd, pgcurd.walks, verbose = F)

plotDim(pgcurd, "tip.clusters", plot.title="Cells in each tip")
saveRDS(pgcurd,file = file.ou.pgcurd.RW)


# Build tree###############################################################
# Load the cells used for each tip into the URD object
pgcurd.tree <- loadTipCells(pgcurd, "tip.clusters")
# Build the tree
pgcurd.tree <- buildTree(pgcurd.tree, pseudotime = "pseudotime", tips.use=1:10, divergence.method = "preference", cells.per.pseudotime.bin = 25, bins.per.pseudotime.window = 8, save.all.breakpoint.info = T, p.thresh=0.001)
# Name the segments based on our previous determination of the identity of tips 1 and 2.
pgcurd.tree <- nameSegments(pgcurd.tree, segments=c(1:10), segment.names = paste('cluster_',c(1:10),sep=''), short.names = paste('cluster_',c(1:10),sep=''))

pdf(file=file.ou.pgcurd.tree.pdf,width = 8,height = 6)
plotTree(pgcurd.tree, "stage", title="Developmental Stage")
plotTree(pgcurd.tree, "NANOS3", title="NANOS3")
plotTree(pgcurd.tree, "NANOG", title="NANOG")
plotTree(pgcurd.tree, "SOX17", title="SOX17")

plotTree(pgcurd.tree, "segment", title="URD tree segment")
plotTree(pgcurd.tree, "segment", title="URD tree segment",label.segments = T)
plotDim(pgcurd.tree, "segment", plot.title="URD tree segment")
dev.off()
saveRDS(pgcurd.tree,file = file.ou.pgcurd.tree)
#### Force-directed layout
# 
# visitation=data.frame(cell=rownames(pgcurd.tree@diff.data),seg=pgcurd.tree@diff.data$segment,stringsAsFactors = F,row.names = rownames(pgcurd.tree@diff.data))
# visitation$visit=log10(apply(visitation,1,function(cr) pgcurd.tree@diff.data[as.character(cr['cell']),paste0('visitfreq.raw.',as.character(cr['seg']))])+1)
# robustly.visited.cells=visitation[visitation$visit>=0.5,'cell']
# final.tips=segTerminal(pgcurd.tree)
# # Generate the force-directed layout
# pgcurd.tree <- treeForceDirectedLayout(pgcurd.tree, num.nn=100, cut.unconnected.segments=2, verbose=T,method = 'fr',cells.to.do = robustly.visited.cells,tips=final.tips)

pgcurd.tree <- treeForceDirectedLayout(pgcurd.tree, num.nn=100, cut.unconnected.segments=2, verbose=T)
saveRDS(pgcurd.tree, file=file.ou.pgcurd.tree.ForceLayout)
# pdf("pgcurd.forceDiretedTree.pdf")
plotTreeForce(pgcurd.tree, "NANOS3", title = "NANOS3", title.cex = 2, title.line=2.5)
plotTreeForceStore3DView(pgcurd.tree, "NANOS3")
plotTreeForce(pgcurd.tree, "NANOG", title = "NANOG", title.cex=2, title.line=2.5)

plotTreeForce(pgcurd.tree, "SOX17", title="SOX17", title.cex=2, title.line=2.5)
# dev.off()




pdf(file='URD.test.pdf',width = 8,height = 7)
plotDim(pgcurd, "stage")
plotDim(pgcurd, "pseudotime")
plotDim(pgcurd, "NANOS3", plot.title="NANOS3")
plotDim(pgcurd, "NANOG", plot.title="NANOG")
plotDim(pgcurd, "SOX17", plot.title="SOX17")


plotDists(pgcurd, "pseudotime", "stage", plot.title="Pseudotime by stage")
dev.off()
# Find tips

pgcurd.d4 <- urdSubset(pgcurd, cells.keep=cellsInCluster(pgcurd, "stage", "U1PGCLCd4"))
pgcurd.d4@var.genes <- var.genes.by.stage[[6]]

# Calculate PCA and tSNE


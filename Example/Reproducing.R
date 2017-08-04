# reproducing Figure 1
library('igraph')
karate <- make_graph("Zachary")
A <- as_adjacency_matrix(karate)
W <- as.matrix(A)
moduleid <- c(9,10,31,34,32,25,23,28,29,15,26,24,30,16,27,19,21,33)

n <- dim(W)[1]
z <- runif(n, min=0, max=1)
z[moduleid] <- runif(length(moduleid), min=0.5, max=1)

## reproducing Figure 1 (a)
V(karate)$color="gray"
V(karate)[moduleid]$color="red"
plot(karate,layout=layout_with_kk(karate))

## reproducing Figure 1 (b)
k <- colSums(W)
D <- diag(k)
L <- D-W
eg <- eigen(L)
u <- eg$vectors[,2]
moduleid <- which(u>0)
V(karate)$color="gray"
V(karate)[moduleid]$color="red"
plot(karate,layout=layout_with_kk(karate))

## reproducing Figure 1 (c)
source('../ModuleExtraction/R/ModuleExtract.R')
pp <- projectedgradient(W,z,lambda=1,maxiter=100)
func <- pp[[1]]
x <- pp[[2]]
predictedid=which(x>0)
V(karate)$color="gray"
V(karate)[predictedid]$color="red"
plot(karate,vertex.size=20*abs((z/max(z))),layout=layout_with_kk(karate))
plot(func,xlab="iteration",ylab="objectiv function")

## rank and select
ps <- sort(z,decreasing=T,index.return=T)
g2 <- ps$ix[1:18]
predictedid=g2
V(karate)$color="gray"
V(karate)[predictedid]$color="red"
plot(karate,vertex.size=20*abs((z/max(z))),layout=layout_with_kk(karate))


# reproducing Figure 2
library(BioNet)
library(DLBCL)
data(dataLym)
data(interactome)
source('../ModuleExtraction/R/ModuleExtract.R')

pvals <- cbind(t=dataLym$t.pval, s=dataLym$s.pval)
rownames(pvals) <- dataLym$label
pval <- aggrPvals(pvals, order=2, plot=FALSE)
subnet <- subNetwork(dataLym$label, interactome)
subnet <- rmSelfLoops(subnet)
subnet

fb <- fitBumModel(pval, plot=FALSE)
scores <- scoreNodes(subnet, fb, fdr=0.001)

module <- runFastHeinz(subnet, scores)
logFC <- dataLym$diff
names(logFC) <- dataLym$label

plotModule(module, scores=scores, diff.expr=logFC)

####################

library(igraph)
subnetg <- igraph.from.graphNEL(subnet)
z <- scores
A <- as_adjacency_matrix(subnetg)
W <- as.matrix(A)

tmpg <- ModuleExtraction(subnetg, z, 100)

moduleg <- igraph.to.graphNEL(tmpg)

############## Plot combined modules Figure
Allnodes <- union(nodes(module),nodes(moduleg))

w =W[Allnodes,Allnodes]
z = scores[Allnodes]
g <- graph.adjacency(w, mode="undirected", weighted=TRUE)

midap=intersect(nodes(module),nodes(moduleg))  # in module and moduleg
leftap=setdiff(nodes(module),midap) #only in module
rightap=setdiff(nodes(moduleg),midap) #only in moduleg
V(g)$color[match(midap,Allnodes)] = 'red'
V(g)$color[match(leftap,Allnodes)] = 'yellow'
V(g)$color[match(rightap,Allnodes)] = 'green'
V(g)$shape[which(z<0)]='square'
V(g)$shape[which(z>0)]='circle'

plot(g,layout=layout_with_kk(g),vertex.size=5,
	vertex.label.dist=0.5,vertex.label.cex = 0.5)

# reproducing Figure 3
library(ggplot2)
alld <- read.csv('stats.csv', header = FALSE)
tmprname <- alld[,1]
alld <- alld[,2:7]
rownames(alld) <- as.character(tmprname)
colnames(alld) <- c('CYC-P','CYC-R','CYC-F','MIPS-P','MIPS-R','MIPS-F')
dat <- data.frame(
    setting = factor(rownames(alld)),
    metric = factor(as.vector(t(replicate(dim(alld)[1],colnames(alld)))),
    	levels=colnames(alld)),
    value = as.vector(as.matrix(alld))
)
ggplot(data=dat, aes(x=metric, y=value, fill=setting))+
geom_bar(stat="identity", position=position_dodge())+
scale_fill_manual(values=rainbow(length(rownames(alld))))
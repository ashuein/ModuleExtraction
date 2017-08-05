# mips_3_100.txt is from www.paccanarolab.org/static_content/clusterone/cl1_gold_standard.zip
rlines=readLines('data/mips_3_100.txt')
MIPS=vector("list", length(rlines))
modsizeMIPS=c()
for (i in 1:length(rlines)) {
    MIPS[[i]]=strsplit(rlines[i],' ')[[1]]
    modsizeMIPS = c(modsizeMIPS,length(MIPS[[i]]))
}

# CYC2008_complex.tab is from http://wodaklab.org/cyc2008/downloads
rlines=read.delim('data/CYC2008_complex.tab')
CYCcomplexlable = unique(rlines[,3])
CYC= vector("list", length(CYCcomplexlable))

modsizeCYC=numeric(length=length(CYCcomplexlable))
for (i in 1:dim(rlines)[1]) {
	j = match(rlines[i,3],CYCcomplexlable)
    CYC[[j]]=c(CYC[[j]],as.character(rlines[i,1]))
    modsizeCYC[j] = modsizeCYC[j]+1
}

# B is the reference set
FScore <- function(A,B){
	m = length(A)
	n = length(B)

	pre = 0
	for (i in 1:m) {
		for (j in 1:n) {
			if ( (length(intersect(A[[i]],B[[j]])))^2/length(union(A[[i]],B[[j]])) > 0.25){
				pre = pre+1
				break
			}

		}
	}
	pre = pre/m

	rec = 0
	for (i in 1:n) {
		for (j in 1:m) {
			if ( (length(intersect(B[[i]],A[[j]])))^2/length(union(B[[i]],A[[j]])) > 0.25){
				rec = rec+1
				break
			}

		}
	}
	rec = rec/n
	f = 2*pre*rec/(pre+rec)
	print(cbind(pre,rec,f))
	return (cbind(pre,rec,f))
}

# F-measure for each time point
for (t in 1:12) {
rlines=readLines(paste('PPI_complexes',t,'.txt',sep=''))
PPIcomplexes=vector("list", length(rlines))
modsize = c()
for (i in 1:length(rlines)) {
    ap=strsplit(rlines[i],'\t')[[1]]
    PPIcomplexes[[i]]=ap
    modsize = c(modsize,length(ap))
}

fp1=FScore(PPIcomplexes,CYC)
fp2=FScore(PPIcomplexes,MIPS)
write.table(cbind(fp1,fp2),'stats.csv',row.names=paste('T',t,sep=''),col.names=F,sep=',',append=T)
}

# F-measure for consensus graph
rlines=readLines('PPI_complexes.txt')
PPIcomplexes=vector("list", length(rlines))
modsize = c()
for (i in 1:length(rlines)) {
    ap=strsplit(rlines[i],'\t')[[1]]
    PPIcomplexes[[i]]=ap
    modsize = c(modsize,length(ap))
}

fp1=FScore(PPIcomplexes,CYC)
fp2=FScore(PPIcomplexes,MIPS)
write.table(cbind(fp1,fp2),'stats.csv',row.names='Consensus',col.names=F,sep=',',append=T)

# GO enrichment analysis
library(STRINGdb)
string_db <- STRINGdb$new( version="10", species=4932,
    score_threshold=0, input_directory="" )
library(xlsx)
rlines=readLines('PPI_complexes.txt')
modulegoscore = matrix(0,nrow=length(rlines),ncol=2)
for (i in 1:length(rlines)) {
    # save the modules gene lists
    gsymbol=strsplit(rlines[i],'\t')[[1]]
    modsize = c(modsize,length(gsymbol))

    # access STRING
    diff_exp_genes <- data.frame(gsymbol,gsymbol)
    colnames(diff_exp_genes)<-c('gene','backup')
    genes_mapped <- string_db$map(diff_exp_genes, "gene", removeUnmappedRows = TRUE )
    hits <- genes_mapped$STRING_id
    
    er <- string_db$get_enrichment(hits)
    er<-er[which(er$pvalue_fdr<0.05),]
    modulegoscore[i,1]=dim(er)[1]
    if (dim(er)[1] > 0){
        write.xlsx(x = er, file = paste(fhead,"/EnrichmentBP.xlsx",sep=''),append = TRUE,
        sheetName = paste('module',i,sep=''), row.names = FALSE, col.names=TRUE)
    }

    er <- string_db$get_enrichment(hits)
    er<-er[which(er$pvalue_fdr<0.05),]
    modulegoscore[i,2]=dim(er)[1]
    if (dim(er)[1] > 0){
        write.xlsx(x = er, file = paste(fhead,"/EnrichmentBP.xlsx",sep=''),append = TRUE,
        sheetName = paste('module',i,sep=''), row.names = FALSE, col.names=TRUE)
    }
}

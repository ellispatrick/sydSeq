countsMi = as.matrix(read.csv('../../Data/pMim/miRNAcounts.csv',row.names = 1))
load("ensembl.gene.RData")
mir = ensembl.gene[ensembl.gene$ensembl_gene_id %in% rownames(countsMi), ]$mirbase_id
names(mir) = ensembl.gene[ensembl.gene$ensembl_gene_id %in% rownames(countsMi), ]$ensembl_gene_id
mir = unlist(lapply(strsplit(mir, "-"), function(x) if (length(x) > 2) paste(x[1:3], 
                                                                             collapse = "-")))
smir = split(names(mir), mir)
lmir = unlist(lapply(smir, length))
cMi = countsMi[unlist(smir[lmir == 1]), ]
rownames(cMi) = names(which(lmir == 1))
y = NULL
for (i in names(which(lmir > 1))) {
  x = (countsMi[smir[[i]], ])
  x = x[which.max(rowSums(x)),]
  y = rbind(y, x)
}
rownames(y) = names(which(lmir > 1))
cMi = rbind(cMi, y)
countsMi = cMi

counts = as.matrix(read.csv('../../Data/pMim/mRNAcounts.csv',row.names = 1))

save(counts,countsMi,file = 'counts.RData')


library(clusterProfiler)
library(org.Hs.eg.db)

#####################################################
data(geneList, package="DOSE")
gene <- names(geneList)[abs(geneList) > 2]
gene.df <- bitr(gene, fromType = "ENTREZID",
                toType = c("ENSEMBL", "SYMBOL"),
                OrgDb = org.Hs.eg.db)
head(gene.df)

# groupGO is designed for gene classification 
# based on GO distribution at a specific level
ggo <- groupGO(gene     = gene,
               OrgDb    = org.Hs.eg.db,
               ont      = "CC",
               level    = 3,
               readable = TRUE)
head(ggo)


######################################################
gmtfile <- system.file("extdata", "c5.cc.v5.0.entrez.gmt", package="clusterProfiler")
c5 <- read.gmt(gmtfile)
gene
head(c5)
egmt <- enricher(gene, TERM2GENE=c5)
head(egmt)

egmt2 <- GSEA(geneList, TERM2GENE=c5, verbose=FALSE)
head(egmt2)


######################################################
datAnno <- data.frame(
    ont=sample(LETTERS[1:5], 100, replace=T),
    gene=1:100)

lst <- split(datAnno$gene, datAnno$ont)
datItem <- unique(c(sample(lst$A, 10), sample(1:100, 20)))

tmp <- enricher(datItem, TERM2GENE=datAnno)
head(tmp)

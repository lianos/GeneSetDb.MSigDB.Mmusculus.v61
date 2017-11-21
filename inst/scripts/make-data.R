# A script describing the steps involved in making the data object(s). Output of
# this script should be files on disk ready to be pushed to S3.

## This scripts requires that the GeneSetDb.MSigDB.Hsapiens data package
## is already built. It will then use biomaRt to map human entrez identifieres
## to mouse.

## Convert GeneSetDb to mouse --------------------------------------------------
## We use the biomaRt::getLDS (get linked datasets) function to map human
## entrez ids to mouse. Look at ?getLDS function in biomaRt
library(multiGSEA)
library(biomaRt)
library(dplyr)

gdb <- readRDS("path/to/GeneSetDb.MSigDb.Hsapiens/inst/extdata/GeneSetDb.MSigDb.Hsapiens.rds")
hdf <- gdb %>%
  as.data.frame %>%
  group_by(collection, name) %>%
  mutate(n=n()) %>%
  ungroup
hids <- hdf %>%
  select(entrezgene=featureId, symbol) %>%
  distinct(entrezgene, .keep_all=TRUE) %>%
  mutate(featureId=entrezgene)

human <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
mouse <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")
xref <- getLDS(attributes = c("hgnc_symbol", "entrezgene", "ensembl_gene_id"),
               filters = "entrezgene", values = hids$entrezgene, mart = human,
               attributesL = c("ensembl_gene_id", "entrezgene", "mgi_symbol"),
               martL = mouse)
xfer <- xref %>%
  select(featureId=NCBI.gene.ID, mm_featureId=NCBI.gene.ID.1, mm_symbol=MGI.symbol) %>%
  distinct(featureId, mm_featureId, .keep_all = TRUE) %>%
  mutate(featureId=as.character(featureId), mm_featureId=as.character(mm_featureId))

mdf <- hdf %>%
  filter(collection != 'c1') %>%
  inner_join(xfer, by="featureId") %>%
  group_by(collection, name) %>%
  mutate(nm=n()) %>%
  ungroup

mdf.go <- mdf %>%
  select(collection, name, featureId=mm_featureId, symbol=mm_symbol) %>%
  filter(!is.na(featureId), !is.na(symbol),
         nchar(featureId) > 0, nchar(symbol) > 0)

mgdb <- GeneSetDb(mdf.go)
## add metadata from gdb@table
take.cols <- c('collection', 'name', setdiff(colnames(gdb@table), colnames(mgdb@table)))
meta <- gdb@table[mgdb@table, take.cols, with=FALSE]
mnew <- mgdb@table[meta]
stopifnot(
  all.equal(mgdb@table[, list(collection, name)], mnew[, list(collection, name)]),
  all.equal(mgdb@table$N, mnew$N))
mgdb@table <- mnew

url.fn <- function(collection, name) {
  url <- "http://www.broadinstitute.org/gsea/msigdb/cards/%s.html"
  sprintf(url, name)
}
for (col in unique(geneSets(mgdb)$collection)) {
  geneSetCollectionURLfunction(mgdb, col) <- url.fn
  featureIdType(mgdb, col) <- EntrezIdentifier()
  mgdb <- addCollectionMetadata(mgdb, col, 'source', 'MSigDB_v6.1')
}

org(mgdb) <- 'Mus_musculus'
# gdb.fn <- sprintf('MSigDB.Mus_musculus.GeneSetDb.rds', species)
gdb.fn <- '../extdata/GeneSetDb.MSigDB.Mmusculus.v61.rds'
saveRDS(mgdb, gdb.fn)

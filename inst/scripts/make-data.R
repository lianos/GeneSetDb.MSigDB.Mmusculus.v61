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
gdb.fn <- '../extdata/GeneSetDb.MSigDB.Mmusculus-entrez.v61.rds'
saveRDS(mgdb, gdb.fn)

# Create Ensembl version
# library(multiGSEA)
# library(biomaRt)
# library(dplyr)
# library(GSEABase)
#
# mgdb <- readRDS("inst/extdata/GeneSetDb.MSigDB.Mmusculus-entrez.v61.rds")
mdf <- mgdb@db
mart <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")
exref <- getBM(
  attributes = c("entrezgene", "ensembl_gene_id", "mgi_symbol"),
  filters = "entrezgene",
  values = unique(mdf$featureId),
  mart = mart) %>%
  transmute(entrezgene = as.character(entrezgene),
            featureId = ensembl_gene_id, symbol = mgi_symbol)

mdf.ens <- mdf %>%
  select(collection, name, entrezgene = featureId) %>%
  inner_join(exref, by = "entrezgene") %>%
  arrange(collection, name, symbol) %>%
  distinct(collection, name, featureId, .keep_all = TRUE) %>%
  select(-entrezgene)

gdb2 <- GeneSetDb(mdf.ens)
take.cols <- c('collection', 'name', setdiff(colnames(mgdb@table), colnames(gdb2@table)))
meta <- mgdb@table[gdb2@table, take.cols, with=FALSE]
mnew <- gdb2@table[meta]
stopifnot(
  all.equal(gdb2@table[, list(collection, name)], mnew[, list(collection, name)]),
  all.equal(gdb2@table$N, mnew$N))
gdb2@table <- mnew

url.fn <- function(collection, name) {
  url <- "http://www.broadinstitute.org/gsea/msigdb/cards/%s.html"
  sprintf(url, name)
}
for (col in unique(geneSets(gdb2)$collection)) {
  geneSetCollectionURLfunction(gdb2, col) <- url.fn
  featureIdType(gdb2, col) <- ENSEMBLIdentifier()
  gdb2 <- addCollectionMetadata(gdb2, col, 'source', 'MSigDB_v6.1')
}

org(gdb2) <- 'Mus_musculus'
gdb.fn <- '../extdata/GeneSetDb.MSigDB.Mmusculus-ensembl.v61.rds'
saveRDS(gdb2, gdb.fn)

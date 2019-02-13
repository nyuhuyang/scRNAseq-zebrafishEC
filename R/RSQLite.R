# make a new EnsDb.Drerio.v95 pacakge from EnsDb.Drerio.v95.sqlite

#https://annotationhub.bioconductor.org/rdataclass/EnsDb search "Danio rerio"
# download Ensembl 95 EnsDb for Danio rerio
file.copy("~/Downloads/EnsDb.Drerio.v95.sqlite", "./data")
#https://datacarpentry.org/R-ecology-lesson/05-r-and-databases.html
library(dplyr)
library(dbplyr)
library(DBI)
Drerio.sqlite <- DBI::dbConnect(RSQLite::SQLite(), "data/EnsDb.Drerio.v95.sqlite")
tbls <- src_dbi(Drerio.sqlite);tbls
tbls <- c("chromosome", "entrezgene", "exon", "gene", "metadata", "protein", 
          "protein_domain", "tx", "tx2exon","uniprot")
df_tbls <- lapply(tbls.names, function(x) tbl(Drerio.sqlite,x))
names(df_tbls) = tbls
df_tbls

Drerio_gene <- dbGetQuery(Drerio.sqlite, 'SELECT * FROM gene')
head(Drerio_gene,2)
table(Drerio_gene$seq_name) %>% tail
#======================
Drerio_metadata <- dbGetQuery(Drerio.sqlite, 'SELECT * FROM metadata')


if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install("MeSHDbi", version = "3.8")

library(MeSHDbi)
data(PAO1)
head(PAO1)
data(metaPAO1)
metaPAO1

MeSHDbi_df <- Drerio_gene[,c("gene_id","gene_name","gene_biotype",
                             "seq_name","description"),]
colnames(MeSHDbi_df) = colnames(PAO1)
colnames(Drerio_metadata) = colnames(metaPAO1)

makeGeneMeSHPackage(pkgname = "EnsDb.Drerio.v95",
                    data = MeSHDbi_df,
                    metadata = Drerio_metadata,
                    organism = "Danio rerio",
                    version = "v95",
                    maintainer = "Yang Hu <yah2014@med.cornell.edu>",
                    author = "Yang",
                    destDir = "./data",
                    license="Artistic-2.0")

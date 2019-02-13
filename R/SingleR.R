library(SingleR)
library(Seurat)
library(reshape2)
library(pheatmap)
library(kableExtra)
library(dplyr)
library(tidyr)
source("../R/Seurat_functions.R")
source("../R/SingleR_functions.R")
path <- paste0(getwd(),"/output/",gsub("-","",Sys.Date()),"/")
if(!dir.exists(path)) dir.create(path, recursive = T)
#====== 3.1 Create Singler Object  ==========================================
(load(file="data/zebrafishEC_1_20190212.Rda"))
(load(file="../SingleR/data/mouse.rnaseq.RData"))
(load(file="../SingleR/data/Blueprint_encode.RData"))
object@scale.data = NULL
object_data <- as.matrix(object@data)
rownames(object_data) = Hmisc::capitalize(tolower(rownames(object_data)))
GC();GC();GC();GC();GC();GC();GC();GC();GC();GC();GC();GC();GC();GC();
singler.m = CreateSinglerObject(object_data, annot = NULL, project.name="zebrafishEC",
                              min.genes = 500,technology = "10X", species = "Mouse", citation = "",
                              ref.list = list(mouse.rnaseq),normalize.gene.length = F, variable.genes = "de",
                              fine.tune = T, do.signatures = F, clusters = NULL,
                              numCores = SingleR.numCores/4)
rownames(object_data) = toupper(rownames(object_data))
GC();GC();GC();GC();GC();GC();GC();GC();GC();GC();GC();GC();GC();GC();
singler.h = CreateSinglerObject(object_data, annot = NULL, project.name="zebrafishEC",
                                min.genes = 500,technology = "10X", species = "Human", citation = "",
                                ref.list = list(Blueprint_encode),normalize.gene.length = F, variable.genes = "de",
                                fine.tune = T, do.signatures = F, clusters = NULL,
                                numCores = SingleR.numCores/4)
# if singler didn't find all cell labels
print(length(singler.h$singler[[1]]$SingleR.single$labels) == ncol(object@data))
remove(object_data)
# if singler didn't find all cell labels
if(length(singler.h$singler[[1]]$SingleR.single$labels) != ncol(object@data)){
        all.cell = object@cell.names;length(all.cell)
        know.cell = rownames(singler.h$singler[[1]]$SingleR.single$labels);length(know.cell)
        object = SubsetData(object, cells.use = know.cell)
}

singler.h$meta.data$orig.ident = object@meta.data$orig.ident # the original identities, if not supplied in 'annot'
singler.h$meta.data$xy = object@dr$tsne@cell.embeddings # the tSNE coordinates
singler.h$meta.data$clusters = object@ident # the Seurat clusters (if 'clusters' not provided)

singler.m$meta.data$orig.ident = object@meta.data$orig.ident # the original identities, if not supplied in 'annot'
singler.m$meta.data$xy = object@dr$tsne@cell.embeddings # the tSNE coordinates
singler.m$meta.data$clusters = object@ident # the Seurat clusters (if 'clusters' not provided)

save(singler.h,singler.m,file="./output/singler_hm_zebrafishEC_1T_20190212.Rda")

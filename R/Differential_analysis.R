########################################################################
#
#  0 setup environment, install libraries if necessary, load libraries
# 
# ######################################################################

library(Seurat)
library(dplyr)
library(tidyr)
library(kableExtra)
library(magrittr)
library(gplots)
source("../R/Seurat_functions.R")
path <- paste0("./output/",gsub("-","",Sys.Date()),"/")
if(!dir.exists(path)) dir.create(path, recursive = T)
#3.1  Compare DE across all major cell types==================
#We would need the data for all clusters, as well the subclusters.
#detect changes in gene expression between young and aged, 
#in the different cell types and subtypes. 
#It will also be interesting to check if there is some subtype enriched in young compared to aged or viceversa. 

# 3.1.1 load data
(load(file="data/zebrafishEC_1_20190212.Rda"))

# cell population
table(object@meta.data$singler1main, object@meta.data$orig.ident) %>% 
        as.data.frame.matrix %>% kable %>%
        kable_styling()

###############################
# Doheatmap
###############################
subset.object <- SubsetData(object, ident.use = c(1,2,5,6,7,8))
gde.markers <- FindAllMarkers.UMI(subset.object, logfc.threshold = 0.25,only.pos = T,
                                  return.thresh = 0.05)
DoHeatmap.1(subset.object, gde.markers, Top_n = 30,
            use.scaled = T,ident.use = "Endothelial and Hematopoietic cells",
            group.label.rot = F,cex.row = 3,remove.key =F,title.size = 12,
            do.print = T)
write.csv(gde.markers,paste(path,"Top_enriched_genes_in_cluster1_2_5_6_7_8.csv"))
all.markers <- FindAllMarkers.UMI(object, logfc.threshold = 0.05,only.pos = F,
                                  return.thresh = 0.05)
write.csv(all.markers,paste(path,"all.markers.csv"))
###############################
# Doheatmap for cluster.1 / cluster.2
###############################
object %<>% SetAllIdent(id = "res.0.6")
#---FindAllMarkers.UMI---- 
pair.markers <- FindPairMarkers(subset.object, 
                               ident.1 = c(1,2,2,5), 
                               ident.2 = c(6,5,6,8),
                               only.pos = F,
                               logfc.threshold = 0.25,min.cells.group =3,
                               min.pct = 0.1,
                               return.thresh = 0.1,save.files = TRUE,
                               save.path = path)

###############################
# Sub-cluster 2 and 6 as 3 groups
###############################
subset.object <- SubsetData(object, ident.use = c(2,6))
system.time({
        subset.object %<>% FindClusters(reduction.type = "pca", resolution = 0.2, dims.use = 1:30,
                             save.SNN = TRUE, n.start = 10, nn.eps = 0.5,
                             force.recalc = TRUE, print.output = FALSE)
})
TSNEPlot.1(subset.object,do.label = T,do.print=T)
subset.object@ident <- plyr::mapvalues(x = subset.object@ident,
                                from = c(0,1,3,2),
                                to = c(0,1,2,3))
subset.object@ident %<>% factor(levels = 0:3)
gde.markers1 <- FindAllMarkers.UMI(subset.object, logfc.threshold = 0.25,only.pos = T,
                                  return.thresh = 0.05)
DoHeatmap.1(subset.object, gde.markers1, Top_n = 40,
            use.scaled = T,ident.use = "Endothelial cells ",
            group.label.rot = F,cex.row = 3,remove.key =F,title.size = 12,
            do.print = T)
write.csv(gde.markers1,paste(path,"Top_enriched_genes_in_EC_cluster4.csv"))

Read10X()
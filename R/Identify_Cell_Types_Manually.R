library(Seurat)
library(dplyr)
library(tidyr)
library(kableExtra)
library(magrittr)
library(gplots)
source("../R/Seurat_functions.R")
source("../R/SingleR_functions.R")
path <- paste0("./output/",gsub("-","",Sys.Date()),"/")
if(!dir.exists(path)) dir.create(path, recursive = T)
#====== 2.1 identify phenotype for each cluster  ==========================================
(load(file="data/zebrafishEC_1_20190212.Rda"))

#blueprint_encode_main = read.csv("../SingleR/output/blueprint_encode_main.csv",row.names =1,header = T,
#                                 stringsAsFactors = F)
df_markers <- readxl::read_excel("doc/zebrafish.markers.xlsx")
colnames(df_markers) = gsub(" ","_",colnames(df_markers))
colnames(df_markers) = gsub(":|\\/","_",colnames(df_markers))
colnames(df_markers) = gsub("\\+","",colnames(df_markers))
markers = df_markers[,-grep("Alias",colnames(df_markers))]
marker.list <- df2list(markers)

marker.list %<>% lapply(function(x) x[1:16]) %>% 
        lapply(function(x) HumanGenes(object,x)) %>% 
        lapply(function(x) x[1:9]) %>%
        lapply(function(x) x[!is.na(x)]) %>%
        .[lapply(.,length) >0]

marker.list %>% list2df %>% t %>% kable() %>% kable_styling()
dev.off()
for(i in 1:length(marker.list)){
        p <- lapply(marker.list[[i]], function(marker) {
                SingleFeaturePlot.1(object = object, feature = marker,pt.size = 0.5,
                                    gradient.use = c("lightblue", "blue3"),threshold=0.1)+
                        ggtitle(paste0(marker))+
                        theme(plot.title = element_text(hjust = 0.5,size = 15, face = "bold"))
                })

        jpeg(paste0(path,names(marker.list)[i],".jpeg"),
             units="in", width=10, height=7,res=600)
        print(do.call(plot_grid, p)+ ggtitle(paste(names(marker.list)[i],"markers"))+
                      theme(plot.title = element_text(hjust = 0.5,size = 20, face = "bold")))
        dev.off()
        print(paste0(i,":",length(marker.list)))
}


# change labels
system.time({
        object %<>% FindClusters(reduction.type = "pca", resolution = 1.2, dims.use = 1:30,
                             save.SNN = TRUE, n.start = 10, nn.eps = 0.5,
                             force.recalc = TRUE, print.output = FALSE)
})

TSNEPlot(object, do.label = T)
idents <- as.data.frame(table(object@ident))
(old.ident.ids <- idents$Var1)
new.cluster.ids <- c(10,2,4,1,9,7,7,6,2,10,3,5,13,8,15,16,17,18, ,, , )
old.ident.ids <-    c(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20)
object@ident <- plyr::mapvalues(x = object@ident,
                                from = old.ident.ids,
                                to = new.cluster.ids)
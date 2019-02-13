library(SingleR)
library(Seurat)
library(reshape2)
library(pheatmap)
library(kableExtra)
library(dplyr)
library(tidyr)
source("../R/Seurat_functions.R")
source("../R/SingleR_functions.R")
path <- paste0("./output/",gsub("-","",Sys.Date()),"/")
if(!dir.exists(path)) dir.create(path, recursive = T)

#====== 3.2 SingleR specifications ==========================================
# Step 1: Spearman coefficient
(load(file="data/zebrafishEC_1_20190212.Rda"))
(load(file="output/singler_hm_zebrafishEC_1T_20190212.Rda"))
# if singler didn't find all cell labels
if(length(singler$singler[[1]]$SingleR.single$labels) != ncol(object@data)){
        all.cell = object@cell.names;length(all.cell)
        know.cell = rownames(singler$singler[[1]]$SingleR.single$labels);length(know.cell)
        object = SubsetData(object, cells.use = know.cell)
}
object
##############################
# add singleR label to Seurat
###############################
singlerDF = data.frame("singler1sub" = singler.m$singler[[1]]$SingleR.single$labels,
                       "singler1main" = singler.m$singler[[1]]$SingleR.single.main$labels,
                       "singler2sub" = singler.h$singler[[1]]$SingleR.single$labels,
                       "singler2main" = singler.h$singler[[1]]$SingleR.single.main$labels,
                       row.names = rownames(singler.h$singler[[1]]$SingleR.single$labels))

table(rownames(singlerDF) %in% object@cell.names)

apply(singlerDF,2,function(x) length(unique(x)))
object <- AddMetaData(object = object,metadata = singlerDF)
object <- SetAllIdent(object = object, id = "singler1main")
TSNEPlot.1(object, do.label = T)
object <- SetAllIdent(object = object, id = "singler1sub")
TSNEPlot.1(object, do.label = T)

##############################
# check the spearman correlation
###############################
#Or by all cell types (showing the top 50 cell types):
jpeg(paste0(path,"DrawHeatmap_sub1.jpeg"), units="in", width=10, height=7,
     res=600)
print(SingleR.DrawHeatmap(singler.m$singler[[1]]$SingleR.single,top.n = 50,normalize = T))
dev.off()
jpeg(paste0(path,"DrawHeatmap_main1.jpeg"), units="in", width=10, height=7,
     res=600)
print(SingleR.DrawHeatmap(singler.m$singler[[1]]$SingleR.single.main, top.n = 50,normalize = F))
dev.off()

#Finally, we can also view the labeling as a table compared to the original identities:
object@meta.data$singler1main %>% table() %>% kable() %>% kable_styling()
object@meta.data$singler1sub %>% table() %>% kable() %>% kable_styling()
object@meta.data$singler2main %>% table() %>% kable() %>% kable_styling()
##############################
# process color scheme
##############################
singler_colors <- readxl::read_excel("./doc/singler.colors.xlsx")
singler.colors <- sapply(paste0("singler.color",1:3),function(x) {
        singler_colors[!is.na(singler_colors[,x]),x]}) %>%
        sapply(function(x) x[!duplicated(x)])
lapply(singler.colors,length)

apply(object@meta.data[,c("singler1main","singler1sub","singler2main")],
      2,function(x) length(unique(x)))
object <- AddMetaColor(object = object, label= "singler1main", colors = singler.colors[[1]][1:14])
object <- AddMetaColor(object = object, label= "singler1sub", colors = singler.colors[[2]][1:21])
object <- AddMetaColor(object = object, label= "singler2main", colors = singler.colors[[3]][1:14])
object <- SetAllIdent(object = object, id = "singler1main")
TSNEPlot.1(object, colors.use = ExtractMetaColor(object),no.legend = F)
##############################
# draw tsne plot
##############################
p3 <- TSNEPlot.1(object = object, do.label = T, group.by = "ident",
                 do.return = TRUE, no.legend = T,
                 colors.use = ExtractMetaColor(object),
                 pt.size = 1,label.size = 4,force = 2)+
  ggtitle("Supervised cell type labeling by mouse.rnaseq")+
  theme(text = element_text(size=10),							
        plot.title = element_text(hjust = 0.5,size = 18, face = "bold")) 

jpeg(paste0(path,"PlotTsne_singler1main.jpeg"), units="in", width=10, height=7,
     res=600)
print(p3)
dev.off()

save(object, file = "data/zebrafishEC_1_20190212.Rda")

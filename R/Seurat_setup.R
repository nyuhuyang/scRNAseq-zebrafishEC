########################################################################
#
#  0 setup environment, install libraries if nLynchessary, load libraries
# 
# ######################################################################

library(Seurat)
library(magrittr)
library(harmony)
library(dplyr)
library(kableExtra)
source("../R/Seurat_functions.R")
path <- paste0("./output/",gsub("-","",Sys.Date()),"/")
if(!dir.exists(path)) dir.create(path, recursive = T)
if(!dir.exists("./data/")) dir.create("data")
########################################################################
#
#  1 Data preprocessing
# 
# ######################################################################
#======1.1 Load the data files and Set up Seurat object =========================
# read sample summary list
df_samples <- readxl::read_excel("doc/190125_scRNAseq_info.xlsx")
colnames(df_samples) <- colnames(df_samples) %>% tolower
sample_n = which(df_samples$tests %in% c("test1"))
df_samples[sample_n,] %>% kable() %>% kable_styling()
table(df_samples$tests);nrow(df_samples)
samples <- df_samples$sample[sample_n]
sample.id <- df_samples$sample.id[sample_n]
conditions <- df_samples$conditions[sample_n]
projects <- df_samples$project[sample_n]
tissues <- df_samples$tissue[sample_n]
tests <- df_samples$tests[sample_n] 

#======1.2 load  SingleCellExperiment =========================
(load(file = "data/sce_1_20190212.Rda"))
names(sce_list)
object_list <- lapply(sce_list, as.seurat) %>%
        lapply(NormalizeData) %>%
        #lapply(ScaleData) %>%
        lapply(FindVariableGenes, do.plot = FALSE)

for(i in 1:length(samples)){
        object_list[[i]]@meta.data$tests <- tests[i]
        object_list[[i]]@meta.data$conditions <- conditions[i]
        #object_list[[i]]@meta.data$projects <- projects[i]
        #object_list[[i]]@meta.data$notes <- notes[i]
        object_list[[i]]@meta.data$tissues <- tissues[i]
        
}
# we will take the union of the top 1k variable genes in each dataset for alignment
genes.use <- object_list %>% 
        lapply(function(object) head(rownames(object@hvg.info), 2000)) %>%
        unlist %>% unique
length(genes.use)

#========1.3 merge ===================================
object <- Reduce(function(x, y) MergeSeurat(x, y, do.normalize = F), object_list)
object@var.genes = genes.use
remove(sce_list,object_list);GC()

object = SetAllIdent(object, id = "orig.ident")
#======1.4 mito, QC, filteration =========================
object@meta.data$percent.mito = object@meta.data$pct_counts_Mito/100

(remove <- which(colnames(object@meta.data) %in%c("is_cell_control",
                                           "pct_counts_in_top_500_features_Mito")))
meta.data = object@meta.data[,-seq(remove[1], remove[2], by=1)]
object@meta.data = meta.data 

(load(file = paste0("output/20190212/g1_1_20190212.Rda")))

object <- FilterCells(object = object, subset.names = c("nGene","nUMI","percent.mito"),
                   low.thresholds = c(500,3000, -Inf), 
                   high.thresholds = c(Inf,Inf, 0.5))

object@ident = factor(object@ident,levels = samples)
g2 <- lapply(c("nGene", "nUMI", "percent.mito"), function(features){
        VlnPlot(object = object, features.plot = features, nCol = 3, 
                point.size.use = 0.2,size.x.use = 10, group.by = "ident",
                x.lab.rot = T, do.return = T)
})
save(g2,file= paste0(path,"g2_1_20190212.Rda"))

jpeg(paste0(path,"/S1_nGene_nUMI.jpeg"), units="in", width=10, height=7,res=600)
print(plot_grid(g2[[1]]+ggtitle("nGene")+ 
                        scale_y_log10(limits = c(500,10000)),
                g2[[2]]+ggtitle("nUMI")+ 
                        scale_y_log10(limits = c(3000,100000))))
dev.off()

#======1.5 FindVariableGenes=======================
object <- NormalizeData(object = object)
jpeg(paste0(path,"/S1_dispersion.jpeg"), units="in", width=10, height=7,res=600)
object <- FindVariableGenes(object = object, mean.function = ExpMean, 
                            dispersion.function = LogVMR, do.plot = T, 
                            x.low.cutoff = 0.025, x.high.cutoff = 8, y.cutoff = 0.25)
dev.off()
length(object@var.genes)
table(object@var.genes %in% genes.use)

#======1.5 Add Cell-cycle score =========================
# Read in a list of cell cycle markers, from Tirosh et al, 2015
cc.genes <- readLines(con = "../R/seurat_resources/regev_lab_cell_cycle_genes.txt")
s.genes <- HumanGenes(object,cc.genes[1:43]) # this function can be used on zebrafish gene...
g2m.genes <- HumanGenes(object,cc.genes[44:97])
object <- CellCycleScoring(object = object, s.genes = s.genes, g2m.genes = g2m.genes, 
                        set.ident = FALSE)
RidgePlot(object = object, features.plot = HumanGenes(object,c("CCND1","CDK4","CCND2","CDK6","CCND3","RB1")), 
          nCol = 2)
object@meta.data$CC.Difference <- object@meta.data$S.Score - object@meta.data$G2M.Score
object@meta.data$S.Score = object@meta.data$S.Score - min(object@meta.data$S.Score)
object@meta.data$G2M.Score = object@meta.data$G2M.Score - min(object@meta.data$G2M.Score)
tail(x = object@meta.data)


#======1.6 PCA =========================
object %<>% ScaleData %>%
        RunPCA(pc.genes = object@var.genes, pcs.compute = 50, do.print = F)

jpeg(paste0(path,"/S1_PCElbowPlot.jpeg"), units="in", width=10, height=7,res=600)
PCElbowPlot(object, num.pc = 50)
dev.off()

jpeg(paste0(path,"/S1_PCHeatmap.jpeg"), units="in", width=10, height=7,res=600)
PCHeatmap(object, pc.use = c(1:3,28:30,38:40), cells.use = 500, do.balanced = TRUE)
dev.off()

GC()

system.time({
        object %<>% RunTSNE(reduction.use = "pca", dims.use = 1:30, do.fast = TRUE) %>%
                FindClusters(reduction.type = "pca", resolution = 0.6, dims.use = 1:30,
                             save.SNN = TRUE, n.start = 10, nn.eps = 0.5,
                             force.recalc = TRUE, print.output = FALSE)
})
p3 <- TSNEPlot(object, do.return = T, pt.size = 0.3, group.by = "orig.ident")
p4 <- TSNEPlot(object, do.label = T, do.return = T, pt.size = 0.3)

jpeg(paste0(path,"/S1_Harmony_TSNEPlot.jpeg"), units="in", width=10, height=7,res=600)
plot_grid(p3+ggtitle("group by samples")+
                  theme(plot.title = element_text(hjust = 0.5,size = 18, face = "bold")),
          p4+ggtitle("group by clusters")+
                  theme(plot.title = element_text(hjust = 0.5,size = 18, face = "bold")))
dev.off()

g <- TSNEPlot.1(object = object, do.label = F, group.by = "ident",
                        do.return = TRUE, no.legend = F, 
                        #colors.use = ExtractMetaColor(object),
                        pt.size = 1,label.size = 6 )+
        ggtitle("Tsne plot of all clusters")+
        theme(plot.title = element_text(hjust = 0.5,size = 18, face = "bold")) 

jpeg(paste0(path,"/TSNEplot~.jpeg"), units="in", width=10, height=7,res=600)
print(g)
dev.off()

save(object, file = "data/zebrafishEC_1_20190212.Rda")

# here, we generate the count tables and meta data which will be loaded by shiny:
library(here)
library(Seurat)
library(data.table)
library(tidyverse)
seurat <- readRDS(here("data", "ctrl_notch_integrated_scent.rds"))
sling <- readRDS(here("data", "seurat_slingshot_object_integrated.rds"))

EEPs <- Cells(sling[,sling$tmp_celltype == "EEP"])

seurat$celltype_manual <- ifelse(seurat$celltype_manual == "unk2","unk",seurat$celltype_manual)
seurat$celltype_manual <- ifelse(Cells(seurat) %in% EEPs, "EEP",seurat$celltype_manual)

seurat$celltype_manual <- factor(seurat$celltype_manual,
                                 levels = c("ISC", "EB", "EEP", "dEC", "daEC" ,"aEC", "mEC" ,"Copper", "LFC", "pEC", "EE", "MT", "unk"))


# remove uncessary data.
reductions <-c("pca","phate","umap")
dim_data <- data.frame(CellID = rownames(seurat@meta.data))

for (reduction in reductions) {
  dim_data <- cbind(dim_data, Embeddings(seurat, reduction)[, 1:2])
}
dim_data <- dim_data %>% select(-CellID)

seurat_small <- DietSeurat(seurat,
                           assays = "RNA", 
                           counts = F)

# remove all genes that are not expressed
seurat_small <- seurat_small[rowSums(seurat_small@assays$RNA@data) != 0,]
seurat_small@meta.data <- seurat_small[[c("celltype_manual", "perturbation")]]

all(Cells(seurat_small) == rownames(dim_data))
seurat_small@meta.data <- cbind(seurat_small@meta.data, dim_data)


# change the gene names
rownames(seurat_small@assays$RNA@data) <- 
  ifelse(rownames(seurat_small@assays$RNA@data) == "CG9650", "Cph", rownames(seurat_small@assays$RNA@data))
rownames(seurat_small@assays$RNA@counts) <- 
  ifelse(rownames(seurat_small@assays$RNA@counts) == "CG9650", "Cph", rownames(seurat_small@assays$RNA@counts))


saveRDS(seurat_small, here("data", "filtered_seurat_object.rds"), compress = T)

sce <- Seurat::as.SingleCellExperiment(seurat)
saveRDS(sce, here::here("data","filtered_sce_object.rds"))
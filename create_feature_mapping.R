# this is returned by the cell ranger output
x <- data.table::fread(here::here("data","features.tsv.gz"), header =F)[,1:2]
colnames(x) <- c("fbid", "symbol")


# available gene names for which we need mapping:
seurat <- readRDS(here::here("data", "filtered_seurat_object.rds"))
seurat_RNAi <- readRDS(here::here("data", "filtered_seurat_object_RNAi.rds"))
available_gene_names <- unique(c(rownames(seurat), rownames(seurat_RNAi)))

# gene to id mapping
mapping <- x %>% filter(symbol %in% available_gene_names)
mapping %>% write.csv(., here("data","fbid_gene_mapping.csv"), row.names = F)






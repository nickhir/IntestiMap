# wget https://ftp.flybase.net/releases/current/precomputed_files/genes/fbgn_annotation_ID.tsv.gz
x <- data.table::fread("C:/Users/Nick/Rprojects/ShinyDECODE/data/fbgn_annotation_ID.tsv", header = T)[,c(1,3,4)]
colnames(x)<- c("symbol","prim_id","second_id")

df_long <- x %>%
  separate_rows(second_id, sep = ",", convert = T) %>% 
  pivot_longer(cols = c(prim_id, second_id), names_to = "Column", values_to = "CombinedID") %>% 
  select(symbol, FBid=CombinedID) %>% 
  distinct(FBid, .keep_all = T) %>% 
  filter(FBid != "") 

df_long %>% write.csv(., file ="C:/Users/Nick/Rprojects/ShinyDECODE/data/features.csv", row.names = F)

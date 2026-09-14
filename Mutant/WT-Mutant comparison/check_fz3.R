library(Seurat)
combined <- readRDS("/home/mattia/Thesis/03Slide4/rfiles/Apc_sc_combined.rds")
DefaultAssay(combined) <- "RNA"
combined <- JoinLayers(combined)
df <- FetchData(combined, vars = c("Apc", "fz3", "Genotype", "orig.ident"), layer = "counts")
cat("\n--- fz3 by sample ---\n")
print(tapply(df$fz3, df$orig.ident, mean))
cat("\n--- Apc by sample ---\n")
print(tapply(df$Apc, df$orig.ident, mean))

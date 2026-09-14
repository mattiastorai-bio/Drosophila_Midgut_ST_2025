library(Seurat)
library(dplyr)
library(pheatmap)

# Helper function
generate_module_heatmap <- function(seurat_obj, genes_list, module_names, group_col, filename, title, img_height = 1800) {
  Idents(seurat_obj) <- group_col
  obj_sub <- subset(seurat_obj, idents = module_names)
  
  all_cells <- colnames(obj_sub)
  cells_to_keep <- unlist(lapply(module_names, function(id) {
    cluster_cells <- colnames(obj_sub)[Idents(obj_sub) == id]
    sample(cluster_cells, min(length(cluster_cells), 200))
  }))
  obj_sub <- obj_sub[, cells_to_keep]
  
  valid_data_genes <- rownames(GetAssayData(obj_sub, assay = 'SCT', layer = 'data'))
  
  gene_to_mod <- data.frame(gene = character(), module = character(), stringsAsFactors = FALSE)
  for(i in 1:length(genes_list)) {
    valid_g <- intersect(genes_list[[i]], valid_data_genes)
    if(length(valid_g) > 0) {
      gene_to_mod <- rbind(gene_to_mod, data.frame(gene = valid_g, module = module_names[i], stringsAsFactors = FALSE))
    }
  }
  
  gene_to_mod <- gene_to_mod[!duplicated(gene_to_mod$gene), ]
  all_genes <- gene_to_mod$gene
  
  mat <- GetAssayData(obj_sub, assay = 'SCT', layer = 'data')[all_genes, ]
  mat_scaled <- t(scale(t(as.matrix(mat))))
  mat_scaled[is.na(mat_scaled)] <- 0
  mat_scaled[mat_scaled > 2.5] <- 2.5
  mat_scaled[mat_scaled < -2.5] <- -2.5
  
  cluster_vals <- as.character(obj_sub@meta.data[[group_col]])
  annotation_col <- data.frame(Cluster = factor(cluster_vals, levels = module_names), row.names = colnames(obj_sub))
  annotation_row <- data.frame(Module = factor(gene_to_mod$module, levels = module_names), row.names = all_genes)
  
  mat_scaled <- mat_scaled[, order(annotation_col$Cluster)]
  annotation_col <- annotation_col[colnames(mat_scaled), , drop = FALSE]

  png(filename, width = 1200, height = img_height, res = 150)
  pheatmap(mat_scaled, 
           cluster_rows = FALSE, 
           cluster_cols = FALSE, 
           show_colnames = FALSE, 
           show_rownames = TRUE, 
           fontsize_row = 6,
           border_color = NA, 
           annotation_col = annotation_col, 
           annotation_row = annotation_row, 
           color = colorRampPalette(c('blue', 'white', 'red'))(100), 
           main = title)
  dev.off()
}

sc <- readRDS('/home/mattia/Thesis/111125ST_def/Random data/sc_lineage_labeled.rds')
markers_df <- read.csv('/home/mattia/Thesis/111125ST_def/Excels/All_Subcluster_Markers.csv')

rename_map <- c(
  'ISC-Sub-1' = 'aISC',
  'ISC-Sub-0' = 'pISC',
  'EB-Sub-1'  = 'aEB',
  'EB-Sub-4'  = 'pLFC',
  'EB-Sub-0'  = 'pEB',
  'EB-Sub-2'  = 'mEB1',
  'EB-Sub-3'  = 'mEB2',
  'dEC-Sub-0' = 'adEC',
  'dEC-Sub-1' = 'pdEC',
  'dEC-Sub-2' = 'mdEC',
  'daEC'      = 'daEC' # keep daEC
)

sc$Combined_Cluster <- gsub('_', '-', as.character(sc$Combined_Cluster))
for(old_name in names(rename_map)) {
  sc$Combined_Cluster[sc$Combined_Cluster == old_name] <- rename_map[old_name]
}

markers_df$cluster <- gsub('_', '-', as.character(markers_df$cluster))
for(old_name in names(rename_map)) {
  markers_df$cluster[markers_df$cluster == old_name] <- rename_map[old_name]
}

save_cluster_components <- function(markers_df, cluster_names, output_csv) {
  comp_df <- markers_df %>% 
    filter(cluster %in% cluster_names) %>%
    group_by(cluster) %>% 
    top_n(n = 50, wt = avg_log2FC)
  write.csv(comp_df, output_csv, row.names = FALSE)
  return(comp_df)
}

# ISC
isc_names <- c('aISC', 'pISC', 'ISC-Sub-2')
isc_comps <- save_cluster_components(markers_df, isc_names, '/home/mattia/Thesis/111125ST_def/Excels/ISC_Subcluster_Components_Renamed.csv')
isc_list <- split(isc_comps$gene, isc_comps$cluster)
isc_list <- isc_list[isc_names]
generate_module_heatmap(sc, isc_list, isc_names, 'Combined_Cluster', '/home/mattia/Thesis/111125ST_def/figures/lineage/ISC_Subcluster_Components_Heatmap_Renamed.png', 'ISC Sub-cluster Marker Components')

# EB
eb_names <- c('aEB', 'pLFC', 'pEB', 'mEB1', 'mEB2')
eb_comps <- save_cluster_components(markers_df, eb_names, '/home/mattia/Thesis/111125ST_def/Excels/EB_Subcluster_Components_Renamed.csv')
eb_list <- split(eb_comps$gene, eb_comps$cluster)
eb_list <- eb_list[eb_names]
generate_module_heatmap(sc, eb_list, eb_names, 'Combined_Cluster', '/home/mattia/Thesis/111125ST_def/figures/lineage/EB_Subcluster_Components_Heatmap_Renamed.png', 'EB Sub-cluster Marker Components')

# dEC
dec_names <- c('adEC', 'pdEC', 'mdEC')
dec_comps <- save_cluster_components(markers_df, dec_names, '/home/mattia/Thesis/111125ST_def/Excels/dEC_Subcluster_Components_Renamed.csv')
dec_list <- split(dec_comps$gene, dec_comps$cluster)
dec_list <- dec_list[dec_names]
generate_module_heatmap(sc, dec_list, dec_names, 'Combined_Cluster', '/home/mattia/Thesis/111125ST_def/figures/lineage/dEC_Subcluster_Components_Heatmap_Renamed.png', 'dEC Sub-cluster Marker Components')

# Terminals
terminal_names <- c('aEC', 'mEC-Zip-neg', 'mEC-Zip-pos', 'Copper', 'LFC', 'pEC')
terminal_comps <- save_cluster_components(markers_df, terminal_names, '/home/mattia/Thesis/111125ST_def/Excels/Terminal_Subcluster_Components_Renamed.csv')
terminal_list <- split(terminal_comps$gene, terminal_comps$cluster)
terminal_list <- terminal_list[terminal_names]
generate_module_heatmap(sc, terminal_list, terminal_names, 'Combined_Cluster', '/home/mattia/Thesis/111125ST_def/figures/lineage/Terminal_Subcluster_Components_Heatmap_Renamed.png', 'Terminal Cell Marker Components', img_height = 2400)

message('Heatmaps regenerated with new names, gene labels, and terminal cells!')

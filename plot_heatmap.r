library(dplyr)
library(pheatmap)
library(RColorBrewer)
library(ComplexHeatmap)
library(grid)

df_change <- read.csv('chemosensory_receptors_change.csv', stringsAsFactors = FALSE, check.names = FALSE)
df_counts <- read.csv('chemosensory_receptors.csv', stringsAsFactors = FALSE, check.names = FALSE)

common_ids <- intersect(df_change$FamilyID, df_counts$FamilyID)
all_change <- df_change %>% filter(FamilyID %in% common_ids) %>% arrange(description)
all_counts <- df_counts %>% filter(FamilyID %in% common_ids) %>% arrange(description)


ordered_cols <- c(
  'Cimex lectularius', 'Rhynocoris fuscipes', 'Rhodnius prolixus', 'Triatoma rubida', 'Triatoma sanguisuga'
)


color_matrix <- as.matrix(all_change[, ordered_cols])
text_matrix <- as.matrix(all_counts[, ordered_cols])


color_matrix[is.na(color_matrix)] <- 0
text_matrix[is.na(text_matrix)] <- 0


rownames(color_matrix) <- paste(all_change$FamilyID, all_change$description, sep = " | ")


vibrant_palette <- colorRampPalette(rev(brewer.pal(11, "RdBu")))(100)
col_fun = colorRamp2(seq(-15, 15, length.out = 100), vibrant_palette)


ht <- Heatmap(
  color_matrix, 
  name = "Gene Family Shifts", 
  col = col_fun,
  column_title = "Top expanded/contracted orthogroups",
  column_title_gp = gpar(fontsize = 16, fontface = "bold"),
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  row_names_side = "left",
  column_names_side = "top",
  column_names_rot = 90,
  row_names_gp = gpar(fontsize = 6, fontface = "bold"),
  column_names_gp = gpar(fontsize = 10, fontface = "bold.italic"),
  rect_gp = gpar(col = "black", lwd = 0.2),
  cell_fun = function(j, i, x, y, width, height, fill) {
    grid.text(
      text_matrix[i, j], x, y, 
      gp = gpar(fontsize = 5, col = "black")
    )
  x, y, width, height, fill) {
    grid.text(
      text_matrix[i, j], x, }
)

plot_height <- nrow(color_matrix) * 0.1 + 5

pdf("All_OGs_Vertical_Heatmap.pdf", width = 10, height = plot_height)
draw(ht, padding = unit(c(10, 20, 10, 7), "mm"))
dev.off()

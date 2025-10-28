library(tidyverse)
library(pheatmap)
library(RColorBrewer)

d1 <- read.csv("CEF_PL_double_bond_mean_heatmap.csv", header = TRUE, row.names = 1)
d2 <- read.csv("TET_PL_double_bond_mean_heatmap.csv", header = TRUE, row.names = 1)

annotation_col = data.frame(
  Group = rep(c("Ctrl", "DTAC", "SDS", "TX-100"), c(1, 1, 1, 1)))
rownames(annotation_col) = colnames(d1)

annotation_row = data.frame(
  Double_bond_number = rep(c("0", "1", "2"), c(2, 2, 2)))
rownames(annotation_row) = rownames(d1)

annotation_colors_1 <- list(
  Group = c(
    Ctrl     = "#C0C0C0", 
    DTAC     = "#bd1e2f", 
    SDS      = "#1f4e9f", 
    "TX-100" = "#7367BE"),
  Double_bond_number =c(
    "0" = "#FFFFFF",
    "1" = "#FFE5CC",
    "2" = "#FFC285"
  ))

heatmap1 <- pheatmap(d1, scale = "row",
                     display_numbers = TRUE, number_format = "%.2f", fontsize_number = 8, number_color = "black",
                     color = colorRampPalette(c("#1A5592", "white", "#B83D3D"))(100), border_color = "grey50",
                     cluster_cols = FALSE, cluster_rows = FALSE,
                     show_rownames = TRUE, show_colnames = FALSE,
                     annotation_col = annotation_col, annotation_row = annotation_row, annotation_colors = annotation_colors_1,
                     fontsize = 13.5, cellwidth = 22, cellheight = 22)

save_pheatmap_pdf <- function(x, filename, width = 6, height = 9) {
  stopifnot(!missing(x))
  stopifnot(!missing(filename))
  pdf(filename, width = width, height = height)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}
save_pheatmap_pdf(heatmap1, "CEF_PL_double_bond_mean_heatmap.pdf")



heatmap2 <- pheatmap(d2, scale = "row",
                     display_numbers = TRUE, number_format = "%.2f", fontsize_number = 8, number_color = "black",
                     color = colorRampPalette(c("#1A5592", "white", "#B83D3D"))(100), border_color = "grey50",
                     cluster_cols = FALSE, cluster_rows = FALSE,
                     show_rownames = TRUE, show_colnames = FALSE,
                     annotation_col = annotation_col, annotation_row = annotation_row, annotation_colors = annotation_colors_1,
                     fontsize = 13.5, cellwidth = 22, cellheight = 22)

save_pheatmap_pdf <- function(x, filename, width = 6, height = 9) {
  stopifnot(!missing(x))
  stopifnot(!missing(filename))
  pdf(filename, width = width, height = height)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}
save_pheatmap_pdf(heatmap2, "TET_PL_double_bond_mean_heatmap.pdf")



d3 <- read.csv("CEF_PL_acyl_carbon_mean_heatmap.csv", header = TRUE, row.names = 1)
d4 <- read.csv("TET_PL_acyl_carbon_mean_heatmap.csv", header = TRUE, row.names = 1)

annotation_col = data.frame(
  Group = rep(c("Ctrl", "DTAC", "SDS", "TX-100"), c(1, 1, 1, 1)))
rownames(annotation_col) = colnames(d3)

annotation_row = data.frame(
  Acyl_carbon_number = rep(c("28", "30", "31", "32", "33", "34", "35", "36", "38"), c(2, 2, 2, 2, 2, 2, 2, 1, 1)))
rownames(annotation_row) = rownames(d3)

annotation_colors_2 <- list(
  Group = c(
    Ctrl     = "#C0C0C0", 
    DTAC     = "#bd1e2f", 
    SDS      = "#1f4e9f", 
    "TX-100" = "#7367BE"),
  Acyl_carbon_number =c(
    "28" = "#FFFFFF",
    "30" = "#F9F6F0",
    "31" = "#EEE6D8",
    "32" = "#D8C9B0",
    "33" = "#BFAE8F",
    "34" = "#A58F6E",
    "35" = "#8A7250",
    "36" = "#6A5738",
    "38" = "#3D2A1D"
  ))

heatmap3 <- pheatmap(d3, scale = "row",
                     display_numbers = TRUE, number_format = "%.2f", fontsize_number = 8, number_color = "black",
                     color = colorRampPalette(c("#1A5592", "white", "#B83D3D"))(100), border_color = "grey50",
                     cluster_cols = FALSE, cluster_rows = FALSE,
                     show_rownames = TRUE, show_colnames = FALSE,
                     annotation_col = annotation_col, annotation_row = annotation_row, annotation_colors = annotation_colors_2,
                     fontsize = 13.5, cellwidth = 22, cellheight = 22)

save_pheatmap_pdf <- function(x, filename, width = 6, height = 9) {
  stopifnot(!missing(x))
  stopifnot(!missing(filename))
  pdf(filename, width = width, height = height)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}
save_pheatmap_pdf(heatmap3, "CEF_PL_acyl_carbon_mean_heatmap.pdf")



normalize_matrix <- function(d4, scale = "row") {
  if (scale == "row") {
    return(t(scale(t(d4))))
  } else if (scale == "column") {
    return(scale(d4))
  } else {
    return(d4)
  }
}

scale_method <- "row"
normalized_matrix <- normalize_matrix(d4, scale = scale_method)

create_normalized_labels <- function(norm_mat, digits = 1) {
  labels <- matrix("", nrow(norm_mat), ncol(norm_mat))
  
  non_na <- !is.na(norm_mat)
  
  labels[non_na] <- sapply(norm_mat[non_na], function(x) {
    if (is.na(x)) return("")
    return(sprintf("%.2f", x))
  })
  
  return(labels)
}

heatmap_labels <- create_normalized_labels(normalized_matrix)

heatmap4 <- pheatmap(d4, scale = "row",
                     display_numbers = heatmap_labels, number_format = "%.2f", fontsize_number = 8, number_color = "black",
                     color = colorRampPalette(c("#1A5592", "white", "#B83D3D"))(100), border_color = "grey50",
                     cluster_cols = FALSE, cluster_rows = FALSE,
                     show_rownames = TRUE, show_colnames = FALSE,
                     annotation_col = annotation_col, annotation_row = annotation_row, annotation_colors = annotation_colors_2,
                     fontsize = 13.5, cellwidth = 22, cellheight = 22)

save_pheatmap_pdf <- function(x, filename, width = 6, height = 9) {
  stopifnot(!missing(x))
  stopifnot(!missing(filename))
  pdf(filename, width = width, height = height)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}
save_pheatmap_pdf(heatmap4, "TET_PL_acyl_carbon_mean_heatmap.pdf")

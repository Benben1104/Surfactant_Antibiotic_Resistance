library(ggplot2)
library(FactoMineR)
library(factoextra)
library(dplyr)
library(tidyr)
library(stringr)
library(tibble)

d1 <- read.csv("CEF_PL_composition.csv")

d1_new <- as.data.frame(t(d1[-1]))
colnames(d1_new) <- d1$type
d1_new <- rownames_to_column(d1_new, var = "group")

d1_new <- d1_new %>%
  mutate(group_new = case_when(
    str_detect(group, "Ctrl\\.") ~ "Ctrl",
    str_detect(group, "DTAC\\.") ~ "DTAC",
    str_detect(group, "SDS\\.") ~ "SDS",
    str_detect(group, "TX\\.100\\.") ~ "TX-100")) %>%
  mutate(group_new = factor(group_new, levels = c("Ctrl", "DTAC", "SDS", "TX-100")))

pca_data_1 <- d1_new[, 2:6] 
pca_result_1 <- prcomp(pca_data_1, scale. = TRUE)
summary(pca_result_1)

scores_1 <- as.data.frame(pca_result_1$x[, 1:2])
scores_1$group_new <- d1_new$group_new

var_explained_1 <- pca_result_1$sdev^2 / sum(pca_result_1$sdev^2)
var_explained_1 <- round(var_explained_1 * 100, 1)
var_explained_1 <- sprintf("%.1f", var_explained_1)

ggplot(scores_1, aes(x = PC1, y = PC2, color = group_new))+
  geom_point(size = 3.25, alpha = 0.85)+
  stat_ellipse(level = 0.95, show.legend = FALSE)+
  labs(
    x = paste0("PC1 (", var_explained_1[1], "%)"),
    y = paste0("PC2 (", var_explained_1[2], "%)"),
    title = "CEF")+
  scale_x_continuous(limits = c(-8, 8))+
  scale_y_continuous(limits = c(-6, 6))+
  theme_bw()+
  theme(panel.border = element_rect(fill = NA, color = "black", size = 2, linetype = "solid"))+
  theme(axis.text.x = element_text(size = 18, face = "bold", angle = 0))+
  theme(axis.text.y = element_text(size = 18, face = "bold"))+
  theme(axis.title.x = element_text(size = 20, vjust = -0.25))+
  theme(axis.title.y = element_text(size = 20, vjust = 2.25))+
  theme(legend.text = element_text(size = 17, face = "bold"))+
  theme(legend.title = element_blank())+
  theme(legend.position = "bottom")+
  theme(title = element_text(size = 18))+
  scale_color_manual(values = c("Ctrl" = "#C0C0C0", "DTAC" = "#bd1e2f", "SDS" = "#1f4e9f", "TX-100" = "#7367BE"))

ggsave("CEF_PL_PCA.png", width = 400/90, height = 485/90, dpi = 600, units = "in")



d2 <- read.csv("TET_PL_composition.csv")

d2_new <- as.data.frame(t(d2[-1]))
colnames(d2_new) <- d2$type
d2_new <- rownames_to_column(d2_new, var = "group")

d2_new <- d2_new %>%
  mutate(group_new = case_when(
    str_detect(group, "Ctrl\\.") ~ "Ctrl",
    str_detect(group, "DTAC\\.") ~ "DTAC",
    str_detect(group, "SDS\\.") ~ "SDS",
    str_detect(group, "TX\\.100\\.") ~ "TX-100")) %>%
  mutate(group_new = factor(group_new, levels = c("Ctrl", "DTAC", "SDS", "TX-100")))

pca_data_2 <- d2_new[, 2:6] 
pca_result_2 <- prcomp(pca_data_2, scale. = TRUE)
summary(pca_result_2)

scores_2 <- as.data.frame(pca_result_2$x[, 1:2])
scores_2$group_new <- d2_new$group_new

var_explained_2 <- pca_result_2$sdev^2 / sum(pca_result_2$sdev^2)
var_explained_2 <- round(var_explained_2 * 100, 1)
var_explained_2 <- sprintf("%.1f", var_explained_1)

ggplot(scores_2, aes(x = PC1, y = PC2, color = group_new))+
  geom_point(size = 3.25, alpha = 0.85)+
  stat_ellipse(level = 0.95, show.legend = FALSE)+
  labs(
    x = paste0("PC1 (", var_explained_2[1], "%)"),
    y = paste0("PC2 (", var_explained_2[2], "%)"),
    title = "TET")+
  scale_x_continuous(limits = c(-8, 8))+
  scale_y_continuous(limits = c(-6, 6))+
  theme_bw()+
  theme(panel.border = element_rect(fill = NA, color = "black", size = 2, linetype = "solid"))+
  theme(axis.text.x = element_text(size = 18, face = "bold", angle = 0))+
  theme(axis.text.y = element_text(size = 18, face = "bold"))+
  theme(axis.title.x = element_text(size = 20, vjust = -0.25))+
  theme(axis.title.y = element_text(size = 20, vjust = 2.25))+
  theme(legend.text = element_text(size = 17, face = "bold"))+
  theme(legend.title = element_blank())+
  theme(legend.position = "bottom")+
  theme(title = element_text(size = 18))+
  scale_color_manual(values = c("Ctrl" = "#C0C0C0", "DTAC" = "#bd1e2f", "SDS" = "#1f4e9f", "TX-100" = "#7367BE"))

ggsave("TET_PL_PCA.png", width = 400/90, height = 485/90, dpi = 600, units = "in")

library(ggplot2)
library(tidyverse)
library(scales)
library(RColorBrewer)

d1 <- read.csv("CEF_KEGG_enrichment.csv", header = TRUE)
head(d1)

d1_top20 <- d1 %>%
  group_by(description) %>%
  summarise(min_p_adj = min(p_adj)) %>%
  arrange(min_p_adj) %>%
  slice_head(n = 20) %>%
  pull(description)

d1_filter <- d1 %>%
  filter(description %in% d1_top20) %>%
  group_by(description) %>%
  mutate(min_p_adj = min(p_adj)) %>%
  ungroup() %>%
  mutate(description = fct_reorder(description, min_p_adj))

color1 <- colorRampPalette(rev(brewer.pal(8, "RdBu")))

ggplot(d1_filter, aes(x = group, 
                      y = reorder(description, min_p_adj),
                      size = gene_number,
                      color = p_adj))+
  geom_point(alpha = 0.8)+
  scale_size_area(max_size = 12.5, labels = number_format(accuracy = 1))+
  labs(size = "Gene number",
       color = expression(bold(paste(bolditalic(p), "-adj", sep = ""))),
       x = NULL,
       y = "Pathway")+
  theme(panel.background = element_rect(fill = NA, color = "black", size = 2),
        strip.background = element_rect(fill = NA , color = NA),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.y = element_text(size = 24, vjust = 2),
        axis.text.x = element_text(vjust = 1, angle = 45, color = "black", face = "bold", size = 18, hjust = 1),
        axis.text.y = element_text(color = "black", face = "bold", size = 18, hjust = 1),
        legend.position = "right",
        legend.title = element_text(size = 18, face = "bold"),
        legend.text = element_text(size = 14, face = "bold"),
        legend.key = element_blank())+
  scale_colour_gradientn(colours = color1(8000), labels = function(x) sprintf("%.2f", x))

ggsave("CEF_KEGG_enrichment.png", width = 910/90, height = 900/90, dpi = 600, units = "in")


d2 <- read.csv("TET_KEGG_enrichment.csv", header = TRUE)
head(d2)

d2_top20 <- d2 %>%
  group_by(description) %>%
  summarise(min_p_adj = min(p_adj)) %>%
  arrange(min_p_adj) %>%
  slice_head(n = 20) %>%
  pull(description)

d2_filter <- d2 %>%
  filter(description %in% d2_top20) %>%
  group_by(description) %>%
  mutate(min_p_adj = min(p_adj)) %>%
  ungroup() %>%
  mutate(description = fct_reorder(description, min_p_adj))

ggplot(d2_filter, aes(x = group, 
                      y = reorder(description, min_p_adj),
                      size = gene_number,
                      color = p_adj))+
  geom_point(alpha = 0.8)+
  scale_size_area(max_size = 12.5, labels = number_format(accuracy = 1))+
  labs(size = "Gene number",
       color = expression(bold(paste(bolditalic(p), "-adj", sep = ""))),
       x = NULL,
       y = "Pathway")+
  theme(panel.background = element_rect(fill = NA, color = "black", size = 2),
        strip.background = element_rect(fill = NA , color = NA),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.y = element_text(size = 24, vjust = 2),
        axis.text.x = element_text(vjust = 1, angle = 45, color = "black", face = "bold", size = 18, hjust = 1),
        axis.text.y = element_text(color = "black", face = "bold", size = 18, hjust = 1),
        legend.position = "right",
        legend.title = element_text(size = 18, face = "bold"),
        legend.text = element_text(size = 14, face = "bold"),
        legend.key = element_blank())+
  scale_colour_gradientn(colours = color1(8000), labels = function(x) sprintf("%.2f", x))

ggsave("TET_KEGG_enrichment.png", width = 893/90, height = 900/90, dpi = 600, units = "in")

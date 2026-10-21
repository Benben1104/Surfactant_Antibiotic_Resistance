library(ggVennDiagram)
library(ggplot2)
library(dplyr)

d_gather1 <- read.csv("CEF_gather.csv")
head(d_gather1)
unique(d_gather1$group)

set1 <- d_gather1 %>%
  filter(group == c("DTAC_vs._Ctrl")) %>%
  distinct(gene) %>%
  pull(gene)
set2 <- d_gather1 %>%
  filter(group == c("SDS_vs._Ctrl")) %>%
  distinct(gene) %>%
  pull(gene)
set3 <- d_gather1 %>%
  filter(group == c("TX-100_vs._Ctrl")) %>%
  distinct(gene) %>%
  pull(gene)

ggVennDiagram(list(set1, set2, set3),
              label_alpha = 0, label_size = 5)+
  theme(legend.position = "none")+
  scale_fill_distiller(palette = "RdBu", direction = -1)+
  scale_color_manual(values = c("black", "black", "black"))

ggsave("CEF_venn.png", width = 520/90, height = 350/90, dpi = 600, unit = "in")


d_gather2 <- read.csv("TET_gather.csv")
head(d_gather2)
unique(d_gather2$group)

set4 <- d_gather2 %>%
  filter(group == c("DTAC_vs._Ctrl")) %>%
  distinct(gene) %>%
  pull(gene)
set5 <- d_gather2 %>%
  filter(group == c("SDS_vs._Ctrl")) %>%
  distinct(gene) %>%
  pull(gene)
set6 <- d_gather2 %>%
  filter(group == c("TX-100_vs._Ctrl")) %>%
  distinct(gene) %>%
  pull(gene)

ggVennDiagram(list(set4, set5, set6),
              label_alpha = 0, label_size = 5)+
  theme(legend.position = "none")+
  scale_fill_distiller(palette = "RdBu", direction = -1)+
  scale_color_manual(values = c("black", "black", "black"))

ggsave("TET_venn.png", width = 520/90, height = 350/90, dpi = 600, unit = "in")

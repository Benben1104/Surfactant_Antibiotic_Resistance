library(ggplot2)
library(sf)
library(tidyverse)
library(RColorBrewer)

d <- read.csv("Gene_knockout_comparision.csv")

d1 <- subset(d, type == "WT")
d2 <- subset(d, type == "acrB")

shapiro.test(d1$FC)
bartlett.test(FC ~ group, data = d1)
t.test(FC ~ group, data = d1, var.equal = TRUE)

ggplot(d1, aes(group, FC, fill = group))+
  theme_bw()+
  theme(panel.border = element_rect(fill = NA, color = "black", size = 2, linetype = "solid"))+
  geom_violin(draw_quantiles = 0.1, aes(fill = group), width = 0.9, size = 1.5, alpha = 0.75)+
  geom_boxplot(width = 0.25, height = 0.5, size = 1, alpha = 0.7, outlier.shape = NA)+
  labs(x = NULL,
       y = "Relative TET accumulation",
       title = "WT")+
  coord_cartesian(ylim = c(0, 2.5))+
  scale_fill_manual(values = c("Ctrl" = "#C0C0C0", "DTAC" = "#bd1e2f"))+
  scale_color_manual(values = c("Ctrl" = "#C0C0C0", "DTAC" = "#bd1e2f"))+
  theme(axis.text.x = element_text(size = 21, face = "bold", angle = 0))+
  theme(axis.text.y = element_text(size = 21, face = "bold"))+
  theme(axis.title.y = element_text(size = 24, vjust = 2))+
  theme(legend.position = "none")+
  theme(title = element_text(size = 21))

ggsave("WT.png", width = 600/90, height = 440/90, dpi = 600, unit = "in")



shapiro.test(d2$FC)
bartlett.test(FC ~ group, data = d2)
t.test(FC ~ group, data = d2, var.equal = TRUE)

ggplot(d2, aes(group, FC, fill = group))+
  theme_bw()+
  theme(panel.border = element_rect(fill = NA, color = "black", size = 2, linetype = "solid"))+
  geom_violin(draw_quantiles = 0.1, aes(fill = group), width = 0.9, size = 1.5, alpha = 0.75)+
  geom_boxplot(width = 0.25, height = 0.5, size = 1, alpha = 0.7, outlier.shape = NA)+
  labs(x = NULL,
       y = "Relative TET accumulation",
       title = "ΔacrB")+
  coord_cartesian(ylim = c(0, 1.5))+
  scale_fill_manual(values = c("Ctrl" = "#C0C0C0", "DTAC" = "#bd1e2f"))+
  scale_color_manual(values = c("Ctrl" = "#C0C0C0", "DTAC" = "#bd1e2f"))+
  theme(axis.text.x = element_text(size = 21, face = "bold", angle = 0))+
  theme(axis.text.y = element_text(size = 21, face = "bold"))+
  theme(axis.title.y = element_text(size = 24, vjust = 2))+
  theme(legend.position = "none")+
  theme(title = element_text(size = 21))

ggsave("ΔacrB.png", width = 600/90, height = 440/90, dpi = 600, unit = "in")

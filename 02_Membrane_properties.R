library(ggplot2)
library(ggalt)
library(RColorBrewer)

d_CEF_1 <-  read.csv("CEF_potential.csv")
head(d_CEF_1)

ggplot(d_CEF_1, aes(x = time, y = mean, color = group))+
  theme_bw()+
  scale_y_continuous(limits = c(0.8, 1.3))+
  geom_point(size = 4, alpha = 0.8)+
  geom_errorbar(aes(ymin = mean-sd, ymax = mean+sd), width = 0)+
  geom_xspline(spline_shape = 1.0)+
  theme(strip.text.x = element_text(size = 16))+
  labs(x = expression("Time (min)"),
       y = expression("Fold change of potential"),
       title = expression("CEF"),
       color = NULL)+
  scale_color_manual(values = c("Ctrl" = "#C0C0C0", "DTAC" = "#bd1e2f", "SDS" = "#1f4e9f", "TX-100" = "#7367BE"),
                     breaks = c("Ctrl", "DTAC", "SDS", "TX-100"))+
  theme(axis.text.x = element_text(size = 22, face = "bold", vjust = 0))+
  theme(axis.text.y = element_text(size = 22, face = "bold"))+
  theme(axis.title.x = element_text(size = 24, face = "bold", vjust = -0.5))+
  theme(axis.title.y = element_text(size = 24, face = "bold", vjust = 2))+
  theme(title = element_text(size = 22, face = "bold"))+
  theme(legend.text = element_text(size = 21, face = "bold"))+
  theme(legend.position = "bottom")+
  theme(panel.border = element_rect(fill = NA, color = "black", size = 2, linetype = "solid"))

ggsave("CEF_potential.png", width = 600/90, height = 500/90, dpi = 600, unit = "in")



d_TET_1 <- read.csv("TET_potential.csv")
head(d_TET_1)

ggplot(d_TET_1, aes(x = time, y = mean, color = group))+
  theme_bw()+
  scale_y_continuous(limits = c(0.6, 1.6))+
  geom_point(size = 4, alpha = 0.8)+
  geom_errorbar(aes(ymin = mean-sd, ymax = mean+sd), width = 0)+
  geom_xspline(spline_shape = 1.0)+
  theme(strip.text.x = element_text(size = 16))+
  labs(x = expression("Time (min)"),
       y = expression("Fold change of potential"),
       title = expression("TET"),
       color = NULL)+
  scale_color_manual(values = c("Ctrl" = "#C0C0C0", "DTAC" = "#bd1e2f", "SDS" = "#1f4e9f", "TX-100" = "#7367BE"),
                     breaks = c("Ctrl", "DTAC", "SDS", "TX-100"))+
  theme(axis.text.x = element_text(size = 22, face = "bold", vjust = 0))+
  theme(axis.text.y = element_text(size = 22, face = "bold"))+
  theme(axis.title.x = element_text(size = 24, face = "bold", vjust = -0.5))+
  theme(axis.title.y = element_text(size = 24, face = "bold", vjust = 2))+
  theme(title = element_text(size = 22, face = "bold"))+
  theme(legend.text = element_text(size = 21, face = "bold"))+
  theme(legend.position = "bottom")+
  theme(panel.border = element_rect(fill = NA, color = "black", size = 2, linetype = "solid"))

ggsave("TET_potential.png", width = 600/90, height = 500/90, dpi = 600, unit = "in")



library(sf)
library(tidyverse)
library(ggpubr)
library(scales)

d_CEF_2 <- read.csv("CEF_membrane.csv")
head(d_CEF_2)

residuals <- (aov(membrane.permeability ~ group, data = d_CEF_2))$residuals
shapiro.test(residuals)
bartlett.test(membrane.permeability ~ group, data = d_CEF_2)
aov.1 <- aov(membrane.permeability ~ group, data = d_CEF_2)
summary(aov.1)
TukeyHSD(aov.1)

ggplot(d_CEF_2, aes(group, membrane.permeability, fill = group))+
  theme_bw()+
  scale_y_continuous(limits = c(0, 1.5))+
  geom_jitter(aes(group, membrane.permeability), width = 0.22, 
              size = 4.5, shape = 21, stroke = 0.01, alpha = 1)+
  stat_summary(geom = "col", fun = "mean",
               position = "dodge", width = 0.625, alpha = 0.7)+
  stat_summary(geom = "errorbar",
               fun.min = function(x) mean(x)-sd(x),
               fun.max = function(x) mean(x)+sd(x),
               width = 0.15)+
  labs(x = NULL,
       y = expression("Fold change of permeability"),
       title = expression("CEF"))+
  scale_fill_manual(values = c("Ctrl" = "#C0C0C0", "DTAC" = "#bd1e2f", "SDS" = "#1f4e9f", "TX-100" = "#7367BE"))+
  theme(axis.text.x = element_text(size = 18, face = "bold", vjust = 0))+
  theme(axis.text.y = element_text(size = 18, face = "bold"))+
  theme(axis.title.y = element_text(size = 20, face = "bold", vjust = 2))+
  theme(legend.position = "none")+
  theme(title = element_text(size = 18, face = "bold"))+
  theme(panel.border = element_rect(fill = NA, color = "black", size = 2, linetype = "solid"))

ggsave("CEF_permeability.png", width = 520/90, height = 357/90, dpi = 600, unit = "in")



residuals <- (aov(membrane.hydrophobicity ~ group, data = d_CEF_2))$residuals
shapiro.test(residuals)
bartlett.test(membrane.hydrophobicity ~ group, data = d_CEF_2)
aov.2 <- aov(membrane.hydrophobicity ~ group, data = d_CEF_2)
summary(aov.2)
TukeyHSD(aov.2)

ggplot(d_CEF_2, aes(group, membrane.hydrophobicity, fill = group))+
  theme_bw()+
  scale_y_continuous(limits = c(0, 1.25))+
  geom_jitter(aes(group, membrane.hydrophobicity), width = 0.22, 
              size = 4.5, shape = 21, stroke = 0.01, alpha = 1)+
  stat_summary(geom = "col", fun = "mean",
               position = "dodge", width = 0.625, alpha = 0.7)+
  stat_summary(geom = "errorbar",
               fun.min = function(x) mean(x)-sd(x),
               fun.max = function(x) mean(x)+sd(x),
               width = 0.15)+
  labs(x = NULL,
       y = expression("Fold change of hydrophobicity"),
       title = expression("CEF"))+
  scale_fill_manual(values = c("Ctrl" = "#C0C0C0", "DTAC" = "#bd1e2f", "SDS" = "#1f4e9f", "TX-100" = "#7367BE"))+
  theme(axis.text.x = element_text(size = 18, face = "bold", vjust = 0))+
  theme(axis.text.y = element_text(size = 18, face = "bold"))+
  theme(axis.title.y = element_text(size = 20, face = "bold", vjust = 2))+
  theme(legend.position = "none")+
  theme(title = element_text(size = 18, face = "bold"))+
  theme(panel.border = element_rect(fill = NA, color = "black", size = 2, linetype = "solid"))

ggsave("CEF_hydrophobicity.png", width = 520/90, height = 357/90, dpi = 600, unit = "in")



d_TET_2 <- read.csv("TET_membrane.csv")
head(d_TET_2)

residuals <- (aov(membrane.permeability ~ group, data = d_TET_2))$residuals
shapiro.test(residuals)
bartlett.test(membrane.permeability ~ group, data = d_TET_2)
aov.3 <- aov(membrane.permeability ~ group, data = d_TET_2)
summary(aov.3)
TukeyHSD(aov.3)

ggplot(d_TET_2, aes(group, membrane.permeability, fill = group))+
  theme_bw()+
  scale_y_continuous(limits = c(0, 3), labels = number_format(accuracy = 0.1))+
  geom_jitter(aes(group, membrane.permeability), width =0.22, 
              size = 4.5, shape = 21, stroke = 0.01, alpha = 1)+
  stat_summary(geom = "col", fun = "mean",
               position = "dodge", width = 0.625, alpha = 0.7)+
  stat_summary(geom = "errorbar",
               fun.min = function(x) mean(x)-sd(x),
               fun.max = function(x) mean(x)+sd(x),
               width = 0.15)+
  labs(x = NULL,
       y = expression("Fold change of permeability"),
       title = expression("TET"))+
  scale_fill_manual(values = c("Ctrl" = "#C0C0C0", "DTAC" = "#bd1e2f", "SDS" = "#1f4e9f", "TX-100" = "#7367BE"))+
  theme(axis.text.x = element_text(size = 18, face = "bold", vjust = 0))+
  theme(axis.text.y = element_text(size = 18, face = "bold"))+
  theme(axis.title.y = element_text(size = 20, face = "bold", vjust = 2))+
  theme(legend.position = "none")+
  theme(title = element_text(size = 18, face = "bold"))+
  theme(panel.border = element_rect(fill = NA, color = "black", size = 2, linetype = "solid"))

ggsave("TET_permeability.png", width = 520/90, height = 357/90, dpi = 600, unit = "in")



residuals <- (aov(membrane.hydrophobicity ~ group, data = d_TET_2))$residuals
shapiro.test(residuals)
bartlett.test(membrane.hydrophobicity ~ group, data = d_TET_2)
aov.4 <- aov(membrane.hydrophobicity ~ group, data = d_TET_2)
summary(aov.4)
TukeyHSD(aov.4)

ggplot(d_TET_2, aes(group, membrane.hydrophobicity, fill = group))+
  theme_bw()+
  scale_y_continuous(limits = c(0, 1.25))+
  geom_jitter(aes(group, membrane.hydrophobicity), width = 0.22, 
              size = 4.5, shape = 21, stroke = 0.01, alpha = 1)+
  stat_summary(geom = "col", fun = "mean",
               position = "dodge", width = 0.625, alpha = 0.7)+
  stat_summary(geom = "errorbar",
               fun.min = function(x) mean(x)-sd(x),
               fun.max = function(x) mean(x)+sd(x),
               width = 0.15)+
  labs(x = NULL,
       y = expression("Fold change of hydrophobicity"),
       title = expression("TET"))+
  scale_fill_manual(values = c("Ctrl" = "#C0C0C0", "DTAC" = "#bd1e2f", "SDS" = "#1f4e9f", "TX-100" = "#7367BE"))+
  theme(axis.text.x = element_text(size = 18, face = "bold", vjust = 0))+
  theme(axis.text.y = element_text(size = 18, face = "bold"))+
  theme(axis.title.y = element_text(size = 20, face = "bold", vjust = 2))+
  theme(legend.position = "none")+
  theme(title = element_text(size = 18, face = "bold"))+
  theme(panel.border = element_rect(fill = NA, color = "black", size = 2, linetype = "solid"))

ggsave("TET_hydrophobicity.png", width = 520/90, height = 357/90, dpi = 600, unit = "in")

library(multcompView)
library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(ggpubr)

d_CEF_DTAC <- read.csv("CEF+DTAC.csv")
head(d_CEF_DTAC)
d_CEF_DTAC$group <- factor(d_CEF_DTAC$group, level = c("Ctrl", "0.02", "0.04", "0.2", "1", "5", "10"))

residuals <- (aov(CEF ~ group, data = d_CEF_DTAC))$residuals
shapiro.test(residuals)
bartlett.test(CEF ~ group, data = d_CEF_DTAC)
aov.1 <- aov(CEF ~ group, data = d_CEF_DTAC)
summary(aov.1)
TukeyHSD(aov.1)

ggplot(d_CEF_DTAC, aes(group, CEF, color = group))+
  geom_point(stat = "summary", fun = "mean", alpha = 1.0, size = 9)+
  geom_errorbar(stat = "summary",
                fun.min = function(x) mean(x) - sd(x),
                fun.max = function(x) mean(x) + sd(x),
                width = 0, size = 3.25)+
  theme_bw()+
  scale_y_continuous(limits = c(1.5*10^-10, 6*10^-10))+
  labs(x = NULL,
       y = expression("CEF accumulation (ng/cell)"),
       title = expression("DTAC"))+
  theme(axis.text.x = element_text(size = 18, face = "bold", vjust = 0.4, angle = 0))+
  theme(axis.text.y = element_text(size = 18, face = "bold"))+
  theme(axis.title.y = element_text(size = 24, face = "bold", vjust = 2))+
  scale_fill_brewer(palette = "Reds")+
  scale_color_brewer(palette = "Reds")+
  theme(legend.position = "none")+
  theme(title = element_text(size = 20, face = "bold"))+
  theme(panel.border = element_rect(fill = NA, color = "black", size = 2, linetype = "solid"))

ggsave("CEF+DTAC.png", width = 640/90, height = 430/90, dpi = 600, unit = "in")



d_CEF_SDS <- read.csv("CEF+SDS.csv")
head(d_CEF_SDS)
d_CEF_SDS$group <- factor(d_CEF_SDS$group, level = c("Ctrl", "0.02", "0.04", "0.2", "1", "5", "10"))

residuals <- (aov(CEF ~ group, data = d_CEF_SDS))$residuals
shapiro.test(residuals)
bartlett.test(CEF ~ group, data = d_CEF_SDS)
aov.2 <- aov(CEF ~ group, data = d_CEF_SDS)
summary(aov.2)
TukeyHSD(aov.2)

ggplot(d_CEF_SDS, aes(group, CEF, color = group))+
  geom_point(stat = "summary", fun = "mean", alpha = 1.0, size = 9)+
  geom_errorbar(stat = "summary",
                fun.min = function(x) mean(x) - sd(x),
                fun.max = function(x) mean(x) + sd(x),
                width = 0, size = 3.25)+
  theme_bw()+
  scale_y_continuous(limits = c(5*10^-11, 5*10^-10))+
  labs(x = NULL,
       y = expression("CEF accumulation (ng/cell)"),
       title = expression("SDS"))+
  theme(axis.text.x = element_text(size = 18, face = "bold", vjust = 0.4, angle = 0))+
  theme(axis.text.y = element_text(size = 18, face = "bold"))+
  theme(axis.title.y = element_text(size = 24, face = "bold", vjust = 2))+
  scale_fill_brewer(palette = "Blues")+
  scale_color_brewer(palette = "Blues")+
  theme(legend.position = "none")+
  theme(title = element_text(size = 20, face = "bold"))+
  theme(panel.border = element_rect(fill = NA, color = "black", size = 2, linetype = "solid"))

ggsave("CEF+SDS.png", width = 640/90, height = 430/90, dpi = 600, unit = "in")



d_CEF_TX <- read.csv("CEF+TX-100.csv")
head(d_CEF_TX)
d_CEF_TX$group <- factor(d_CEF_TX$group, level = c("Ctrl", "0.02", "0.04", "0.2", "1", "5", "10"))

residuals <- (aov(CEF ~ group, data = d_CEF_TX))$residuals
shapiro.test(residuals)
bartlett.test(CEF ~ group, data = d_CEF_TX)
aov.3 <- aov(CEF ~ group, data = d_CEF_TX)
summary(aov.3)
TukeyHSD(aov.3)

ggplot(d_CEF_TX, aes(group, CEF, color = group))+
  geom_point(stat = "summary", fun = "mean", alpha = 1.0, size = 9)+
  geom_errorbar(stat = "summary",
                fun.min = function(x) mean(x) - sd(x),
                fun.max = function(x) mean(x) + sd(x),
                width = 0, size = 3.25)+
  theme_bw()+
  scale_y_continuous(limits = c(5*10^-11, 5*10^-10))+
  labs(x = NULL,
       y = expression("CEF accumulation (ng/cell)"),
       title = expression("TX-100"))+
  theme(axis.text.x = element_text(size = 18, face = "bold", vjust = 0.4, angle = 0))+
  theme(axis.text.y = element_text(size = 18, face = "bold"))+
  theme(axis.title.y = element_text(size = 24, face = "bold", vjust = 2))+
  scale_fill_brewer(palette = "Purples")+
  scale_color_brewer(palette = "Purples")+
  theme(legend.position = "none")+
  theme(title = element_text(size = 20, face = "bold"))+
  theme(panel.border = element_rect(fill = NA, color = "black", size = 2, linetype = "solid"))

ggsave("CEF+TX-100.png", width = 640/90, height = 430/90, dpi = 600, unit = "in")



d_TET_DTAC <- read.csv("TET+DTAC.csv")
head(d_TET_DTAC)
d_TET_DTAC$group <- factor(d_TET_DTAC$group, level = c("Ctrl", "0.02", "0.04", "0.2", "1", "5", "10"))

residuals <- (aov(TET ~ group, data = d_TET_DTAC))$residuals
shapiro.test(residuals)
bartlett.test(TET ~ group, data = d_TET_DTAC)
aov.4 <- aov(TET ~ group, data = d_TET_DTAC)
summary(aov.4)
TukeyHSD(aov.4)

ggplot(d_TET_DTAC, aes(group, TET, color = group))+
  geom_point(stat = "summary", fun = "mean", alpha = 1.0, size = 9)+
  geom_errorbar(stat = "summary",
                fun.min = function(x) mean(x) - sd(x),
                fun.max = function(x) mean(x) + sd(x),
                width = 0, size = 3.25)+
  theme_bw()+
  scale_y_continuous(limits = c(1*10^-8, 7*10^-8))+
  labs(x = NULL,
       y = expression("TET accumualtion (ng/cell)"),
       title = expression("DTAC"))+
  theme(axis.text.x = element_text(size = 18, face = "bold", vjust = 0.4, angle = 0))+
  theme(axis.text.y = element_text(size = 18, face = "bold"))+
  theme(axis.title.y = element_text(size = 24, face = "bold", vjust = 2))+
  scale_fill_brewer(palette = "Reds")+
  scale_color_brewer(palette = "Reds")+
  theme(legend.position = "none")+
  theme(title = element_text(size = 20, face = "bold"))+
  theme(panel.border = element_rect(fill = NA, color = "black", size = 2, linetype = "solid"))

ggsave("TET+DTAC.png", width = 640/90, height = 430/90, dpi = 600, unit = "in")



d_TET_SDS <- read.csv("TET+SDS.csv")
head(d_TET_SDS)
d_TET_SDS$group <- factor(d_TET_SDS$group, level = c("Ctrl", "0.02", "0.04", "0.2", "1", "5", "10"))

residuals <- (aov(TET ~ group, data = d_TET_SDS))$residuals
shapiro.test(residuals)
bartlett.test(TET ~ group, data = d_TET_SDS)
aov.5 <- aov(TET ~ group, data = d_TET_SDS)
summary(aov.5)
TukeyHSD(aov.5)

ggplot(d_TET_SDS, aes(group, TET, color = group))+
  geom_point(stat = "summary", fun = "mean", alpha = 1.0, size = 9)+
  geom_errorbar(stat = "summary",
                fun.min = function(x) mean(x) - sd(x),
                fun.max = function(x) mean(x) + sd(x),
                width = 0, size = 3.25)+
  theme_bw()+
  scale_y_continuous(limits = c(3.0*10^-8, 9*10^-8))+
  labs(x = NULL,
       y = expression("TET accumulation (ng/cell)"),
       title = expression("SDS"))+
  theme(axis.text.x = element_text(size = 18, face = "bold", vjust = 0.4, angle = 0))+
  theme(axis.text.y = element_text(size = 18, face = "bold"))+
  theme(axis.title.y = element_text(size = 24, face = "bold", vjust = 2))+
  scale_fill_brewer(palette = "Blues")+
  scale_color_brewer(palette = "Blues")+
  theme(legend.position = "none")+
  theme(title = element_text(size = 20, face = "bold"))+
  theme(panel.border = element_rect(fill = NA, color = "black", size = 2, linetype = "solid"))

ggsave("TET+SDS.png", width = 640/90, height = 430/90, dpi = 600, unit = "in")



d_TET_TX <- read.csv("TET+TX-100.csv")
head(d_TET_TX)
d_TET_TX$group <- factor(d_TET_TX$group, level = c("Ctrl", "0.02", "0.04", "0.2", "1", "5", "10"))

residuals <- (aov(TET ~ group, data = d_TET_TX))$residuals
shapiro.test(residuals)
bartlett.test(TET ~ group, data = d_TET_TX)
aov.6 <- aov(TET ~ group, data = d_TET_TX)
summary(aov.6)
TukeyHSD(aov.6)

ggplot(d_TET_TX, aes(group, TET, color = group))+
  geom_point(stat = "summary", fun = "mean", alpha = 1.0, size = 9)+
  geom_errorbar(stat = "summary",
                fun.min = function(x) mean(x) - sd(x),
                fun.max = function(x) mean(x) + sd(x),
                width = 0, size = 3.25)+
  theme_bw()+
  scale_y_continuous(limits = c(1*10^-8, 5*10^-8))+
  labs(x = NULL,
       y = expression("TET accumulation (ng/cell)"),
       title = expression("TX-100"))+
  theme(axis.text.x = element_text(size = 18, face = "bold", vjust = 0.4, angle = 0))+
  theme(axis.text.y = element_text(size = 18, face = "bold"))+
  theme(axis.title.y = element_text(size = 24, face = "bold", vjust = 2))+
  scale_fill_brewer(palette = "Purples")+
  scale_color_brewer(palette = "Purples")+
  theme(legend.position = "none")+
  theme(title = element_text(size = 20, face = "bold"))+
  theme(panel.border = element_rect(fill = NA, color = "black", size = 2, linetype = "solid"))

ggsave("TET+TX-100.png", width = 640/90, height = 430/90, dpi = 600, unit = "in")

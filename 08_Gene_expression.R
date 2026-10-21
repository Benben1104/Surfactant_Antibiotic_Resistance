library(ggplot2)
library(dplyr)
library(tidyr)
library(car)
library(RColorBrewer)
library(scales)

d1 <- read.csv("TET_gene_expression_1.csv")
head(d1)
d1$gene <- factor(d1$gene, level = c("marA", "acrA", "acrB", "ompF"))

d_marA_TET <- subset(d1, gene == "marA")
residuals <- (aov(log2FC ~ group, data = d_marA_TET))$residuals
shapiro.test(residuals)
leveneTest(log2FC ~ group, data = d_marA_TET)
aov.1 <- aov(log2FC ~ group, data = d_marA_TET)
summary(aov.1)
TukeyHSD(aov.1)

d_acrA_TET <- subset(d1, gene == "acrA")
residuals <- (aov(log2FC ~ group, data = d_acrA_TET))$residuals
shapiro.test(residuals)
leveneTest(log2FC ~ group, data = d_acrA_TET)
aov.2 <- aov(log2FC ~ group, data = d_acrA_TET)
summary(aov.2)
TukeyHSD(aov.2)

d_acrB_TET <- subset(d1, gene == "acrB")
residuals <- (aov(log2FC ~ group, data = d_acrB_TET))$residuals
shapiro.test(residuals)
leveneTest(log2FC ~ group, data = d_acrB_TET)
aov.3 <- aov(log2FC ~ group, data = d_acrB_TET)
summary(aov.3)
TukeyHSD(aov.3)

d_ompF_TET <- subset(d1, gene == "ompF")
residuals <- (aov(log2FC ~ group, data = d_ompF_TET))$residuals
shapiro.test(residuals)
leveneTest(log2FC ~ group, data = d_ompF_TET)
aov.4 <- aov(log2FC ~ group, data = d_ompF_TET)
summary(aov.4)
TukeyHSD(aov.4)

ggplot(d1, aes(gene, log2FC, group = group, color = group))+
  theme_bw()+
  theme(panel.border = element_rect(fill = NA, color = "black", size = 2, linetype = "solid"))+
  geom_point(stat = "summary", fun = "mean", alpha = 1.0, size = 6,
             position = position_dodge(width = 0.75))+
  geom_errorbar(stat = "summary",
                fun.min = function(x) mean(x) - sd(x),
                fun.max = function(x) mean(x) + sd(x),
                width = 0, size = 2.25,
                position = position_dodge(width = 0.75))+
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", size = 1)+
  labs(x = NULL,
       y = expression(Log[2]("FC")),
       title = expression("TET"))+
    scale_y_continuous(limits = c(-2, 6), labels = label_number(accuracy = 0.1))+
    scale_fill_manual(values = c("Ctrl" = "#C0C0C0", "DTAC" = "#bd1e2f", "SDS" = "#1f4e9f", "TX-100" = "#7367BE"))+
    scale_color_manual(values = c("Ctrl" = "#C0C0C0", "DTAC" = "#bd1e2f", "SDS" = "#1f4e9f", "TX-100" = "#7367BE"))+
    theme(axis.text.x = element_text(size = 21, face = "bold.italic", angle = 0))+
    theme(axis.text.y = element_text(size = 21, face = "bold"))+
    theme(axis.title.y = element_text(size = 24, vjust = 2))+
    theme(legend.text = element_text(size = 21, face = "bold"))+
    theme(legend.title = element_blank())+
    theme(legend.position = "bottom")+
    theme(title = element_text(size = 20))

ggsave("TET_gene_expression_1.png", width = 900/90, height = 480/90, dpi = 600, units = "in")



d2 <- read.csv("TET_gene_expression_2.csv")
head(d2)
d2$gene <- factor(d2$gene, level = c("cdsA", "pgsA", "pgpA", "pssA", "psd", "pldA"))

d_cdsA_TET <- subset(d2, gene == "cdsA")
residuals <- (aov(log2FC ~ group, data = d_cdsA_TET))$residuals
shapiro.test(residuals)
leveneTest(log2FC ~ group, data = d_cdsA_TET)
aov.5 <- aov(log2FC ~ group, data = d_cdsA_TET)
summary(aov.5)
TukeyHSD(aov.5)

d_pgsA_TET <- subset(d2, gene == "pgsA")
residuals <- (aov(log2FC ~ group, data = d_pgsA_TET))$residuals
shapiro.test(residuals)
leveneTest(log2FC ~ group, data = d_pgsA_TET)
aov.6 <- aov(log2FC ~ group, data = d_pgsA_TET)
summary(aov.6)
TukeyHSD(aov.6)

d_pgpA_TET <- subset(d2, gene == "pgpA")
residuals <- (aov(log2FC ~ group, data = d_pgpA_TET))$residuals
shapiro.test(residuals)
leveneTest(log2FC ~ group, data = d_pgpA_TET)
aov.7 <- aov(log2FC ~ group, data = d_pgpA_TET)
summary(aov.7)
TukeyHSD(aov.7)

d_pssA_TET <- subset(d2, gene == "pssA")
residuals <- (aov(log2FC ~ group, data = d_pssA_TET))$residuals
shapiro.test(residuals)
leveneTest(log2FC ~ group, data = d_pssA_TET)
aov.8 <- aov(log2FC ~ group, data = d_pssA_TET)
summary(aov.8)
TukeyHSD(aov.8)

d_psd_TET <- subset(d2, gene == "psd")
residuals <- (aov(log2FC ~ group, data = d_psd_TET))$residuals
shapiro.test(residuals)
leveneTest(log2FC ~ group, data = d_psd_TET)
aov.9<- aov(log2FC ~ group, data = d_psd_TET)
summary(aov.9)
TukeyHSD(aov.9)

d_pldA_TET <- subset(d2, gene == "pldA")
residuals <- (aov(log2FC ~ group, data = d_pldA_TET))$residuals
shapiro.test(residuals)
leveneTest(log2FC ~ group, data = d_pldA_TET)
aov.10 <- aov(log2FC ~ group, data = d_pldA_TET)
summary(aov.10)
TukeyHSD(aov.10)

ggplot(d2, aes(gene, log2FC, group = group, color = group))+
  theme_bw()+
  theme(panel.border = element_rect(fill = NA, color = "black", size = 2, linetype = "solid"))+
  geom_point(stat = "summary", fun = "mean", alpha = 1.0, size = 6,
             position = position_dodge(width = 0.75))+
  geom_errorbar(stat = "summary",
                fun.min = function(x) mean(x) - sd(x),
                fun.max = function(x) mean(x) + sd(x),
                width = 0, size = 2.25,
                position = position_dodge(width = 0.75))+
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", size = 1)+
  labs(x = NULL,
       y = expression(Log[2]("FC")),
       title = expression("TET"))+
  scale_y_continuous(limits = c(-5, 1), labels = label_number(accuracy = 0.1))+
  scale_fill_manual(values = c("Ctrl" = "#C0C0C0", "DTAC" = "#bd1e2f", "SDS" = "#1f4e9f", "TX-100" = "#7367BE"))+
  scale_color_manual(values = c("Ctrl" = "#C0C0C0", "DTAC" = "#bd1e2f", "SDS" = "#1f4e9f", "TX-100" = "#7367BE"))+
  theme(axis.text.x = element_text(size = 21, face = "bold.italic", angle = 0))+
  theme(axis.text.y = element_text(size = 21, face = "bold"))+
  theme(axis.title.y = element_text(size = 24, vjust = 2))+
  theme(legend.text = element_text(size = 21, face = "bold"))+
  theme(legend.title = element_blank())+
  theme(legend.position = "bottom")+
  theme(title = element_text(size = 20))

ggsave("TET_gene_expression_2.png", width = 1350/90, height = 480/90, dpi = 600, units = "in")



d3 <- read.csv("TET_gene_expression_3.csv")
head(d3)
d3$group <- factor(d3$group, level = c("Ctrl", "0.02", "0.04", "0.2", "1", "5", "10"))

residuals <- (aov(log2FC~ group, data = d3))$residuals
shapiro.test(residuals)
leveneTest(log2FC ~ group, data = d3)
aov.11 <- aov(log2FC ~ group, data = d3)
summary(aov.11)
TukeyHSD(aov.11)

ggplot(d3, aes(x = group, y = log2FC, color = group))+
  theme_bw()+
  geom_point(stat = "summary", fun = "mean", alpha = 1.0, size = 6)+
  geom_errorbar(stat = "summary",
                fun.min = function(x) mean(x) - sd(x),
                fun.max = function(x) mean(x) + sd(x),
                width = 0, size = 2.25)+
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", size = 1)+
  labs(x = expression("DTAC treatment (μmol/L)"),
       y = expression(Log[2]("FC")~of~italic(acrB)),
       title = expression("TET"))+
  scale_y_continuous(limits = c(-2, 6), labels = label_number(accuracy = 0.1))+
  scale_fill_brewer(palette = "Reds")+
  scale_color_brewer(palette = "Reds")+
  theme(axis.text.x = element_text(size = 21, face = "bold", vjust = 0.4, angle = 0),
        axis.text.y = element_text(size = 21, face = "bold"),
        axis.title.x = element_text(size = 24, face = "bold", vjust = -1),
        axis.title.y = element_text(size = 24, face = "bold", vjust = 2),
        legend.position = "none",
        title = element_text(size = 20, face = "bold"),
        panel.border = element_rect(fill = NA, color = "black", size = 2, linetype = "solid"),
        plot.title = element_text(margin = margin(t = 0, b = 1.5, unit = "pt")))

ggsave("TET_gene_expression_3.png", width = 660/90, height = 480/90, dpi = 600, units = "in")



d4 <- read.csv("CEF_gene_expression_1.csv")
head(d4)
d4$gene <- factor(d4$gene, level = c("marA", "ampC", "acrB", "ompF"))

d_marA_CEF <- subset(d4, gene == "marA")
residuals <- (aov(log2FC ~ group, data = d_marA_CEF))$residuals
shapiro.test(residuals)
leveneTest(log2FC ~ group, data = d_marA_CEF)
aov.12 <- aov(log2FC ~ group, data = d_marA_CEF)
summary(aov.12)
TukeyHSD(aov.12)

d_ampC_CEF <- subset(d4, gene == "ampC")
residuals <- (aov(log2FC ~ group, data = d_ampC_CEF))$residuals
shapiro.test(residuals)
leveneTest(log2FC ~ group, data = d_ampC_CEF)
aov.13 <- aov(log2FC ~ group, data = d_ampC_CEF)
summary(aov.13)
TukeyHSD(aov.13)

d_acrB_CEF <- subset(d4, gene == "acrB")
residuals <- (aov(log2FC ~ group, data = d_acrB_CEF))$residuals
shapiro.test(residuals)
leveneTest(log2FC ~ group, data = d_acrB_CEF)
aov.14 <- aov(log2FC ~ group, data = d_acrB_CEF)
summary(aov.14)
TukeyHSD(aov.14)

d_ompF_CEF <- subset(d4, gene == "ompF")
residuals <- (aov(log2FC ~ group, data = d_ompF_CEF))$residuals
shapiro.test(residuals)
leveneTest(log2FC ~ group, data = d_ompF_CEF)
aov.15 <- aov(log2FC ~ group, data = d_ompF_CEF)
summary(aov.15)
TukeyHSD(aov.15)

ggplot(d4, aes(gene, log2FC, group = group, color = group))+
  theme_bw()+
  theme(panel.border = element_rect(fill = NA, color = "black", size = 2, linetype = "solid"))+
  geom_point(stat = "summary", fun = "mean", alpha = 1.0, size = 6,
             position = position_dodge(width = 0.75))+
  geom_errorbar(stat = "summary",
                fun.min = function(x) mean(x) - sd(x),
                fun.max = function(x) mean(x) + sd(x),
                width = 0, size = 2.25,
                position = position_dodge(width = 0.75))+
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", size = 1)+
  labs(x = NULL,
       y = expression(Log[2]("FC")),
       title = expression("CEF"))+
  scale_y_continuous(limits = c(-2.5, 7.5))+
  scale_fill_manual(values = c("Ctrl" = "#C0C0C0", "DTAC" = "#bd1e2f", "SDS" = "#1f4e9f", "TX-100" = "#7367BE"))+
  scale_color_manual(values = c("Ctrl" = "#C0C0C0", "DTAC" = "#bd1e2f", "SDS" = "#1f4e9f", "TX-100" = "#7367BE"))+
  theme(axis.text.x = element_text(size = 21, face = "bold.italic", angle = 0))+
  theme(axis.text.y = element_text(size = 21, face = "bold"))+
  theme(axis.title.y = element_text(size = 24, vjust = 2))+
  theme(legend.text = element_text(size = 21, face = "bold"))+
  theme(legend.title = element_blank())+
  theme(legend.position = "bottom")+
  theme(title = element_text(size = 20))

ggsave("CEF_gene_expression_1.png", width = 900/90, height = 480/90, dpi = 600, units = "in")



d5 <- read.csv("CEF_gene_expression_2.csv")
head(d5)
d5$group <- factor(d5$group, level = c("Ctrl", "0.02", "0.04", "0.2", "1", "5", "10"))

residuals <- (aov(log2FC~ group, data = d5))$residuals
shapiro.test(residuals)
leveneTest(log2FC ~ group, data = d5)
aov.16 <- aov(log2FC ~ group, data = d5)
summary(aov.16)
TukeyHSD(aov.16)

ggplot(d5, aes(x = group, y = log2FC, color = group))+
  theme_bw()+
  geom_point(stat = "summary", fun = "mean", alpha = 1.0, size = 6)+
  geom_errorbar(stat = "summary",
                fun.min = function(x) mean(x) - sd(x),
                fun.max = function(x) mean(x) + sd(x),
                width = 0, size = 2.25)+
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", size = 1)+
  labs(x = expression("TX-100 treatment (μmol/L)"),
       y = expression(Log[2]("FC")~of~italic(ampC)),
       title = expression("CEF"))+
  scale_y_continuous(limits = c(-2.5, 7.5))+
  scale_fill_brewer(palette = "Purples")+
  scale_color_brewer(palette = "Purples")+
  theme(axis.text.x = element_text(size = 21, face = "bold", vjust = 0.4, angle = 0),
        axis.text.y = element_text(size = 21, face = "bold"),
        axis.title.x = element_text(size = 24, face = "bold", vjust = -1),
        axis.title.y = element_text(size = 24, face = "bold", vjust = 2),
        legend.position = "none",
        title = element_text(size = 20, face = "bold"),
        panel.border = element_rect(fill = NA, color = "black", size = 2, linetype = "solid"),
        plot.title = element_text(margin = margin(t = 0, b = 1.5, unit = "pt")))

ggsave("CEF_gene_expression_2.png", width = 660/90, height = 480/90, dpi = 600, units = "in")

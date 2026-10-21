library(ggplot2)
library(ggrepel)

d1 <- read.csv("CEF_DTAC_vs_Ctrl.csv")
head(d1)

d1$diff_expression <- "NS"
d1$diff_expression[d1$log2FC >= 1 & d1$p_adj < 0.05] <- "Up"
d1$diff_expression[d1$log2FC <= -1 & d1$p_adj < 0.05] <- "Down" 

d1$diff_expression <- factor(d1$diff_expression, levels = c("Down", "NS", "Up"))

sig_genes_1 <- d1[order(d1$p_adj), ][1:10, ]

ggplot(data = d1, aes(x = log2FC, y = -log10(p_adj)))+
  theme_bw()+
  theme(panel.border = element_rect(fill = NA, color = "black", size = 2, linetype = "solid"))+
  geom_point(aes(color = diff_expression), alpha = 0.8, size = 3)+
  scale_color_manual(values = c("#1f4e9f", "#a9a9a9", "#bd1e2f"))+
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black", size = 1)+
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black", size = 1)+
  scale_x_continuous(limits = c(-7.75, 7.75))+
  labs(title = "DTAC vs. Ctrl")+
  xlab(expression(Log[2]("FC")))+
  ylab(expression(-Log[10](italic(p)*"-adj")))+
  theme(axis.text.x = element_text(size = 18, face = "bold", vjust = 0.4, angle = 0))+
  theme(axis.text.y = element_text(size = 18, face = "bold"))+
  theme(axis.title.x = element_text(size = 22, vjust = 1))+
  theme(axis.title.y = element_text(size = 22, vjust = 2))+
  theme(title = element_text(size = 18))+
  theme(legend.text = element_text(size = 15, face = "bold"))+
  theme(legend.title = element_blank())+
  geom_text_repel(
    data = sig_genes_1,
    aes(label = name),
    size = 3.5,
    fontface = "italic",
    box.padding = 0.5, 
    segment.color = "grey50",
    segment.alpha = 0.5,
    max.overlaps = 20)

ggsave("CEF_DTAC_vs_Ctrl.png", width = 680/90, height = 435/90, dpi = 600, unit = "in")


d2 <- read.csv("CEF_SDS_vs_Ctrl.csv")
head(d2)

d2$diff_expression <- "NS"
d2$diff_expression[d2$log2FC >= 1 & d2$p_adj < 0.05] <- "Up"
d2$diff_expression[d2$log2FC <= -1 & d2$p_adj < 0.05] <- "Down" 

d2$diff_expression <- factor(d2$diff_expression, levels = c("Down", "NS", "Up"))

sig_genes_2 <- d2[order(d2$p_adj), ][1:10, ]

ggplot(data = d2, aes(x = log2FC, y = -log10(p_adj)))+
  theme_bw()+
  theme(panel.border = element_rect(fill = NA, color = "black", size = 2, linetype = "solid"))+
  geom_point(aes(color = diff_expression), alpha = 0.8, size = 3)+
  scale_color_manual(values = c("#1f4e9f", "#a9a9a9", "#bd1e2f"))+
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black", size = 1)+
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black", size = 1)+
  scale_x_continuous(limits = c(-6, 6))+
  scale_y_continuous(limits = c(0, 115))+
  labs(title = "SDS vs. Ctrl")+
  xlab(expression(Log[2]("FC")))+
  ylab(expression(-Log[10](italic(p)*"-adj")))+
  theme(axis.text.x = element_text(size = 18, face = "bold", vjust = 0.4, angle = 0))+
  theme(axis.text.y = element_text(size = 18, face = "bold"))+
  theme(axis.title.x = element_text(size = 22, vjust = 1))+
  theme(axis.title.y = element_text(size = 22, vjust = 2))+
  theme(title = element_text(size = 18))+
  theme(legend.text = element_text(size = 15, face = "bold"))+
  theme(legend.title = element_blank())+
  geom_text_repel(
    data = sig_genes_2,
    aes(label = name),
    size = 3.5,
    fontface = "italic",
    box.padding = 0.5, 
    segment.color = "grey50",
    segment.alpha = 0.5,
    max.overlaps = 20)

ggsave("CEF_SDS_vs_Ctrl.png", width = 680/90, height = 435/90, dpi = 600, unit = "in")



d3 <- read.csv("CEF_TX_vs_Ctrl.csv")
head(d3)

d3$diff_expression <- "NS"
d3$diff_expression[d3$log2FC >= 1 & d3$p_adj < 0.05] <- "Up"
d3$diff_expression[d3$log2FC <= -1 & d3$p_adj < 0.05] <- "Down" 

d3$diff_expression <- factor(d3$diff_expression, levels = c("Down", "NS", "Up"))

sig_genes_3 <- d3[order(d3$p_adj), ][1:10, ]

ggplot(data = d3, aes(x = log2FC, y = -log10(p_adj)))+
  theme_bw()+
  theme(panel.border = element_rect(fill = NA, color = "black", size = 2, linetype = "solid"))+
  geom_point(aes(color = diff_expression), alpha = 0.8, size = 3)+
  scale_color_manual(values = c("#1f4e9f", "#a9a9a9", "#bd1e2f"))+
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black", size = 1)+
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black", size = 1)+
  scale_x_continuous(limits = c(-6, 6))+
  scale_y_continuous(limits = c(0, 115))+
  labs(title = "TX-100 vs. Ctrl")+
  xlab(expression(Log[2]("FC")))+
  ylab(expression(-Log[10](italic(p)*"-adj")))+
  theme(axis.text.x = element_text(size = 18, face = "bold", vjust = 0.4, angle = 0))+
  theme(axis.text.y = element_text(size = 18, face = "bold"))+
  theme(axis.title.x = element_text(size = 22, vjust = 1))+
  theme(axis.title.y = element_text(size = 22, vjust = 2))+
  theme(title = element_text(size = 18))+
  theme(legend.text = element_text(size = 15, face = "bold"))+
  theme(legend.title = element_blank())+
  geom_text_repel(
    data = sig_genes_3,
    aes(label = name),
    size = 3.5,
    fontface = "italic",
    box.padding = 0.5, 
    segment.color = "grey50",
    segment.alpha = 0.5,
    max.overlaps = 20)

ggsave("CEF_TX-100_vs_Ctrl.png", width = 680/90, height = 435/90, dpi = 600, unit = "in")



d4 <- read.csv("TET_DTAC_vs_Ctrl.csv")
head(d4)

d4$diff_expression <- "NS"
d4$diff_expression[d4$log2FC >= 1 & d4$p_adj < 0.05] <- "Up"
d4$diff_expression[d4$log2FC <= -1 & d4$p_adj < 0.05] <- "Down" 

d4$diff_expression <- factor(d4$diff_expression, levels = c("Down", "NS", "Up"))

sig_genes_4 <- d4[order(d4$p_adj), ][1:10, ]

ggplot(data = d4, aes(x = log2FC, y = -log10(p_adj)))+
  theme_bw()+
  theme(panel.border = element_rect(fill = NA, color = "black", size = 2, linetype = "solid"))+
  geom_point(aes(color = diff_expression), alpha = 0.8, size = 3)+
  scale_color_manual(values = c("#1f4e9f", "#a9a9a9", "#bd1e2f"))+
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black", size = 1)+
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black", size = 1)+
  scale_x_continuous(limits = c(-5, 5))+
  scale_y_continuous(limits = c(0, 25))+
  labs(title = "DTAC vs. Ctrl")+
  xlab(expression(Log[2]("FC")))+
  ylab(expression(-Log[10](italic(p)*"-adj")))+
  theme(axis.text.x = element_text(size = 18, face = "bold", vjust = 0.4, angle = 0))+
  theme(axis.text.y = element_text(size = 18, face = "bold"))+
  theme(axis.title.x = element_text(size = 22, vjust = 1))+
  theme(axis.title.y = element_text(size = 22, vjust = 2))+
  theme(title = element_text(size = 18))+
  theme(legend.text = element_text(size = 15, face = "bold"))+
  theme(legend.title = element_blank())+
  geom_text_repel(
    data = sig_genes_4,
    aes(label = name),
    size = 3.5,
    fontface = "italic",
    box.padding = 0.5, 
    segment.color = "grey50",
    segment.alpha = 0.5,
    max.overlaps = 20)

ggsave("TET_DTAC_vs_Ctrl.png", width = 680/90, height = 435/90, dpi = 600, unit = "in")


d5 <- read.csv("TET_SDS_vs_Ctrl.csv")
head(d5)

d5$diff_expression <- "NS"
d5$diff_expression[d5$log2FC >= 1 & d5$p_adj < 0.05] <- "Up"
d5$diff_expression[d5$log2FC <= -1 & d5$p_adj < 0.05] <- "Down" 

d5$diff_expression <- factor(d5$diff_expression, levels = c("Down", "NS", "Up"))

sig_genes_5 <- d5[order(d5$p_adj), ][1:10, ]

ggplot(data = d5, aes(x = log2FC, y = -log10(p_adj)))+
  theme_bw()+
  theme(panel.border = element_rect(fill = NA, color = "black", size = 2, linetype = "solid"))+
  geom_point(aes(color = diff_expression), alpha = 0.8, size = 3)+
  scale_color_manual(values = c("#1f4e9f", "#a9a9a9", "#bd1e2f"))+
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black", size = 1)+
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black", size = 1)+
  scale_x_continuous(limits = c(-5, 5))+
  scale_y_continuous(limits = c(0, 10))+
  labs(title = "SDS vs. Ctrl")+
  xlab(expression(Log[2]("FC")))+
  ylab(expression(-Log[10](italic(p)*"-adj")))+
  theme(axis.text.x = element_text(size = 18, face = "bold", vjust = 0.4, angle = 0))+
  theme(axis.text.y = element_text(size = 18, face = "bold"))+
  theme(axis.title.x = element_text(size = 22, vjust = 1))+
  theme(axis.title.y = element_text(size = 22, vjust = 2))+
  theme(title = element_text(size = 18))+
  theme(legend.text = element_text(size = 15, face = "bold"))+
  theme(legend.title = element_blank())+
  geom_text_repel(
    data = sig_genes_5,
    aes(label = name),
    size = 3.5,
    fontface = "italic",
    box.padding = 0.5, 
    segment.color = "grey50",
    segment.alpha = 0.5,
    max.overlaps = 20)

ggsave("TET_SDS_vs_Ctrl.png", width = 680/90, height = 435/90, dpi = 600, unit = "in")



d6 <- read.csv("TET_TX_vs_Ctrl.csv")
head(d6)

d6$diff_expression <- "NS"
d6$diff_expression[d6$log2FC >= 1 & d6$p_adj < 0.05] <- "Up"
d6$diff_expression[d6$log2FC <= -1 & d6$p_adj < 0.05] <- "Down" 

d6$diff_expression <- factor(d6$diff_expression, levels = c("Down", "NS", "Up"))

sig_genes_6 <- d6[order(d6$p_adj), ][1:10, ]

ggplot(data = d6, aes(x = log2FC, y = -log10(p_adj)))+
  theme_bw()+
  theme(panel.border = element_rect(fill = NA, color = "black", size = 2, linetype = "solid"))+
  geom_point(aes(color = diff_expression), alpha = 0.8, size = 3)+
  scale_color_manual(values = c("#1f4e9f", "#a9a9a9", "#bd1e2f"))+
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black", size = 1)+
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black", size = 1)+
  scale_x_continuous(limits = c(-5, 5))+
  scale_y_continuous(limits = c(0, 10))+
  labs(title = "TX-100 vs. Ctrl")+
  xlab(expression(Log[2]("FC")))+
  ylab(expression(-Log[10](italic(p)*"-adj")))+
  theme(axis.text.x = element_text(size = 18, face = "bold", vjust = 0.4, angle = 0))+
  theme(axis.text.y = element_text(size = 18, face = "bold"))+
  theme(axis.title.x = element_text(size = 22, vjust = 1))+
  theme(axis.title.y = element_text(size = 22, vjust = 2))+
  theme(title = element_text(size = 18))+
  theme(legend.text = element_text(size = 15, face = "bold"))+
  theme(legend.title = element_blank())+
  geom_text_repel(
    data = sig_genes_6,
    aes(label = name),
    size = 3.5,
    fontface = "italic",
    box.padding = 0.5, 
    segment.color = "grey50",
    segment.alpha = 0.5,
    max.overlaps = 20)

ggsave("TET_TX-100_vs_Ctrl.png", width = 680/90, height = 435/90, dpi = 600, unit = "in")

library(ggplot2)

d_gene_1 <- read.csv("CEF_mechanism_heatmap.csv")
d_gene_1$gene_new <- factor(d_gene_1$gene, levels = rev(c("csgD", "rfbC", "waaB", "waaO", "rfaJ", "rcsA", "glf", "galE", "cydX",
                                                          "lsrB", "lsrF", "rrf", "gabD", "nrdE", "guaC", "carA", "argB", "argC")))
                                                             
d_gene_1$signif <- ifelse(d_gene_1$p_adj < 0.05, "*", "")

ggplot(d_gene_1, aes(group, gene_new))+
  geom_tile(aes(fill = log2FC), color="black", size = 1)+
  theme_minimal()+
  geom_text(aes(label = signif),
            nudge_y = -0.2,
            color = "black",
            size = 4.5,             
            fontface = "bold")+
  theme(axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_text(size = 12, face = "bold", color = "black", angle = 90, hjust = 0.75, vjust = 0.5),
        axis.text.y = element_text(size = 12, face = "italic", color = "black"),
        legend.position = "top")+
  scale_y_discrete(position="right")+
  labs(fill = expression(bold(Log[2]("FC"))))+
  theme(legend.title = element_text(size = 12, face = "bold"),
        legend.text = element_text(size = 10, face = "bold"),
        legend.position = "right")+
  scale_fill_gradient2(limits = c(-3, 5.5),
                       breaks = c(-2.5, 0, 2.5, 5),
                       low = "#1A5592", mid = "white", high = "#B83D3D")

ggsave("CEF_mechanism_heatmap.png", width = 205/90, height = 373/90, dpi = 600, unit = "in")



d_gene_2 <- read.csv("TET_mechanism_heatmap.csv")
d_gene_2$gene_new <- factor(d_gene_2$gene, levels = rev(c("srlE", "gatA", "gatB", "lpdA", "aceE", "aceF", "mqo", "sdhA", "cyoE", "rplB", "rpmC",
                                                          "gcd", "gntK", "ppsA", "dhaL", "ndk", "hcp", "phsC", "carB", "putA", "tdcB", "tdcG")))
d_gene_2$signif <- ifelse(d_gene_2$p_adj < 0.05, "*", "")

ggplot(d_gene_2, aes(group, gene_new))+
  geom_tile(aes(fill = log2FC), color = "black", size = 1)+
  theme_minimal()+
  geom_text(aes(label = signif),
            nudge_y = -0.2,
            color = "black",
            size = 4.5,             
            fontface = "bold")+
  theme(axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_text(size = 12, face = "bold", color = "black", angle = 90, hjust = 0.75, vjust = 0.5),
        axis.text.y = element_text(size = 12, face = "italic", color = "black"),
        legend.position = "top")+
  scale_y_discrete(position = "right")+
  labs(fill = expression(bold(Log[2]("FC"))))+
  theme(legend.title = element_text(size = 12, face = "bold"),
        legend.text = element_text(size = 10, face = "bold"),
        legend.position = "right")+
  scale_fill_gradient2(limits = c(-3.25, 1.75),
                       breaks = c(-3, -1.5, 0, 1.5),
                       low = "#1A5592", mid = "white", high = "#B83D3D")

ggsave("TET_mechanism_heatmap.png", width = 205/90, height = 440/90, dpi = 600, unit="in")

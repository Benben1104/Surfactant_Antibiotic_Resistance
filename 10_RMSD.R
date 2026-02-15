library(ggplot2)
library(RColorBrewer)

d <- read.csv("RMSD.csv")
head(d)

ggplot(d, aes(x = time, y = RMSD, color = group))+
  theme_bw()+
  theme(panel.border = element_rect(fill = NA, color = "black", size = 2, linetype = "solid"))+
  scale_y_continuous(limits = c(0, 1))+
  geom_point(size = 1.5, alpha = 0.8)+
  theme(strip.text.x = element_text(size = 16))+
  labs(x = expression("Time (ns)"),
       y = expression("RMSD (nm)"),
       color = NULL)+
  scale_color_manual(values = c("DTAC-AcrB" = "#bd1e2f", "SDS-AcrB" = "#1f4e9f", "TX-100-AcrB" = "#7367BE"),
                     breaks = c("DTAC-AcrB", "SDS-AcrB", "TX-100-AcrB"))+
  theme(axis.text.x = element_text(size = 18, face = "bold", vjust = 0))+
  theme(axis.text.y = element_text(size = 18, face = "bold"))+
  theme(axis.title.x = element_text(size = 20, face = "bold", vjust = -0.5))+
  theme(axis.title.y = element_text(size = 20, face = "bold", vjust = 2))+
  theme(title = element_text(size = 18, face = "bold"))+
  theme(legend.text = element_text(size = 18, face = "bold"))+
  theme(legend.position = "bottom")
 
ggsave("RMSD.png", width = 580/90, height = 450/90, dpi = 600, unit = "in")

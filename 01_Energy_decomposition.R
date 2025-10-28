library(ggplot2)

d <- read.csv("Energy_decomposition.csv")
head(d)

d_CEF <- subset(d, antibiotic == "CEF")

d_CEF$sign <- ifelse(d_CEF$value >= 0, "Repulsive", "Attractive")
d_CEF$group <- factor(d_CEF$group, level = c("TX-100", "SDS", "DTAC"))
d_CEF$type <- factor(d_CEF$type, levels = c("electrostatic", "repulsive", "dispersion"))

ggplot(d_CEF, aes(x = group, y = value, fill = sign,
                  group = interaction(group, type)))+
  theme_bw()+
  theme(panel.border = element_rect(fill = NA, color = "black", size = 2, linetype = "solid"))+
  geom_col(position = position_dodge(width = 0.75), alpha = 0.85)+
  coord_flip()+
  scale_fill_manual(values = c("Attractive" = "#1f4e9f",
                               "Repulsive" = "#bd1e2f"))+
  scale_y_continuous(limits = c(-300, 100))+
  labs(title = "CEF",
       x = NULL,
       y = "Energy (kJ"~mol^"-1"*")",
       fill = NULL)+
  theme(axis.text.x = element_text(size = 18, face = "bold", vjust = 0.4, angle = 0))+
  theme(axis.text.y = element_text(size = 18, face = "bold"))+
  theme(axis.title.x = element_text(size = 20, vjust = 1))+
  theme(axis.title.y = element_text(size = 20, vjust = 2))+
  theme(title = element_text(size = 18))+
  theme(legend.text = element_text(size = 15, face = "bold"))

ggsave("CEF_energy_decomposition.png", width = 590/90, height = 360/90, dpi = 600, unit = "in")


d_TET <- subset(d, antibiotic == "TET")

d_TET$sign <- ifelse(d_TET$value >= 0, "Repulsive", "Attractive")
d_TET$group <- factor(d_TET$group, level = c("TX-100", "SDS", "DTAC"))
d_TET$type <- factor(d_TET$type, levels = c("electrostatic", "repulsive", "dispersion"))

ggplot(d_TET, aes(x = group, y = value, fill = sign,
                  group = interaction(group, type)))+
  theme_bw()+
  theme(panel.border = element_rect(fill = NA, color= "black", size = 2, linetype = "solid"))+
  geom_col(position = position_dodge(width = 0.75), alpha = 0.85)+
  coord_flip()+
  scale_fill_manual(values = c("Attractive" = "#1f4e9f",
                               "Repulsive" = "#bd1e2f"))+
  scale_y_continuous(limits = c(-100, 10))+
  labs(title = "TET",
       x = NULL,
       y = "Energy (kJ"~mol^"-1"*")",
       fill = NULL)+
  theme(axis.text.x=element_text(size = 18, face = "bold", vjust = 0.4, angle = 0))+
  theme(axis.text.y = element_text(size = 18, face = "bold"))+
  theme(axis.title.x = element_text(size = 20, vjust = 1))+
  theme(axis.title.y = element_text(size = 20, vjust = 2))+
  theme(title = element_text(size = 18))+
  theme(legend.text = element_text(size = 15, face = "bold"))

ggsave("TET_energy_decomposition.png", width = 590/90, height = 360/90, dpi = 600, unit = "in")

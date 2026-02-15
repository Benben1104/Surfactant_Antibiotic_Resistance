library(ggplot2)
library(dplyr)
library(tidyr)
library(stringr)

d1_wide <- read.csv("CEF_PL_composition.csv")

d1_long <- gather(d1_wide, key = "group", value = "percentage", -'type')
head(d1_long)

d1_long <- d1_long %>%
  mutate(group_new = case_when(
    str_detect(group, "Ctrl\\.") ~ "Ctrl",
    str_detect(group, "DTAC\\.") ~ "DTAC",
    str_detect(group, "SDS\\.") ~ "SDS",
    str_detect(group, "TX\\.100\\.") ~ "TX-100")) %>%
  mutate(group_new = factor(group_new, levels = c("Ctrl", "DTAC", "SDS", "TX-100")),
         type_new = factor(d1_long$type, levels = c("PE", "PG", "PS", "PA", "LPE")))

ggplot(d1_long, aes(type_new, percentage, fill = group_new))+
  theme_bw()+
  theme(panel.border = element_rect(fill = NA, color = "black", size = 2, linetype = "solid"))+
  stat_summary(geom = "col", fun = "mean",
               position = "dodge", width = 0.75, alpha = 0.8)+
  coord_polar(start = 0)+
  labs(x = NULL,
       y = "Relative abundance (%)",
       title = "CEF")+
  scale_y_continuous(limits = c(0, 80))+
  scale_fill_manual(values = c("Ctrl" = "#C0C0C0", "DTAC" = "#bd1e2f", "SDS" = "#1f4e9f", "TX-100" = "#7367BE"))+
  scale_color_manual(values = c("Ctrl" = "#C0C0C0", "DTAC" = "#bd1e2f", "SDS" = "#1f4e9f", "TX-100" = "#7367BE"))+
  theme(axis.text.y = element_text(margin = margin(r = -157)))+
  theme(axis.line.y = element_blank(), axis.ticks.y = element_blank())+
  theme(axis.text.x = element_text(size = 18, face = "bold", angle = 0))+
  theme(axis.text.y = element_text(size = 10, face = "bold"))+
  theme(axis.title.y = element_text(size = 11, vjust = 2, hjust = 0.86))+
  theme(legend.text = element_text(size = 18, face = "bold"))+
  theme(legend.title = element_blank())+
  theme(legend.position = "bottom")+
  theme(title = element_text(size = 20))

ggsave("CEF_PL_composition.png", width = 695/90, height = 485/90, dpi = 600, units = "in")



d2_wide <- read.csv("TET_PL_composition.csv")

d2_long <- gather(d2_wide, key = "group", value = "percentage", -'type')
head(d2_long)

d2_long <- d2_long %>%
  mutate(group_new = case_when(
    str_detect(group, "Ctrl\\.") ~ "Ctrl",
    str_detect(group, "DTAC\\.") ~ "DTAC",
    str_detect(group, "SDS\\.") ~ "SDS",
    str_detect(group, "TX\\.100\\.") ~ "TX-100")) %>%
  mutate(group_new = factor(group_new, levels = c("Ctrl", "DTAC", "SDS", "TX-100")),
         type_new = factor(d1_long$type, levels = c("PE", "PG", "PS", "PA", "LPE")))

ggplot(d2_long, aes(type_new, percentage, fill = group_new))+
  theme_bw()+
  theme(panel.border = element_rect(fill = NA, color = "black", size = 2, linetype = "solid"))+
  stat_summary(geom = "col", fun = "mean",
               position = "dodge", width = 0.75, alpha = 0.8)+
  coord_polar(start = 0)+
  labs(x = NULL,
       y = "Relative abundance (%)",
       title = "TET")+
  scale_y_continuous(limits = c(0, 80))+
  scale_fill_manual(values = c("Ctrl" = "#C0C0C0", "DTAC" = "#bd1e2f", "SDS" = "#1f4e9f", "TX-100" = "#7367BE"))+
  scale_color_manual(values = c("Ctrl" = "#C0C0C0", "DTAC" = "#bd1e2f", "SDS" = "#1f4e9f", "TX-100" = "#7367BE"))+
  theme(axis.text.y = element_text(margin = margin(r = -157)))+
  theme(axis.line.y = element_blank(), axis.ticks.y = element_blank())+
  theme(axis.text.x = element_text(size = 18, face = "bold", angle = 0))+
  theme(axis.text.y = element_text(size = 10, face = "bold"))+
  theme(axis.title.y = element_text(size = 11, vjust = 2, hjust = 0.86))+
  theme(legend.text = element_text(size = 18, face = "bold"))+
  theme(legend.title = element_blank())+
  theme(legend.position = "bottom")+
  theme(title = element_text(size = 20))

ggsave("TET_PL_composition.png", width = 695/90, height = 485/90, dpi = 600, units = "in")

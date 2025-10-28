library(ggplot2)
library(dplyr)
library(tidyr)

d1 <- read.csv("DTAC-AcrB_energy.csv")

d1_labeled <- d1 %>%
  mutate(original_order = row_number()) %>%
  mutate(rank = min_rank(energy)) %>%
  mutate(
    group = case_when(
      rank <= 10 ~ "attractive",
      rank >= max(rank, na.rm = TRUE) - 9 ~ "repulsive",
      TRUE ~ "other"
    )
  ) %>%
  arrange(original_order) %>%
  select(-rank, -original_order)

ggplot(d1_labeled, aes(amino.acid.residue, energy, color = group))+
  theme_bw()+
  theme(panel.border = element_rect(fill = NA, color = "black", size = 2, linetype = "solid"))+
  theme(panel.grid.major.x  = element_blank(), panel.grid.minor.x = element_blank(),
        axis.text.x = element_blank())+
  coord_polar(start = 0)+
  geom_point(size = 2.75)+
  labs(x = NULL,
       y = "Energy (kJ"~mol^"-1"*")",
       title = "DTAC-AcrB")+
  scale_y_continuous(limits = c(-34, 42))+
  scale_color_manual(values = c("other" = "#E89E6D", "attractive" = "#1f4e9f", "repulsive" = "#bd1e2f"))+
  theme(axis.text.y = element_text(margin = margin(r = -165)))+
  theme(axis.line.y = element_blank(), axis.ticks.y = element_blank())+
  theme(axis.text.y = element_text(size = 10, face = "bold"))+
  theme(axis.title.y = element_text(size = 11, vjust = 0, hjust = 0.86))+
  theme(legend.position = "none")+
  theme(title = element_text(size = 20))

ggsave("DTAC-AcrB.png", width = 695/90, height = 455/90, dpi = 600, units = "in")



d2 <- read.csv("SDS-AcrB_energy.csv")

d2_labeled <- d2 %>%
  mutate(original_order = row_number()) %>%
  mutate(rank = min_rank(energy)) %>%
  mutate(
    group = case_when(
      rank <= 10 ~ "attractive",
      rank >= max(rank, na.rm = TRUE) - 9 ~ "repulsive",
      TRUE ~ "other"
    )
  ) %>%
  arrange(original_order) %>%
  select(-rank, -original_order)

ggplot(d2_labeled, aes(amino.acid.residue, energy, color = group))+
  theme_bw()+
  theme(panel.border = element_rect(fill = NA, color = "black", size = 2, linetype = "solid"))+
  theme(panel.grid.major.x  = element_blank(), panel.grid.minor.x = element_blank(),
        axis.text.x = element_blank())+
  coord_polar(start = 0)+
  geom_point(size = 2.75)+
  labs(x = NULL,
       y = "Energy (kJ"~mol^"-1"*")",
       title = "SDS-AcrB")+
  scale_y_continuous(limits = c(-26, 22))+
  scale_color_manual(values = c("other" = "#8CAE7E", "attractive" = "#1f4e9f", "repulsive" = "#bd1e2f"))+
  theme(axis.text.y = element_text(margin = margin(r = -165)))+
  theme(axis.line.y = element_blank(), axis.ticks.y = element_blank())+
  theme(axis.text.y = element_text(size = 10, face = "bold"))+
  theme(axis.title.y = element_text(size = 11, vjust = 0, hjust = 0.86))+
  theme(legend.position = "none")+
  theme(title = element_text(size = 20))

ggsave("SDS-AcrB.png", width = 695/90, height = 455/90, dpi = 600, units = "in")



d3 <- read.csv("TX-100-AcrB_energy.csv")

d3_labeled <- d3 %>%
  mutate(original_order = row_number()) %>%
  mutate(rank = min_rank(energy)) %>%
  mutate(
    group = case_when(
      rank <= 10 ~ "attractive",
      rank >= max(rank, na.rm = TRUE) - 9 ~ "repulsive",
      TRUE ~ "other"
    )
  ) %>%
  arrange(original_order) %>%
  select(-rank, -original_order)

ggplot(d3_labeled, aes(amino.acid.residue, energy, color = group))+
  theme_bw()+
  theme(panel.border = element_rect(fill = NA, color = "black", size = 2, linetype = "solid"))+
  theme(panel.grid.major.x  = element_blank(), panel.grid.minor.x = element_blank(),
        axis.text.x = element_blank())+
  coord_polar(start = 0)+
  geom_point(size = 2.75)+
  labs(x = NULL,
       y = "Energy (kJ"~mol^"-1"*")",
       title = "TX-100-AcrB")+
  scale_y_continuous(limits = c(-2.5, 1.2))+
  scale_color_manual(values = c("other" = "#B29AC9", "attractive" = "#1f4e9f", "repulsive" = "#bd1e2f"))+
  theme(axis.text.y = element_text(margin = margin(r = -165)))+
  theme(axis.line.y = element_blank(), axis.ticks.y = element_blank())+
  theme(axis.text.y = element_text(size = 10, face = "bold"))+
  theme(axis.title.y = element_text(size = 11, vjust = 0, hjust = 0.86))+
  theme(legend.position = "none")+
  theme(title = element_text(size = 20))

ggsave("TX-100-AcrB.png", width = 695/90, height = 455/90, dpi = 600, units = "in")

library(ggplot2)
library(dplyr)
library(tidyr)
library(stringr)
library(car)
library(scales)

d1_wide <- read.csv("CEF_PL_composition.csv")

d1_long <- gather(d1_wide, key = "group", value = "percentage", -'type')
head(d1_long)

d1_long <- d1_long %>%
  mutate(group_new = case_when(
    str_detect(group, "Ctrl\\.") ~ "Ctrl",
    str_detect(group, "DTAC\\.") ~ "DTAC",
    str_detect(group, "SDS\\.") ~ "SDS",
    str_detect(group, "TX\\.100\\.") ~ "TX-100")) %>%
  mutate(group_new = factor(group_new, levels = c("Ctrl", "DTAC", "SDS", "TX-100")))

d1_PE <- subset(d1_long, type == "PE")
residuals <- (aov(percentage ~ group_new, data = d1_PE))$residuals
shapiro.test(residuals)
leveneTest(percentage ~ group_new, data = d1_PE)
aov_PE <- aov(percentage ~ group_new, data = d1_PE)
summary(aov_PE)
TukeyHSD(aov_PE)

d1_PG <- subset(d1_long, type == "PG")
residuals <- (aov(percentage ~ group_new, data = d1_PG))$residuals
shapiro.test(residuals)
leveneTest(percentage ~ group_new, data = d1_PG)
aov_PG <- aov(percentage ~ group_new, data = d1_PG)
summary(aov_PG)
TukeyHSD(aov_PG)

d1_PS <- subset(d1_long, type == "PS")
residuals <- (aov(percentage ~ group_new, data = d1_PS))$residuals
shapiro.test(residuals)
leveneTest(percentage ~ group_new, data = d1_PS)
aov_PS <- aov(percentage ~ group_new, data = d1_PS)
summary(aov_PS)
TukeyHSD(aov_PS)

d1_PA <- subset(d1_long, type == "PA")
residuals <- (aov(percentage ~ group_new, data = d1_PA))$residuals
shapiro.test(residuals)
leveneTest(percentage ~ group_new, data = d1_PA)
aov_PA <- aov(percentage ~ group_new, data = d1_PA)
summary(aov_PA)
TukeyHSD(aov_PA)

ggplot(d1_PE, aes(group_new, percentage, fill = group_new))+
  theme_bw()+
  theme(panel.border = element_rect(fill = NA, color = "black", size = 2, linetype = "solid"))+
  geom_boxplot(width = 0.6, size = 1.25, alpha = 0.7)+
  geom_jitter(aes(group_new, percentage),
              width = 0.15, size = 4.5, shape = 21, stroke = 0.01, alpha = 1)+
  labs(x = NULL,
       y = "%PE",
       title = "CEF")+
  scale_y_continuous(limits = c(62, 67), labels = scales::number_format(accuracy = 0.1))+
  scale_fill_manual(values = c("Ctrl" = "#C0C0C0", "DTAC" = "#bd1e2f", "SDS" = "#1f4e9f", "TX-100" = "#7367BE"))+
  scale_color_manual(values = c("Ctrl" = "#C0C0C0", "DTAC" = "#bd1e2f", "SDS" = "#1f4e9f", "TX-100" = "#7367BE"))+
  theme(axis.text.x = element_text(size = 18, face = "bold", angle = 0))+
  theme(axis.text.y = element_text(size = 18, face = "bold"))+
  theme(axis.title.y = element_text(size = 20, vjust = 2.25))+
  theme(legend.position = "none")+
  theme(title = element_text(size = 18))

ggsave("CEF_PE.png", width = 497/90, height = 375/90, dpi = 600, unit = "in")


ggplot(d1_PG, aes(group_new, percentage, fill = group_new))+
  theme_bw()+
  theme(panel.border = element_rect(fill = NA, color = "black", size = 2, linetype = "solid"))+
  geom_boxplot(width = 0.6, size = 1.25, alpha = 0.7)+
  geom_jitter(aes(group_new, percentage),
              width = 0.15, size = 4.5, shape = 21, stroke = 0.01, alpha = 1)+
  labs(x = NULL,
       y = "%PG",
       title = "CEF")+
  scale_y_continuous(limits = c(15, 17.5), labels = scales::number_format(accuracy = 0.1))+
  scale_fill_manual(values = c("Ctrl" = "#C0C0C0", "DTAC" = "#bd1e2f", "SDS" = "#1f4e9f", "TX-100" = "#7367BE"))+
  scale_color_manual(values = c("Ctrl" = "#C0C0C0", "DTAC" = "#bd1e2f", "SDS" = "#1f4e9f", "TX-100" = "#7367BE"))+
  theme(axis.text.x = element_text(size = 18, face = "bold", angle = 0))+
  theme(axis.text.y = element_text(size = 18, face = "bold"))+
  theme(axis.title.y = element_text(size = 20, vjust = 2.25))+
  theme(legend.position = "none")+
  theme(title = element_text(size = 18))

ggsave("CEF_PG.png", width = 497/90, height = 375/90, dpi = 600, unit = "in")


ggplot(d1_PS, aes(group_new, percentage, fill = group_new))+
  theme_bw()+
  theme(panel.border = element_rect(fill = NA, color = "black", size = 2, linetype = "solid"))+
  geom_boxplot(width = 0.6, size = 1.25, alpha = 0.7)+
  geom_jitter(aes(group_new, percentage),
              width = 0.15, size = 4.5, shape = 21, stroke = 0.01, alpha = 1)+
  labs(x = NULL,
       y = "%PS",
       title = "CEF")+
  scale_y_continuous(limits = c(12, 16), labels = scales::number_format(accuracy = 0.1))+
  scale_fill_manual(values = c("Ctrl" = "#C0C0C0", "DTAC" = "#bd1e2f", "SDS" = "#1f4e9f", "TX-100" = "#7367BE"))+
  scale_color_manual(values = c("Ctrl" = "#C0C0C0", "DTAC" = "#bd1e2f", "SDS" = "#1f4e9f", "TX-100" = "#7367BE"))+
  theme(axis.text.x = element_text(size = 18, face = "bold", angle = 0))+
  theme(axis.text.y = element_text(size = 18, face = "bold"))+
  theme(axis.title.y = element_text(size = 20, vjust = 2.25))+
  theme(legend.position = "none")+
  theme(title = element_text(size = 18))

ggsave("CEF_PS.png", width = 497/90, height = 375/90, dpi = 600, unit = "in")


ggplot(d1_PA, aes(group_new, percentage, fill = group_new))+
  theme_bw()+
  theme(panel.border = element_rect(fill = NA, color = "black", size = 2, linetype = "solid"))+
  geom_boxplot(width = 0.6, size = 1.25, alpha = 0.7)+
  geom_jitter(aes(group_new, percentage),
              width = 0.15, size = 4.5, shape = 21, stroke = 0.01, alpha = 1)+
  labs(x = NULL,
       y = "%PA",
       title = "CEF")+
  scale_y_continuous(limits = c(4, 6), labels = scales::number_format(accuracy = 0.1))+
  scale_fill_manual(values = c("Ctrl" = "#C0C0C0", "DTAC" = "#bd1e2f", "SDS" = "#1f4e9f", "TX-100" = "#7367BE"))+
  scale_color_manual(values = c("Ctrl" = "#C0C0C0", "DTAC" = "#bd1e2f", "SDS" = "#1f4e9f", "TX-100" = "#7367BE"))+
  theme(axis.text.x = element_text(size = 18, face = "bold", angle = 0))+
  theme(axis.text.y = element_text(size = 18, face = "bold"))+
  theme(axis.title.y = element_text(size = 20, vjust = 2.25))+
  theme(legend.position = "none")+
  theme(title = element_text(size = 18))

ggsave("CEF_PA.png", width = 485/90, height = 375/90, dpi = 600, unit = "in")





d2_wide <- read.csv("TET_PL_composition.csv")

d2_long <- gather(d2_wide, key = "group", value = "percentage", -'type')
head(d2_long)

d2_long <- d2_long %>%
  mutate(group_new = case_when(
    str_detect(group, "Ctrl\\.") ~ "Ctrl",
    str_detect(group, "DTAC\\.") ~ "DTAC",
    str_detect(group, "SDS\\.") ~ "SDS",
    str_detect(group, "TX\\.100\\.") ~ "TX-100")) %>%
  mutate(group_new = factor(group_new, levels = c("Ctrl", "DTAC", "SDS", "TX-100")))

d2_PE <- subset(d2_long, type == "PE")
residuals <- (aov(percentage ~ group_new, data = d2_PE))$residuals
shapiro.test(residuals)
leveneTest(percentage ~ group_new, data = d2_PE)
aov_PE <- aov(percentage ~ group_new, data = d2_PE)
summary(aov_PE)
TukeyHSD(aov_PE)

d2_PG <- subset(d2_long, type == "PG")
residuals <- (aov(percentage ~ group_new, data = d2_PG))$residuals
shapiro.test(residuals)
leveneTest(percentage ~ group_new, data = d2_PG)
aov_PG <- aov(percentage ~ group_new, data = d2_PG)
summary(aov_PG)
TukeyHSD(aov_PG)

d2_PS <- subset(d2_long, type == "PS")
residuals <- (aov(percentage ~ group_new, data = d2_PS))$residuals
shapiro.test(residuals)
leveneTest(percentage ~ group_new, data = d2_PS)
aov_PS <- aov(percentage ~ group_new, data = d2_PS)
summary(aov_PS)
TukeyHSD(aov_PS)

d2_PA <- subset(d2_long, type == "PA")
residuals <- (aov(percentage ~ group_new, data = d2_PA))$residuals
shapiro.test(residuals)
leveneTest(percentage ~ group_new, data = d2_PA)
aov_PA <- aov(percentage ~ group_new, data = d2_PA)
summary(aov_PA)
TukeyHSD(aov_PA)

ggplot(d2_PE, aes(group_new, percentage, fill = group_new))+
  theme_bw()+
  theme(panel.border = element_rect(fill = NA, color = "black", size = 2, linetype = "solid"))+
  geom_boxplot(width = 0.6, size = 1.25, alpha = 0.7)+
  geom_jitter(aes(group_new, percentage),
              width = 0.15, size = 4.5, shape = 21, stroke = 0.01, alpha = 1)+
  labs(x = NULL,
       y = "%PE",
       title = "TET")+
  scale_y_continuous(limits = c(75, 79), labels = scales::number_format(accuracy = 0.1))+
  scale_fill_manual(values = c("Ctrl" = "#C0C0C0", "DTAC" = "#bd1e2f", "SDS" = "#1f4e9f", "TX-100" = "#7367BE"))+
  scale_color_manual(values = c("Ctrl" = "#C0C0C0", "DTAC" = "#bd1e2f", "SDS" = "#1f4e9f", "TX-100" = "#7367BE"))+
  theme(axis.text.x = element_text(size = 18, face = "bold", angle = 0))+
  theme(axis.text.y = element_text(size = 18, face = "bold"))+
  theme(axis.title.y = element_text(size = 20, vjust = 2.25))+
  theme(legend.position = "none")+
  theme(title = element_text(size = 18))

ggsave("TET_PE.png", width = 497/90, height = 375/90, dpi = 600, unit = "in")


ggplot(d2_PG, aes(group_new, percentage, fill = group_new))+
  theme_bw()+
  theme(panel.border = element_rect(fill = NA, color = "black", size = 2, linetype = "solid"))+
  geom_boxplot(width = 0.6, size = 1.25, alpha = 0.7)+
  geom_jitter(aes(group_new, percentage),
              width = 0.15, size = 4.5, shape = 21, stroke = 0.01, alpha = 1)+
  labs(x = NULL,
       y = "%PG",
       title = "TET")+
  scale_y_continuous(limits = c(13, 18), labels = scales::number_format(accuracy = 0.1))+
  scale_fill_manual(values = c("Ctrl" = "#C0C0C0", "DTAC" = "#bd1e2f", "SDS" = "#1f4e9f", "TX-100" = "#7367BE"))+
  scale_color_manual(values = c("Ctrl" = "#C0C0C0", "DTAC" = "#bd1e2f", "SDS" = "#1f4e9f", "TX-100" = "#7367BE"))+
  theme(axis.text.x = element_text(size = 18, face = "bold", angle = 0))+
  theme(axis.text.y = element_text(size = 18, face = "bold"))+
  theme(axis.title.y = element_text(size = 20, vjust = 2.25))+
  theme(legend.position = "none")+
  theme(title = element_text(size = 18))

ggsave("TET_PG.png", width = 497/90, height = 375/90, dpi = 600, unit = "in")


ggplot(d2_PS, aes(group_new, percentage, fill = group_new))+
  theme_bw()+
  theme(panel.border = element_rect(fill = NA, color = "black", size = 2, linetype = "solid"))+
  geom_boxplot(width = 0.6, size = 1.25, alpha = 0.7)+
  geom_jitter(aes(group_new, percentage),
              width = 0.15, size = 4.5, shape = 21, stroke = 0.01, alpha = 1)+
  labs(x = NULL,
       y = "%PS",
       title = "TET")+
  scale_y_continuous(limits = c(3, 7), labels = scales::number_format(accuracy = 0.1))+
  scale_fill_manual(values = c("Ctrl" = "#C0C0C0", "DTAC" = "#bd1e2f", "SDS" = "#1f4e9f", "TX-100" = "#7367BE"))+
  scale_color_manual(values = c("Ctrl" = "#C0C0C0", "DTAC" = "#bd1e2f", "SDS" = "#1f4e9f", "TX-100" = "#7367BE"))+
  theme(axis.text.x = element_text(size = 18, face = "bold", angle = 0))+
  theme(axis.text.y = element_text(size = 18, face = "bold"))+
  theme(axis.title.y = element_text(size = 20, vjust = 2.25))+
  theme(legend.position = "none")+
  theme(title = element_text(size = 18))

ggsave("TET_PS.png", width = 485/90, height = 375/90, dpi = 600, unit = "in")


ggplot(d2_PA, aes(group_new, percentage, fill = group_new))+
  theme_bw()+
  theme(panel.border = element_rect(fill = NA, color = "black", size = 2, linetype = "solid"))+
  geom_boxplot(width = 0.6, size = 1.25, alpha = 0.7)+
  geom_jitter(aes(group_new, percentage),
              width = 0.15, size = 4.5, shape = 21, stroke = 0.01, alpha = 1)+
  labs(x = NULL,
       y = "%PA",
       title = "TET")+
  scale_y_continuous(limits = c(1, 3), labels = scales::number_format(accuracy = 0.1))+
  scale_fill_manual(values = c("Ctrl" = "#C0C0C0", "DTAC" = "#bd1e2f", "SDS" = "#1f4e9f", "TX-100" = "#7367BE"))+
  scale_color_manual(values = c("Ctrl" = "#C0C0C0", "DTAC" = "#bd1e2f", "SDS" = "#1f4e9f", "TX-100" = "#7367BE"))+
  theme(axis.text.x = element_text(size = 18, face = "bold", angle = 0))+
  theme(axis.text.y = element_text(size = 18, face = "bold"))+
  theme(axis.title.y = element_text(size = 20, vjust = 2.25))+
  theme(legend.position = "none")+
  theme(title = element_text(size = 18))

ggsave("TET_PA.png", width = 485/90, height = 375/90, dpi = 600, unit = "in")

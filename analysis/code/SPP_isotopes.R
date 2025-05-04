library(tidyverse)
library(readxl)
library(stringr)
isotopes_2023 <- read_excel(here::here("analysis","data","raw_data","SPP-manual-qunat-Jun-2023.xls"))
SBP_ceramics <- read_excel(here::here("analysis","data","raw_data","SPP5_ceramic_samples.xlsx"), sheet =2)

# tidy metadata
SBP_ceramics <-
  SBP_ceramics %>%
  mutate(Sam_no = as.character(`流水號`))

# samples with low yields
lowyield <- c("SPP006", "SPP017", "SPP032", "SPP010", "SPP007", "SPP-008", "SPP-028")

# datasets for analysis
SBP_isotopes_rm <-
  isotopes_2023 %>%
  filter(str_detect(Filename, "SPP")) %>%
  filter(!str_detect(Filename, paste(lowyield, collapse = "|"))) %>% # remove low yields
  filter(`Ampl  44 (mV)` > 250) # background noise

# tidy up isotope data
tidy_SBP_isotopes <-
  SBP_isotopes_rm %>%
  select("Filename", "Rt (s)", "d13C_PDB (permil)") %>%
  mutate(Filename = sapply(str_extract_all(
    Filename, "SPP[0-9]+|SPP-[0-9]+|SPP-D[0-9]+-[0-9]|SPP-[:alpha:]+[0-9]"), toString)) %>%
  mutate(FA = case_when(`Rt (s)`<1260&`Rt (s)`>1240  ~ "C16",
                        `Rt (s)`<1369&`Rt (s)`>1352  ~ "C18")) %>%
  group_by(Filename, FA) %>%
  summarize(Ave = mean(`d13C_PDB (permil)`)) %>%
  pivot_wider(names_from = FA, values_from = Ave) %>%
  dplyr::mutate(Delta = C18 -C16) %>%
  mutate(Sam_no = sapply(str_extract_all(Filename, "(?<=SPP).*"), toString)) %>%
  mutate(Sam_no = str_remove_all(Sam_no, "^-")) %>% #"[0-9]+|D[0-9]+-[0-9]|[A-Z]+[0-9]"
  mutate(Sam_no = sapply(str_remove(Sam_no, "^0+"), toString)) %>%
  left_join(SBP_ceramics) %>%
  mutate(group = case_when(`所屬文化` == "大湖烏山頭期" ~ "Neolithic-Wushantou",
                            `所屬文化` == "蔦松早期" ~ "Iron Age",
                            `所屬文化` == "疑似大湖晚期" ~ "Neolithic-Final", TRUE ~ "Bone")) %>%
  mutate(context = case_when(str_detect(`出土脈絡`, "文化") ~ "cultural layer",
                             str_detect(`出土脈絡`, "灰坑") ~ "midden",
                             str_detect(`出土脈絡`, "墓葬") ~ "burial"))

# plot isotopes of 16 and C18 with reference range
library(png)
library(ggpubr)
tem <- readPNG("tem.png")
tem2 <- readPNG(here::here("analysis", "figures","iso_ellipse.png"))

isotope_period <-
  tidy_SBP_isotopes %>%
  filter(!group == "Bone") %>%
  ggplot(aes(C16,C18)) +
  geom_point(size = 2, alpha = 0.9, aes(color = group))+ # period
  ggrepel::geom_text_repel(aes(label = Sam_no), size = 2) +
  labs(x = bquote(delta*{}^13*"C 16:0 \u2030"),
       y = bquote(delta*{}^13*"C 18:0 \u2030")) +
  xlim(-40,-10) +
  ylim(-40,-10) +
  coord_fixed(ratio=1) +
  geom_abline(intercept = -1, linetype = "dashed") +
  geom_abline(intercept = -3.1, linetype = "dashed") +
  annotate("text", x = -12.5, y = -11, angle = 45, vjust = 1,
           label = bquote(Delta*{}^13*"C = -1 \u2030")) +
  annotate("text", x = -12, y = -15.5, angle = 45, vjust = 1,
           label = bquote(Delta*{}^13*"C = -3.3 \u2030")) +
  annotation_raster(tem2, ymin = -40, ymax= -10 ,
                    xmin = -40.1 , xmax = -10) +
  annotate("text", size = 2.5, x = -35, y = -30.5, label = "FW") +
  annotate("text", size = 2.5, x = -28, y = -33, label = "R") +
  annotate("text", size = 2.5, x = -28.5, y = -25, label = "P") +
  annotate("text", size = 2.5, x = -22, y = -17, label = "M") +
  theme_minimal(base_size = 14)

ggsave(here::here("SBP_isotope_period.png"), width = 8, height = , unit = "in")

# plot isotopes of 16 and C18 by context
isotope_context <-
  tidy_SBP_isotopes %>%
  filter(!is.na(context)) %>%
  ggplot(aes(C16, C18)) +
  geom_point(size = 2, alpha = 0.9, aes(color = context)) +
  ggrepel::geom_text_repel(aes(label = Sam_no)) +
  labs(x = bquote(delta*{}^13*"C 16:0 \u2030"),
       y = bquote(delta*{}^13*"C 18:0 \u2030")) +
  scale_x_continuous(limits = c(-35, -20),
                     breaks = seq(-35, -20, 2)) +
  scale_y_continuous(limits = c(-35, -20),
                     breaks = seq(-35, -20, 2)) +
  xlim(-40,-10) +
  ylim(-40,-10) +
  coord_fixed(ratio=1) +
  geom_abline(intercept = -1, linetype = "dashed") + #Suryanarayan et al.2021
  geom_abline(intercept = -3.1, linetype = "dashed") +
  annotation_raster(tem2, ymin = -40, ymax= -10 ,
                    xmin = -40.1 , xmax = -10) +
  annotate("text", size = 2.5, x = -35, y = -30.5, label = "FW") +
  annotate("text", size = 2.5, x = -28, y = -33, label = "R") +
  annotate("text", size = 2.5, x = -28.5, y = -25, label = "P") +
  annotate("text", size = 2.5, x = -22, y = -17, label = "M") +
  theme_minimal()

# big delta
Delta <-
  ggplot(tidy_SBP_isotopes, # delta_meth_cor or delta_blank_cor
       aes(C16, Delta)) +
  geom_point(size = 2, alpha = 0.9, aes(color = group))+
  ggrepel::geom_text_repel(aes(label = Sam_no)) +
  labs(x = bquote(delta*{}^13*"C 16:0 \u2030"),
       y = bquote(Delta*{}^13*"C")) +
  scale_x_continuous(limits = c(-33, -21),
                     breaks = seq(-33, -21, 2)) +
  scale_y_continuous(limits = c(-5, 2.5),
                     breaks = seq(-4, 2, 2)) +
  geom_hline(yintercept = c(-1, -3.1), linetype = "dashed") +
  annotate("text", x = -32.5, y = -0.9,
           label = bquote(Delta*{}^13*"C = -1 \u2030")) +
  annotate("text", x = -32.5, y = -3,
           label = bquote(Delta*{}^13*"C = -3.1 \u2030")) +
  annotate("text", x = c(-32, -21.5), y = c(2.3, 2.3),
           label = c("C3-diet", "C4-diet")) +
  annotate("segment", x = -29, y = 2.3, xend = -24, yend = 2.3,
           arrow = arrow(type = "closed", length = unit(0.02, "npc"))) +
  annotate("text", x = c(-21.5,-21.5,-21.5), y = c(-0.5, -1.5, -3.6),
           label = c("Non-ruminant", "Ruminant", "Ruminant dairy")) +
  theme_minimal()


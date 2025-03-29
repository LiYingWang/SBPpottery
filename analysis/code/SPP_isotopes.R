library(tidyverse)
library(readxl)
library(stringr)
isotopes_2023 <- read_excel(here::here("analysis","data","raw_data","SPP-manual-qunat-Jun-2023.xls"))
SBP_ceramics <- read_excel(here::here("analysis","data","raw_data","SPP5脂肪酸分析陶片挑選清單_1110.xlsx"))

# tidy metadata
SBP_ceramics <-
  SBP_ceramics %>%
  mutate(Sam_no = as.character(`流水號`))

# join two datasets
SBP_isotopes <-
  rbind(isotopes_2023, isotopes_2024) %>%
  filter(str_detect(Filename, "SPP"))

# remove low yields
SPP_isotopes_rm <-
  rbind(isotopes_2023, isotopes_2024) %>%
  filter(`Ampl  44 (mV)` > 150 & str_detect(Filename, "SPP"))

# tidy up
tidy_SPP_isotopes <-
  SPP_isotopes_rm %>% # SPP_isotopes
  select("Filename", "Rt (s)","d13C_PDB (permil)") %>%
  mutate(Filename = sapply(str_extract_all(Filename,
                                           "SPP[0-9]+|SPP-[0-9]+|SPP-D[0-9]+-[0-9]|SPP-[:alpha:]+[0-9]"), toString)) %>%
  mutate(FA = case_when(`Rt (s)`<1260&`Rt (s)`>1240  ~ "C16",
                        `Rt (s)`<1370&`Rt (s)`>1352  ~ "C18")) %>%
  group_by(Filename, FA) %>%
  filter(!is.na(FA)) %>%
  summarize(Ave = mean(`d13C_PDB (permil)`)) %>%
  pivot_wider(names_from = FA, values_from = Ave) %>%
  dplyr::mutate(Delta = C18 -C16) %>%
  mutate(Sam_no = sapply(str_extract_all(Filename, "(?<=SPP).*"), toString)) %>%
  mutate(Sam_no = str_remove_all(Sam_no, "^-")) %>% #"[0-9]+|D[0-9]+-[0-9]|[A-Z]+[0-9]"
  mutate(Sam_no = sapply(str_remove(Sam_no, "^0+"), toString)) %>%
  left_join(SPP_ceramics) %>%
  mutate(period = case_when(`所屬文化` == "大湖烏山頭期" ~ "Late Neolithic",
                            `所屬文化` == "蔦松早期" ~ "Iron Age",
                            `所屬文化` == "疑似大湖晚期" ~ "Final Neolithic")) %>%
  mutate(period = ifelse(str_detect(Sam_no, "D"), "Iron Age", period)) %>%
  mutate(group = ifelse(is.na(period), "Bone", period))

# plot isotopes of 16 and C18 with reference range
library(png)
library(ggpubr)
tem <- readPNG("tem.png")

ggplot(tidy_SPP_isotopes,
       aes(C16,C18)) +
  geom_point(size = 2, alpha = 0.9, aes(color = group))+ #period
  ggrepel::geom_text_repel(aes(label = Sam_no), size = 2) +
  theme_minimal(base_size = 14) +
  labs(x = bquote(delta*{}^13*"C 16:0 \u2030"),
       y = bquote(delta*{}^13*"C 18:0 \u2030")) +
  xlim(-40,-10) +
  ylim(-40,-10) +
  coord_fixed(ratio=1) +
  annotation_raster(tem, ymin = -42.8, ymax= -8 ,
                    xmin = -44.5 , xmax = -5.65)

# plot isotopes of 16 and C18 using a different reference range
tem2 <- readPNG("tem2.png")
tidy_SPP_isotopes %>%
  ggplot(aes(C16, C18)) +
  geom_point(size = 2, alpha = 0.9, aes(color = group))+
  ggrepel::geom_text_repel(aes(label = Sam_no), size = 2) +
  theme_minimal(base_size = 14) +
  labs(x = bquote(delta*{}^13*"C 16:0 \u2030"),
       y = bquote(delta*{}^13*"C 18:0 \u2030")) +
  xlim(-40,-20) +
  ylim(-40,-20) +
  coord_fixed(ratio=1) +
  #scale_colour_viridis_d(direction = -1) +
  annotation_raster(tem2, ymin = -45.1, ymax= -15.1 ,
                    xmin = -45.55, xmax = -15.55)

# plot isotopes of 16 and C18 with reference line only
ggplot(tidy_SPP_isotopes,
       aes(C16, C18)) +
  geom_point(size = 2, alpha = 0.9, aes(color = group)) +
  ggrepel::geom_text_repel(aes(label = Sam_no)) +
  theme_minimal() +
  labs(x = bquote(delta*{}^13*"C 16:0 \u2030"),
       y = bquote(delta*{}^13*"C 18:0 \u2030")) +
  scale_x_continuous(limits = c(-35, -20),
                     breaks = seq(-35, -20, 2)) +
  scale_y_continuous(limits = c(-35, -20),
                     breaks = seq(-35, -20, 2)) +
  #coord_fixed(ratio = 1) +
  geom_abline(intercept = -1, linetype = "dashed") + #Suryanarayan et al.2021
  geom_abline(intercept = -3.1, linetype = "dashed") +
  annotate("text", x = -21.5, y = -21.5, angle = 45, vjust = 1.5,
           label = bquote(Delta*{}^13*"C = -1 \u2030")) +
  annotate("text", x = -21.5, y = -24.5, angle = 45, vjust = 1.5,
           label = bquote(Delta*{}^13*"C = -3.1 \u2030"))

# big delta
ggplot(tidy_SPP_isotopes, # delta_meth_cor or delta_blank_cor
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


library(data.table)
library(ggplot2)
library(tidyverse)
library(rstatix)
library(ggpubr)

#############################################################################
#### file reading
#### the major file used for analysis was generated modkit, mijamp, snappy
setwd("")
motif_combo <- fread("modkit_snappy_mijamp_combo.txt", header = T)

#############################################################################
#### file clean up
motif_combo <- motif_combo %>% filter(!(replicate %in% c("replicate_4", "replicate_5")))
motif_combo$motif[motif_combo$motif %in% c("GATCNTm4C","GATCNTm4CV")] <- "GATCNTm4C(V)"

#############################################################################
#### statistical analsysis
motif_combo_replicate <- motif_combo %>%
  group_by(tool, genus, phase, replicate) %>%
  summarise(motif_num=mean(length(unique(motif))), .groups = "drop")

motif_combo_summary <- motif_combo_replicate %>%
  group_by(tool, genus, phase) %>%
  summarise(mean_motif_num=mean(motif_num),
            sd_motif_num=sd(motif_num), .groups = "drop")

motif_combo_compare <- motif_combo_replicate %>%
  group_by(genus, phase) %>%
  t_test(data = ., motif_num ~ tool, paired = T) %>%
  add_significance() %>%
  add_xy_position(x="tool")

motif_combo_validate <- motif_combo %>%
  group_by(motif, genus, phase) %>%
  summarise(motif_occcurence=n(), .groups = "drop")

#############################################################################
#### visualization
ggplot(motif_combo_summary[motif_combo_summary$genus=="Escherichia",], 
       aes(x=tool, y=mean_motif_num, color = phase)) +
  facet_wrap(vars(phase)) +
  geom_point(alpha=0.9, size=5) +
  geom_errorbar(aes(ymin=mean_motif_num-sd_motif_num,
                    ymax =mean_motif_num+sd_motif_num),
                color="black", alpha=0.7, linewidth=0.5, width=0.1) +
  geom_linerange(aes(x = tool, ymin = 0, ymax = mean_motif_num),
                 alpha=0.5) +
  stat_pvalue_manual(motif_combo_compare[motif_combo_compare$genus=="Escherichia",], 
                     hide.ns = TRUE,
                     label = "{p.adj.signif}", size=6, bracket.size=0.8) +
  scale_color_manual("Phase", breaks=c("exponential", "stationary"),
                     values = c("#444574", "#cd534cff")) +
  labs(y="Unique motif numbers") +
  theme_bw() +
  theme(
    axis.title.y = element_text(face = "bold"),
    axis.title.x = element_blank(),
    legend.title = element_text(face = "bold")) +
  theme(text = element_text(size = 19),
        legend.position = "none",
        strip.background = element_blank())

ggplot(motif_combo_summary[motif_combo_summary$genus=="Lactobacillus",], 
       aes(x=tool, y=mean_motif_num, color = phase)) +
  facet_wrap(vars(phase)) +
  geom_point(alpha=0.9, size=5) +
  geom_errorbar(aes(ymin=mean_motif_num-sd_motif_num,
                    ymax =mean_motif_num+sd_motif_num),
                color="black", alpha=0.7, linewidth=0.5, width=0.1) +
  geom_linerange(aes(x = tool, ymin = 0, ymax = mean_motif_num),
                 alpha=0.5) +
  scale_color_manual("Phase", breaks=c("exponential", "stationary"),
                     values = c("#444574", "#cd534cff")) +
  labs(y="Unique motif numbers") +
  theme_bw() +
  theme(
    axis.title.y = element_text(face = "bold"),
    axis.title.x = element_blank(),
    legend.title = element_text(face = "bold")) +
  theme(text = element_text(size = 19),
        legend.position = "none",
        strip.background = element_blank())

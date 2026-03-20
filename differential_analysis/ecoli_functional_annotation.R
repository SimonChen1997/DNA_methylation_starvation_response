library(ggplot2)
library(tidyverse)
library(data.table)
library(grid)
library(janitor)

#############################################################################
### file reading
### these are the files generated from NCBI DAVID
setwd("")
ecoli_annotation_cluster_dna_rna_seq <- fread("annotation_cluster_dna_rna.tsv",sep="\t", header=TRUE)
ecoli_annotation_cluster_dna <- fread("annotation_cluster_methyl.tsv",sep="\t", header=TRUE)
ecoli_annotation_cluster_rna <- fread("annotation_cluster_rna.tsv",sep="\t", header=TRUE)

options(scipen=999)

#############################################################################
### file cleaning
ecoli_annotation_cluster_dna_rna_seq$Term <- gsub("^.*:|^.*~", "", ecoli_annotation_cluster_dna_rna_seq$Term)
ecoli_annotation_cluster_dna$Term <- gsub("^.*:|^.*~", "", ecoli_annotation_cluster_dna$Term)
ecoli_annotation_cluster_rna$Term <- gsub("^.*:|^.*~", "", ecoli_annotation_cluster_rna$Term)

ecoli_annotation_cluster_dna_rna_seq$tpye <- "DMGs_DEGs"
ecoli_annotation_cluster_dna$tpye <- "DMGs"
ecoli_annotation_cluster_rna$tpye <- "DEGs"

ecoli_annotation_cluster_combine <- rbind(ecoli_annotation_cluster_dna_rna_seq, ecoli_annotation_cluster_dna,
                                          ecoli_annotation_cluster_rna)

ecoli_annotation_cluster_combine_enriched <- ecoli_annotation_cluster_combine %>%
  mutate(cluster = paste(tpye, Annotation_Cluster, sep = "_")) %>%
  filter(Enrichment_Score >= 2)

ecoli_annotation_cluster_combine_enriched_sig <- ecoli_annotation_cluster_combine %>%
  mutate(cluster = paste(tpye, Annotation_Cluster, sep = "_")) %>%
  filter(Enrichment_Score >= 2 & FDR <= .05)

#############################################################################
### visualization
ggplot(ecoli_annotation_cluster_combine_enriched, aes(x=Count, y=Term, color=tpye)) +
  facet_grid(vars(Annotation_Cluster), scales = "free") +
  geom_point(alpha=0.9, size=4, position = position_dodge(width = 0.5), shape=16) +
  scale_color_manual("Data", breaks = c("DMGs", "DMGs_DEGs", "DEGs"),
                     values = c("#444574", "#2b8cbe", "#cd534cff")) +
  theme_light() +
  labs(x="Gene number", y= "Funtional annotation") +
  theme(axis.title.x = element_text(face = "bold"),
        axis.title.y = element_text(face = "bold"),
        axis.text.x = element_text(angle = 45, hjust=1),
        legend.position = "bottom",
        legend.title = element_text(face="bold"),
        strip.background = element_blank(),
        strip.text = element_text(size=20, color = "black"),
  ) +
  theme(text = element_text(size=26)) +
  theme(plot.margin = margin(t = 5, r = 10, b = 5, l = 5))

ggplot(ecoli_annotation_cluster_combine_enriched_sig, aes(x=Count, y=Term, color=tpye)) +
  facet_grid(vars(Annotation_Cluster), scales = "free") +
  geom_point(alpha=0.9, size=4, position = position_dodge(width = 0.5), shape=16) +
  scale_color_manual("Data", breaks = c("DMGs", "DMGs_DEGs", "DEGs"),
                     values = c("#444574", "#2b8cbe", "#cd534cff")) +
  theme_light() +
  labs(x="Gene number", y= "Funtional annotation") +
  theme(axis.title.x = element_text(face = "bold"),
        axis.title.y = element_text(face = "bold"),
        axis.text.x = element_text(angle = 45, hjust=1),
        legend.position = "bottom",
        legend.title = element_text(face="bold"),
        strip.background = element_blank(),
        strip.text = element_text(size=20, color = "black"),
  ) +
  theme(text = element_text(size=26)) +
  theme(plot.margin = margin(t = 5, r = 10, b = 5, l = 5))

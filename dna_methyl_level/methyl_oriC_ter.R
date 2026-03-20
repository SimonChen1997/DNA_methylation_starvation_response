library(data.table)
library(ggplot2)
library(tidyverse)
library(rstatix)
library(ggpubr)


#############################################################################
#### file reading
#### the major file used for analysis was generated modkit pileup
#### e. coli
setwd("")
ecoli_m4C_position <- fread("ecoli_primary_pileup_combine_motif_m4c.bed",sep="\t", header=TRUE)
ecoli_m5C_position <- fread("ecoli_primary_pileup_combine_motif_m5c.bed",sep="\t", header=TRUE)
ecoli_m6A_position <- fread("ecoli_primary_pileup_combine_motif_m6a.bed",sep="\t", header=TRUE)

ecoli_oriC <- fread("ecoli_oriC_extend.bed", header=F)
ecoli_ter <- fread("ecoli_ter_extend.bed", header=F)
ecoli_other_regions <- fread("ecoli_other_regions_extend.bed", header=F)

colnames(ecoli_oriC) <- c("chromosome", "start", "end")
colnames(ecoli_ter) <- c("chromosome", "start", "end")
colnames(ecoli_other_regions) <- c("chromosome", "start", "end")

options(scipen=999)

#### l. acidophilus
lacidophilus_m4C_position <- fread("lacidophilus_primary_pileup_combine_motif_m4c.bed",sep="\t", header=TRUE)
lacidophilus_m5C_position <- fread("lacidophilus_primary_pileup_combine_motif_m5c.bed",sep="\t", header=TRUE)
lacidophilus_m6A_position <- fread("lacidophilus_primary_pileup_combine_motif_m6a.bed",sep="\t", header=TRUE)

lacidophilus_oriC <- fread("lacidophilus_oriC_extend.bed", header=F)
lacidophilus_ter <- fread("lacidophilus_ter_extend.bed", header=F)
lacidophilus_other_regions <- fread("lacidophilus_other_regions_extend.bed", header=F)

colnames(lacidophilus_oriC) <- c("chromosome", "start", "end")
colnames(lacidophilus_ter) <- c("chromosome", "start", "end")
colnames(lacidophilus_other_regions) <- c("chromosome", "start", "end")

options(scipen=999)

#############################################################################
#### file claen up
ecoli_m6A_position_filter <- ecoli_m6A_position[ecoli_m6A_position$valid_coverage>=20,]
ecoli_m4C_position_filter <- ecoli_m4C_position[ecoli_m4C_position$valid_coverage>=20,]
ecoli_m5C_position_filter <- ecoli_m5C_position[ecoli_m5C_position$valid_coverage>=20,]

ecoli_m6A_position_filter$modified_base_code[ecoli_m6A_position_filter$modified_base_code=="a"] <- "m6A" 
ecoli_m6A_position_filter <- ecoli_m6A_position_filter%>%
  subset(modified_base_code=="m6A")
ecoli_m4C_position_filter$modified_base_code[ecoli_m4C_position_filter$modified_base_code=="21839"] <- "m4C" 
ecoli_m4C_position_filter <- ecoli_m4C_position_filter%>%
  subset(modified_base_code=="m4C")
ecoli_m5C_position_filter$modified_base_code[ecoli_m5C_position_filter$modified_base_code=="m"] <- "m5C" 
ecoli_m5C_position_filter <- ecoli_m5C_position_filter%>%
  subset(modified_base_code=="m5C")


lacidophilus_m6A_position_filter <- lacidophilus_m6A_position[lacidophilus_m6A_position$valid_coverage>=20,]
lacidophilus_m4C_position_filter <- lacidophilus_m4C_position[lacidophilus_m4C_position$valid_coverage>=20,]
lacidophilus_m5C_position_filter <- lacidophilus_m5C_position[lacidophilus_m5C_position$valid_coverage>=20,]

lacidophilus_m6A_position_filter$modified_base_code[lacidophilus_m6A_position_filter$modified_base_code=="a"] <- "m6A" 
lacidophilus_m6A_position_filter <- lacidophilus_m6A_position_filter%>%
  subset(modified_base_code=="m6A")
lacidophilus_m4C_position_filter$modified_base_code[lacidophilus_m4C_position_filter$modified_base_code=="21839"] <- "m4C" 
lacidophilus_m4C_position_filter <- lacidophilus_m4C_position_filter%>%
  subset(modified_base_code=="m4C")
lacidophilus_m5C_position_filter$modified_base_code[lacidophilus_m5C_position_filter$modified_base_code=="m"] <- "m5C" 
lacidophilus_m5C_position_filter <- lacidophilus_m5C_position_filter%>%
  subset(modified_base_code=="m5C")

#############################################################################
#### functions to calculate DNA methylation level at oriC, ter, and other regions
windows_methyl_form <- function(genome_slide_bed, filter_pileup_file, region_name) {
  fraction_modified_windows <- as.data.frame(matrix(nrow=0, ncol = 10))
  colnames(fraction_modified_windows) <- c("chrom", "start", "end", "modified_base_code", "phase", "replicate",
                                           "sum_valid_coverage", "sum_modified_number", "fraction_modified", "region")
  
  replicate_name_list <- unique(filter_pileup_file$replicate)
  inter_addition_form <- as.data.frame(matrix(nrow=0, ncol = 10))
  colnames(inter_addition_form) <- c("chrom", "start", "end", "modified_base_code","phase", "replicate",
                                     "sum_valid_coverage", "sum_modified_number", "fraction_modified", "region")
  
  n_iter <- length(seq_len(nrow(genome_slide_bed)))
  
  pb <- txtProgressBar(min = 0,      # Minimum value of the progress bar
                       max = n_iter, # Maximum value of the progress bar
                       style = 3,    # Progress bar style (also available style = 1 and style = 2)
                       width = 50,   # Progress bar width. Defaults to getOption("width")
                       char = "=")
  
  for (i in seq_len(nrow(genome_slide_bed))){
    start_label <- as.integer(genome_slide_bed[i, "start"])
    end_label <- as.integer(genome_slide_bed[i, "end"])
    chromosome <- as.character(genome_slide_bed[i, "chromosome"])
    modified_base_code <- as.character(unique(filter_pileup_file$modified_base_code))
    
    inter_form <- subset(filter_pileup_file,
                         filter_pileup_file$start >= start_label & filter_pileup_file$end < end_label & filter_pileup_file$chrom == chromosome)
    
    inter_form_summary <- inter_form %>%
      group_by(chrom, modified_base_code, phase, replicate) %>%
      summarise(sum_valid_coverage = sum(valid_coverage),
                sum_modified_number = sum(modified_number), .groups = "drop")
    inter_form_summary <- inter_form_summary %>%
      mutate(fraction_modified = 100*sum_modified_number/sum_valid_coverage)
    
    inter_form_summary$start <- start_label
    inter_form_summary$end <- end_label
    inter_form_summary$region <- region_name
    
    fraction_modified_windows <- rbind(fraction_modified_windows, inter_form_summary)
    
    
    setTxtProgressBar(pb, i)
  }
  cat("\n done!")
  close(pb) # Close the connection
  return(fraction_modified_windows)
}
#############################################################################
#### execute functions to calculate DNA methylation level
#### ecoli
ecoli_m6a_oric <- windows_methyl_form(ecoli_oriC, ecoli_m6A_position_filter, "oriC")
ecoli_m4c_oric <- windows_methyl_form(ecoli_oriC, ecoli_m4C_position_filter, "oriC")
ecoli_m5c_oric <- windows_methyl_form(ecoli_oriC, ecoli_m5C_position_filter, "oriC")

ecoli_m6a_ter <- windows_methyl_form(ecoli_ter, ecoli_m6A_position_filter, "ter")
ecoli_m4c_ter <- windows_methyl_form(ecoli_ter, ecoli_m4C_position_filter, "ter")
ecoli_m5c_ter <- windows_methyl_form(ecoli_ter, ecoli_m5C_position_filter, "ter")

ecoli_m6a_others <- windows_methyl_form(ecoli_other_regions, ecoli_m6A_position_filter, "others")
ecoli_m4c_others <- windows_methyl_form(ecoli_other_regions, ecoli_m4C_position_filter, "others")
ecoli_m5c_others <- windows_methyl_form(ecoli_other_regions, ecoli_m5C_position_filter, "others")

ecoli_metyhl_region <- rbind(ecoli_m6a_oric, ecoli_m4c_oric, ecoli_m5c_oric, 
                             ecoli_m6a_ter, ecoli_m4c_ter, ecoli_m5c_ter,
                             ecoli_m6a_others, ecoli_m4c_others, ecoli_m5c_others)

ecoli_metyhl_region$region <- factor(ecoli_metyhl_region$region, levels = c("oriC","ter","others"))
ecoli_metyhl_region$genus <- "Escherichia"

#### lacidophilus
lacidophilus_m6a_oric <- windows_methyl_form(lacidophilus_oriC, lacidophilus_m6A_position_filter, "oriC")
lacidophilus_m4c_oric <- windows_methyl_form(lacidophilus_oriC, lacidophilus_m4C_position_filter, "oriC")
lacidophilus_m5c_oric <- windows_methyl_form(lacidophilus_oriC, lacidophilus_m5C_position_filter, "oriC")

lacidophilus_m6a_ter <- windows_methyl_form(lacidophilus_ter, lacidophilus_m6A_position_filter, "ter")
lacidophilus_m4c_ter <- windows_methyl_form(lacidophilus_ter, lacidophilus_m4C_position_filter, "ter")
lacidophilus_m5c_ter <- windows_methyl_form(lacidophilus_ter, lacidophilus_m5C_position_filter, "ter")

lacidophilus_m6a_others <- windows_methyl_form(lacidophilus_other_regions, lacidophilus_m6A_position_filter, "others")
lacidophilus_m4c_others <- windows_methyl_form(lacidophilus_other_regions, lacidophilus_m4C_position_filter, "others")
lacidophilus_m5c_others <- windows_methyl_form(lacidophilus_other_regions, lacidophilus_m5C_position_filter, "others")

lacidophilus_metyhl_region <- rbind(lacidophilus_m6a_oric, lacidophilus_m4c_oric, lacidophilus_m5c_oric, 
                                    lacidophilus_m6a_ter, lacidophilus_m4c_ter, lacidophilus_m5c_ter,
                                    lacidophilus_m6a_others, lacidophilus_m4c_others, lacidophilus_m5c_others)

lacidophilus_metyhl_region$region <- factor(lacidophilus_metyhl_region$region, levels = c("oriC","ter","others"))
lacidophilus_metyhl_region$genus <- "Lactobacillus"

ecoli_lacidophilus_metyhl_region <- rbind(ecoli_metyhl_region, lacidophilus_metyhl_region)
ecoli_lacidophilus_metyhl_region$region <- factor(ecoli_lacidophilus_metyhl_region$region, levels = c("oriC","ter","others"))

#############################################################################
#### statistical analysis
#### e. coli
ecoli_lacidophilus_ter_oric_other_phase_compare <- ecoli_lacidophilus_metyhl_region %>%
  group_by(genus, modified_base_code, region) %>%
  t_test(data = ., fraction_modified ~ phase, paired = T) %>%
  add_significance()

ecoli_lacidophilus_ter_oric_other_phase_compare <- ecoli_lacidophilus_ter_oric_other_phase_compare %>%
  add_xy_position(x="region", dodge = 0.8)

ecoli_acidophilus_region_modified_phase_summary <- ecoli_lacidophilus_metyhl_region %>%
  group_by(genus, modified_base_code, phase, region) %>%
  summarise(mean_fraction_modified = mean(fraction_modified),
            sd_fraction_modified = sd(fraction_modified), .groups = "drop")

#### l. acidophilus
lacidophilus_ter_oric_other_phase_compare <- lacidophilus_metyhl_region %>%
  group_by(modified_base_code, region) %>%
  t_test(data = ., fraction_modified ~ phase, paired = T) %>%
  add_significance()

lacidophilus_ter_oric_other_phase_compare <- lacidophilus_ter_oric_other_phase_compare %>%
  add_xy_position(x="region", dodge = 0.8)

lacidophilus_region_modified_phase_summary <- lacidophilus_metyhl_region %>%
  group_by(modified_base_code, phase, region) %>%
  summarise(mean_fraction_modified = mean(fraction_modified),
            sd_fraction_modified = sd(fraction_modified), .groups = "drop")

#############################################################################
#### visualization
ggplot(ecoli_lacidophilus_metyhl_region[ecoli_lacidophilus_metyhl_region$genus=="Escherichia", ], aes(x=region,y=fraction_modified)) +
  facet_wrap(vars(modified_base_code)) +
  geom_boxplot(aes(color=phase), alpha=0.8) +
  scale_color_manual("Phase", breaks=c("exponential", "stationary"),
                     values = c("#444574", "#cd534cff")) +
  stat_pvalue_manual(ecoli_lacidophilus_ter_oric_other_phase_compare[ecoli_lacidophilus_ter_oric_other_phase_compare$genus=="Escherichia", ], hide.ns = T,
                     label = "{p.signif}", size=9, bracket.size=1) +
  coord_trans(y="log1p") +
  labs(y="Modified fraction (%)") +
  theme_bw() +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(face = "bold"),
        legend.title = element_text(face = "bold")) +
  theme(text = element_text(size = 20),
        axis.title.x = element_blank(),
        legend.position = "bottom",
        strip.background = element_blank())

ggplot(ecoli_lacidophilus_metyhl_region[ecoli_lacidophilus_metyhl_region$genus=="Lactobacillus", ], aes(x=region,y=fraction_modified)) +
  facet_wrap(vars(modified_base_code)) +
  geom_boxplot(aes(color=phase), alpha=0.8) +
  scale_color_manual("Phase", breaks=c("exponential", "stationary"),
                     values = c("#444574", "#cd534cff")) +
  stat_pvalue_manual(ecoli_lacidophilus_ter_oric_other_phase_compare[ecoli_lacidophilus_ter_oric_other_phase_compare$genus=="Lactobacillus", ],, hide.ns = T,
                     label = "{p.signif}", size=9, bracket.size=1) +
  coord_trans(y="log1p") +
  labs(y="Modified fraction (%)") +
  theme_bw() +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(face = "bold"),
        legend.title = element_text(face = "bold")) +
  theme(text = element_text(size = 20),
        axis.title.x = element_blank(),
        legend.position = "bottom",
        strip.background = element_blank())

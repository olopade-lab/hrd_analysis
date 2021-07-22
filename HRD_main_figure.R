library(ggplot2)
library(readxl)
library(tidyr)
library(grid)
library(cowplot)
library(dplyr)
library(rstatix)
library(tidyverse)
library(flextable)
library(RColorBrewer)
library(facetscales)

# Setting up the data
setwd("/Users/jbreynier/Desktop/Research/Olopade Lab/hrd")
sample_info <- read.table("NIG_BC_metadata_clean_173.txt", sep = "\t", header = 1)
sample_info <- sample_info[, c("ID", "Subtype", "RACE_geno")]
colnames(sample_info) <- c("sample", "subtype", "race")
mutation_info <- read_xlsx("mutation_info.xlsx")
sv_counts <- read.csv("sv_counts.csv")
chord_pred_nigerian <- read_xlsx("chord_final/nigerian_chord_pred_delly_manta_lumpy_2vote.xlsx")
chord_columns_keep <- c("sample", "p_BRCA1", "p_BRCA2", "p_hrd", "hr_status", "hrd_type")
chord_pred_nigerian <- chord_pred_nigerian[, chord_columns_keep]
chord_pred_nigerian$group <- "Nigerian"
chord_pred_TCGA <- read_xlsx("chord_final/TCGA_chord_pred_delly_manta_lumpy_2vote.xlsx")
chord_pred_TCGA <- chord_pred_TCGA[, chord_columns_keep]
chord_pred_TCGA$group <- "TCGA"
chord_pred <- rbind(chord_pred_nigerian, chord_pred_TCGA)
sv_signatures <- read.table("SV7SignaturesContributions-173Samples.tsv", sep = "\t", header = 1)
indel_signatures <- read.table("HRD_CHORD_SBS_INDEL_signatures_173_NAP.txt", sep = "\t", header = 1)
indel_signatures$ID <- indel_signatures$ID6 + indel_signatures$ID8
indel_signatures <- indel_signatures[, c("sample", "SBS3", "ID")]
total_results <- merge(chord_pred, sample_info, all.x = TRUE)
total_results <- merge(total_results, mutation_info, all.x = TRUE)
total_results <- merge(total_results, sv_counts, all.x = TRUE)
total_results <- merge(indel_signatures, total_results, all.x = TRUE)
total_results <- merge(total_results, sv_signatures, all.x = TRUE)
total_results$hrd_type[total_results$hrd_type == "none"] <- NA
total_results <- drop_na(total_results, "subtype")

total_results_brca1 <- total_results[which(total_results$hrd_type == "BRCA1_type"), ]
total_results_brca1 <- total_results_brca1[order(total_results_brca1$p_BRCA1), ]
total_results_brca2 <- total_results[which(total_results$hrd_type == "BRCA2_type"), ]
total_results_brca2 <- total_results_brca2[order(total_results_brca2$p_BRCA2), ]
total_results_undetermined <- total_results[which(total_results$hrd_type == "cannot_be_determined"), ]
total_results_undetermined <- total_results_undetermined[order(total_results_undetermined$p_hrd), ]
total_results_hrp <- total_results[which(total_results$hr_status == "HR_proficient"), ]
total_results_hrp <- total_results_hrp[order(total_results_hrp$p_hrd), ]
total_results <- bind_rows(total_results_hrp, total_results_undetermined, total_results_brca2, total_results_brca1)
total_results <- mutate(group_by(total_results, race), order = row_number())
total_results$race <- factor(total_results$race, levels=c('Nigerian','Black','White'))
colbrew=brewer.pal(7,"Dark2")
colors_signatures=c("violet",colbrew[5],"red2", "turquoise1","orange",colbrew[c(2:3)])
scales_x <- list(
  `Nigerian` = scale_x_continuous(limits = c(0.5, sum(total_results$race == "Nigerian")+0.5),
                                  breaks = seq(0, sum(total_results$race == "Nigerian")+0.5, 15),
                                  expand = c(0, 0)),
  `Black` = scale_x_continuous(limits = c(0.5, sum(total_results$race == "Black")+0.5),
                               breaks = seq(0, sum(total_results$race == "Black")+0.5, 15),
                               expand = c(0, 0)),
  `White` = scale_x_continuous(limits = c(0.5, sum(total_results$race == "White")+0.5),
                               breaks = seq(0, sum(total_results$race == "White")+0.5, 15),
                               expand = c(0, 0))
)

make_chord_score_axis <- function(df_results) {
  
  df_results_stackedformat_hrd <- pivot_longer(df_results, cols=c(p_BRCA1, p_BRCA2),
                                               names_to = "p_mutation", values_to = "p_mutation_score")
  df_results_stackedformat_hrd$p_mutation <- gsub('p_', '',
                                                  df_results_stackedformat_hrd$p_mutation)
  df_results_stackedformat_hrd$p_mutation <- factor(df_results_stackedformat_hrd$p_mutation,
                                                    levels=c('BRCA2','BRCA1'))
  chord_score_axis <- ggplot(df_results_stackedformat_hrd, aes(y=p_mutation, x=order)) + 
    geom_tile(aes(fill=p_mutation_score), show.legend=TRUE) + theme_bw() + 
    theme(rect = element_rect(fill = "transparent"), plot.background = element_blank()) +
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          panel.grid.minor.x=element_blank(),
          panel.grid.major.x=element_blank(),
          strip.switch.pad.grid = unit(0, "cm")) +
    theme(axis.title.y = element_text(size = 8,
                                      margin = margin(t = 0, r = 10, b = 0, l = 0)),
          axis.text.y = element_text(size = 7),
          panel.grid.minor.y=element_blank(),
          panel.grid.major.y=element_blank()) +
    theme(legend.key.size = unit(0.03, "npc"),
          legend.title = element_blank(),
          legend.text = element_text(size =7),
          legend.direction = "horizontal") +
    scale_fill_gradient(low="white", high="darkviolet", limits = c(0, 1)) +
    scale_y_discrete(breaks = c("BRCA2", "BRCA1"),
                     labels = c(expression(italic("BRCA2")), expression(italic("BRCA1"))), expand = c(0, 0)) + 
    ylab("CHORD HRD\nprobability") + 
    facet_grid_sc(cols = vars(race), scales = list(x = scales_x), space = "free_x") +
    theme(plot.margin=unit(c(0.05,0,0,0.025),"npc")) +
    theme(strip.text.x = element_text(size = 9, face = "bold"))
  return(chord_score_axis)
}


make_chord_prediction_axis <- function(df_results) {
  
  chord_prediction_axis <- ggplot(df_results, aes(y=1, x=order)) + 
    geom_tile(aes(fill=hrd_type), show.legend=TRUE) + theme_bw() + 
    theme(rect = element_rect(fill = "transparent"), plot.background = element_blank()) +
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          panel.grid.minor.x=element_blank(),
          panel.grid.major.x=element_blank()) +
    theme(axis.title.y = element_text(size = 8,
                                      margin = margin(t = 0, r = 10, b = 0, l = 0)),
          axis.text.y=element_blank(),
          axis.ticks.y = element_blank(),
          panel.grid.minor.y=element_blank(),
          panel.grid.major.y=element_blank()) +
    theme(legend.key.size = unit(0.03, "npc"),
          legend.title = element_blank(),
          legend.text = element_text(size =7)) + 
    scale_fill_manual(breaks = c("BRCA1_type", 'BRCA2_type', "cannot_be_determined"),
                      labels = c(expression(paste(italic("BRCA1"), "-like HR deficiency")),
                                 expression(paste(italic("BRCA2"), "-like HR deficiency")),
                                 "Unspecified HR deficiency"),
                      values = c("#F8766D", "#619CFF", "gray"),
                      na.value = "white") + 
    scale_y_discrete(expand = c(0, 0)) +
    facet_grid_sc(cols = vars(race), scales = list(x = scales_x), space = "free_x") +
    theme(strip.background = element_blank(),
          strip.text.x = element_blank()) +
    theme(plot.margin=unit(c(0,0,0,0.025),"npc")) + 
    ylab("CHORD HRD\ntype prediction")
  return(chord_prediction_axis)
}

make_mutation_axis <- function(df_results) {
  mutation_axis <- ggplot(df_results, aes(y=1, x=order)) + 
    geom_tile(aes(fill=mutation), show.legend=TRUE) + theme_bw() + 
    theme(rect = element_rect(fill = "transparent"), plot.background = element_blank()) +
    scale_x_continuous(limits = c(0, nrow(df_results) + 1), expand = c(0, 0)) +
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          panel.grid.minor.x=element_blank(),
          panel.grid.major.x=element_blank()) +
    theme(axis.title.y = element_text(size = 8,
                                      margin = margin(t = 0, r = 10, b = 0, l = 0)),
          axis.text.y=element_blank(),
          axis.ticks.y = element_blank(),
          panel.grid.minor.y=element_blank(),
          panel.grid.major.y=element_blank()) +
    theme(legend.key.size = unit(0.03, "npc"),
          legend.title = element_blank(),
          legend.text = element_text(size =7)) + 
    scale_fill_manual(breaks = c("BRCA1", "BRCA2", "PALB2"),
                      labels = c(expression(italic("BRCA1")), expression(italic("BRCA2")),
                                 expression(italic("PALB2"))),
                      name = "Gene mutated", values = c("#F8766D", "#619CFF", "#00BA38"),
                      na.value = "white") +
    scale_y_discrete(expand = c(0, 0)) +
    facet_grid_sc(cols = vars(race), scales = list(x = scales_x), space = "free_x") +
    theme(strip.background = element_blank(),
          strip.text.x = element_blank()) +
    theme(plot.margin=unit(c(0,0,0,0.025),"npc")) + 
    ylab("Gene mutated")
}

make_stacked_barplot_sv_signatures <- function(df_results) {
  
  df_results_stackedformat_sv <- pivot_longer(df_results, cols=c(S_1, S_2, S_3, S_4, S_5, S_6, S_7),
                                              names_to = "sv_signature", values_to = "sv_signature_score")
  df_results_stackedformat_sv$sv_signature <- gsub('S_', 'Signature ',
                                                   df_results_stackedformat_sv$sv_signature)
  stacked_barplot_sv_signatures <- ggplot(df_results_stackedformat_sv,
                                          aes(fill=sv_signature, y=sv_signature_score, x=order)) + 
    geom_bar(position="fill", stat="identity", show.legend = TRUE) + theme_bw() + 
    theme(rect = element_rect(fill = "transparent"), plot.background = element_blank()) +
    scale_y_continuous(limits = c(0,1), expand = c(0, 0)) +
    theme(axis.title.x=element_blank(),
          axis.text.x = element_blank(),
          panel.grid.minor.x=element_blank(),
          panel.grid.major.x=element_blank(),
          strip.switch.pad.grid = unit(0, "cm")) +
    theme(axis.title.y = element_text(size = 8,
                                      margin = margin(t = 0, r = 10, b = 0, l = 0)),
          axis.text.y = element_text(size = 7)) +
    theme(legend.key.size = unit(0.03, "npc"),
          legend.title = element_blank(),
          legend.text = element_text(size =7)) + 
    scale_fill_manual(values = colors_signatures) +
    guides(fill=guide_legend(title = "SV signatures")) + ylab("Proportion") +
    facet_grid_sc(cols = vars(race), scales = list(x = scales_x), space = "free_x") +
    theme(strip.background = element_blank(),
          strip.text.x = element_blank()) +
    theme(plot.margin=unit(c(0,0,0,0.025),"npc")) + 
    ylab("SV signatures\nproportion")
  return(stacked_barplot_sv_signatures)
}

make_sbs_axis <- function(df_results) {
  sbs_axis <- ggplot(df_results, aes(y=1, x=order)) + 
    geom_tile(aes(fill=SBS3), show.legend=TRUE) + theme_bw() + 
    theme(rect = element_rect(fill = "transparent"), plot.background = element_blank()) +
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          panel.grid.minor.x=element_blank(),
          panel.grid.major.x=element_blank()) +
    theme(axis.title.y = element_text(size = 8,
                                      margin = margin(t = 0, r = 10, b = 0, l = 0)),
          axis.text.y=element_blank(),
          axis.ticks.y = element_blank(),
          panel.grid.minor.y=element_blank(),
          panel.grid.major.y=element_blank()) +
    theme(legend.key.size = unit(0.03, "npc"),
          legend.title = element_blank(),
          legend.text = element_text(size =7),
          legend.direction = "horizontal") + 
    scale_fill_gradient(low="white", high="#EA21F7", limits = c(0, 1), name = "SBS3 proportion") +
    scale_y_discrete(expand = c(0, 0)) +
    facet_grid_sc(cols = vars(race), scales = list(x = scales_x), space = "free_x") +
    theme(strip.background = element_blank(),
          strip.text.x = element_blank()) +
    theme(plot.margin=unit(c(0,0,0,0.025),"npc")) + 
    ylab("SBS3\nproportion")
}

make_indel_axis <- function(df_results) {
  indel_axis <- ggplot(df_results, aes(y=1, x=order)) + 
    geom_tile(aes(fill=ID), show.legend=TRUE) + theme_bw() + 
    theme(rect = element_rect(fill = "transparent"), plot.background = element_blank()) +
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          panel.grid.minor.x=element_blank(),
          panel.grid.major.x=element_blank()) +
    theme(axis.title.y = element_text(size = 8,
                                      margin = margin(t = 0, r = 10, b = 0, l = 0)),
          axis.text.y=element_blank(),
          axis.ticks.y = element_blank(),
          panel.grid.minor.y=element_blank(),
          panel.grid.major.y=element_blank()) +
    theme(legend.key.size = unit(0.03, "npc"),
          legend.title = element_blank(),
          legend.text = element_text(size =7),
          legend.direction = "horizontal") + 
    scale_fill_gradient(low="white", high="orangered3", limits = c(0, 1), name = "ID6 + ID8 proportion") +
    scale_y_discrete(expand = c(0, 0)) +
    facet_grid_sc(cols = vars(race), scales = list(x = scales_x), space = "free_x") +
    theme(strip.background = element_blank(),
          strip.text.x = element_blank()) +
    theme(plot.margin=unit(c(0,0,0,0.025),"npc")) + 
    ylab("ID6 + ID8\nproportion")
}

make_subtype_axis <- function(df_results) {
  subtype_axis <- ggplot(df_results, aes(y=1, x=order)) + 
    geom_tile(aes(fill=subtype), show.legend=TRUE) + theme_bw() + 
    theme(rect = element_rect(fill = "transparent"), plot.background = element_blank()) +
    theme(axis.title.x=element_blank(),
          axis.text.x=element_text(size = 8),
          panel.grid.minor.x=element_blank(),
          panel.grid.major.x=element_blank()) +
    theme(axis.title.y = element_text(size = 8,
                                      margin = margin(t = 0, r = 10, b = 0, l = 0)),
          axis.ticks.y = element_blank(),
          axis.text.y=element_blank(),
          panel.grid.minor.y=element_blank(),
          panel.grid.major.y=element_blank()) +
    theme(legend.key.size = unit(0.03, "npc"),
          legend.title = element_blank(),
          legend.text = element_text(size =7)) + 
    scale_fill_manual(breaks = c('HER2+', 'HR-/HER2-', 'HR+/HER2-'),
                      name = "Subtype", values = c("#F5C3CB", "#EB4CAF", "#7F170E")) +
    scale_y_discrete(expand = c(0, 0)) +
    facet_grid_sc(cols = vars(race), scales = list(x = scales_x), space = "free_x") +
    theme(strip.background = element_blank(),
          strip.text.x = element_blank()) +
    theme(plot.margin=unit(c(0,0,0,0.025),"npc")) + 
    ylab("Subtype")
  return(subtype_axis)
}

make_x_labels <- function(df_results) {
  add_counts <- c(`Nigerian` = paste(sum(total_results$race == "Nigerian"), "samples"),
                       `Black` = paste(sum(total_results$race == "Black"), "samples"),
                       `White` = paste(sum(total_results$race == "White"), "samples"))
  df_results$ignore <- "1"
  x_labels <- ggplot(df_results, aes(y=1, x=order)) + ylim(1, 1) +
    geom_tile(aes(fill=ignore), show.legend=TRUE) + theme_minimal() +
    facet_grid_sc(cols = vars(race), scales = list(x = scales_x), space = "free_x",
                  labeller = as_labeller(add_counts)) +
    scale_fill_manual(values = c("transparent")) + 
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          panel.grid.minor.x=element_blank(),
          panel.grid.major.x=element_blank()) +
    theme(axis.title.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.text.y=element_blank(),
          panel.grid.minor.y=element_blank(),
          panel.grid.major.y=element_blank()) +
    theme(legend.text = element_blank(),
          legend.title = element_blank()) +
    theme(plot.margin=unit(c(0,0,0,0.025),"npc")) +
    theme(strip.text.x = element_text(size = 9))
  return(x_labels)
}

combined_signatures <- plot_grid(
  make_chord_score_axis(total_results),
  NULL,
  make_chord_prediction_axis(total_results),
  NULL,
  make_mutation_axis(total_results),
  NULL,
  make_stacked_barplot_sv_signatures(total_results),
  NULL,
  make_sbs_axis(total_results),
  NULL,
  make_indel_axis(total_results),
  NULL,
  make_subtype_axis(total_results),
  NULL,
  make_x_labels(total_results),
  ncol = 1, align = "hv", rel_heights = c(2, -0.65, 2, -0.65, 2, -0.75,
                                          4, -0.65, 2, -0.65, 2, -0.65, 2, 0, 0.5),
  labels = c('A', '', 'B', '', 'C', '', 'D', '', 'E', '', 'F', '', 'G', '', ''),
  label_size = 12,
  vjust = 3.5
)

ggsave("HRD_analysis_figure_JB.pdf", combined_signatures, dpi = 300, width = 10, height = 7.25)

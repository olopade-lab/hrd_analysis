library(ggplot2)
library(readxl)
library(tidyr)
library(grid)
library(cowplot)
library(dplyr)
library(knitr)
library(fmsb)
library(ggpubr)
library(tibble)
library(rstatix)
library(scales)

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
total_results <- merge(chord_pred, sample_info, all.x = TRUE)
total_results <- merge(total_results, mutation_info, all.x = TRUE)
total_results <- merge(total_results, sv_counts, all.x = TRUE)
df_colors <- data.frame(c("#F8766D", "#619CFF", "#00BA38"),
                        c("Nigerian", "Black", "White"))
colnames(df_colors) <- c("color", "race")

dotplot_hrd <- function(df_results, subtype) {
  df_results <- df_results[complete.cases(df_results[ , c("race", "subtype")]),]
  df_results <- df_results[order(df_results$p_hrd),]
  df_results$race <- factor(df_results$race, levels=c("Nigerian", "Black", "White"))
  df_results_subtype <- df_results[df_results$subtype == subtype, ]
  df_frequency <- as.data.frame(table(df_results_subtype$race))
  colnames(df_frequency) <- c("race", "count")
  df_frequency <- df_frequency[df_frequency$count > 0,]
  colors_to_use <- df_colors[df_colors$race %in% df_frequency$race,"color"]
  dotplot_HRD <- ggplot(df_results_subtype, aes(x=race, y=p_hrd, fill=race))  + theme_bw() +
    geom_dotplot(binaxis= "y", stackdir = "center", dotsize = 0.52,
                 stroke = 0.3, position = "dodge", binpositions="all", method = "histodot") +
    geom_hline(yintercept=0.5, linetype="dashed", color = "#CC0000") + 
    ylab("HRD probability score") + xlab("Group") + scale_y_continuous(limits = c(0,1)) +
    scale_x_discrete(labels=paste0(df_frequency$race, "\n(N=", df_frequency$count_total, ")")) +
    scale_fill_manual(values = colors_to_use) +
    theme(plot.title = element_text(hjust = 0.5)) + ggtitle(subtype)
  if (subtype == "HR-/HER2-") {
    dotplot_HRD <- dotplot_HRD + annotate("text", x = 3.05, y = 0.45,
                                          label = "HRD threshold",
                                          size = 3.5, color = "red", hjust = 0.5)
  }
  return(dotplot_HRD)
}

barplot_hrd <- function(df_results, subtype) {
  df_results <- df_results[complete.cases(df_results[ , c("race", "subtype")]),]
  df_results <- df_results[order(df_results$p_hrd),]
  df_results$race <- factor(df_results$race, levels=c("Nigerian", "Black", "White"))
  df_results_subtype <- df_results[df_results$subtype == subtype, ]
  df_frequency <- as.data.frame(table(df_results_subtype$race))
  colnames(df_frequency) <- c("race", "count_total")
  df_frequency <- df_frequency[df_frequency$count > 0,]
  colors_to_use <- df_colors[df_colors$race %in% df_frequency$race, "color"]
  df_count_hrd <- count(df_results_subtype, race, hr_status)
  df_frequency$count_hrd <- df_count_hrd[df_count_hrd$hr_status == "HR_deficient", "n"]
  df_frequency$prop_hrd <- df_frequency$count_hrd / df_frequency$count_total
  df_frequency$se <- sqrt(df_frequency$prop_hrd*(1-df_frequency$prop_hrd)/df_frequency$count_total)
  
  combinations_race <- t(combn(df_frequency$race, 2))
  df_stat_test <- data.frame(group1=character(),
                   group2=character(), 
                   p.adjust=character()) 
  pairwise_fisher <- pairwise.fisher.test(df_frequency$count_hrd, df_frequency$count_total)
  p_values_fisher <- as.list(pairwise_fisher$p.value)
  p_values_fisher <- p_values_fisher[!is.na(p_values_fisher)]
  
  for (index in 1:nrow(combinations_race)) {
    if (p_values_fisher[index] < 0.2) {
      if (p_values_fisher[index] < 0.05) {
        df_stat_test <- add_row(df_stat_test, group1 = combinations_race[index, 1],
                             group2 = combinations_race[index, 2],
                             p.adjust = character(round(p_values_fisher[[index]], digits = 3))
        )
      } else {
        
        df_stat_test <- add_row(df_stat_test, group1 = combinations_race[index, 1],
                             group2 = combinations_race[index, 2],
                             p.adjust = paste(as.character(round(p_values_fisher[[index]], digits = 3)), "(n.s.)")
        )
      }
    }
  }
  
  barplot_HRD <- ggplot(df_frequency, aes(x = race, y = prop_hrd, fill = race)) + 
    geom_bar(stat="identity") +
    scale_fill_manual(values = colors_to_use) + theme_bw() + 
    ylab("Proportion of HRD samples") + xlab("Group") + 
    scale_x_discrete(labels=paste0(df_frequency$race, "\n(N=", df_frequency$count_total, ")")) +
    geom_errorbar(aes(ymin=prop_hrd-se, ymax=prop_hrd+se), width=0.2, size = 0.3) +
    scale_y_continuous(limits = c(0,1))
  if (nrow(df_stat_test) > 0) {
    barplot_HRD <- barplot_HRD + stat_pvalue_manual(df_stat_test,
                                                    y.position = max(df_frequency$prop_hrd+df_frequency$se) + 0.1,
                                                    step.increase = 0.1, label = "p.adjust",
                                                    inherit.aes = FALSE, size = 3.5)
  }
  return(barplot_HRD)
}

combined_dotplot_plots <- plot_grid(
  dotplot_hrd(total_results, "HR+/HER2-") +
    theme(legend.position="none", axis.text.x=element_blank()) + xlab(NULL),
  dotplot_hrd(total_results, "HER2+") +
    theme(legend.position="none", axis.text.x=element_blank(), axis.text.y=element_blank()) +
    ylab(NULL) + xlab(NULL),
  dotplot_hrd(total_results, "HR-/HER2-") +
    theme(legend.position="none", axis.text.x=element_blank(), axis.text.y=element_blank()) +
    ylab(NULL) + xlab(NULL),
  nrow = 1, align='h', axis = 'l'
)

combined_barplot_plots <- plot_grid(
  barplot_hrd(total_results, "HR+/HER2-") + theme(legend.position="none") + xlab(NULL),
  barplot_hrd(total_results, "HER2+") + theme(legend.position="none", axis.text.y=element_blank()) +
    ylab(NULL),
  barplot_hrd(total_results, "HR-/HER2-") + theme(legend.position="none", axis.text.y=element_blank())+
    ylab(NULL) + xlab(NULL),
  nrow = 1, align = 'h', axis = 'l'
)

combined_all <- plot_grid(
  combined_dotplot_plots,
  combined_barplot_plots,
  labels = c("A", "B"),
  ncol = 1, align = 'v'
)

grid.draw(combined_all)

ggsave("supplementary_hrd_JB.pdf", combined_all, width = 10, height = 7.25)

library(ggplot2)
library(readxl)
library(tidyr)
library(grid)
library(cowplot)
library(dplyr)

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
total_results <- merge(chord_pred, sample_info, all.x = TRUE)
total_results <- merge(total_results, mutation_info, all.x = TRUE)
total_results <- merge(total_results, sv_counts, all.x = TRUE)
df_colors <- data.frame(c("#999999", "#DCA037", "#6CB2E3"),
                        c("Nigerian", "Black", "White"))
colnames(df_colors) <- c("color", "race")


# Stacked barplot
hrd_genotype_graph <- function(df_results, group) {
  
  df_results_stackedformat <- pivot_longer(df_results, cols=c(p_BRCA1, p_BRCA2),
                                              names_to = "p_mutation", values_to = "p_mutation_score")
  df_results_stackedformat$p_mutation <- gsub('p_', '', df_results_stackedformat$p_mutation)
  df_results_stackedformat_nigerian <- df_results_stackedformat[df_results_stackedformat$group == group,]
  df_results_stackedformat_nigerian$sample <- with(df_results_stackedformat_nigerian, reorder(sample, p_hrd))
  stacked_barplot <- ggplot(df_results_stackedformat_nigerian, aes(fill=p_mutation, y=p_mutation_score, x=sample)) + 
    geom_bar(position="stack", stat="identity", show.legend = TRUE) + theme_minimal() + scale_y_continuous(limits = c(0,1), expand = c(0, 0)) +
    theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
    scale_fill_manual(values = c("#F8766D", "#619CFF")) + ylab("HRD Signature Score")
  
  mutation_axis <- ggplot(df_results_stackedformat_nigerian, aes(y=1, x=sample)) + 
    geom_tile(aes(fill=mutation), show.legend=TRUE) + 
    theme(axis.title.x=element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
    theme(axis.text.y = element_blank(), axis.ticks.y = element_blank()) +
    scale_fill_manual(na.value = rgb(241/255, 241/255, 241/255, 1), breaks = c('BRCA1', 'BRCA2', 'PALB2'),
                      name = "Mutation", values = c("#F8766D", "#619CFF", "#00BA38")) + ylab("Genotype")
  
  legend_mutation_axis <- get_legend(mutation_axis + 
                                       guides(color = guide_legend(nrow = 1)) +
                                       theme(legend.position = "bottom", legend.title=element_blank()))
  
  combined_plots <- plot_grid(
    stacked_barplot + theme(legend.position="none"),
    mutation_axis + theme(legend.position="none"),
    ncol = 1, align='v', rel_heights = c(10, 1)
  )
  
  combined_plots_withlegend <- plot_grid(
    combined_plots, legend_mutation_axis, ncol = 1, rel_heights = c(8, 1)
  )
  
  final_plot <- grid.draw(combined_plots_withlegend)
  return(final_plot)
}

# Testing continuous

df_results <- total_results
group <- "Nigerian"
df_results_group <- df_results[df_results$group == group,]
df_results_group <- df_results_group[order(df_results_group$p_hrd), ]
df_results_group$order <- 1:nrow(df_results_group)
# df_results_group$sample <- with(df_results_group, reorder(sample, p_hrd))
df_results_group_stackedformat <- pivot_longer(df_results_group, cols=c(p_BRCA1, p_BRCA2),
                                         names_to = "p_mutation", values_to = "p_mutation_score")
df_results_group_stackedformat$p_mutation <- gsub('p_BRCA1', 'BRCA1-type HRD', df_results_group_stackedformat$p_mutation)
df_results_group_stackedformat$p_mutation <- gsub('p_BRCA2', 'BRCA2-type HRD', df_results_group_stackedformat$p_mutation)

ggplot(df_results_group_stackedformat, aes(fill=p_mutation, y=p_mutation_score, x=order)) + 
  geom_bar(position="stack", stat="identity", show.legend = TRUE) + theme_bw() + 
  scale_y_continuous(limits = c(0,1), expand = c(0, 0.03)) +
  theme(axis.title.x=element_blank(),
        panel.grid.minor.x=element_blank(),
        panel.grid.major.x=element_blank()) +
  scale_fill_manual(values = c("#F8766D", "#619CFF")) + ylab("HRD Signature Score")






hrd_genotype_graph(total_results, "Nigerian")
hrd_genotype_graph(total_results, "TCGA")

table_chord <- function(df_results, group) {
  df_results <- df_results[df_results$group == group, c("sample", "subtype", "mutation",
                                                        "p_hrd", "hr_status", "hrd_type")]
  df_results$hr_status <- gsub('HR_proficient', 'Proficient', 
                               gsub('HR_deficient', 'Deficient', df_results$hr_status))
  df_results$hrd_type <- gsub('cannot_be_determined', 'Cannot be determined',
                              gsub('_type', '', df_results$hrd_type))
  df_results$hrd_type <- gsub("none", "", df_results$hrd_type)
  df_results <- df_results[order(df_results$p_hrd, decreasing = TRUE), ]
  df_results[is.na(df_results)] <- ""
  columns <- c("Sample", "Subtype", "Mutation", "CHORD HRD Score",
               "CHORD HR Status", "CHORD HRD Type")
  colnames(df_results) <- columns
  kable(df_results, "simple")
}

table_chord(total_results, "Nigerian")

df_colors <- data.frame(c("#999999", "#DCA037", "#6CB2E3"),
                        c("Nigerian", "Black", "White"))
colnames(df_colors) <- c("color", "race")
                        

dotplot_hrd <- function(df_results, subtype) {
  df_results <- df_results[complete.cases(df_results[ , c("race", "subtype")]),]
  df_results <- df_results[order(df_results$p_hrd),]
  df_results$race <- factor(df_results$race, levels=c("Nigerian", "Black", "White"))
  df_results_subtype <- df_results[df_results$subtype == subtype, ]
  df_frequency <- as.data.frame(table(df_results_subtype$race))
  colnames(df_frequency) <- c("race", "frequency")
  df_frequency <- df_frequency[df_frequency$frequency > 0,]
  colors_to_use <- df_colors[df_colors$race %in% df_frequency$race,"color"]
  print(df_results_subtype)
  # df_results_subtype <- df_results_subtype[df_results_subtype$p_hrd > 0,]
  plot_HRD <- ggplot(df_results_subtype, aes(x=race, y=p_hrd, fill=race))  + theme_minimal() +
    geom_dotplot(binaxis= "y", stackdir = "center", dotsize = 0.5, stroke = 0.3, position = "dodge", binpositions="all", method = "histodot") +
    ylab("HRD probability") + xlab("Number of samples") + scale_y_continuous(limits = c(0,1)) +
    scale_x_discrete(labels=df_frequency$frequency) + scale_fill_manual(values = colors_to_use) +
    theme(plot.title = element_text(hjust = 0.5)) + ggtitle(subtype)
  return(plot_HRD)
}

combined_dotplot_plots <- plot_grid(
  dotplot_hrd(total_results, "HR+/HER2-") + theme(legend.position="none") + xlab(NULL),
  dotplot_hrd(total_results, "HER2+") + theme(legend.position="none", axis.text.y=element_blank()) + ylab(NULL),
  dotplot_hrd(total_results, "HR-/HER2-") + theme(legend.position="none", axis.text.y=element_blank())+ ylab(NULL) + xlab(NULL),
  nrow = 1, align='h', axis = 'l', ncol = 3
)

grid.draw(combined_dotplot_plots)

barplot_hrd <- function(df_results, subtype) {
  df_results <- df_results[complete.cases(df_results[ , c("race", "subtype")]),]
  df_results <- df_results[order(df_results$p_hrd),]
  df_results$race <- factor(df_results$race, levels=c("Nigerian", "Black", "White"))
  df_results_subtype <- df_results[df_results$subtype == subtype, ]
  df_frequency <- as.data.frame(table(df_results_subtype$race))
  colnames(df_frequency) <- c("race", "count")
  df_frequency <- df_frequency[df_frequency$count > 0,]
  colors_to_use <- df_colors[df_colors$race %in% df_frequency$race,"color"]
  df_count_hrd <- count(df_results_subtype, race, hr_status)
  df_frequency$prop_hrd <- df_count_hrd[df_count_hrd$hr_status == "HR_deficient", "n"] / df_frequency$count
  barplot_HRD <- ggplot(df_frequency, aes(x = race, y = prop_hrd, fill = race)) + geom_bar(stat="identity") + 
    scale_fill_manual(values = colors_to_use) + theme_bw() + 
    ylab("Proportion HRD") + xlab("Number of samples") + scale_y_continuous(limits = c(0,1)) +
    scale_x_discrete(labels=df_frequency$count)
  return(barplot_HRD)
}

combined_barplot_plots <- plot_grid(
  barplot_hrd(total_results, "HR+/HER2-") + theme(legend.position="none") + xlab(NULL),
  barplot_hrd(total_results, "HER2+") + theme(legend.position="none", axis.text.y=element_blank()) + ylab(NULL),
  barplot_hrd(total_results, "HR-/HER2-") + theme(legend.position="none", axis.text.y=element_blank())+ ylab(NULL) + xlab(NULL),
  nrow = 1, align='h', axis = 'l', ncol = 3
)

x_axis_barplot <- get_x_axis(barplot_hrd(total_results, "HR+/HER2-"))
ggdraw(x_axis_barplot)
grid.draw(combined_barplot_plots)

# Testing boxplot

df_results <- total_results
subtype <- "HER2+"
df_results <- df_results[complete.cases(df_results[ , c("race", "subtype")]),]
df_results <- df_results[order(df_results$p_hrd),]
df_results$race <- factor(df_results$race, levels=c("Nigerian", "Black", "White"))
df_results_subtype <- df_results[df_results$subtype == subtype, ]
df_frequency <- as.data.frame(table(df_results_subtype$race))
colnames(df_frequency) <- c("race", "count")
df_frequency <- df_frequency[df_frequency$count > 0,]
colors_to_use <- df_colors[df_colors$race %in% df_frequency$race,"color"]
df_count_hrd <- count(df_results_subtype, race, hr_status)
df_frequency$prop_hrd <- df_count_hrd[df_count_hrd$hr_status == "HR_deficient", "n"] / df_frequency$count
ggplot(df_frequency, aes(x = race, y = prop_hrd, fill = race)) + geom_bar(stat="identity") + 
  scale_fill_manual(values = colors_to_use) + theme_bw() + 
  ylab("Proportion HRD") + xlab("Number of samples") + scale_y_continuous(limits = c(0,1)) +
  scale_x_discrete(labels=df_frequency$count)

df_frequency
df_results_subtype
count(df_results_subtype, race, hr_status)[]
colSums()



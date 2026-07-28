# Script to calculate synergy scores and generate finalized figures for the manuscript
# This script generates Figure 3 A,B,C and Supplementary Figure 4 A & B.

PROJECT.DIR <- Sys.getenv("MANUSCRIPT_FIGURE_CODE_DIR", unset = getwd())

# Keep your original workflow (everything relative to project dir)
setwd(PROJECT.DIR)

DATA.DIR    <- "data"
PLOTS.DIR   <- "plots"
SCRIPTS.DIR <- "scripts"

if (!dir.exists(PLOTS.DIR)) dir.create(PLOTS.DIR, recursive = TRUE)

library(synergyfinder)
library(tidyverse)
library(patchwork)
library(RColorBrewer)


# CALCULATE SYNERGY SCORES ----------------------
# ...RUN synergy Finder (each replicate separately) ------------------------------------
load(file.path(DATA.DIR, "all_cellLines_drug_combo_wrangle_replicate_DMSO_normalized_batch1.RData"))
batch1 <- obj.list
load(file.path(DATA.DIR, "2024-01-08_all_cellLines_drug_combo_wrangle_replicate_DMSO_normalized_batch2.RData"))
batch2 <- obj.list

# batch 1
# Use lapply to convert elements to data frames
lapply(names(batch1), function(x) assign(x, batch1[[x]], envir = .GlobalEnv))
# renaming objects from batch 1 (i.e. replicate 1)
# disregard IMR05 and LAN5
rm(list=ls(pattern="^IMR05_"))
rm(list=ls(pattern="^LAN5_"))

# discard some more combinations (2023-12-27)
rm(list=ls(pattern="^KELLY_cyclophosphamide"))
rm(list=ls(pattern="^SKNAS_cyclophosphamide"))
rm(list=ls(pattern="^SKNSH_cyclophosphamide"))
rm(list=ls(pattern="^SKNSH_irinotecan"))

# duplicate SKNSH-ETOPO
SKNSH_etoposide_1_duplicated <- SKNSH_etoposide
SKNSH_etoposide_2_duplicated <- SKNSH_etoposide

all_obj <- ls(pattern = '_')
for (old_name in all_obj) {
  new_name <- paste0(old_name, '_batch1')
  if (new_name != old_name) {
    assign(new_name, get(old_name))
    rm(old_name)
  }
}

# batch 2
lapply(names(batch2), function(x) assign(x, batch2[[x]], envir = .GlobalEnv))

# running for batch1_2
all_obj <- ls(pattern = '_batch[1,2]')
# 90 = 6 cl * 5 chemo * 3 reps

# to store all synergy scores in one data frame
synergy.scores.all <- data.frame()


# batch 2
lapply(names(batch2), function(x) assign(x, batch2[[x]], envir = .GlobalEnv))

# running for batch1_2
all_obj <- ls(pattern = '_batch[1,2]')
# 90 = 6 cl * 5 chemo * 3 reps

# to store all synergy scores in one data frame
synergy.scores.all <- data.frame()

# ....plotting function ---------------

generateViz <- function(obj, title){
  
  res <- obj
  
  # plot heatmap for synergyFinder
  p1 <- Plot2DrugHeatmap(
    data = res,
    plot_block = 1,
    drugs = c(1, 2),
    plot_value = "ZIP_synergy",
    dynamic = FALSE,
    summary_statistic = c("mean", "median")) +
    ggtitle(title)
  
  ggsave(p1,
         filename = file.path(PLOTS.DIR, paste0(Sys.Date(), "_", title, "_heatmap.pdf")),
         width = 10, height = 8)
  
}


# run synergy finder and create visualizations
for(i in 1:length(all_obj)){
  
  print(paste0(i,': ',all_obj[i]))
  
  
  # grepping cell line name
  if(stringr::str_detect(pattern = fixed('imr05', ignore_case = TRUE), all_obj[i])){
    cl <- 'IMR05'
    
  } else if(stringr::str_detect(pattern = fixed('imr5', ignore_case = TRUE), all_obj[i])){
    cl <- 'IMR05'
    
  } else if(stringr::str_detect(pattern = fixed('lan5', ignore_case = TRUE), all_obj[i])){
    cl <- 'LAN5'
    
  } else if(stringr::str_detect(pattern = fixed('shep', ignore_case = TRUE), all_obj[i])){
    cl <- 'SHEP'
    
  } else if(stringr::str_detect(pattern = fixed('sknas', ignore_case = TRUE), all_obj[i])){
    cl <- 'SKNAS'
    
  } else if(stringr::str_detect(pattern = fixed('kelly', ignore_case = TRUE), all_obj[i])){
    cl <- 'KELLY'
    
  } else if(stringr::str_detect(pattern = fixed('sknsh', ignore_case = TRUE), all_obj[i])){
    cl <- 'SKNSH'
    
  }
  
  # grepping chemo drug name
  if(stringr::str_detect(pattern = fixed('cyclo', ignore_case = TRUE), all_obj[i])){
    drug <- 'Cyclophosphamide'
    
  } else if(stringr::str_detect(pattern = fixed('etop', ignore_case = TRUE), all_obj[i])){
    drug <- 'Etoposide'
    
  } else if(stringr::str_detect(pattern = fixed('topo', ignore_case = TRUE), all_obj[i])){
    drug <- 'Topotecan'
    
  } else if(stringr::str_detect(pattern = fixed('irino', ignore_case = TRUE), all_obj[i])){
    drug <- 'Irinotecan'
    
  } else if(stringr::str_detect(pattern = fixed('temo', ignore_case = TRUE), all_obj[i])){
    drug <- 'Temozolomide'
    
  } else if(stringr::str_detect(pattern = fixed('tmz', ignore_case = TRUE), all_obj[i])){
    drug <- 'Temozolomide'
    
  }
  
  
  # grepping replicate number
  if(grepl('_1_|_1\\s', all_obj[i])){
    if(grepl('duplicated', all_obj[i])){
      replicate <- 3
    } else {
      replicate <- 1
    }
  } else if(grepl('_2_|_2xlsx|_2\\s', all_obj[i])) {
    if(grepl('duplicated', all_obj[i])){
      replicate <- 3
    } else {
      replicate <- 2
    }
  } else {
    replicate <- 1
  }
  
  # grepping batch number i.e. replicates
  if(grepl('_batch1', all_obj[i])){
    batch <- 'batch1'
  } else if(grepl('_batch2', all_obj[i])) {
    batch <- 'batch2'
  } else {
    batch <-  'batch1'
  }
  
  df_name <- paste0(drug,'_',cl,'_',batch)
  title <- df_name
  
  # get the data from env
  df <- get(all_obj[i])
  
  
  
  # reshapeData using synergyFinder function
  res <- ReshapeData(
    data = df,
    data_type = "viability",
    impute = TRUE,
    impute_method = NULL,
    noise = TRUE,
    seed = 1)
  
  
  # calculate synergy scores
  res <- CalculateSynergy(
    data = res,
    method = c("ZIP", "HSA", "Bliss", "Loewe"),
    correct_baseline = "non")
  
  
  # call function to generate plots
  if(title %in% c('Topotecan_LAN5_batch2','Topotecan_SHEP_batch2','Topotecan_KELLY_batch1','Topotecan_IMR05_batch2','Topotecan_SKNAS_batch2','Topotecan_SKNSH_batch2','Cyclophosphamide_LAN5_batch2','Cyclophosphamide_SHEP_batch2','Cyclophosphamide_KELLY_batch2','Cyclophosphamide_IMR05_batch2','Cyclophosphamide_SKNAS_batch1','Cyclophosphamide_SKNSH_batch2','Irinotecan_LAN5_batch2','Irinotecan_SHEP_batch2','Irinotecan_KELLY_batch2','Irinotecan_IMR05_batch2','Irinotecan_SKNAS_batch2','Irinotecan_SKNSH_batch1','Etoposide_LAN5_batch2','Etoposide_SHEP_batch2','Etoposide_KELLY_batch2','Etoposide_IMR05_batch2','Etoposide_SKNAS_batch2','Etoposide_SKNSH_batch1','Temozolomide_LAN5_batch2','Temozolomide_SHEP_batch2','Temozolomide_KELLY_batch2','Temozolomide_IMR05_batch2','Temozolomide_SKNAS_batch2','Temozolomide_SKNSH_batch2')){
    generateViz(obj = res, title = title)
  }
  
  # add cell line names
  res$synergy_scores$cellLine <- cl
  
  # add df_name
  res$synergy_scores$filename <- all_obj[i]
  
  # add drug name
  res$synergy_scores$combination <- df_name
  
  assign(df_name, res)
  
  # 2023-12-26 to generate final synergy plots using all replicates (from repeated replicates experiment as well)
  synergy.scores.all <- rbind(synergy.scores.all, res$synergy_scores)
  
  
}

length(unique(synergy.scores.all$filename))
# 90

length(unique(synergy.scores.all$cellLine))
# 6


length(unique(synergy.scores.all$combination))
# 45


backup <- synergy.scores.all

synergy.scores.all <- synergy.scores.all |>
  group_by(combination) |>
  mutate(replicate = paste0('replicate_',dense_rank(filename))) |>
  mutate(combination = paste0(combination, '_', replicate))


# VISUALIZATIONS -------------------

# added synergistic concenrtations in excel
synergistic.combos <- read.delim(file.path(DATA.DIR, "final_synergistic_combinations_2.txt"), header = T)
synergistic.combos$conc1 <- as.numeric(synergistic.combos$conc1)
synergistic.combos$conc2 <- as.numeric(synergistic.combos$conc2)


head(synergy.scores.all)

# filter synergistic combinations
most.synergistic <- synergy.scores.all |>
  mutate(combination_name = gsub('_batch.*_replicate.*','', combination)) |>
  inner_join(synergistic.combos, by = c('combination_name' = 'combination', 'conc1', 'conc2')) |>
  select(2:3,6,15,17) |>
  distinct()

# 360 rows = 3 reps x 4 conc x 5 chemo x 6 cls


# plot - one representative replicate barplot for each chemo:DT:cl combo

# to determine which replicate out of the 3 to plot

tmp <- most.synergistic |>
  mutate(synergy_combo = paste(conc1,':',conc2)) |>
  select(-conc1, -conc2) |>
  group_by(combination) |>
  mutate(sum_ZIP = sum(ZIP_synergy))


tmp <- tmp |>
  group_by(combination_name) |>
  slice(which.max(sum_ZIP))


tmp[duplicated(tmp$combination),]


# set ordering of cell lines
most.synergistic$combination_name2 <- most.synergistic$combination_name
most.synergistic <- separate(most.synergistic, col = 'combination_name2', into = c('chemo','cl'), sep = '_')
most.synergistic$chemo <- factor(most.synergistic$chemo, levels = c('Topotecan', 'Cyclophosphamide', 'Irinotecan', 'Temozolomide', 'Etoposide'))



specific.order.vec <- most.synergistic |>
  arrange(chemo, cl) |>
  ungroup() |>
  select(combination_name) |>
  distinct()


# Reordering the factor levels within combination_name based on a specific order
most.synergistic$combination_name <- fct_relevel(most.synergistic$combination_name,
                                                 specific.order.vec$combination_name)




#plotted.replicates <-  # add sd values to these replicates

# calculate sd for error bars
sd.zip <- most.synergistic |>
  mutate(synergy_combo = paste(conc1,':',conc2)) |>
  #filter(grepl('Temozolomide', combination)) |>
  #filter(grepl('KELLY', combination)) |>
  mutate(ZIP_synergy = as.numeric(ZIP_synergy)) |>
  group_by(combination_name, synergy_combo) |>
  mutate(sd_ZIP_synergy_combo = sd(ZIP_synergy),
         mean_ZIP_synergy_combo = mean(ZIP_synergy)) |>
  distinct()

sd.zip <- sd.zip[,-c(1:3)]
sd.zip <- sd.zip |>
  distinct()

# individual plots put together in a grid
most.synergistic$cl <- factor(most.synergistic$cl, levels = c("IMR05", "LAN5", "SKNSH","SHEP","KELLY", "SKNAS"))

col_vec <- c("IMR05" = "brown", "LAN5" = "cornflowerblue", "SKNSH" = "darkorange3", "SHEP" = "lightgoldenrod3", "KELLY" = "black","SKNAS" = "aquamarine4")

for(i in unique(most.synergistic$chemo)){
  for(j in unique(most.synergistic$cl)){
    print(paste0(i,':',j))
    
    plt_obj <- paste0(i,'_',j,'_plot')
    
    assign(plt_obj, most.synergistic |>
             mutate(synergy_combo = paste(conc1,':',conc2)) |>
             filter(combination %in% tmp$combination) |>
             left_join(sd.zip, by = c('combination', 'synergy_combo', 'combination_name', 'cl', 'chemo')) |>
             filter(chemo == i) |>
             filter(cl == j) |>
             ggplot(aes(synergy_combo, ZIP_synergy)) +
             geom_bar(stat = 'identity', fill = col_vec[j]) +
             geom_errorbar(aes(ymin=ifelse(ZIP_synergy-sd_ZIP_synergy_combo < 0, 0, ZIP_synergy-sd_ZIP_synergy_combo), ymax=ifelse(ZIP_synergy+sd_ZIP_synergy_combo > 50, 50, ZIP_synergy+sd_ZIP_synergy_combo)), width=.2,
                           position=position_dodge(.9)) +
             labs(x = '', y = 'ZIP Score', title = j) +
             theme_minimal() +
             geom_hline(yintercept = 10, linetype="dotted",
                        color = "red", size=1, alpha = 0.7) +
             ylim(0,51) +
             theme(plot.title = element_text(hjust = 0.5)))
    
  }
}


# put it all together in a grid plot
all_plots <- grep('_plot', ls(), value = TRUE)
topo_plots <- ls()[grepl("Topotecan", ls()) & grepl("_plot", ls())]
most.synergistic |>
  filter(chemo == 'Topotecan') |>
  arrange(chemo, cl) |>
  select(combination_name) |>
  distinct() |>
  pull(combination_name) -> topo.vec

topo.vec <- paste0(topo.vec, '_plot')
#order by topo.vec
topo_plots <- topo_plots[order(match(topo_plots, topo.vec))]



cyclo_plots <- ls()[grepl("Cyclophosphamide", ls()) & grepl("_plot", ls())]
most.synergistic |>
  filter(chemo == 'Cyclophosphamide') |>
  arrange(chemo, cl) |>
  select(combination_name) |>
  distinct() |>
  pull(combination_name) -> cyclo.vec

cyclo.vec <- paste0(cyclo.vec, '_plot')
#order by cyclo.vec
cyclo_plots <- cyclo_plots[order(match(cyclo_plots, cyclo.vec))]



irino_plots <- ls()[grepl("Irinotecan", ls()) & grepl("_plot", ls())]
most.synergistic |>
  filter(chemo == 'Irinotecan') |>
  arrange(chemo, cl) |>
  select(combination_name) |>
  distinct() |>
  pull(combination_name) -> irino.vec

irino.vec <- paste0(irino.vec, '_plot')
#order by irino.vec
irino_plots <- irino_plots[order(match(irino_plots,irino.vec))]


etop_plots <- ls()[grepl("Etoposide", ls()) & grepl("_plot", ls())]
most.synergistic |>
  filter(chemo == 'Etoposide') |>
  arrange(chemo, cl) |>
  select(combination_name) |>
  distinct() |>
  pull(combination_name) -> etop.vec

etop.vec <- paste0(etop.vec, '_plot')
#order by cyclo.vec
etop_plots <- etop_plots[order(match(etop_plots,etop.vec))]


tmz_plots <- ls()[grepl("Temozolomide", ls()) & grepl("_plot", ls())]
most.synergistic |>
  filter(chemo == 'Temozolomide') |>
  arrange(chemo, cl) |>
  select(combination_name) |>
  distinct() |>
  pull(combination_name) -> tmz.vec

tmz.vec <- paste0(tmz.vec, '_plot')
#order by tmz.vec
tmz_plots <- tmz_plots[order(match(tmz_plots,tmz.vec))]



# Create a list of plots using get() to retrieve the plot objects by name
topo_plots_to_arrange <- lapply(topo_plots, get)
cyclo_plots_to_arrange <- lapply(cyclo_plots, get)
irino_plots_to_arrange <- lapply(irino_plots, get)
etop_plots_to_arrange <- lapply(etop_plots, get)
tmz_plots_to_arrange <- lapply(tmz_plots, get)

# Arrange plots side by side using patchwork
a1 <- wrap_plots(topo_plots_to_arrange, nrow = 1, ncol = length(topo_plots) + 1) + grid::textGrob('Topotecan:DT (nM)')
a2 <- wrap_plots(cyclo_plots_to_arrange, nrow = 1, ncol = length(topo_plots) + 1) + grid::textGrob('Cyclophosphamide:DT (uM)')
a3 <- wrap_plots(irino_plots_to_arrange, nrow = 1, ncol = length(topo_plots) + 1) + grid::textGrob('Irinotecan:DT (nM)')
a4 <- wrap_plots(tmz_plots_to_arrange, nrow = 1, ncol = length(topo_plots) + 1) + grid::textGrob('Temozolomide:DT (uM)')
a5 <- wrap_plots(etop_plots_to_arrange, nrow = 1, ncol = length(topo_plots) + 1) + grid::textGrob('Etoposide:DT (uM)')


grid.plot <- a1/a2/a3/a4/a5
grid.plot

ggsave(grid.plot,
       filename = file.path(PLOTS.DIR, paste0(Sys.Date(), "_SynergyScores_gridPlot_all_cl_chemo_most_synergistic_combo.pdf")),
       width = 23, height = 10)

# to add error bars to the mean ZIP score bar plot
# get mean sd across all 4 synergistic combinations and across all 3 replicates
sd.zip.mean <- sd.zip |>
  group_by(combination_name) |>
  mutate(mean_SD = mean(sd_ZIP_synergy_combo)) |>
  select(2,8) |>
  distinct()


# this plot type separates topo-irino-etopo and cyclo-temo into separate plots
plot.type.1a <- most.synergistic |>
  filter(chemo %in% c('Topotecan', 'Irinotecan', 'Etoposide')) |>
  group_by(combination_name) |>
  mutate(mean_ZIP_all_reps_combos = mean(ZIP_synergy)) |>
  select(5:8) |>
  distinct() |>
  left_join(sd.zip.mean, by = 'combination_name') |>
  ggplot(aes(chemo, mean_ZIP_all_reps_combos, fill = cl)) +
  geom_col(position = 'dodge', width = 0.5) +
  geom_errorbar(aes(ymin=ifelse(mean_ZIP_all_reps_combos-mean_SD < 0, 0, mean_ZIP_all_reps_combos-mean_SD), ymax=ifelse(mean_ZIP_all_reps_combos+mean_SD > 30, 30, mean_ZIP_all_reps_combos+mean_SD)), width=.2,
                position=position_dodge(0.5)) +
  theme_linedraw() +
  labs(x = '', y = 'Mean ZIP Score across all replicates', fill = 'Cell Line', title = "Topoisomerase Inhibitors")  +
  scale_fill_manual(values = c("brown", "cornflowerblue", "darkorange3","lightgoldenrod3", "black","aquamarine4")) +
  geom_hline(yintercept = 10, linetype="dotted",
             color = "red", size=0.5, alpha = 1) +
  ylim(0,30) +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 15),
        plot.title = element_text(hjust = 0.5, size = 15),
        legend.position = 'none')


plot.type.1b <- most.synergistic |>
  filter(chemo %in% c('Cyclophosphamide', 'Temozolomide')) |>
  group_by(combination_name) |>
  mutate(mean_ZIP_all_reps_combos = mean(ZIP_synergy)) |>
  select(5:8) |>
  distinct() |>
  left_join(sd.zip.mean, by = 'combination_name') |>
  ggplot(aes(chemo, mean_ZIP_all_reps_combos, fill = cl)) +
  geom_col(position = 'dodge', width = 0.5) +
  geom_errorbar(aes(ymin=ifelse(mean_ZIP_all_reps_combos-mean_SD < 0, 0, mean_ZIP_all_reps_combos-mean_SD), ymax=ifelse(mean_ZIP_all_reps_combos+mean_SD > 30, 30, mean_ZIP_all_reps_combos+mean_SD)), width=.2,
                position=position_dodge(0.5)) +
  theme_linedraw() +
  labs(x = '', y = '', fill = 'Cell line', caption = "NOTE: Synergistic concentrations vary between the chemotherapeutics", title = "Alkylating Agents")  +
  scale_fill_manual(values = c("brown", "cornflowerblue", "darkorange3","lightgoldenrod3", "black","aquamarine4")) +
  geom_hline(yintercept = 10, linetype="dotted",
             color = "red", size=0.5, alpha = 1) +
  ylim(0,30) +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 15),
        plot.title = element_text(hjust = 0.5, size = 15))


plot.type.1 <- plot.type.1a|plot.type.1b

ggsave(plot.type.1,
       filename = file.path(PLOTS.DIR, paste0(Sys.Date(), "_barplot_SynergyScores_between_chemotherapeutics.pdf")),
       width = 23, height = 10)
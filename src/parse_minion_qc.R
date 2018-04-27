#!/usr/bin/env Rscript

library(data.table)
library(ggplot2)
library(yaml)

#############
# FUNCTIONS #
#############

ParseMinionQCYaml <- function(x){
    my_dt <- data.table(t(as.data.frame(x)),
                        keep.rownames = TRUE)
    my_dt[, rn := gsub("[^[:alnum:]]+", "_", rn)]
    my_dt[startsWith(rn, "All_reads"),
          c("category",
            "variable") := .("all_reads",
                             sub("All_reads_", "", rn)) ]
    my_dt[startsWith(rn, "Q_7"),
          c("category",
            "variable") := .("q7",
                             sub("Q_7_", "", rn)) ]
    my_dt[!is.na(category)][, .(category,
                                variable,
                                value = as.numeric(V1))]
}

###########
# GLOBALS #
###########

yaml_files <- snakemake@input[["yaml"]]
log_file <- snakemake@log[["log"]]

# dev
# qc_dir <- "output/010_minion-qc"
# yaml_files <- list.files(path = qc_dir,
#                          pattern = "summary.yaml",
#                          recursive = TRUE,
#                          full.names = TRUE)

########
# MAIN #
########

# set log
log <- file(log_file, open = "wt")
sink(log, type = "message")
sink(log, append = TRUE, type = "output")

# read yaml
names(yaml_files) <- gsub(".+/", "", dirname(yaml_files))
yaml_list <- lapply(yaml_files, read_yaml)

# convert to data table
minion_qc_list <- lapply(yaml_list, ParseMinionQCYaml)
plot_data_all <- rbindlist(minion_qc_list, idcol = "sample")
sample_order <- c("blue_cod" = "Blue cod",
                  "sample_3" = "Sample 3",
                  "sample_5" = "Sample 5",
                  "sample_7" = "Sample 7",
                  "sample_8" = "Sample 8",
                  "sample_10" = "Sample 10")
category_order <- c("all_reads" = "All reads",
                    "q7" = "Q7")
plot_data_all[, sample := factor(plyr::revalue(sample, sample_order),
                                 levels = sample_order)]
plot_data_all[, category := factor(plyr::revalue(category, category_order),
                                   levels = category_order)]

# what do we actually want to plot
kept_vars <- c("total_reads" = "Total reads\n(thousands)",
               "total_gigabases" = "Total\ngigabases",
               "gigabases_10kb" = "Gigabases\n> 10 Kb",
               "N50_length" = "N50 length\n(Kb)", 
               "mean_length" = "Mean length\n(Kb)")
plot_data <- plot_data_all[variable %in% names(kept_vars)]
plot_data[variable == "total_reads", value := value/1e3]
plot_data[variable == "N50_length", value := value/1e3]
plot_data[variable == "mean_length", value := value/1e3]
plot_data[, variable := factor(plyr::revalue(variable, kept_vars),
                               levels = kept_vars)]

# basic plots
gp <- ggplot(plot_data,
       aes(x = sample, y = value)) +
    theme_minimal(base_size = 14) +
    theme(strip.placement = "outside",
          strip.background = element_blank(),
          axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1)) +
    facet_grid(variable ~ category, scales = "free_y", switch = "y") +
    xlab(NULL) + ylab(NULL) +
    scale_fill_brewer(palette = "Set1",
                      guide = FALSE) +
    geom_col(width = 0.8)

# website png
dpi <- 144
Set1 <- RColorBrewer::brewer.pal(9, "Set1")
ggsave(snakemake@output[["png"]],
       gp + geom_col(fill = Set1[1]),
       width = 1200/dpi,
       height = 675*1.5/dpi,
       dpi = dpi,
       units = "in")

# write output
ggsave(filename = snakemake@output[["gp"]],
       plot = gp,
       device = "pdf",
       width = 10,
       height = 7.5,
       units = "in")
saveRDS(plot_data_all, snakemake@output[["parsed_data"]])

# write log
sessionInfo()


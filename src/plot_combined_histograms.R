#!/usr/bin/env Rscript

library(data.table)
library(ggplot2)

###########
# GLOBALS #
###########

hist_files <- snakemake@input[["lhist"]]
log_file <- snakemake@log[["log"]]

# dev
# hist_files <- c("output/020_stats/asw_11/all_hist.txt",
#                 "output/020_stats/asw_11/pass_hist.txt",
#                 "output/020_stats/asw_12a/all_hist.txt",
#                 "output/020_stats/asw_12a/pass_hist.txt",
#                 "output/020_stats/asw_12b/all_hist.txt",
#                 "output/020_stats/asw_12b/pass_hist.txt",
#                 "output/020_stats/blue_cod/all_hist.txt",
#                 "output/020_stats/blue_cod/pass_hist.txt",
#                 "output/020_stats/lambda_qc/all_hist.txt",
#                 "output/020_stats/lambda_qc/pass_hist.txt",
#                 "output/020_stats/mh_78/all_hist.txt",
#                 "output/020_stats/mh_78/pass_hist.txt",
#                 "output/020_stats/sample_10/all_hist.txt",
#                 "output/020_stats/sample_10/pass_hist.txt",
#                 "output/020_stats/sample_3/all_hist.txt",
#                 "output/020_stats/sample_3/pass_hist.txt",
#                 "output/020_stats/sample_5/all_hist.txt",
#                 "output/020_stats/sample_5/pass_hist.txt",
#                 "output/020_stats/sample_7/all_hist.txt",
#                 "output/020_stats/sample_7/pass_hist.txt",
#                 "output/020_stats/sample_8/all_hist.txt",
#                 "output/020_stats/sample_8/pass_hist.txt")

########
# MAIN #
########

# set log
log <- file(log_file, open = "wt")
sink(log, type = "message")
sink(log, append = TRUE, type = "output")

# read hists
names(hist_files) <- paste(gsub(".+/", "", dirname(hist_files)),
                           gsub("_.+", "", basename(hist_files)),
                           sep = "|")
hist_data_list <- lapply(hist_files, fread)
hist_data <- rbindlist(hist_data_list, idcol = "filename")
hist_data[, c("sample_name", "type") := tstrsplit(filename, "\\|")]
hist_data[, filename := NULL]

# set sample orders
sample_order <- c("lambda_qc" = "Lambda QC",
                  "asw_11" = "ASW 11",
                  "asw_12a" = "ASW 12a",
                  "asw_12b" = "ASW 12b",
                  "asw_14" = "ASW 14",
                  "asw_31" = "ASW 31",
                  "asw_33" = "ASW 33",
                  "asw_34" = "ASW 34",
                  "mh_78" = "M. hyp 7+8",
                  'bm_wt' = "BM WT",
                  'bm_c' = "BM mutant",
                  'lug_8' = "Stonefly",
                  "blue_cod" = "Blue cod",
                  "sample_3" = "Sample 3",
                  "sample_5" = "Sample 5",
                  "sample_7" = "Sample 7",
                  "sample_8" = "Sample 8",
                  "sample_10" = "Sample 10")
category_order <- c("all" = "All reads",
                    "pass" = "Q7")
hist_data[, sample_name := factor(plyr::revalue(sample_name, sample_order),
                                  levels = sample_order)]
hist_data[, type := factor(plyr::revalue(type, category_order),
                           levels = category_order)]

# transform into log4
hist_data[, l4_length := log(`#Length`, 4)]

# bin manually
n_bins <- 50
qa <- hist_data[, seq(min(l4_length), max(l4_length), length.out = n_bins)]
bin_labels <- sapply(1:length(qa), function(i)
    mean(c(qa[i], qa[i + 1])))[c(1:n_bins - 1)]
hist_data[, lbin := as.numeric(
    as.character(
        cut(l4_length, breaks = qa,
            labels = bin_labels,
            include.lowest = TRUE)))]
lhist_log4 <- hist_data[, .(count = sum(Count)), by = .(lbin, sample_name, type)]

# weighted by number of reads
hist_data[, total_bases := `#Length` * Count]
wlhist <- hist_data[, .(total_bases = sum(total_bases)), by = .(lbin, sample_name, type)]

# plots
lh <- ggplot(lhist_log4, aes(x = lbin, y = count)) +
    facet_grid(sample_name ~ type, scales = "free_y") +
    scale_x_continuous(labels = function(x) 4^x) +
    ylab("Count") + xlab("Read length") +
    geom_col()

wlh <- ggplot(wlhist , aes(x = lbin, y = total_bases/1e6)) +
    facet_grid(sample_name ~ type, scales = "free_y") +
    scale_x_continuous(labels = function(x) 4^x) +
    ylab("Total megabases") + xlab("Read length") +
    geom_col()

# website png
dpi <- 144
Set1 <- RColorBrewer::brewer.pal(9, "Set1")
ggsave(snakemake@output[["png"]],
       wlh + theme_minimal(base_size = 14) + geom_col(fill = Set1[1]),
       width = 1200/dpi,
       height = 675*2/dpi,
       dpi = dpi,
       units = "in")

# write output
ggsave(filename = snakemake@output[["lh"]],
       plot = lh,
       device = "pdf",
       width = 10,
       height = 7.5,
       units = "in")

ggsave(filename = snakemake@output[["wh"]],
       plot = wlh,
       device = "pdf",
       width = 10,
       height = 7.5,
       units = "in")

# write session info
sessionInfo()

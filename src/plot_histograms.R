#!/usr/bin/env Rscript

library(data.table)
library(ggplot2)

###########
# GLOBALS #
###########

lhist_file <- snakemake@input[["lhist"]]
log_file <- snakemake@log[["log"]]

########
# MAIN #
########


# set log
log <- file(log_file, open = "wt")
sink(log, type = "message")
sink(log, append = TRUE, type = "output")

# read hist
lhist <- fread(lhist_file)

# transform into log4
lhist[, l4_length := log(`#Length`, 4)]

# bin manually
n_bins <- 50
qa <- lhist[, seq(min(l4_length), max(l4_length), length.out = n_bins)]
bin_labels <- sapply(1:length(qa), function(i)
    mean(c(qa[i], qa[i + 1])))[c(1:n_bins - 1)]
lhist[, lbin := as.numeric(
    as.character(
        cut(l4_length, breaks = qa,
            labels = bin_labels,
            include.lowest = TRUE)))]
lhist_log4 <- lhist[, .(count = sum(Count)), by = lbin]

# weighted by number of reads
lhist[, total_bases := `#Length` * Count]
wlhist <- lhist[, .(total_bases = sum(total_bases)), by = lbin]

# plot
lh <- ggplot(lhist_log4, aes(x = lbin, y = count)) +
    scale_x_continuous(labels = function(x) 4^x) +
    ylab("Count") + xlab("Read length") +
    geom_col()

wlh <- ggplot(wlhist , aes(x = lbin, y = total_bases/1e6)) +
    scale_x_continuous(labels = function(x) 4^x) +
    ylab("Total megabases") + xlab("Read length") +
    geom_col()

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




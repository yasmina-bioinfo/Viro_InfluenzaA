#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(readr)
})

samples_csv <- "data/metadata/samples.csv"
salmon_dir  <- "results/salmon"
out_file    <- "results/viral_load/viral_load_reads.tsv"

dir.create("results/viral_load", showWarnings = FALSE, recursive = TRUE)

meta <- read_csv(samples_csv, show_col_types = FALSE) %>%
  mutate(sample = as.character(sample_id),
         condition = as.factor(condition))

# Define how to detect viral transcripts in Name column of quant.sf:
# Option A: viral FASTA headers have a prefix like "IAV_" or "WSN33_"
# Option B: put viral segments names that appear in quant.sf Name
viral_pattern <- "(WSN|IAV|Influenza|segment|PB2|PB1|PA|HA|NP|NA|M|NS)"

get_viral_reads <- function(sample) {
  qf <- file.path(salmon_dir, sample, "quant.sf")
  stopifnot(file.exists(qf))
  q <- fread(qf)
  qv <- q %>% filter(grepl(viral_pattern, Name, ignore.case = TRUE))
  tibble(sample = sample,
         viral_NumReads = sum(qv$NumReads))
}

viral <- bind_rows(lapply(meta$sample, get_viral_reads)) %>%
  left_join(meta, by = "sample") %>%
  arrange(condition, sample)

fwrite(viral, out_file, sep = "\t")
message("Saved: ", out_file)

#!/usr/bin/env Rscript

# Input files are either comma or tab delimited
# Output file is tab delimited
library(argparse)
parser <- ArgumentParser(description = "Rearrange read depths from a file")
parser$add_argument("readdepths", help = "Input file (delimited)")
parser$add_argument("metadata", help = "Sample metadata (delimited)")
parser$add_argument("output", help = "Output file (tsv)")
args <- parser$parse_args()

if (!file.exists(args$readdepths)) {
    stop(sprintf("Input file does not exist: %s", args$readdepths))
}

if (!file.exists(args$metadata)) {
    stop(sprintf("Metadata file does not exist: %s", args$metadata))
}

library(data.table)
readdepths <- fread(args$readdepths, header = TRUE, sep = "\t")
metadata <- fread(args$metadata, header = TRUE, sep = ",")

# Rearrange read depths
matched_samples <- metadata[!is.na(hostname), samplename]

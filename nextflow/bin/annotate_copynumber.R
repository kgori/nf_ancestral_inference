#!/usr/bin/env Rscript

library(argparse)
parser <- ArgumentParser()
parser$add_argument("cnv_dataset")
parser$add_argument("metadata")
parser$add_argument("variants_file")
parser$add_argument("output")
args <- parser$parse_args()

if (!dir.exists(args$cnv_dataset)) {
    stop(sprintf("Directory does not exist: %s", args$cnv_dataset))
}

if (!file.exists(args$metadata)) {
    stop(sprintf("File does not exist: %s", args$metadata))
}

if (!file.exists(args$variants_file)) {
    stop(sprintf("File does not exist: %s", args$variants_file))
}

suppressPackageStartupMessages({
    library(arrow)
    library(dplyr)
    library(data.table)
})

cnv_ds <- open_dataset(args$cnv_dataset)
metadata <- fread(args$metadata, na.strings = "")
variants <- fread(args$variants_file, na.strings = "")
variants <- variants[!is.na(hostname)]

# Annotate with purity and ploidy information
variants[metadata, purity := i.purity, on = "samplename"]
variants[metadata, ploidy := i.estimatedPloidy, on = "samplename"]

# Transfer over the logR from the cnv dataset
cnv_ds %>%
    filter(CHROM %in% c("1", "7", "8", "9", "21", "34")) %>%
    select(CHROM, START, END, samplename, T_logr, total_cn,
        dataset_segment_id, pipeline_cn) %>%
    collect() -> cnv_dt
setkey(cnv_dt, samplename, CHROM, START, END)

variants[, CHROM := as.character(CHROM)]
variants[, START := POS]
variants[, END := POS]
setkey(variants, samplename, CHROM, START, END)

variants <- foverlaps(variants, cnv_dt, nomatch = NULL)
variants[, c("START", "END", "i.START", "i.END") := NULL]
setcolorder(variants, c("CHROM", "POS", "REF", "ALT",
                        "samplename", "T_total_reads",
                        "T_alt_reads", "H_total_reads",
                        "H_alt_reads", "T_vaf", "H_vaf",
                        "T_logr", "total_cn", "dataset_segment_id",
                        "pipeline_cn", "purity", "ploidy"))

variants[, tumour_prob := cnpipe:::prob_read_came_from_tumour(T_logr,
    purity, ploidy, 2)]
variants[pipeline_cn == 0, tumour_prob := 0]
variants[, T_corrected_vaf := cnpipe:::fast_estimate_tumour_vaf(
                     T_total_reads,
                     T_alt_reads,
                     T_logr,
                     pmax(1, H_total_reads), # protect against divide by 0
                     H_alt_reads,
                     purity,
                     ploidy,
                     2)]
variants[pipeline_cn == 0, T_corrected_vaf := 0]

variants[, expected_t_reads := tumour_prob * T_total_reads]
variants[T_alt_reads > 3 * expected_t_reads & pipeline_cn == 0,
    T_corrected_vaf := 0]
variants[, is_present := T_alt_reads >= 3 & T_corrected_vaf > 0.1]
variants[pipeline_cn == 1,
    is_present := T_alt_reads >= 3 & T_corrected_vaf > 0.4]

fwrite(variants, args$output, sep = "\t")

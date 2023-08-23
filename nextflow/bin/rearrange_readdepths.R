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

# Open the readdepths file, and do some tidying
# Rename the samples to match the metadata
readdepths <- fread(args$readdepths, header = TRUE, sep = "\t")
setorder(readdepths, SAMPLE, CHROM, POS, REF, ALT)
readdepths[, SAMPLE := gsub("-Dog$", "", SAMPLE)]
readdepths[SAMPLE == "DCFAM.HKM1", SAMPLE := "1380H"]
readdepths[SAMPLE == "DCFAM.HKM2", SAMPLE := "1381H"]
readdepths[SAMPLE == "DT.KM1", SAMPLE := "1380T"]
readdepths[SAMPLE == "DT.KM2", SAMPLE := "1381T"]
readdepths[SAMPLE == "399H", SAMPLE := "1322_399H"]
readdepths[SAMPLE == "399T", SAMPLE := "1322_399T"]

metadata <- fread(args$metadata, header = TRUE, sep = ",",
    na.strings = "")

# Rearrange read depths
matched_samples <- metadata[!is.na(host), tumour]
all_samples <- intersect(readdepths[, unique(SAMPLE)], metadata[, tumour])
dt <- rbindlist(lapply(all_samples, function(samplename_) {
    if (samplename_ %in% matched_samples) {
        matched_host <- metadata[tumour == samplename_, host]
        tumour_df <- readdepths[SAMPLE == samplename_,
            .(CHROM, POS, REF, ALT, samplename = samplename_,
              hostname = matched_host,
              T_total_reads = NR, T_alt_reads = NV,
              T_vaf = NV / NR)]
        host_df <- readdepths[SAMPLE == matched_host,
            .(CHROM, POS, REF, ALT,
              H_total_reads = NR, H_alt_reads = NV,
              H_vaf = NV / NR)]
        merged_df <- merge(tumour_df, host_df,
            by = c("CHROM", "POS", "REF", "ALT"))
    } else {
        merged_df <- readdepths[SAMPLE == samplename_,
            .(CHROM, POS, REF, ALT, samplename = samplename_,
              hostname = NA_character_,
              T_total_reads = NR, T_alt_reads = NV, T_vaf = NV / NR,
              H_total_reads = NA_integer_, H_alt_reads = NA_integer_,
              H_vaf = NA_real_)]
    }
    return(merged_df)
}))

fwrite(dt, args$output, sep = "\t")

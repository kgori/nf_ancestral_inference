#!/usr/bin/env Rscript

library(argparse)
parser <- ArgumentParser()
parser$add_argument("variants_file")
parser$add_argument("outfile")
args <- parser$parse_args()

if (!file.exists(args$variants_file)) {
  stop(sprintf("Variants file does not exist: %s", args$variants_file))
}

getmode <- function(v) {
    uniqv <- unique(v)
    uniqv[which.max(tabulate(match(v, uniqv)))]
}

library(data.table)
variants <- fread(args$variants_file, na.strings = "")
variants[, is_ctvt_a := samplename %in% c("982T",
    "401T", "1322_399T", "1382_399T")]

variants[pipeline_cn > 5 & T_corrected_vaf <= 0.08,
     c("refcount", "altcount") := .(pipeline_cn, 0)]
variants[pipeline_cn > 5 & T_corrected_vaf > 0.08 & T_corrected_vaf <= 0.5,
     c("refcount", "altcount") := .(pipeline_cn - 1, 1)]
variants[pipeline_cn > 5 & T_corrected_vaf > 0.5 & T_corrected_vaf <= 0.91,
     c("refcount", "altcount") := .(1, pipeline_cn - 1)]
variants[pipeline_cn > 5 & T_corrected_vaf > 0.91,
     c("refcount", "altcount") := .(0, pipeline_cn)]

variants[pipeline_cn == 5 & T_corrected_vaf <= 0.09,
     c("refcount", "altcount") := .(5, 0)]
variants[pipeline_cn == 5 & T_corrected_vaf > 0.09 & T_corrected_vaf <= 0.33,
     c("refcount", "altcount") := .(4, 1)]
variants[pipeline_cn == 5 & T_corrected_vaf > 0.33 & T_corrected_vaf <= 0.5,
     c("refcount", "altcount") := .(3, 2)]
variants[pipeline_cn == 5 & T_corrected_vaf > 0.5 & T_corrected_vaf <= 0.67,
     c("refcount", "altcount") := .(2, 3)]
variants[pipeline_cn == 5 & T_corrected_vaf > 0.67 & T_corrected_vaf <= 0.91,
     c("refcount", "altcount") := .(1, 4)]
variants[pipeline_cn == 5 & T_corrected_vaf > 0.91,
     c("refcount", "altcount") := .(0, 5)]

variants[pipeline_cn == 4 & T_corrected_vaf <= 0.12,
     c("refcount", "altcount") := .(4, 0)]
variants[pipeline_cn == 4 & T_corrected_vaf > 0.12 & T_corrected_vaf <= 0.39,
     c("refcount", "altcount") := .(3, 1)]
variants[pipeline_cn == 4 & T_corrected_vaf > 0.39 & T_corrected_vaf <= 0.64,
     c("refcount", "altcount") := .(2, 2)]
variants[pipeline_cn == 4 & T_corrected_vaf > 0.64 & T_corrected_vaf <= 0.86,
     c("refcount", "altcount") := .(1, 3)]
variants[pipeline_cn == 4 & T_corrected_vaf > 0.86,
     c("refcount", "altcount") := .(0, 4)]

variants[pipeline_cn == 3 & T_corrected_vaf <= 0.15,
     c("refcount", "altcount") := .(3, 0)]
variants[pipeline_cn == 3 & T_corrected_vaf > 0.15 & T_corrected_vaf <= 0.5,
     c("refcount", "altcount") := .(2, 1)]
variants[pipeline_cn == 3 & T_corrected_vaf > 0.5 & T_corrected_vaf <= 0.82,
     c("refcount", "altcount") := .(1, 2)]
variants[pipeline_cn == 3 & T_corrected_vaf > 0.82,
     c("refcount", "altcount") := .(0, 3)]

variants[pipeline_cn == 2 & T_corrected_vaf <= 0.25,
     c("refcount", "altcount") := .(2, 0)]
variants[pipeline_cn == 2 & T_corrected_vaf > 0.25 & T_corrected_vaf <= 0.75,
     c("refcount", "altcount") := .(1, 1)]
variants[pipeline_cn == 2 & T_corrected_vaf > 0.75,
     c("refcount", "altcount") := .(0, 2)]

variants[pipeline_cn == 1 & T_corrected_vaf <= 0.5,
     c("refcount", "altcount") := .(1, 0)]
variants[pipeline_cn == 1 & T_corrected_vaf > 0.5,
     c("refcount", "altcount") := .(0, 1)]

variants[pipeline_cn == 0, c("refcount", "altcount") := .(0, 0)]

variants[, is_het := refcount > 0 & altcount > 0]
variants[, is_hom_ref := FALSE]
variants[, is_hom_alt := FALSE]
variants[is_het == FALSE, is_hom_ref := T_corrected_vaf < 0.5]
variants[is_het == FALSE, is_hom_alt := T_corrected_vaf > 0.5]

variants[, genotype := ifelse(is_het, "HET",
                              ifelse(is_hom_ref, "HOMREF", "HOMALT"))]

variants[, modeCN := getmode(pipeline_cn),
    by = .(is_ctvt_a, CHROM, dataset_segment_id)]
variants[unique(variants[is_ctvt_a == TRUE, .(CHROM, POS, REF, ALT, modeCN)]),
         modeCN_ctvtA := i.modeCN, on = c("CHROM", "POS", "REF", "ALT")]
variants[unique(variants[is_ctvt_a == FALSE, .(CHROM, POS, REF, ALT, modeCN)]),
         modeCN_ctvtBG := i.modeCN, on = c("CHROM", "POS", "REF", "ALT")]
variants[, step_change := modeCN_ctvtA - modeCN_ctvtBG]

anc_reconstr <- variants[(is_ctvt_a == FALSE &
                          step_change == 1 &
                          pipeline_cn == modeCN),
                         .(n_homref = sum(genotype == "HOMREF"),
                           n_het = sum(genotype == "HET"),
                           n_homalt = sum(genotype == "HOMALT"),
                           mode_alleleA = getmode(refcount),
                           mode_alleleB = getmode(altcount),
                           modeCN = getmode(pipeline_cn)),
                         by = .(CHROM, POS, REF, ALT)]
setkey(anc_reconstr, CHROM, POS, REF, ALT)

anc_reconstr[, mode_value := apply(.SD, 1, max),
             .SDcols = c("n_homref", "n_het", "n_homalt")]

# How many categories match the max value?
anc_reconstr[, num_modes := anc_reconstr[,
    .(n_homref == mode_value, n_het == mode_value, n_homalt == mode_value)][,
        rowSums(.SD)]]

# Initially all variants have unknown genotype, '?'. There should be no '?'
# left after these next steps.
anc_reconstr[, genotype := "?"]
anc_reconstr[modeCN > 1 & num_modes == 1,
             genotype := ifelse(n_homref == mode_value, "0/0",
                                ifelse(n_het == mode_value,
                                       "0/1",
                                       "1/1"))]
anc_reconstr[modeCN == 1 & num_modes == 1,
             genotype := ifelse(n_homref == mode_value, "0",
                                ifelse(n_het == mode_value,
                                       ".",
                                       "1"))]
anc_reconstr[modeCN == 0,
             genotype := "./."]

anc_reconstr[num_modes > 1, genotype := ifelse(modeCN == 1,
                                               "NA",
                                               "NA")]

if (any(anc_reconstr$genotype == "?")) {
  warning("There are still unknown genotypes in the ancestral reconstruction.")
}

breakdown <- variants[(is_ctvt_a == FALSE &
    step_change == 1 &
    pipeline_cn == modeCN),
  .(mode_refcount = getmode(refcount), mode_altcount = getmode(altcount)),
  by = .(CHROM, POS, REF, ALT, total_cn = modeCN)]
breakdown[anc_reconstr[num_modes > 1],
    c("mode_refcount", "mode_altcount") := -1,
    on = .(CHROM, POS, REF, ALT)]
rbind(breakdown[mode_refcount != -1, .N,
    by = .(total_cn, mode_refcount, mode_altcount)][order(total_cn,
        mode_altcount)],
    breakdown[mode_refcount == -1, .N,
    by = .(total_cn, mode_refcount, mode_altcount)][order(total_cn,
        mode_altcount)])

anc_reconstr <- anc_reconstr[genotype != "NA"]

ht_inference <- variants[(is_ctvt_a == TRUE &
                          step_change == 1 &
                          pipeline_cn == modeCN),
                         .(mode_alleleA_ctvtA = getmode(refcount),
                           mode_alleleB_ctvtA = getmode(altcount),
                           modeCN_ctvtA = getmode(pipeline_cn)),
                         by = .(CHROM, POS, REF, ALT)]
setkey(ht_inference, CHROM, POS, REF, ALT)
setkey(variants, CHROM, POS, REF, ALT)

# Look for mode ambiguity
test <- variants[(is_ctvt_a == TRUE &
    step_change == 1 &
    pipeline_cn == modeCN)][ht_inference]
ambiguous_positions <- test[,
  .(mode_match_A = sum(refcount == mode_alleleA_ctvtA),
    mode_match_B = sum(altcount == mode_alleleB_ctvtA),
    N = .N),
    by = .(CHROM, POS, REF, ALT)]

# Uncallable positions are those where all samples disagree, so there is
# no consensus.
uncallable <- ambiguous_positions[mode_match_A == 1 & mode_match_B == 1]
ht_inference <- ht_inference[!uncallable]

# Recoverable positions are those where there are two possible values
# for the mode. In these cases, we choose the value supported by the
# higher purity samples (not 401T).
# To find the ambiguous mode sites, first find the sites where 2 out of 4
# samples match the mode.
recoverable <- ambiguous_positions[mode_match_A == 2 & N == 4]
setkey(recoverable, CHROM, POS, REF, ALT)
setkey(test, CHROM, POS, REF, ALT)

# Then, find the sites where the remaining 2 out of 4 have the same value
# (there is only one unique other value).
to_keep <- test[recoverable][refcount != mode_alleleA_ctvtA,
    .(N_other = uniqueN(refcount)),
    by = .(CHROM, POS, REF, ALT)][N_other == 1]
recoverable <- recoverable[to_keep]

# Now find the mode allele A and B values at the Recoverable sites, dropping
# sample 401T
recovery_table <- variants[recoverable][(
    is_ctvt_a == TRUE &
    step_change == 1 &
    pipeline_cn == modeCN &
    samplename != "401T"),
    .(mode_alleleA_ctvtA_new = getmode(refcount),
      mode_alleleB_ctvtA_new = getmode(altcount),
      modeCN_ctvtA_new = getmode(pipeline_cn)),
    by = .(CHROM, POS, REF, ALT)]

max_recover <- recovery_table[, .N]
changes <- ht_inference[recovery_table,
    sum(mode_alleleA_ctvtA != mode_alleleA_ctvtA_new)]

ht_inference[recovery_table, c("mode_alleleA_ctvtA", "mode_alleleB_ctvtA",
                               "modeCN_ctvtA") :=
    .(mode_alleleA_ctvtA_new, mode_alleleB_ctvtA_new, modeCN_ctvtA_new)]

ht_inference <- ht_inference[anc_reconstr[, .(CHROM, POS, REF, ALT,
                                              mode_alleleA, mode_alleleB,
                                              modeCN)], ,
                             on = c("CHROM", "POS", "REF", "ALT"),
                             nomatch = NULL]

ht_inference[, diffA := mode_alleleA_ctvtA - mode_alleleA]
ht_inference[, diffB := mode_alleleB_ctvtA - mode_alleleB]
ht_inference[, diffCN := modeCN_ctvtA - modeCN]

# Initially all variants have unknown haplotype, '?'. There should be no '?'
# left after these next steps.
ht_inference[, haplotype := "?"]
ht_inference[diffA == 1 & diffB == 0, haplotype := "0"]
ht_inference[diffA == 0 & diffB == 1, haplotype := "1"]
ht_inference[haplotype == "?" & diffA == 2 & diffB == -1, haplotype := "."]
ht_inference[haplotype == "?" & diffA == -1 & diffB == 2, haplotype := "."]
ht_inference[haplotype == "?" & diffA == 3 & diffB == -2, haplotype := "."]
ht_inference[haplotype == "?" & diffA == -2 & diffB == 3, haplotype := "."]

if (any(ht_inference$haplotype == "?")) {
  warning("There are still unknown haplotypes in the haplotype inference.")
}

conn <- file(args$outfile)
writeLines(c(
    "##fileformat=VCFv4.2",
    "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\"",
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tCTVT\tHT"
), con = conn)

vcf <- ht_inference[anc_reconstr, , nomatch = NULL][,
    .(CHROM, POS, ID = ".", REF, ALT, QUAL = ".",
      FILTER = ".", INFO = ".", FORMAT = "GT",
      CTVT = genotype, HT = haplotype)][HT != "."]
fwrite(vcf, args$outfile, sep = "\t", append = TRUE)

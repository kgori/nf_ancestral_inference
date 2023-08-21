#!/bin/bash -ue

nextflow run nextflow/main.nf \
    --genotypedVcf=data/platypus_result.vcf.gz \
    --metadata=data/metadata.tsv \
    -resume


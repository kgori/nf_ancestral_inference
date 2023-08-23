#!/bin/bash -ue

nextflow run nextflow/main.nf \
    --genotypedVcf=data/platypus_result.vcf.gz \
    --metadata=data/dog_metadata.csv \
    --csvDataSet=data/copynumber_arrow_dataset \
    --outDir=results \
    -resume



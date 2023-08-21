nextflow.enable.dsl=2

params.genotypedVcf
params.metadata

process extract_readdepths {

    input:
    path vcfFile

    output:
    path "${vcfFile.getBaseName(2)}.readdepths.tsv.gz"

    script:
    """
    echo "CHROM\\tPOS\\tREF\\tALT\\tSAMPLE\\tNR\\tNV\\n" > ${vcfFile.getBaseName(2)}.readdepths.tsv
    bcftools view -m2 -M2 ${vcfFile} | \
      bcftools filter -i '(((REF="A"|REF="G") & (ALT="C"|ALT="T")) | ((REF="C"|REF="T") & (ALT="A"|ALT="G")))' | \
      bcftools filter -e '(REF="A" & ALT="G") | (REF="G" & ALT="A") | (REF="C" & ALT="T") | (REF="T" & ALT="C")' | \
      bcftools query -f '[%CHROM\\t%POS\\t%REF\\t%ALT\\t%SAMPLE\\t%NR\\t%NV\\n]' ${vcfFile} >> ${vcfFile.getBaseName(2)}.readdepths.tsv

    gzip ${vcfFile.getBaseName(2)}.readdepths.tsv
    """
}

process rearrange_readdepths {

    input:
    path readdepthsFile
    path metadata

    output:
    path "${readdepthsFile}.rearranged.tsv.gz"

    script:
    """
    rearrange_readdepths.R ${readdepthsFile} ${metadata} ${readdepthsFile}.rearranged.tsv
    gzip ${readdepthsFile}.rearranged.tsv
    """
}

workflow {
    
    genotyped_vcf = Channel.fromPath(params.genotypedVcf, checkIfExists: true)
    read_depths = extract_readdepths(genotyped_vcf)
    read_depths.view()
}

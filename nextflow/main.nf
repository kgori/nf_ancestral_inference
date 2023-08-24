nextflow.enable.dsl=2

params.genotypedVcf
params.panelVcf
params.metadata
params.csvDataSet
params.outDir

process extract_readdepths {

    input:
    path vcfFile

    output:
    path "${vcfFile.getBaseName(2)}.readdepths.tsv.gz"

    script:
    """
    echo "CHROM\tPOS\tREF\tALT\tSAMPLE\tNR\tNV" > ${vcfFile.getBaseName(2)}.readdepths.tsv
    bcftools view -m2 -M2 ${vcfFile} | \
      bcftools filter -i '(((REF="A"|REF="G") & (ALT="C"|ALT="T")) | ((REF="C"|REF="T") & (ALT="A"|ALT="G")))' | \
      bcftools filter -e '(REF="A" & ALT="G") | (REF="G" & ALT="A") | (REF="C" & ALT="T") | (REF="T" & ALT="C")' | \
      bcftools query -f '[%CHROM\\t%POS\\t%REF\\t%ALT\\t%SAMPLE\\t%NR\\t%NV\\n]' >> ${vcfFile.getBaseName(2)}.readdepths.tsv

    gzip ${vcfFile.getBaseName(2)}.readdepths.tsv
    """
}

process rearrange_readdepths {

    input:
    path readdepthsFile
    path metadata

    output:
    path "${readdepthsFile.getBaseName(2)}.rearranged.tsv.gz"

    publishDir "${params.outDir}/readdepths"

    script:
    """
    rearrange_readdepths.R ${readdepthsFile} ${metadata} ${readdepthsFile.getBaseName(2)}.rearranged.tsv
    gzip ${readdepthsFile.getBaseName(2)}.rearranged.tsv
    """
}

process annotate_copynumber {

    input:
    path csvs
    path variants
    path metadata

    output:
    path "${variants.getBaseName(2)}.annotated.tsv.gz"

    publishDir "${params.outDir}/annotated"

    script:
    """
    annotate_copynumber.R \
      ${csvs} \
      ${metadata} \
      ${variants} \
      ${variants.getBaseName(2)}.annotated.tsv
    gzip ${variants.getBaseName(2)}.annotated.tsv
    """
}

process ancestral_reconstruction {

    input:
    path variants

    output:
    tuple path("${variants.getBaseName(2)}.ancestral.vcf.gz"), path("${variants.getBaseName(2)}.ancestral.vcf.gz.csi")

    publishDir "${params.outDir}/ancestral"

    script:
    """
    ancestral_reconstruction.R ${variants} ${variants.getBaseName(2)}.ancestral.vcf
    bgzip ${variants.getBaseName(2)}.ancestral.vcf
    bcftools index ${variants.getBaseName(2)}.ancestral.vcf.gz
    """
}

process intersect_vcfs {
    cpus 4
    executor 'lsf'
    queue 'normal'
    memory 1.GB
    time 2.h

    input:
    tuple path(ancestralVcf), path(ancestralVcfIndex), path(otherVcf), path(otherVcfIndex)

    output:
    tuple path("${otherVcf.getBaseName(2)}.intersected.vcf.gz"), path("${otherVcf.getBaseName(2)}.intersected.vcf.gz.csi")

    script:
    """
    bcftools isec -c none -p . -Oz -n2 -w2 \
      "${ancestralVcf}" "${otherVcf}"
    mv "0001.vcf.gz" "${otherVcf.getBaseName(2)}.intersected.vcf.gz"
    bcftools index "${otherVcf.getBaseName(2)}.intersected.vcf.gz"
    """
}

process merge_vcfs {

    cpus 4
    executor 'lsf'
    queue 'normal'
    memory 8.GB
    time 2.h

    publishDir "${params.outDir}", mode: 'copy'

    input:
    tuple path(vcf1), path(vcf1Index), path(vcf2), path(vcf2Index), path(vcf3), path(vcf3Index)

    output:
    tuple path("final_merged.vcf.gz"), path("final_merged.vcf.gz.csi")

    script:
    """
    bcftools merge -Oz -o final_merged.vcf.gz ${vcf1} ${vcf2} ${vcf3}
    bcftools index final_merged.vcf.gz
    """
}

workflow {
    
    genotyped_vcf = Channel.fromPath(params.genotypedVcf, checkIfExists: true)
    panel_vcf = Channel.fromPath(params.panelVcf, checkIfExists: true)
    metadata = Channel.fromPath(params.metadata, checkIfExists: true)
    csv_dataset = Channel.fromPath(params.csvDataSet, checkIfExists: true)
    read_depths = extract_readdepths(genotyped_vcf)
    rearranged_read_depths = rearrange_readdepths(read_depths, metadata)
    annotated_copynumber = annotate_copynumber(csv_dataset, rearranged_read_depths, metadata)
    ancestral = ancestral_reconstruction(annotated_copynumber)

    other_vcfs = genotyped_vcf.map { it -> tuple(it, it + ".csi") } |
        concat(panel_vcf.map { it -> tuple(it, it + ".csi") })
    intersected = intersect_vcfs(ancestral.combine(other_vcfs))
    intersected.collect().combine(ancestral) | merge_vcfs
}

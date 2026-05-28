process FILTER_VARIANTS {
    tag "${meta.id}"
    label 'filter'
    conda 'bioconda::bcftools=1.19'
    publishDir "${params.outdir}/variants_filtered", mode: 'copy', pattern: '*.vcf.gz'

    input:
    tuple val(meta), path(vcf), path(tbi)

    output:
    tuple val(meta), path("${meta.id}.filtered.vcf.gz"), path("${meta.id}.filtered.vcf.gz.tbi"), emit: vcf

    script:
    """
    bcftools view -f PASS -Oz -o ${meta.id}.filtered.vcf.gz ${vcf} \\
        || bcftools view -Oz -o ${meta.id}.filtered.vcf.gz ${vcf}
    tabix -p vcf ${meta.id}.filtered.vcf.gz
    """

    stub:
    """
    echo '##fileformat=VCFv4.2' | gzip -c > ${meta.id}.filtered.vcf.gz
    touch ${meta.id}.filtered.vcf.gz.tbi
    """
}

workflow filter_variants {
    take:
        vcf

    main:
        FILTER_VARIANTS(vcf)

    emit:
        vcf = FILTER_VARIANTS.out.vcf
}

process BCFTOOLS_CALL {
    tag "$meta.id"
    label 'variant'
    conda 'bioconda::bcftools=1.19 bioconda::samtools=1.19'

    publishDir "${params.outdir}/variants", mode: 'copy', pattern: '*.{vcf.gz,tbi}'

    input:
    tuple val(meta), path(bam), path(bai), path(ref_fasta), val(index_fasta)

    output:
    tuple val(meta), path("${meta.id}.vcf.gz"), path("${meta.id}.vcf.gz.tbi"), emit: vcf

    script:
    """
    if [ "${index_fasta}" = "true" ]; then
        samtools faidx ${ref_fasta}
    fi

    bcftools mpileup -Ou -f ${ref_fasta} ${bam} \\
    | bcftools call -mv -Oz -o ${meta.id}.vcf.gz

    tabix -p vcf ${meta.id}.vcf.gz
    """
}

workflow BCFTOOLS_MPILEUP {
    take:
        bam
        ref
        index_fasta

    main:
        ch = bam.combine(ref).combine(index_fasta)
        BCFTOOLS_CALL(ch)

    emit:
        vcf = BCFTOOLS_CALL.out.vcf
}

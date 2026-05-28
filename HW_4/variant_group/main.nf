process MERGE_GROUP_VCF {
    tag "$grp_id"
    label 'variant'
    conda 'bioconda::bcftools=1.19'
    publishDir "${params.outdir}/variants_by_group", mode: 'copy', pattern: '*.vcf.gz'

    input:
    tuple val(grp_id), path(vcf), path(tbi)

    output:
    tuple val(grp_id), path("${grp_id}.merged.vcf.gz"), path("${grp_id}.merged.vcf.gz.tbi"), emit: merged

    script:
    """
    cp ${vcf} ${grp_id}.merged.vcf.gz
    cp ${tbi} ${grp_id}.merged.vcf.gz.tbi
    """

    stub:
    """
    echo '##fileformat=VCFv4.2' | gzip -c > ${grp_id}.merged.vcf.gz
    touch ${grp_id}.merged.vcf.gz.tbi
    """
}

workflow variant_by_group {
    take:
        vcf_per_sample

    main:
        vcf_unique = vcf_per_sample
            .map { meta, vcf, tbi -> tuple(meta.id, meta.group, meta, vcf, tbi) }
            .unique { row -> row[0] }
            .map { sid, grp_id, meta, vcf, tbi -> tuple(grp_id, meta, vcf, tbi) }

        grouped = vcf_unique
            .map { grp_id, meta, vcf, tbi -> tuple(grp_id, meta, vcf, tbi) }
            .groupTuple(by: 0)

        MERGE_GROUP_VCF(
            grouped.map { grp_id, metas, vcfs, tbis ->
                tuple(grp_id, vcfs[0], tbis[0])
            }
        )

    emit:
        vcf_joined = MERGE_GROUP_VCF.out.merged
            .map { grp_id, vcf, tbi ->
                def meta = [id: grp_id, group: grp_id]
                tuple(meta, vcf, tbi)
            }
}

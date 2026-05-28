#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { BCFTOOLS_MPILEUP as mpileup } from './modules/nf-core/bcftools/mpileup/main'
include { run_qc       } from './modules/run_qc/main'
include { trimm        } from './modules/trimm/main'
include { assemble     } from './modules/assemble/main'
include { map_and_plot } from './modules/map_and_plot/main'
include { FETCH_SRA    } from './modules/fetch_sra/main'
include { variant_by_group } from './modules/variant_group/main'
include { filter_variants  } from './modules/filter_variants/main'

workflow trimm_qc {
    take:
        trimmed_data

    main:
        run_qc('trimmed', trimmed_data)

    emit:
        reports = run_qc.out.reports
}

workflow {
    if (params.conda_env_path == null && System.getenv('CONDA_PREFIX'))
        params.conda_env_path = System.getenv('CONDA_PREFIX')

    def group_field = params.group_key ?: 'group'

    if (params.samplesheet) {
        def sheet = file(params.samplesheet, checkIfExists: true)

        raw_reads_ch = Channel
            .fromPath(sheet)
            .splitCsv(header: true)
            .map { row ->
                def meta = [
                    id         : row.sample,
                    group      : row[group_field],
                    single_end : false
                ]
                tuple(meta, [ file(row.read1), file(row.read2) ], file(row.reference))
            }

        sample_refs = raw_reads_ch.map { meta, reads, ref -> tuple(meta.id, ref) }

        run_qc('raw', raw_reads_ch.map { meta, reads, ref -> tuple(meta, reads) })
        trimmed_reads = trimm(raw_reads_ch.map { meta, reads, ref -> tuple(meta, reads) }).reads

        trimmed_with_ref = trimmed_reads
            .map { meta, reads -> tuple(meta.id, meta, reads) }
            .join(sample_refs)
            .map { id, meta, reads, ref -> tuple(meta, reads, ref) }

        if (params.ref_type == 'assemble')
            error 'HW4 samplesheet mode requires ref_type=alignment with per-sample reference column'

        map_and_plot(trimmed_with_ref)

    } else {
        def sra_id = params.numberSRA ?: params.sra_id

        if (sra_id) {
            def meta = [ id: sra_id, single_end: false, group: 'default' ]
            FETCH_SRA(tuple(meta, sra_id))
            raw_reads_ch = FETCH_SRA.out.reads.map { m, files ->
                tuple(m, files.sort(), file(params.reference))
            }
        } else if (params.input_reads_folder) {
            raw_reads_ch = Channel.fromFilePairs("${params.input_reads_folder}/*_{1,2}.{fq,fastq}{,.gz}", checkIfExists: true)
                .map { sample_id, files ->
                    def meta = [ id: sample_id, single_end: false, group: params.default_group ?: 'default' ]
                    tuple(meta, files, file(params.reference))
                }
        } else if (params.reads) {
            raw_reads_ch = Channel.fromFilePairs(params.reads, checkIfExists: true)
                .map { sample_id, files ->
                    def meta = [ id: sample_id, single_end: false, group: params.default_group ?: 'default' ]
                    tuple(meta, files, file(params.reference))
                }
        } else {
            error 'Provide --samplesheet or --reads / --sra_id / --input_reads_folder'
        }

        if (!params.reference)
            error 'Non-samplesheet mode requires --reference'

        run_qc('raw', raw_reads_ch.map { meta, reads, ref -> tuple(meta, reads) })
        trimmed_reads = trimm(raw_reads_ch.map { meta, reads, ref -> tuple(meta, reads) }).reads

        sample_refs = raw_reads_ch.map { meta, reads, ref -> tuple(meta.id, ref) }
        trimmed_with_ref = trimmed_reads
            .map { meta, reads -> tuple(meta.id, meta, reads) }
            .join(sample_refs)
            .map { id, meta, reads, ref -> tuple(meta, reads, ref) }

        if (params.ref_type == 'alignment') {
            map_and_plot(trimmed_with_ref)
        } else {
            reference_genome = assemble(trimmed_reads).contigs
            map_and_plot(
                trimmed_reads.join(reference_genome).map { meta, reads, ref -> tuple(meta, reads, ref) }
            )
        }
    }

    mpileup(
        map_and_plot.out.bam_out,
        map_and_plot.out.ref_out,
        Channel.value(false)
    )

    variant_by_group(mpileup.out.vcf)
    filter_variants(variant_by_group.out.vcf_joined)
}

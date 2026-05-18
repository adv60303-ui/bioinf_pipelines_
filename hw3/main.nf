#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { BCFTOOLS_MPILEUP as mpileup } from './modules/nf-core/bcftools/mpileup/main'
include { run_qc       } from './modules/run_qc/main'
include { trimm        } from './modules/trimm/main'
include { assemble     } from './modules/assemble/main'
include { map_and_plot } from './modules/map_and_plot/main'
include { FETCH_SRA    } from './modules/fetch_sra/main'

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

    def sra_id = params.numberSRA ?: params.sra_id

    if (sra_id) {
        def meta = [ id: sra_id, single_end: false ]
        FETCH_SRA(tuple(meta, sra_id))
        raw_reads = FETCH_SRA.out.reads.map { m, files ->
            [ m, files.sort() ]
        }
    } else if (params.input_reads_folder) {
        raw_reads = Channel.fromFilePairs("${params.input_reads_folder}/*_{1,2}.{fq,fastq}{,.gz}", checkIfExists: true)
            .map { sample_id, files ->
                def meta = [ id: sample_id, single_end: false ]
                [ meta, files ]
            }
    } else if (params.reads) {
        raw_reads = Channel.fromFilePairs(params.reads, checkIfExists: true)
            .map { sample_id, files ->
                def meta = [ id: sample_id, single_end: false ]
                [ meta, files ]
            }
    } else {
        error 'Provide --numberSRA / --sra_id, --input_reads_folder, or --reads'
    }

    run_qc('raw', raw_reads)
    trimmed_reads = trimm(raw_reads)
    trimm_qc(trimmed_reads)

    if (params.ref_type == 'alignment') {
        if (!params.reference)
            error 'ref_type=alignment requires --reference'
        reference_genome = Channel.fromPath(params.reference, checkIfExists: true).first()
        map_and_plot(reference_genome, trimmed_reads)
    } else {
        reference_genome = assemble(trimmed_reads).contigs
        map_and_plot(reference_genome, trimmed_reads)
    }

    mpileup(
        map_and_plot.out.bam_out,
        map_and_plot.out.ref_out,
        Channel.value(false)
    )
}

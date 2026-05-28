process FASTQC_RAW {
    tag "${meta.id} (raw)"
    label 'qc'
    conda 'bioconda::fastqc=0.12.1'
    publishDir "${params.outdir}/qc_raw", mode: 'copy'

    input:
    tuple val(meta), path(reads)

    output:
    path '*_fastqc.*', emit: reports

    script:
    """
    fastqc -t ${task.cpus} -o . ${reads.join(' ')}
    """
}

process FASTQC_TRIM {
    tag "${meta.id} (trimmed)"
    label 'qc'
    conda 'bioconda::fastqc=0.12.1'
    publishDir "${params.outdir}/qc_trimmed", mode: 'copy'

    input:
    tuple val(meta), path(reads)

    output:
    path '*_fastqc.*', emit: reports

    script:
    """
    fastqc -t ${task.cpus} -o . ${reads.join(' ')}
    """
}

workflow run_qc {
    take:
        reads_type
        reads

    main:
        if (reads_type == 'raw')
            FASTQC_RAW(reads)
        else
            FASTQC_TRIM(reads)

    emit:
        reports = reads_type == 'raw' ? FASTQC_RAW.out.reports : FASTQC_TRIM.out.reports
}

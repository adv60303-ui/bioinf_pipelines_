process FETCH_SRA {
    tag "${meta.id}"
    label 'download'
    conda 'bioconda::sra-tools=3.1.0'

    publishDir "${params.outdir}/raw_reads", mode: 'copy'

    input:
    tuple val(meta), val(sra_id)

    output:
    tuple val(meta), path('*.fastq.gz'), emit: reads

    script:
    """
    fasterq-dump --threads ${task.cpus} --split-files --outdir . ${sra_id}
    gzip -f *.fastq
    """
}

process TRIMMOMATIC {
    tag "${meta.id}"
    label 'trim'
    conda 'bioconda::trimmomatic=0.39'

    publishDir "${params.outdir}/trimmed_reads", mode: 'copy', pattern: '*_p.fastq.gz'

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("${meta.id}_R{1,2}_p.fastq.gz"), emit: reads

    script:
    """
    ADP="${params.adapters ?: ''}"
    if [ -z "\$ADP" ]; then
        ADP=\$(find \${CONDA_PREFIX:-/opt/conda} -name TruSeq3-PE.fa | head -n 1)
    fi
    test -n "\$ADP"

    trimmomatic PE -threads ${task.cpus} ${reads[0]} ${reads[1]} \\
        ${meta.id}_R1_p.fastq.gz /dev/null ${meta.id}_R2_p.fastq.gz /dev/null \\
        ILLUMINACLIP:\$ADP:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:36
    """
}

workflow trimm {
    take:
        reads

    main:
        TRIMMOMATIC(reads)

    emit:
        reads = TRIMMOMATIC.out.reads
}

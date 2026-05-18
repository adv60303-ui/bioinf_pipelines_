process MAP_READS {
    tag "${meta.id}"
    label 'map'
    conda 'bioconda::bwa=0.7.17 bioconda::samtools=1.19'

    publishDir "${params.outdir}/mapping", mode: 'copy'

    input:
    tuple val(meta), path(reads), path(reference)

    output:
    tuple val(meta), path("${meta.id}.sorted.bam"), path("${meta.id}.sorted.bam.bai"), path('ref.fa'), emit: bam

    script:
    """
    cp ${reference} ref.fa
    bwa index ref.fa
    bwa mem -t ${task.cpus} ref.fa ${reads[0]} ${reads[1]} | samtools sort -@ ${task.cpus} -o ${meta.id}.sorted.bam
    samtools index ${meta.id}.sorted.bam
    """
}

process PLOT_COVERAGE {
    tag "${meta.id}"
    label 'plot'
    conda 'bioconda::samtools=1.19 conda-forge::matplotlib=3.8.0'

    publishDir "${params.outdir}/coverage", mode: 'copy'

    input:
    tuple val(meta), path(bam), path(bai)

    output:
    path "${meta.id}_depth.txt", emit: depth
    path "${meta.id}_coverage.png", emit: plot

    script:
    """
    samtools depth -aa ${bam} > ${meta.id}_depth.txt
    python3 -c "import matplotlib.pyplot as plt;d=[int(l.split()[2]) for l in open('${meta.id}_depth.txt')];plt.figure(figsize=(12,4));plt.plot(d,c='darkred',lw=.5);plt.title('${meta.id}');plt.xlabel('Position');plt.ylabel('Depth');plt.grid(True,alpha=.6);plt.tight_layout();plt.savefig('${meta.id}_coverage.png',dpi=300)"
    """
}

workflow map_and_plot {
    take:
        reference_genome
        trimmed_reads

    main:
        map_in = (params.ref_type == 'alignment')
            ? trimmed_reads.combine(reference_genome)
            : trimmed_reads.join(reference_genome)
        alignments = MAP_READS(map_in)
        PLOT_COVERAGE(alignments.map { meta, bam, bai, ref -> tuple(meta, bam, bai) })

    emit:
        bam_out = alignments.map { meta, bam, bai, ref -> tuple(meta, bam, bai) }
        ref_out = alignments.map { meta, bam, bai, ref -> ref }
}

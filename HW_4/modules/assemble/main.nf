process MEGAHIT_ASSEMBLE {
    tag "${meta.id}"
    label 'assemble'
    conda 'bioconda::megahit=1.2.9 bioconda::seqtk=1.5'

    publishDir "${params.outdir}/assembly", mode: 'copy'

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("${meta.id}_contigs.fasta"), emit: contigs

    script:
    if (params.assembly_fraction < 1.0) {
        """
        BIN=\$(dirname \$(which megahit))
        ln -sf megahit_core "\$BIN/megahit_core_popcnt"

        seqtk sample -s42 ${reads[0]} ${params.assembly_fraction} | gzip > sub_R1.fq.gz
        seqtk sample -s42 ${reads[1]} ${params.assembly_fraction} | gzip > sub_R2.fq.gz
        megahit -1 sub_R1.fq.gz -2 sub_R2.fq.gz -o out -t ${task.cpus}
        cp out/final.contigs.fa ${meta.id}_contigs.fasta
        """
    } else {
        """
        BIN=\$(dirname \$(which megahit))
        ln -sf megahit_core "\$BIN/megahit_core_popcnt"

        megahit -1 ${reads[0]} -2 ${reads[1]} -o out -t ${task.cpus}
        cp out/final.contigs.fa ${meta.id}_contigs.fasta
        """
    }
}

workflow assemble {
    take:
        trimmed_reads

    main:
        MEGAHIT_ASSEMBLE(trimmed_reads)

    emit:
        contigs = MEGAHIT_ASSEMBLE.out.contigs
}

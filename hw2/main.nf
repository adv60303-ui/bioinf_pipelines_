#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

params.reads               = null
params.sra_id              = null
params.reference           = null
params.outdir              = 'results'
params.assembly_fraction   = 0.05

process FETCH_SRA {
    conda 'bioconda::sra-tools=3.1.0'
    publishDir "${params.outdir}/raw_reads", mode: 'copy'
    input:  val sra_id
    output: tuple val(sra_id), path('*.fastq.gz')
    script:
    """
    fasterq-dump --threads ${task.cpus} --split-files --outdir . ${sra_id}
    gzip -f *.fastq
    """
}

process FASTQC_RAW {
    conda 'bioconda::fastqc=0.12.1'
    publishDir "${params.outdir}/qc_raw", mode: 'copy'
    input:  tuple val(sample_id), path(reads)
    output: path '*_fastqc.*'
    script:
    """
    fastqc -t ${task.cpus} -o . ${reads.join(' ')}
    """
}

process FASTQC_TRIM {
    conda 'bioconda::fastqc=0.12.1'
    publishDir "${params.outdir}/qc_trimmed", mode: 'copy'
    input:  tuple val(sample_id), path(reads)
    output: path '*_fastqc.*'
    script:
    """
    fastqc -t ${task.cpus} -o . ${reads.join(' ')}
    """
}

process TRIMMOMATIC {
    conda 'bioconda::trimmomatic=0.39'
    publishDir "${params.outdir}/trimmed_reads", mode: 'copy', pattern: '*_p.fastq.gz'
    input:  tuple val(sample_id), path(reads)
    output: tuple val(sample_id), path("${sample_id}_R{1,2}_p.fastq.gz")
    script:
    """
    ADP=\$(find \${CONDA_PREFIX:-/opt/conda} -name TruSeq3-PE.fa | head -1)
    trimmomatic PE -threads ${task.cpus} ${reads[0]} ${reads[1]} \\
        ${sample_id}_R1_p.fastq.gz /dev/null ${sample_id}_R2_p.fastq.gz /dev/null \\
        ILLUMINACLIP:\$ADP:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:36
    """
}

process SPADES {
    label 'assemble'
    conda 'bioconda::megahit=1.2.9 bioconda::seqtk=1.5'
    publishDir "${params.outdir}/assembly", mode: 'copy'
    input:  tuple val(sample_id), path(reads)
    output: tuple val(sample_id), path("${sample_id}_contigs.fasta")
    script:
    if (params.assembly_fraction < 1.0) {
        """
        BIN=\$(dirname \$(which megahit))
        ln -sf megahit_core "\$BIN/megahit_core_popcnt"

        seqtk sample -s42 ${reads[0]} ${params.assembly_fraction} | gzip > sub_R1.fq.gz
        seqtk sample -s42 ${reads[1]} ${params.assembly_fraction} | gzip > sub_R2.fq.gz
        megahit -1 sub_R1.fq.gz -2 sub_R2.fq.gz -o out -t ${task.cpus}
        cp out/final.contigs.fa ${sample_id}_contigs.fasta
        """
    } else {
        """
        BIN=\$(dirname \$(which megahit))
        ln -sf megahit_core "\$BIN/megahit_core_popcnt"

        megahit -1 ${reads[0]} -2 ${reads[1]} -o out -t ${task.cpus}
        cp out/final.contigs.fa ${sample_id}_contigs.fasta
        """
    }
}

process MAP_READS {
    conda 'bioconda::bwa=0.7.17 bioconda::samtools=1.19'
    publishDir "${params.outdir}/mapping", mode: 'copy'
    input:  tuple val(sample_id), path(reads), path(reference)
    output: tuple val(sample_id), path("${sample_id}.sorted.bam"), path("${sample_id}.sorted.bam.bai")
    script:
    """
    cp ${reference} ref.fa && bwa index ref.fa
    bwa mem -t ${task.cpus} ref.fa ${reads[0]} ${reads[1]} | samtools sort -@ ${task.cpus} -o ${sample_id}.sorted.bam
    samtools index ${sample_id}.sorted.bam
    """
}

process PLOT_COVERAGE {
    conda 'bioconda::samtools=1.19 conda-forge::matplotlib=3.8.0'
    publishDir "${params.outdir}/coverage", mode: 'copy'
    input:  tuple val(sample_id), path(bam), path(bai)
    output:
        path "${sample_id}_depth.txt"
        path "${sample_id}_coverage.png"
    script:
    """
    samtools depth -aa ${bam} > ${sample_id}_depth.txt
    python3 -c "import matplotlib.pyplot as plt;d=[int(l.split()[2]) for l in open('${sample_id}_depth.txt')];plt.figure(figsize=(12,4));plt.plot(d,c='darkred',lw=.5);plt.title('${sample_id}');plt.xlabel('Position');plt.ylabel('Depth');plt.grid(True,alpha=.6);plt.tight_layout();plt.savefig('${sample_id}_coverage.png',dpi=300)"
    """
}

workflow {
    if (params.reads && params.sra_id)
        error 'Use either --reads or --sra_id, not both'
    if (!params.reads && !params.sra_id)
        error "Provide --reads 'data/*_{1,2}.fastq.gz' or --sra_id DRR030302"

    reads_ch = params.sra_id
        ? FETCH_SRA(Channel.of(params.sra_id))
        : Channel.fromFilePairs(params.reads, checkIfExists: true)

    FASTQC_RAW(reads_ch)
    trimmed = TRIMMOMATIC(reads_ch)
    FASTQC_TRIM(trimmed)

    map_in = params.reference
        ? trimmed.combine(Channel.fromPath(params.reference, checkIfExists: true).first())
        : trimmed.join(SPADES(trimmed))

    PLOT_COVERAGE(MAP_READS(map_in))
}

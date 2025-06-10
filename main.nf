#!/usr/bin/env nextflow

// Workflow to 
// 1) Generate genome index (optional)
// 2) Align reads to a reference genome (index)
// Quick start: 
// Start with "nextflow run main.nf" (using default parameters and example input files) (nextflow needs to be installed)
// Attention: Docker and nextflow need to installed
// Custom start:
// Adjust parameters in nextflow.config!

nextflow.enable.dsl=2

// Process to download the reference genome if a URL is provided
process DownloadGenome {
    container "curlimages/curl:8.13.0"

    input:
    val refGenomePathOrUrl // Path to local file or URL

    output:
    path "refGenome.fa.gz"

    script:
    """
    echo "Downloading reference genome from: curl -v -k -L -o refGenome.fa.gz '${refGenomePathOrUrl}"
    curl -v -k -L -o refGenome.fa.gz '${refGenomePathOrUrl}'
    """
}

process UnzipReference {
    container "actions/gzip"

    input:
    path refGenome // Downloaded genome file

    output:
    path "refGenome.fa"

    script:
    """
    gunzip -c ${refGenome} > refGenome.fa
    """
}

// Process to generate the STAR genome index
process STARIndex {
    container 'quay.io/biocontainers/star:2.6.1d--0'
    maxForks 1 // STAR index generation is resource-intensive, limit to one fork
    cpus 8

    input:
    path refGenome // Downloaded or local reference genome

    output:
    path "genome_index"

    script:
    """
    mkdir -p genome_index
    STAR --runThreadN 8 --runMode genomeGenerate --genomeDir genome_index --genomeFastaFiles ${refGenome}
    """
}

// Process to align reads using STAR
process AlignReads {
    container 'quay.io/biocontainers/star:2.6.1d--0'
    publishDir "${params.outdir}/STAR_Alignments", mode: 'copy'

    input:
    path genome_index
    tuple val(sampleID), path(samples)

    output:
    path "${sampleID}_Aligned.out.bam"

    script:
    """
    STAR --genomeDir genome_index \
         --readFilesIn ${samples[0]} ${samples[1]} \
         --readFilesCommand zcat \
         --outFileNamePrefix ${sampleID}_ \
         --outSAMtype BAM Unsorted
    """
}

// Workflow definition
workflow {
    // Load sample information
    samples = Channel.fromPath(params.sampleSheet)
                 .splitCsv(header: true)
                 .map { row -> 
                     tuple(
                         row.sampleID,
                         [file(row.fastq1), file(row.fastq2)] // Collect R1 and R2 files in a list
                     )
                 }
    samples.view()

    // Define workflow order
    refGenome = DownloadGenome(params.refGenomePathOrUrl)
    refGenome = UnzipReference(refGenome)
    genomeIndex = STARIndex(refGenome)
    AlignReads(genomeIndex, samples)
}

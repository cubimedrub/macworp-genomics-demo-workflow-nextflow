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
    path "refGenome.fa"

    script:
    """
    if [[ ${refGenomePathOrUrl} =~ ^http ]]; then
        echo "Downloading reference genome from URL: ${refGenomePathOrUrl}"
        curl -L -o refGenome.fa.gz ${refGenomePathOrUrl}
        echo "Decompressing the downloaded file..."
        gunzip -c refGenome.fa.gz > refGenome.fa
    else
        echo "Using local reference genome file: ${refGenomePathOrUrl}"
        if [[ ${refGenomePathOrUrl} == *.gz ]]; then
            echo "Decompressing the local gzipped file..."
            gunzip -c ${refGenomePathOrUrl} > refGenome.fa
        else
            echo "Copying the local reference genome file..."
            cp ${refGenomePathOrUrl} refGenome.fa
        fi
    fi
    """
}

// Process to generate the STAR genome index
process STARIndex {
    input:
    path refGenome // Downloaded or local reference genome

    output:
    path "genome_index"

    script:
    """
    mkdir -p genome_index
    STAR --runMode genomeGenerate --genomeDir genome_index --genomeFastaFiles ${refGenome}
    """
}

// Process to align reads using STAR
process AlignReads {
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
    genomeIndex = STARIndex(refGenome)
    AlignReads(genomeIndex, samples)
}

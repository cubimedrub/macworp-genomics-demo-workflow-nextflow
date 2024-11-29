#!/usr/bin/env nextflow

//Workflow to 
// 1) generate  genome index (optional)
// 2) align reads to a reference genome (index)
// Quick start: 
// Start with "nextflow run main.nf" (using default parameters and example input files) (nextflow needs to be installed)
// 
// Custom start:
// Adjust parameters in nexflow.config!




nextflow.enable.dsl=2




// Process for indexing the genome with STAR
process STARIndex {
    input:
    path refGenome

    output:
    path "genome_index"

    script:
    """
    mkdir -p genome_index
    STAR --runMode genomeGenerate --genomeDir genome_index --genomeFastaFiles $refGenome
    """
}

// Process for aligning FASTQ reads with STAR
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

// Process for quality control with MultiQC
/*process MultiQC {
    publishDir "${params.outdir}/MultiQC", mode: 'copy'

    input:
    path "${params.outdir}/STAR_Alignments/*_Aligned.out.bam"

    output:
    path "multiqc_report.html"

    script:
    """
    multiqc ${params.outdir}/STAR_Alignments -o .
    """
}*/

// Read sample sheet and start workflow
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

    // Define the workflow order
    STARIndex(params.refGenome) 
    AlignReads(STARIndex.out,samples)
    //MultiQC(AlignReads.out)

    
}
# macworp-genomics-demo-workflow-nextflow
## Workflow
Tasks included in this workflow:

1) Generate genome index (optional)
Input: local or download link for reference genome fasta file
Tool: Curl (local installation), STAR (workflow downloads dockerized STAR)
Output: reference genome index
2) Align reads to a reference genome
Input: reference genome index, fastq files, sample sheet
Tool: STAR (workflow downloads dockerized STAR)
Output: Aligned reads from fastq files as .bam files
## How to run
Quick start: 
1) Navigate to folder with "main.nf"
2) Start with "nextflow run main.nf" (using default parameters and example input files) 
// Attention: nextflow, docker and curl need to be installed
Custom start:
1) Adjust parameters in nextflow.config!
2) Start with "nextflow run main.nf" (using default parameters and example input files) 

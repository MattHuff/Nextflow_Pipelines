# Nextflow_Pipelines
A set of pipelines with which you are able to run multiple processes at a time.

This repository will contain key Nextflow pipelines I have developed.
When you enter a directory for a specific pipeline, there will be a subdirectory for the server the pipeline was developed under.
Within the server subdirectory will be subdirectories for each version of the pipeline.
Within the version subdirectory, you may find further subdirectories depending on if you are running code on paired-end reads or single-end reads.

## RNA Sequencing Pipeline
This pipeline will run trimming, alignment, and counting for all fastq files in your input directory.
You will need to index your reference genome prior to running STAR.
The RNA Sequencing Pipeline is available on both Dr. Staton's server and the ACF. It may be run either as a SE pipeline or a PE pipeline.
Note that this version of the pipeline is currently designed for fastq files downloaded from SRA using SRAToolKits.

## BLAST/IPS Gene Annotation
This pipeline will schedule two BLAST jobs - one comparing input genes to the Swissprot database, the other to TrEMBL's plant database - and one InterProScan job.
This is currently present only on the ACF.

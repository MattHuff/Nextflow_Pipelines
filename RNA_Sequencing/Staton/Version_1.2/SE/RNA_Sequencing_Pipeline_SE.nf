#!/usr/bin/env nextflow

import java.nio.channels.FileLock
import java.nio.channels.FileChannel
import java.nio.channels.OverlappingFileLockException;

params.input = null
params.output = null
params.ref_genome = null

query_path = "${params.input}/*_1.fastq*"
if ( query_path.isEmpty() )
{
 error "ERROR: No valid fastq files are present in the specified input directory! Please ensure you entered the correct directory.\
 \nMissing input files"
}
query_data = Channel.fromPath("${query_path}")
     .map { file -> tuple(file.baseName.replaceAll(/_1/,""), file) }

//Step 1: Trim Input

process TrimInput {
 publishDir "$params.output/1_skewer", mode: 'copy'

 input:
 set sampleID, file(R) from query_data
 
 output:
 set sampleID, file("${sampleID}-trimmed.fastq") into trimmed_reads

 script:
 """
 /staton/software/skewer/skewer -x /staton/software/Trimmomatic-0.39/adapters/all.fa -l 30 ${R} -o ${sampleID} >& ${sampleID}.trim_output
 """
}

ref_genome_path = "$params.ref_genome/*.fa"
ref_genome_data = Channel.fromPath("${ref_genome_path}")

gff_path = "$params.ref_genome/*.gff*"
gff_data = Channel.fromPath("${gff_path}")

//Step 2: Align Reads to Genome (STAR)

process STAR {
 publishDir "$params.output/2_star", mode: 'copy'
 
 input:
 set sampleID, file(R) from trimmed_reads

 output:
 set sampleID, file("${sampleID}.Aligned.sortedByCoord.out.bam") into alignments

 script:
 """
 /staton/software/STAR-2.7.3a/STAR --genomeDir ${params.ref_genome}/genome_index --readFilesIn ${R} --runThreadN 2 --outFileNamePrefix ${sampleID}. --outSAMtype BAM SortedByCoordinate
 """
}

//Step 3: Get HTSeq counts

if ( gff_path.isEmpty() )
{
 error "ERROR: A valid GFF3 file is missing! Please make sure a valid GFF3 file exists in '${params.ref_genome}' ()"
}

gff_htseq =  Channel.fromPath("${gff_path}")

process htseq {
 publishDir "$params.output/3_htseq", mode: 'copy'
	
 input:
 set sampleID, file(BAM) from alignments
	
 output:
 set sampleID, file("${sampleID}.counts.txt") into raw_counts
	
 script:
 """ 
 /staton/software/htseq/scripts/htseq-count --format=bam --stranded=no --type=mRNA --idattr=ID ${BAM} ${gff_path} > ${sampleID}.counts.txt
 """
}

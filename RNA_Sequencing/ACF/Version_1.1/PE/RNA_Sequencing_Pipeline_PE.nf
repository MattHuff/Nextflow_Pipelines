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
     .map { file -> tuple(file.baseName.replaceAll(/_1/,""), file, file.getParent()/file.getName().replaceAll(/_1/,"_2")) }
 
//Step 1: Trim Input

process TrimInput {
 publishDir "${params.output}/1_skewer", mode: 'copy' 

 input:
 set sampleID, file(R1), file(R2) from query_data
 
 output:
 set sampleID, file("${sampleID}-trimmed-pair1.fastq"), file("${sampleID}-trimmed-pair2.fastq") into trimmed_reads

 script:
 """
 set +eu
 module load skewer
 skewer -x /lustre/haven/proj/UTK0032/software/Trimmomatic-0.36/adapters/adapters.fa -l 30 ${R1} ${R2} -o ${sampleID} >& ${sampleID}.trim_output
 """
}

ref_genome_path = "${params.ref_genome}/*.fa*"
ref_genome_data = Channel.fromPath("${ref_genome_path}")

gff_path = "${params.ref_genome}/*.gff*"
gff_data = Channel.fromPath("${gff_path}")

index_path = file("${params.ref_genome}/genome_index/")

if ( !index_path.isDirectory() )
{
 error "ERROR: A valid index folder does not exist for your genome of interest! Please make sure the following directory exists:\
 \n'${index_path}'"
}

//Step 2: Align Reads to Genome (STAR)

process STAR {
 publishDir "${params.output}/2_star", mode: 'copy'

 input:
 set sampleID, file(R1), file(R2) from trimmed_reads 

 output:
 set sampleID, file("${sampleID}.Aligned.sortedByCoord.out.bam") into alignments

 script:
 """
 set +eu
 module load star 
 STAR --genomeDir ${params.ref_genome}/genome_index --readFilesIn ${R1} ${R2} --runThreadN 2 --outFileNamePrefix ${sampleID}. --outSAMtype BAM SortedByCoordinate --limitBAMsortRAM 17365764332
 """
}

//Step 3: Get HTSeq counts

if ( gff_path.isEmpty() )
{
 error "ERROR: A valid GFF3 file is missing! Please make sure a valid GFF3 file exists in '${params.ref_genome}' ()"
}

gff_htseq =  Channel.fromPath("${gff_path}")

process htseq {
 publishDir "${params.output}/3_htseq", mode: 'copy'
	
 input:
 set sampleID, file(BAM) from alignments
	
 output:
 set sampleID, file("${sampleID}.counts.txt") into raw_counts
	
 script:
 """ 
 set +eu
 module load anaconda3/4.4.0
 . /sw/acf/anaconda3/4.4.0/centos7.3_gnu6.3.0/anaconda3-4.4.0/etc/profile.d/conda.sh
 conda activate /lustre/haven/proj/UTK0032/projects/undergrads/gene_expression_projects/upython 
 htseq-count --format=bam --stranded=no --type=mRNA --idattr=ID ${BAM} ${gff_path} > ${sampleID}.counts.txt
 """
}

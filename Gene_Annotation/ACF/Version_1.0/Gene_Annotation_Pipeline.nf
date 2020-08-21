#!/usr/bin/env nextflow

import java.nio.channels.FileLock
import java.nio.channels.FileChannel
import java.nio.channels.OverlappingFileLockException;

params.input = null
params.output = null

input_path = "${params.input}/*.fa*"
if ( input_path.isEmpty() )
{
 error "ERROR: No valid fasta files are present in the specified input directory! Please ensure you entered the correct directory.\
 \nMissing input files"
}

sprot_data = Channel.fromPath("${input_path}")
     .map { file -> tuple(file.baseName, file) }

trembl_data = Channel.fromPath("${input_path}")
     .map { file -> tuple(file.baseName, file) }

IPS_data = Channel.fromPath("${input_path}")
     .map { file -> tuple(file.baseName, file) }

// Step 1: BLAST - Swissprot

process swissprot {
 publishDir "$params.output/blast/swissprot", mode: 'copy'

 input:
 set fastaID, file(inputFasta) from sprot_data

 output:
 set fastaID, file("${fastaID}_swissprot.tsv") into swissprot_results

 script:
 """
 module load blast
 blastp -query ${inputFasta} -db /lustre/haven/proj/UTK0032/library/uniprot/uniprot_sprot.fasta -out ${fastaID}_swissprot.tsv -evalue 1e-5 -outfmt 6 
 """
}

// Step 2: BLAST - TrEMBL

process trembl {
 publishDir "$params.output/blast/trembl", mode: 'copy'

 input:
 set fastaID, file(inputFasta) from trembl_data

 output:
 set fastaID, file("${fastaID}_trembl.tsv") into trembl_results

 script:
 """
 module load blast
 blastp -query ${inputFasta} -db /lustre/haven/proj/UTK0032/library/uniprot/uniprot_trembl_plants_July_2018.fasta -out ${fastaID}_trembl.tsv -evalue 1e-5 -outfmt 6
 """
}

// Step 3: InterProScan

process IPS {
 publishDir "$params.output/IPS/tsvs", mode: 'copy'

 input:
 set fastaID, file(inputFasta) from IPS_data

 //output:
 //set fastaID, file("${fastaID}_IPS.tsv") into IPS_results

 script:
 """
 module load python3
 /lustre/haven/proj/UTK0032/software/interproscan-5.34-73.0/interproscan.sh -i ${inputFasta} -f TSV -o ${fastaID}_IPS.tsv --disable-precalc --iprlookup --goterms --pathways -T ${params.output}/IPS/TMP
 """
}

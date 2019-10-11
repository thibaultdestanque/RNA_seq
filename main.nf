#!/usr/bin/env nextflow




def helpMessage() {
  log.info"""

  =============================
    RNA-seq pipeline nextflow 
  =============================
  Usage:
  The typical command for running the pipeline is as follows:
  nextflow run main.nf

  Tree folders should look like this: 
    DATA_TO_ANALYSE_FOLDER                  # Folder
        '=> main.nf                         # Script nextflow
        '=> nextflow.config                 # config nextflow
        '=> Ref_Genome                      # Folder
            '=> genome.fa                   # Reference genome in fasta
            '=> genomeIndex                 # Folder which will store index reference genome
                '=> genome index files      # Files
            '=> genome_annotation.gff3      # Annotation of genome in gff3 
        '=> Fastq                           # Folder
            '=> fastq files to analyse      # Files to analysis in fastq
        '=> Other_files                     # Folder
            '=> adapter file                # adapter file



  """.stripIndent()
}

/*
 * SET UP CONFIGURATION VARIABLES
 */
 
// Show help message
params.help = false
if (params.help){
  helpMessage()
  exit 0
}




/*
 * Start RNA seq pipeline nextflow
 */
log.info ""
log.info ""
log.info "        ================================================="
log.info "         R N A  S E Q   P I P E L I N E  N E X T F L O W "
log.info "        ================================================="
log.info ""
log.info ""
log.info "        Parameters used:"
log.info "            - fastq             : ${params.reads}"
log.info "            - Genome path       : ${params.genome_path}"
log.info "            - Genome fasta      : ${params.genome_fa}"
log.info "            - Genome Index      : ${params.genomeIndex}"
log.info "            - Genome Annotation : ${params.annot}"
log.info ""
log.info "        Results will be store in:"
log.info "            ${params.outdir}"
log.info ""
log.info "        Process are stored in:"
log.info "            work/"
log.info ""
//log.info "        Scratch: ${params.scratch_path}"
log.info ""
log.info ""
log.info "        ==================="
log.info "          Launch pipeline  "
log.info "        ==================="
log.info ""


/*
 * The reference genome file
 */
genome_fa_file  = file(params.genome_fa)
annotation_file = file(params.annot)
adapter_file    = file(params.adapter_file)
reads           = file(params.reads)



/*
 * Create the `read_pairs` channel that emits tuples containing three elements:
 * the pair ID, the first read-pair file and the second read-pair file 
 */
Channel
    .fromFilePairs( params.reads )
    .ifEmpty { exit 1, "Cannot find any reads matching: ${params.reads}" }
    .set { read_pairs } 


 
/*
 * Step 1. Trim sequence with trimmomatic
 */
process Trim {
	tag "$pair_id"
	publishDir params.outdir, mode: 'copy'

	time'2h'
    cpus 8
    queue 'omp'
    memory '60 GB'
    echo true
    scratch '/home1/scratch/tdestanq/'
    conda 'bioconda::trimmomatic=0.36'

    input:
    set pair_id, file(reads) from read_pairs
    file(adapter_file) from adapter_file
     
    output:
    file "*_paired.fastq.gz" into paired_read_trimed
    file "*_unpaired.fastq.gz" into unpaired_read_trimed


    shell:
    """
    trimmomatic PE -threads 8 -phred33 ${reads} \
    	${pair_id}_trim_R1_paired.fastq.gz \
    	${pair_id}_trim_R1_unpaired.fastq.gz \
    	${pair_id}_trim_R2_paired.fastq.gz \
    	${pair_id}_trim_R2_unpaired.fastq.gz \
    	ILLUMINACLIP:${adapter_file}:${params.illuminaclip_1}:${params.illuminaclip_2}:${params.illuminaclip_3} \
    	LEADING:${params.leading} TRAILING:${params.trailing} SLIDINGWINDOW:${params.slidingwindows_1}:${params.slidingwindows_2} MINLEN:${params.minlen} \
        >& /home1/scratch/tdestanq/Trim.log 2>&1
    """
}

/*
 * Step 2. Index reference genome
 */
process Index_Genome {

    time'24h'
    cpus 1
    queue 'sequentiel'
    memory '30 GB'
    echo true
    scratch '/home1/scratch/tdestanq/'
    conda 'bioconda::gmap=2018.07.04'
    
    input:
    file genome_fa from genome_fa_file
    params.genome_path
    params.genomeIndex

    output:
    file "genome.fa" into genome

    shell:
    """
    gmap_build --dir=${params.genome_path} ${genome_fa} -d ${params.genomeIndex} >& /home1/scratch/tdestanq/Index_Genome.log 2>&1
    """
}

/*
 * Step 3. Align reads on reference genome
 */
process Alignment {
	publishDir params.outdir, mode: 'copy'

	time'23h'
    cpus 16
    queue 'omp'
    memory '60 GB'
    echo true
    scratch '/home1/scratch/tdestanq/'
    conda 'bioconda::gmap=2018.07.04'
    
    input:
    file genome from genome
    file read_trimed from paired_read_trimed

    output:
    file "*" into all_file_alignment
    file "*.concordant_uniq" into concordant_uniq_files

    shell:
    """
    bash
    #source activate Gmap
    # Retrieve base name from samples
    ls *_trim_R1_paired.fastq.gz | sed 's/_trim_R1_paired.fastq.gz//g' > Name_tmp.txt;
    base_name=`cat Name_tmp.txt`;
    gsnap --gunzip -t 8 -A sam --min-coverage=${params.min_coverage} \
        --dir=${params.genome_path} -d ${params.genomeIndex}  \
        --max-mismatches=${params.max_mismatches} --novelsplicing=${params.novel_splicing} \
        --split-output=\$base_name \
        --read-group-id=\$base_name \
        --read-group-platform=Illumina \
        ${read_trimed} >& /home1/scratch/tdestanq/Alignment.log 2>&1
    rm Name_tmp.txt
    """
}

/*
 * Step 4. Format data
 */

process Format_Data {
    publishDir params.outdir, mode: 'copy'

    time'5h'
    cpus 1
    queue 'sequentiel'
    memory '60 GB'
    scratch '/home1/scratch/tdestanq/'
    //conda 'samtools=1.9'
    
    input:
    file concordant_uniq from concordant_uniq_files

    output:
    file "*" into all_format_data
    file "*.sorted.bam" into read_mapped_sort_bam

    shell:
    """
    . /appli/bioinfo/samtools/1.9/env.sh
    samtools view -b ${concordant_uniq} >& ${concordant_uniq}.bam 2> /home1/scratch/tdestanq/Format_Data_sam_view.log ;
    samtools sort ${concordant_uniq}.bam >& ${concordant_uniq}.sorted.bam 2> /home1/scratch/tdestanq/Format_Data_sam_sort.log ;
    samtools index ${concordant_uniq}.sorted.bam >& /home1/scratch/tdestanq/Format_Data_sam_sort_index.log 2>&1 ;
    """
}



/*
 * Step 5. Count read with Htseq
 */

process Htseq_count {
    publishDir params.outdir, mode: 'copy'

    time'4h'
    cpus 1
    queue 'sequentiel'
    memory '50 GB'
    echo true
    scratch '/home1/scratch/tdestanq/'
    conda 'bioconda::htseq=0.6.1'
    
    input:
    file sorted_bam from read_mapped_sort_bam
    file GFF3_annotation from annotation_file

    output:
    file "*_htseq_count.txt" into htseq_count

    shell:
    """
    #. /appli/bioinfo/htseq-count/latest/env.sh
    htseq-count -f "bam" -s "no" -r "pos" -t "gene" -i "Name" --mode "union" ${sorted_bam} ${GFF3_annotation} >& ${sorted_bam}_htseq_count.txt 2> /home1/scratch/tdestanq/HtseqCount.log
    """
    // .${params.add_header} >& /home1/scratch/tdestanq/Format_Data_for_R_add_header.log 2>&1

}


process Format_data_for_R_1 {
    publishDir params.outdir, mode: 'copy'

    time'4h'
    cpus 1
    queue 'sequentiel'
    memory '50 GB'
    echo true
    //conda 'bioconda::samtools=1.9'
    scratch '/home1/scratch/tdestanq/'
    
    input:
    file htseq_count from htseq_count

    output:
    file "*.trim" into format_data_for_R_1

    shell:
    """
    #echo "${htseq_count}"
    echo -e "genes\t\$(basename ${htseq_count})" | cat - ${htseq_count} >& ${htseq_count}.temp 2> /home1/scratch/tdestanq/Format_Data_for_R_add_header.log
    mv ${htseq_count}.temp ${htseq_count}.trim
    """
}

process Format_data_for_R_2 {
    publishDir params.outdir, mode: 'copy'

    time'4h'
    cpus 1
    queue 'sequentiel'
    memory '50 GB'
    echo true
    //conda 'bioconda::samtools=1.9'

    scratch '/home1/scratch/tdestanq/'
    
    input:
    file htseq_count_header from format_data_for_R_1.collect()

    output:
    file "join_devlarve.txt" into all_file

    shell:
    """
    join ${htseq_count_header} >& join_devlarve.txt 2> /home1/scratch/tdestanq/Format_Data_for_R_join_multiplefile.log
    """
}





workflow.onComplete {
	log.info ( workflow.success ? "\n\nDone! Results files are stocked in --> $params.outdir ;)\n" : "Oops... something went wrong" )
}

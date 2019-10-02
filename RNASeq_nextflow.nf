#!/usr/bin/env nextflow





/*
 * Defines some parameters in order to specify the refence genomes
 * and read pairs by using the command line options
 */
params.reads = "*_{R1,R2}.fastq.gz"
params.annot = "/home1/datawork/tdestanq/Nextflow_workflow/genome_Pmarg_v2.gff3"

params.genome_path = "/home1/datawork/tdestanq/Nextflow_workflow/"
params.genome_fa = "/home1/datawork/tdestanq/Nextflow_workflow/genome.fa"
params.genomeIndex = "genomeIndex"

params.adapter_file = "/home1/datawork/tdestanq/00_ressources/univec/univec.fa"
params.outdir = 'results'
reads = file(params.reads)


/*
 * Start RNA seq pipeline nextflow
 */
log.info """\n
        \n
         ===============================================
         R N A  S E Q   P I P E L I N E  N E X T F L O W  
         ===============================================
         Genome     : ${params.genome_fa}
         Annotation : ${params.annot}
         Results    : ${params.outdir}
         \n
         """
         .stripIndent()


/*
 * the reference genome file
 */
genome_fa_file = file(params.genome_fa)
annotation_file = file(params.annot)
adapter_file = file(params.adapter_file)



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

    input:
    set pair_id, file(reads) from read_pairs
    file(adapter_file) from adapter_file
     
    output:
    file "*_paired.fastq.gz" into paired_read_trimed
    file "*_unpaired.fastq.gz" into unpaired_read_trimed


    shell:
    """
    . /appli/bioinfo/trimmomatic/latest/env.sh
    trimmomatic PE -threads 8 -phred33 ${reads} \\
    	${pair_id}_trim_R1_paired.fastq.gz \\
    	${pair_id}_trim_R1_unpaired.fastq.gz \\
    	${pair_id}_trim_R2_paired.fastq.gz \\
    	${pair_id}_trim_R2_unpaired.fastq.gz \\
    	ILLUMINACLIP:${adapter_file}:2:20:7 \\
    	LEADING:20 TRAILING:20 SLIDINGWINDOW:30:30 MINLEN:60 \\
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
    //conda 'Gmap'
    echo true
    scratch '/home1/scratch/tdestanq/'
    
    input:
    file genome_fa from genome_fa_file
    params.genome_path
    params.genomeIndex

    output:
    file "genome.fa" into genome

    shell:
    """
    source activate Gmap
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
    //conda 'Gmap'
    echo true
    scratch '/home1/scratch/tdestanq/'
    
    input:
    file genome from genome
    file read_trimed from paired_read_trimed

    output:
    file "*" into all_file_alignment
    file "*.concordant_uniq" into concordant_uniq_files

    shell:
    """
    bash
    source activate Gmap
    # Retrieve base name from samples
    ls *_trim_R1_paired.fastq.gz | sed 's/_trim_R1_paired.fastq.gz//g' > Name_tmp.txt;
    base_name=`cat Name_tmp.txt`;
    echo \$base_name;
    gsnap --gunzip -t 8 -A sam --min-coverage=0.9 \
        --dir=${params.genome_path} -d ${params.genomeIndex}  \
        --max-mismatches=2 --novelsplicing=1 \
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

//process Format_Data_index {
//    publishDir params.outdir, mode: 'copy'
//
//    time'5h'
//    cpus 1
//    queue 'sequentiel'
//    memory '60 GB'
//    scratch '/home1/scratch/tdestanq/'
//    
//    input:
//    file sort_bam from read_mapped_sort_bam
//
//    output:
//    file "*" into bam_index
//
//    shell:
//    """
//    . /appli/bioinfo/samtools/1.9/env.sh
//    """
//}


/*
 * Step 5. Count read with Htseq
 */

process Htseq_count {
    publishDir params.outdir, mode: 'copy'

    time'4h'
    cpus 1
    queue 'sequentiel'
    memory '50 GB'
    //conda 'Gmap'
    echo true
    scratch '/home1/scratch/tdestanq/'
    
    input:
    file sorted_bam from read_mapped_sort_bam
    file GFF3_annotation from annotation_file

    output:
    file "*_htseq_count.txt" into all_file

    shell:
    """
    . /appli/bioinfo/htseq-count/latest/env.sh
    htseq-count -f "bam" -s "no" -r "pos" -t "gene" -i "Name" --mode "union" ${sorted_bam} ${GFF3_annotation} > ${sorted_bam}_htseq_count.txt >& /home1/scratch/tdestanq/HtseqCount.log 2>&1
    """
}





workflow.onComplete {
	log.info ( workflow.success ? "\n\nDone! Results files are stocked in --> $params.outdir ;)\n" : "Oops... something went wrong" )
}

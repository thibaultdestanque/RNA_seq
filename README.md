# RNAseq pipeline analysis


# Prerequisites

## Clone git hub directory
```
git clone https://github.com/thibaultdestanque/RNA_seq_Nextflow.git
```

## Import univec
```
wget Other_files/univec.fasta ftp://ftp.ncbi.nlm.nih.gov/pub/UniVec/UniVec
```

## Import Genome
```
cp Path/to/your/genome/genome.fa Folder/where/you/clone/github/directory/Ref_Genome/genome.fa
```

## Import Fastq
```
cp Path/to/your/fastq/files/*fastq Folder/where/you/clone/github/directory/Fastq/*fastq
```

## Launch pipeline
```
nextflow run main.nf
```

# Documentation

## Steps:
  - Trimmomatic
  - Gsnap index genome
  - Gsnap alignment
  - Samtools format data
  - Htseq count





# Code
```
My code
```

# highligt
`test'

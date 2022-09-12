

# Assignment 4 [*30 marks*]

**Due before 5pm, Friday 16th September**

## Assignment overview [*5 marks for general format*]

This assignment is based on the `Alignmen/NGS` practicals. 

Your answers to all questions should be submitted to myUni as a `.zip` file containing three bash scripts (see below for details), and a single file containing the answers to the questions 3,4 and 5 (can be in any readable format, but plain text or PDF are preferred). [*1 mark*]

The `.zip` filename must start with your student number [*1 mark*] and *your bash script must be able to run without errors*. This means it should run regardless of where it is located in the directory tree.
Meaningful comments are strongly advised [*2 marks*]

For all scripts, please use the directory `~/Assignment4` as the parent directory for all downloads and analysis.
**You will be expected to hard code this into your script.**
All other folders must be created as subdirectories of this. It will be simplest for you to stick with the structure in NGS_skeleton [*1 mark*]


## Practical questions [*20 marks*] 

### Scripts and associated questions

1. Script 1 must do the following:
    + Download the genomic sequence (i.e. fasta file) and annotation (i.e. gff file) of the model plant _Arabidopsis thaliana_ to the subdirectory `Refs/Athaliana` from the Ensembl ftp directory (link below) [*2 marks*]
    + Identify how many chromosomes are in the genome, and write this information to standard output (`stdout`) [*2 marks*]
    + Identify how many unique genes are located in the genome, and send this information to `stdout`? [*2 marks*]
    + Create a genome index for use with `bwa`[*1 marks*] 
2. Script 2 must do the following: download the sequencing data contained at the following links (using `curl` or `wget`) to the directory `01_rawData/fastq` (note that the filename in the Box link is not informative; **make sure you rename your files as below**()): [*2 marks*]
    - [SRR5882797_10M_1.fastq.gz](https://universityofadelaide.box.com/shared/static/egl3n16r0ziaxlvbs9074xqd1liktnuz.gz)
    - [SRR5882797_10M_2.fastq.gz](https://universityofadelaide.box.com/shared/static/g2ly4kzz1blus5juy426i37zl45o38pu.gz)
    + Trim your data for poor quality bases and remove adapters using cutadapt. Write your output to the directory `02_trimmedData/fastq` [*3 marks*]
    + Align paired-end reads to the genome index using `bwa mem`, resulting in a single `.bam` file in the directory `03_alignedData/bam`. Ensure there are no intermediary `.sam` files saved. [*3 marks*]
    + Sort and index your bam file [*2 marks*] 
3. Script 3 must export the following information as a tab-delimited file to the parent directory. Include the category as the first column, and the numbers as the second column
    + How many reads were mapped in total? [*1 marks*]
    + How many reads were mapped as a pair? [*1 marks*]
    + How many reads were mapped as a "proper" pair? [*1 marks*]


## Theory questions [*5 marks*]

4. What is the difference between an index and a barcode? [*2 marks*]
5. What is the difference between a SAM file and a BAM file? [*1 marks*]
6. What does the CIGAR string "26M2I78M" mean? [*2 marks*]




### Resources

- The genome (fasta) for the model plant _Arabidopsis thaliana_ can be found [here: http://ftp.ebi.ac.uk/ensemblgenomes/pub/release-51/plants/fasta/arabidopsis_thaliana/dna/](http://ftp.ebi.ac.uk/ensemblgenomes/pub/release-51/plants/fasta/arabidopsis_thaliana/dna/) and the annotation (gff3) files can be found [here: http://ftp.ebi.ac.uk/ensemblgenomes/pub/release-51/plants/gff3/arabidopsis_thaliana/](http://ftp.ebi.ac.uk/ensemblgenomes/pub/release-51/plants/gff3/arabidopsis_thaliana/). There are multiple files in these directories, so take care to choose the full genome sequence and full genome gff files.
- The genome file should be "Arabidopsis_thaliana.TAIR10.dna.toplevel.fa.gz" and not contain an "rm" or "sm" in the filename. These refer to the genome being "repeat-masked" i.e. all the repetitive elements have been ignored, or "soft-masked" were the repeats and low-complexity regions are in lower-case letters (acgt not ACGT)

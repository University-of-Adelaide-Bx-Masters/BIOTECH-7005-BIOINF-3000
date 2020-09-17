# Assignment 5 - [*n marks*]

**Due before Tuesday 20th October**

Your answers to all questions should be submitted to myUni as a `.zip` file containing three bash scripts, and a single file containing the answers to the questions 3,4 and 5 (can be in any readable format). [*1 marks*]

The `.zip` filename must start with your student number [*1 marks*] and your bash script must be able to run without errors.
Meaningful comments are strongly advised [*1 mark*]

For all scripts, please use the directory `~/Assignment5` as the parent directory for all downloads and analysis.
**You will be expected to hard code this into your script.**
All other folders must be created as subdirectories of this. [*1 mark*]

## Practical questions [*25 marks*]

This practical will use human exome sequence data from [BioProject_659398]. These are 200nt single end reads

1. Install the SRA Toolkit using conda `conda install -c bioconda sra-tools` this will give you the `fastq-dump` tool that you will use to get the data. (no marks for this)
    
2. Write a script to: 
    + download the fastqc read data using sra-tools to download the dataset you will use from the Short Read Archive (SRA) at NCBI with the following command: `fastq-dump --gzip SRR12538042 -O ~/Assignment5/0_rawData/fastq/` [1 marks]
    + download the genome reference and gff3 files for human chromosome 1 to ~/Assignment5/Ref/ (see `https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.28_GRCh38.p13/GCA_000001405.28_GRCh38.p13_assembly_structure/Primary_Assembly/assembled_chromosomes/FASTA/chr1.fna.gz` and `https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.39_GRCh38.p13/GCF_000001405.39_GRCh38.p13_genomic.gff.gz` )
    + run fastqc and trim the reads as for previous practicals/assignment (note that this dataset is NOT paired end, so no need to use the `-t 2` flag). You will only need one adapter, so the command for trimming will be: `cutadapt -m 35 -q 20 -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA  -o 1_trimmedData/fastq/SRR12538042.fastq.gz  0_rawData/fastq/SRR12538042.fastq.gz > 1_trimmedData/log/cutadapt.log` [1 mark]
    + run fastqc on the trimmed reads [1 mark]
    + 
## Theoretical questions [*5 marks*]


[BioProject_659398]:https://www.ncbi.nlm.nih.gov//bioproject/659398

### Genes For Question 3
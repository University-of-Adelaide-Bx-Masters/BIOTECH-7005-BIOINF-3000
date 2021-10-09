# Assignment 6 - Transcriptome analysis [*35 marks*]

**Due before 5pm, Tuesday 19th October - extended deadline**

Your answers to all questions should be submitted to myUni as a `.zip` file containing **four** bash scripts, and a `5_final_assembly` folder including files required in **step 5**. [*1 mark*]

The `.zip` filename must start with your student number and your bash script must be able to run without errors.
Meaningful comments are strongly advised [*1 mark*]

For all scripts, please use the directory `~/Assignment6` as the parent directory for all downloads and analysis. 
**You will be expected to hard code this into your scripts.**
Use an organised folder structure to store files generated in this assignment. [*1 mark*]

```
Assignment6/
├── data
├── results
│   ├── 1_QC
│   ├── 2_clean_data
│   ├── 3_denovo_assembly
│   ├── 4_genome_guided_assembly
│   └── 5_final_assembly
└── scripts
```

## Practical questions [*28 marks* in total]

__Step 1. Write the first script to:__

+ Download the genomic sequence (i.e. fasta file) and annotation (i.e. gff file) of the model plant _Arabidopsis thaliana_ to the subdirectory `DB` from the Ensembl `ftp` directory as you did in assignment 4 [*1 mark*]
+ Download the BUSCO "brassicales" lineage dataset using this [link](https://busco-data.ezlab.org/v5/data/lineages/brassicales_odb10.2020-08-05.tar.gz). Store the lineage dataset in the `DB` subdirectory and decompress and untar it [*1 mark*] 
+ Build the genome index for GMAP [*2 marks*]
+ Build the genome index for STAR [*2 marks*]

__Step 2. Write the second script to copy the sequencing data to your `data` directory  and carry out QC__. The data are in `~/data/Transcriptomics_data_2021/Assignment/`. [*1 mark*] Include the following steps:

+ Get QC report for raw reads [*1 mark*]
+ Trim adaptor and low-quality sequences using cutadapt [*1 mark*]
+ Get QC report for clean reads [*1 mark*]

__Step 3. Write the third script to:__

+ Do _de novo_ transcriptome assembly for clean reads from chromosome 4 using `Trinity` [*3 marks*]
+ Get the genome coordinates file (GFF3) of _de novo_ assembled transcripts using `GMAP` [*1 mark*]
+ Assess the _de novo_ assembly quality using `BUSCO` [*2 marks*]

__Step 4. Write the fourth script to:__

+ Do genome guided transcriptome assembly for clean reads from chromosome 4 with STAR and StringTie [*3 marks*]
+ Get a `fasta` file of genome guided assembled transcripts using `gffread` [*1 mark*]
+ Assess the genome guided assembly quality using `BUSCO` [*2 marks*]

__Step 5. Organise the final results:__

+ Copy QC reports for raw reads (two html files) and clean reads (two html files) into `5_final_assembly` [*1 mark*]
+ Copy _de novo_ assembly results (including one fasta file and one GFF3 file) into `5_final_assembly` [*1 mark*]
+ Copy genome guided assembly results (including one GTF file and one fasta file) into `5_final_assembly` [*1 mark*]

__Step 6. Identify one assembled transcript including alternative splicing events with read evidence__ (Take a screenshot from IGV and put it in the `5_final_assembly` subdirectory) [*3 marks*].

## Theoretical questions [*4 marks* in total]

1. Can we get all genes assembled from one transcriptome sequenced from a sample collected from one particular tissue at one particular developmental stage? [*1 mark*]
2. Give the reason(s) supporting your answer in quesion 1 [*3 marks*]

### Resources (also in assignment 4)

+ The genome (fasta) for the model plant _Arabidopsis thaliana_ can be found [here: http://ftp.ebi.ac.uk/ensemblgenomes/pub/release-51/plants/fasta/arabidopsis_thaliana/dna/](http://ftp.ebi.ac.uk/ensemblgenomes/pub/release-51/plants/fasta/arabidopsis_thaliana/dna/) and the annotation (gff3) files can be found [here: http://ftp.ebi.ac.uk/ensemblgenomes/pub/release-51/plants/gff3/arabidopsis_thaliana/](http://ftp.ebi.ac.uk/ensemblgenomes/pub/release-51/plants/gff3/arabidopsis_thaliana/). There are multiple files in these directories, so take care to choose the full genome sequence and full genome gff files.
+ The genome file should be "Arabidopsis_thaliana.TAIR10.dna.toplevel.fa.gz" and not contain an "rm" or "sm" in the filename. These refer to the genome being "repeat-masked" i.e. all the repetitive elements have been ignored, or "soft-masked" were the repeats and low-complexity regions are in lower-case letters (acgt not ACGT)

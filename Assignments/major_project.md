# Major Project (20%)

## Introduction

In this course, the following next-generation sequencing (NGS) datasets/protocols are described in detail:

- Transcriptome Sequencing (RNA-Seq) - We provide a dataset for this as a Major Project. 

Analysis of the major project dataset requires quality control (quality and sequencing adapter trimming), genome alignment, feature counts, differential expression and downstream visualisation. The data also aim to address a particular scientific question and investigate a hypothesis. For the major project, you will be assigned a dataset (from a published project) and complete all the analysis tasks (from raw data to final results) and write up a report. This report must be structured like a small journal article, with abstract (summarising the project), introduction (background on the study and identification of the research hypothesis), methods (analysis steps and programs used), results (what you found) and discussion (how the results relate to the research question) sections. Marks will also be awarded to the bash/R or RMarkdown scripts that you use.

|Section                    |Mark |
|:--------------------------|:----|
|Abstract                   |5%   |
|Introduction + hypothesis  |10%  |
|Methods                    |20%  |
|Results & Discussion       |30%  |
|References                 |5%   |
|Analysis scripts           |30%  |

## Transcriptome/RNAseq data set

Your major project RNA-Seq data is from [a gene expression study](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE155018). The study (Kumar et al., 2021) was trying to identify candidate genes that contribute to a long circadian period in the model organism *Drosophila melanogaster* by performing transcriptional profiling of whole fly heads from two genotypes: DGRP_892 (DGRP892), a line with a long 31-hour circadian period and Canton-S B (CantonSB), a line with a normal 24-hour circadian period. If you want to learn more about the study, you can read their published paper [1](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8390771/). This RNA-Seq data was sequenced on Illumina HiSeq4000 (PE76). 


For your major project, you will compare a pair of transcriptome data between male or female long circadian line DGRP892 and normal circadian line CantonSB at a specific circadian time (every 2 hours of a 24 hours circadian time period, from CT02 to CT24). See details in the **Materials and Methods** section of the reference paper [1](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8390771/) if you want to get more information about the experimental design. Your research goal will be to investigate the transcriptional difference between your designated pair of Drasophila samples and interpret your results. 

### Description of your unique transcriptome dataset

Each student will be provided with a unique pair of transcriptome dataset (12 zipped `fastq` files), including 6 files representing the RNA-Seq of three biological replicates from the long circadian line DGRP892 and 6 files representing the RNA-seq of corresponding normal circadian line CantonSB. For example, the following is a typical list of files for one student:

| Files (paired end)                                                                 | File format | Line (genotype) | Tissue | Circadian time | Sex  | Description of files                                                                               |
|------------------------------------------------------------------------------------|-------------|-----------------|--------|----------------|------|----------------------------------------------------------------------------------------------------|
| DGRP892_Head_CT02_male_rep1_R1.fastq.gz, DGRP892_Head_CT02_male_rep1_R2.fastq.gz   | fastq       | DGRP_892        | Head   | CT02           | male | Biological replicate 1 (rep1) of long circadian line DGRP_892 at 2 hours circadian time (CT02)     |
| DGRP892_Head_CT02_male_rep2_R1.fastq.gz, DGRP892_Head_CT02_male_rep2_R2.fastq.gz   | fastq       | DGRP_892        | Head   | CT02           | male | Biological replicate 2 (rep2) of long circadian line DGRP_892 at 2 hours circadian time (CT02)     |
| DGRP892_Head_CT02_male_rep3_R1.fastq.gz, DGRP892_Head_CT02_male_rep3_R2.fastq.gz   | fastq       | DGRP_892        | Head   | CT02           | male | Biological replicate 3 (rep3) of long circadian line DGRP_892 at 2 hours circadian time (CT02)     |
| CantonSB_Head_CT02_male_rep1_R1.fastq.gz, CantonSB_Head_CT02_male_rep1_R2.fastq.gz | fastq       | Canton-S B      | Head   | CT02           | male | Biological replicate 1 (rep1) of normal circadian line Canton-S B at 2 hours circadian time (CT02) |
| CantonSB_Head_CT02_male_rep2_R1.fastq.gz, CantonSB_Head_CT02_male_rep2_R2.fastq.gz | fastq       | Canton-S B      | Head   | CT02           | male | Biological replicate 2 (rep2) of normal circadian line Canton-S B at 2 hours circadian time (CT02) |
| CantonSB_Head_CT02_male_rep3_R1.fastq.gz, CantonSB_Head_CT02_male_rep3_R2.fastq.gz | fastq       | Canton-S B      | Head   | CT02           | male | Biological replicate 3 (rep3) of normal circadian line Canton-S B at 2 hours circadian time (CT02) |

### [**IMPORTANT, PLEASE READ**] Access/use your unique transcriptome dataset 

Your unique transcriptome dataset can be accessed through VM path `/shared/data/major_project_2024/01_raw_data/axxxxxxx`, in which `axxxxxxx` is your university `a` number. You have two options to access/use your unique transcriptome dataset:

- Option 1, create soft links pointing to the original raw data files

For example, if your major project top level directory is `~/major_project`, and you want to access your raw fastq data in directory `~/major_project/01_raw_data` in your analysis, you can run the following commands (replace `axxxxxxx` with your `a` number) to create soft links of your unique transcriptome raw data files:

```
mkdir -p ~/major_project/01_raw_data
cd ~/major_project/01_raw_data
ln -s /shared/data/major_project_2024/01_raw_data/axxxxxxx/*.fastq.gz ./
```

- Option 2, directly access the raw data

You can hard code your unique transcriptome dataset path (`/shared/data/major_project_2024/01_raw_data/axxxxxxx`, replace `axxxxxxx` with your `a` number) in your bash script.

**[Important] Please DO NOT copy the transcriptome raw data files into your home directory** If everyone did this, we would run out of storage.

### Additional information
The major project is mainly based on the practicals of `Week 6, Alignment/NGS` and `Week 11, Transciptomics/Gene Expression`, and also requires skills that you have learned from other practicals (R, bash).

### Some useful tips

1. Set up your project folder structure in an organised way.

2. We have provided the Drosophila reference genome files in shared folder `/shared/data/major_project_2024/00_DB`. There are three files in the folder:
+ `Drosophila_melanogaster.BDGP6.46.dna.primary_assembly.no_scaffolds.fa`, this is the reference genome sequences file in `fasta` format.
+ `Drosophila_melanogaster.BDGP6.46.112.chr.gtf`, this is the reference gene annotation file in `gtf` format.
+ `Drosophila_melanogaster.ensDB.112.gene_annotation.txt`, this is a text file including additional information about reference genes.

3. You may want to include figures/tables in your report, however, please don't try to include every figure/table from your analysis in your report. The general rule for figures is no more than 6-8 figures (you can have multiple sub-panels in one figure) in a scientific article.

4. There is no hard word limit for your report. However, please remember sometimes **Less is more**.

### Reference
[1] Kumar S, Tunc I, Tansey TR, Pirooznia M, Harbison ST. Identification of Genes Contributing to a Long Circadian Period in Drosophila Melanogaster. J Biol Rhythms. 2021 Jun;36(3):239-253. doi: 10.1177/0748730420975946. Epub 2020 Dec 4. PMID: 33274675; PMCID: PMC8390771. 

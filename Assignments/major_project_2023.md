# Major Project (15%)

## Introduction

In this course, the following next-generation sequencing (NGS) datasets/protocols are described in detail:

- Transcriptome Sequencing (RNA-Seq) - We provide a dataset for this as a Major Project. 

Analysis of the major project dataset requires quality control (quality and sequencing adapter trimming), genome alignment, feature counts, differential expression and downstream visualisation. The data also aim to address a particular scientific question and investigate a hypothesis. For the major project, you will be assigned a dataset (from a published project) and complete all the analysis tasks (from raw data to final results) and write up a report. This report must be structured like a small journal article, with abstract (summarising the project), introduction (background on the study and identification of the research hypothesis), methods (analysis steps and programs used), results (what you found) and discusson (how the results relate to the research question) sections. Marks will also be awarded to the bash/R or RMarkdown scripts that you use.

|Section                    |Mark |
|:--------------------------|:----|
|Abstract                   |5%   |
|Introduction + hypothesis  |10%  |
|Methods                    |20%  |
|Results & Discussion       |30%  |
|References                 |5%   |
|Analysis scripts           |30%  |

## Transcriptome/RNAseq data set

Your major project RNA-Seq data is from [a gene expression study](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE182617). The study (Yin et al., 2023) was trying to understand the roles of plant hormones on the regulation of plant growth and development using Arabidopsis as a plant model. If you want to learn more about the study, you can read their preprint paper [1](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC10028877/). This RNA-Seq data was sequenced on Illumina HiSeq2500 (single end 95 bp). 

For your major project, you will compare Arabidopsis transcriptome data under the treatment of plant hormone (Salicylic Acid (SA), Brassinosteroid (BR) or Strigolactone/karrikin (SL_KAR)) for a specific time period (0.25, 0.5, 1, 2, 4, 8, 12 or 24 hours) with the control transcriptome data (treated with plant hormone for 0 hour). Your research goal will be to understand what kind of genes/pathways are potentially regulated by the plant hormone. 

### Description of your unique transcriptome dataset

Each student will be provided a unique pair of transcriptome dataset (6 zipped `fastq` files), including 3 files representing the RNA-Seq of three biological replicates under hormone treatment and 3 files representing the RNA-seq of three corresponding control biological replicates. For example, the following is a typical list of files for one student:

| File name                          | File format | Group     | Description of file                                                                                                     |
|------------------------------------|-------------|-----------|-------------------------------------------------------------------------------------------------------------------------|
| Col0_SA_treated_025h_rep1.fastq.gz | fastq       | Treatment | Biological replicate 1 (rep1) of Arabidopsis Col-0 ecotype (Col0) treated with Salicylic Acid (SA) for 0.25 hours (025h) |
| Col0_SA_treated_025h_rep2.fastq.gz | fastq       | Treatment | Biological replicate 2 (rep2) of Arabidopsis Col-0 ecotype (Col0) treated with Salicylic Acid (SA) for 0.25 hours (025h) |
| Col0_SA_treated_025h_rep3.fastq.gz | fastq       | Treatment | Biological replicate 3 (rep3) of Arabidopsis Col-0 ecotype (Col0) treated with Salicylic Acid (SA) for 0.25 hours (025h) |
| Col0_SA_treated_0h_rep1.fastq.gz   | fastq       | Control   | Biological replicate 1 (rep1) of Arabidopsis Col-0 ecotype (Col0) treated with Salicylic Acid (SA) for 0 hours (0h)    |
| Col0_SA_treated_0h_rep2.fastq.gz   | fastq       | Control   | Biological replicate 2 (rep2) of Arabidopsis Col-0 ecotype (Col0) treated with Salicylic Acid (SA) for 0 hours (0h)    |
| Col0_SA_treated_0h_rep3.fastq.gz   | fastq       | Control   | Biological replicate 3 (rep3) of Arabidopsis Col-0 ecotype (Col0) treated with Salicylic Acid (SA) for 0 hours (0h)    |

### [**IMPORTANT, PLEASE READ**] Access/use your unique transcriptome dataset 

Your unique transcriptome dataset can be accessed through VM path `/shared/data/major_project_2023/axxxxxxx`, in which `axxxxxxx` is your university `a` number. You have two options to access/use your unique transcriptome dataset:

- Option 1, create soft links pointing to the original raw data files

For example, if your major project top level directory is `~/major_project`, and you want to access your raw fastq data in directory `~/major_project/01_raw_data` in your analysis, you can run the following commands (replace `axxxxxxx` with your `a` number) to create soft links of your unique transcriptome raw data files:

```
mkdir -p ~/major_project/01_raw_data
cd ~/major_project/01_raw_data
ln -s /shared/data/major_project_2023/axxxxxxx/*.fastq.gz ./
```

- Option 2, directly access the raw data

You can hard code your unique transcriptome dataset path (`/shared/data/major_project_2023/axxxxxxx`, replace `axxxxxxx` with your `a` number) in your bash script.

**[Important] Please DO NOT copy the transcriptome raw data files into your home directory** If everyone did this, we would run out of storage again.

### Additional information
The major project is mainly based on the practicals of `Week 6, Alignment/NGS` and `Week 11, Transciptomics/Gene Expression`, and also requires skills that you have learned from other practicals (R, bash).

### Reference
[1] Yin L, Zander M, Huang SC, Xie M, Song L, Guzm N JPS, Hann E, Shanbhag BK, Ng S, Jain S, Janssen BJ, Clark NM, Walley JW, Beddoe T, Bar-Joseph Z, Lewsey MG, Ecker JR. Transcription Factor Dynamics in Cross-Regulation of Plant Hormone Signaling Pathways. bioRxiv [Preprint]. 2023 Mar 9:2023.03.07.531630. doi: 10.1101/2023.03.07.531630. PMID: 36945593; PMCID: PMC10028877.

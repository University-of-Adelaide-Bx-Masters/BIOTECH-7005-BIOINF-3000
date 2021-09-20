# Major Project (15%)

In this course, the following next-generation sequencing (NGS) datasets/protocols are described in detail:

- Whole genome sequencing/Resequencing - We provide a dataset for this as a Major Project.
- Transcriptome Sequencing (RNAseq) - We provide a dataset for this as a Major Project. 

Each of these NGS approaches uses similar programs and analysis approaches, such as quality control (quality and sequencing adapter trimming), genome alignment, and downstream visualisation and statistical methods. They also aim to address a particular scientific question and investigate a scientific question. For the major project, you will take a published dataset and complete all the analysis tasks (from raw data to final results) and write up a report. This report must be structured like a small journal article, with abstract (summarising the project), introduction (background on the study and identification of the research hypothesis), methods (analysis steps and programs used), results (what you found) and discusson (how the results relate to the research question) sections. Marks will also be awarded to the bash/R or RMarkdown scripts that you use.

|Section                    |Mark |
|:--------------------------|:----|
|Abstract                   |5%   |
|Introduction + hypothesis  |10%  |
|Methods                    |20%  |
|Results & Discussion       |30%  |
|References                 |5%   |
|Analysis scripts           |30%  |

**You have the freedom to choose any dataset from any research article you would like**, however you need to let Dave know beforehand.

If you cannot find a suitable dataset, we have provided the following datasets:

- Resequencing/variant calling data set. This is from a whole genome resequencing of a female [Belgian Tervuren Dog](./images/sage.jpg) and for a female [Shiba Inu dog](./images/Taka_Shiba.jpg). You can see how these breeds are related on [this tree](https://research.nhgri.nih.gov/dog_genome/downloads/studies-figure1_032017.pdf). We have restricted the data to reads coming from Chromosomes 1 and X in order to ensure that your VMs will be able to  process the data. **You will need to request the dataset from Dave**. You will need to carry out QC on the sequencing data as you have done before in the practicals and then map the QC'ed reads to the supplied reference sequences.  Once you have mapped the reads, you will need to call variants  with respect to the reference genome as you did in the practicals. It is worth bearing in mind that the reference genome is from a female Boxer dog called Tasha. For the analysis, you will need to identify the common/shared genomic positions for variants in both breeds and for specific genomic positions of variants in each breed. Your research question is *What proportion of variant nucleotide positions are breed specific and what proportion are shared between breeds*. 

- Transcriptome/RNAseq data set. This data is from [a gene expression study](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA503894). We will supply you with fastq files that you can use to generate assembled RNA transcripts as in the transcriptomic tutorial and practical. **You will need to request the dataset from Dave**. You are expected to carry QC on the raw fastq files in order to generate clean data for transcript assembly. You will assemble the transcripts with [trinity](https://github.com/trinityrnaseq/trinityrnaseq/wiki) and then determine how many complete protein coding sequences you have assembled by comparison to *Swissprot* and the quality of your transcript assemblies using [BUSCO](https://pubmed.ncbi.nlm.nih.gov/26059717/). You will compare the assemblies you will generate using different amounts of input data to answer your research questions *What is the effect of increasing read numbers on the quality of Trinity assembly?*. 


**Determine the adapter sequences needed for trimming for the transcriptome data set. It looks like these should work:
IDT for Illumina–TruSeq DNA and RNA UD Indexes
These unique dual (UD) index adapters are arranged in the plate to enforce the recommended pairing  strategy.

Adapter Trimming
The following sequences are used for adapter trimming.Read 1
AGATCGGAAGAGCACACGTCTGAACTCCAGTCA
Read 2
AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT
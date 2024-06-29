# Major Project (20%)

In this course, the following next-generation sequencing (NGS) datasets/protocols are described in detail:

- Transcriptome Sequencing (RNAseq) - We provide a dataset for this as a Major Project. 

Analysis of the major project dataset requires quality control (quality and sequencing adapter trimming), genome alignment, feature counts, differential expression and downstream visualisation and statistical methods. The data also aim to address a particular scientific question and investigate a hypothesis. For the major project, you will be assigned a dataset (from a published project of ours) and complete all the analysis tasks (from raw data to final results) and write up a report. This report must be structured like a small journal article, with abstract (summarising the project), introduction (background on the study and identification of the research hypothesis), methods (analysis steps and programs used), results (what you found) and discusson (how the results relate to the research question) sections. Marks will also be awarded to the bash/R or RMarkdown scripts that you use.

|Section                    |Mark |
|:--------------------------|:----|
|Abstract                   |5%   |
|Introduction + hypothesis  |10%  |
|Methods                    |20%  |
|Results & Discussion       |30%  |
|References                 |5%   |
|Analysis scripts           |30%  |

## Transcriptome/RNAseq data set

This data is from [a gene expression study](). We will supply you with `fastq` files that you can use to generate assembled RNA transcripts as in the transcriptomic tutorial and practical. **You will need to request the dataset from Dave**. You are expected to carry out QC on the raw fastq files in order to generate clean data for transcript mapping. You will map the transcripts using [STAR](https://github.com/alexdobin/STAR) and then generate raw counts. You will then normalise the counts and remove noise (minimally expressed transcripts) and carry out differential gene expression analysis. 



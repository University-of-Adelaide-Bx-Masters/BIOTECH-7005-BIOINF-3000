# Major Project (15%)

In this course, the following next-generation sequencing (NGS) datasets/protocols are described in detail:

- Whole genome sequencing/Resequencing
- Transcriptome Sequencing (RNAseq) - We provide a dataset for this as a Major Project, see below. 

Each of these NGS approaches uses similar programs and analysis approaches, such as quality control (quality and sequencing adapter trimming), genome alignment, and downstream visualisation and statistical methods. They also aim to address a particular scientific question and investigate a scientific hypothesis. For the major project, you will take a published dataset and complete all the analysis tasks (from raw data to final results) and write up a report. This report must be structured like a small journal article, with abstract (summarising the project), introduction (background on the study and identification of the research hypothesis), methods (analysis steps and programs used), results (what you found) and discusson (how the results relate to the research hypothesis) sections. Marks will also be awarded to the bash/R or RMarkdown scripts that you use.

|Section                    |Mark |
|:--------------------------|:----|
|Abstract                   |5%   |
|Introduction + hypothesis  |10%  |
|Methods                    |20%  |
|Results & Discussion       |30%  |
|References                 |5%   |
|Analysis scripts           |30%  |

**You have the freedom to choose any dataset from any research article you would like**, however you need to let Dave know beforehand.

If you cannot find a suitable dataset, we have provided datasets from two plant RNAseq profiling studies, one has multiple mutants of the histone deacetylase gene [`hda`](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4848314/), which is involved in the regulation of flowering time in *Arabidopsis thaliana*. The details of the sequencing experiment [are found at this GEO link](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE78946). The other study looks at fatty acid (FA) synthesis and has [three mutations that affect FA synthesis](https://www.pnas.org/content/111/3/1204). The information for the sequencing experiment [are found at this GEO link](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE53952). To ensure that everyone is not working on the same data, each student should work on a separate sample group:

|ID      |Group |
|:-------|:-----|
|a1792295|hda9-1 vs wild-type (arabidopsis dataset)|
|a1731521|fae1 at 7-8 vs 9-10 (fatty acid synthesis dataset)|
|a1794835|PDAT at 7-8 vs 9-10 (fatty acid synthesis dataset)|
|a1627307|fae1 at 11-12 vs 9-10 (fatty acid synthesis dataset)|
|a1704339|CL37 vs fae1 at 11-12 (fatty acid synthesis dataset)|
|a1776797|PDAT at 11-12 vs 9-10 (fatty acid synthesis dataset)|
|a1819113|PDAT at 7-8 vs 9-10 (fatty acid synthesis dataset)|
|a1812900|fae1 at 7-8 vs 9-10 (fatty acid synthesis dataset)|
|a1724529|CL37 vs fae1 at 9-10 (fatty acid synthesis dataset)|
|a1805535|PDAT at 11-12 vs 7-8 (fatty acid synthesis dataset)|
|a1811637|PDAT at 11-12 vs 7-8 (fatty acid synthesis dataset)|
|a1820124|CL37 vs PDAT at 11-12 (fatty acid synthesis dataset)|
|a1783527|PDAT at 11-12 vs 7-8 (fatty acid synthesis dataset)|
|a1778718|CL37 vs fae1 at 11-12 (fatty acid synthesis dataset)|
|a1795646|CL37 vs PDAT at 11-12 (fatty acid synthesis dataset)|
|a1681781|hda9-1 vs wild-type (arabidopsis dataset)|
|a1686683|CL37 vs PDAT at 11-12 (fatty acid synthesis dataset)|
|a1787867|CL37 vs PDAT at 7-8 (fatty acid synthesis dataset)|
|a1789951|hda5-1 vs wild-type (arabidopsis dataset)|
|a1754994|hda5-1 vs wild-type (arabidopsis dataset)|
|a1720273|CL37 at 11-12 vs 9-10 (fatty acid synthesis dataset)|
|a1786706|CL37 vs fae1 at 11-12 (fatty acid synthesis dataset)|
|a1792442|hda9-1 vs wild-type (arabidopsis dataset)|
|a1758223|CL37 at 7-8 vs 9-10 (fatty acid synthesis dataset)|
|a1791929|CL37 vs PDAT at 11-12 (fatty acid synthesis dataset)|
|a1786514|CL37 vs PDAT at 7-8 (fatty acid synthesis dataset)|
|a1786921|CL37 at 7-8 vs 9-10 (fatty acid synthesis dataset)|
|a1680671|CL37 vs PDAT at 9-10 (fatty acid synthesis dataset)|
|a1817954|hda6-6 vs wild-type (arabidopsis dataset)|
|a1805954|CL37 vs fae1 at 9-10 (fatty acid synthesis dataset)|
|a1811380|PDAT vs fae1 at 11-12 (fatty acid synthesis dataset)|


For this particular dataset, we expect you to run a differential expression analysis between the replicates of your sample group against the Col wildtype control sample. You will report differentially expressed genes in your results section and discuss how this relates to the study's experimental hypothesis.

# Updated instructions

You will have seen how to proceed from the point where you have count files in the [DE tutorial](https://university-of-adelaide-bx-masters.github.io/BIOTECH-7005/DE_gene_tutorial/Tutorial_DE_Genes.html). However you will first have to run QC on the reads (including adapter trimming).

Note that these datasets are all single read as opposed to paired end read. For adapter trimming please see  [Cutadapt documentation](https://cutadapt.readthedocs.io/en/stable/guide).html#basic-usage .  For the adapter sequence that you will need, please see [Illumina adapter documentation](https://support.illumina.com/content/dam/illumina-support/documents/documentation/chemistry_documentation/experiment-design/illumina-adapter-sequences-1000000002694-14.pdf).

After you carry out adapter trimming, you will need to map the reads as you have done in Assignments 4 and 5. Hint you will use BWA. 

Once you have mapped the reads, you will need to summarise or quantify the reads and create count files from the `.bam` files. For this you will need the `Rsubread` package for `R`. 

You should read the [Rsubread user guide](https://bioconductor.org/packages/release/bioc/vignettes/Rsubread/inst/doc/SubreadUsersGuide.pdf) to determine the syntax and arguments needed to create a counts file. Hint, see Chapter 10 of the user guide for installation instructions and use case examples. 

Once you have `count` files, you can proceed as shown in the `DE tutorial`. 


Raw FASTQ files will be provided via a data link provided in a myuni announcement/email and in the MyUni Major Project assignment page.

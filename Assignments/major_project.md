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

## Resequencing/variant calling data set 

This is from the [dog 1000 genome project](https://www.ncbi.nlm.nih.gov/bioproject/188158) for a number of dog breeds. We will provide sequence reads for a female [Belgian Tervuren Dog](./images/sage.jpg) and for a female [Akita Inu dog](https://en.wikipedia.org/wiki/Akita_(dog)#/media/File:Akita_Inu_dog.jpg) and for other breeds depending on demand. You can see how these breeds are related on [this tree](https://research.nhgri.nih.gov/dog_genome/downloads/studies-figure1_032017.pdf). We have restricted the data to reads coming from Chromosomes 1 and X in order to ensure that your VMs will be able to  process the data. **You will need to request the dataset from Dave**. We will provide you with download links for the data from two dog breeds for your analysis. You will need to carry out QC on the sequencing data as you have done before in the practicals and then map the QC'ed reads to the supplied reference sequences.  Once you have mapped the reads, you will need to call variants  with respect to the reference genome as you did in the practicals. It is worth bearing in mind that the reference genome is from a female Boxer dog called [Tasha](https://www.broadinstitute.org/files/news/stories/full/tasha-12072005.jpg). For the analysis, you will need to identify the common/shared genomic positions for variants in both breeds and for specific genomic positions of variants in each breed. Your principal research question is *What proportion of variant nucleotide positions are breed specific and what proportion are shared between breeds*. 

Additional information regarding annotating variants using `SnpEff`. `SnpEff` is a variant annotation and effect prediction tool. It annotates and predicts the effects of genetic variants (such as amino acid changes). The detailed documentation can be found [here](http://pcingola.github.io/SnpEff/se_introduction/). `SnpEff` can be installed via `conda` in VM (Assume you want to install it under `3000` conda environment):

```bash
conda activate 3000
conda install -c bioconda snpeff
```

`SnpEff` requires databases to produce the annotations. Dog database has been pre-built in `SnpEff`. To find available dog databases in `SnpEff`, run following command:

```bash
snpEff databases | grep -i familiaris
```

You will get a list of dog genome assemblies and annotations as follows:

```
Basenji_breed-1.1.99                                            Canis_lupus_familiarisbasenji                                                                                   https://snpeff.blob.core.windows.net/databases/v5_0/snpEff_v5_0_Basenji_breed-1.1.99.zip
CanFam3.1.75                                                    Canis_familiaris                                                                                                https://snpeff.blob.core.windows.net/databases/v5_0/snpEff_v5_0_CanFam3.1.75.zip
CanFam3.1.99                                                    Canis_familiaris                                                                                                https://snpeff.blob.core.windows.net/databases/v5_0/snpEff_v5_0_CanFam3.1.99.zip
UMICH_Zoey_3.1.99                                               Canis_lupus_familiarisgreatdane                                                                                 https://snpeff.blob.core.windows.net/databases/v5_0/snpEff_v5_0_UMICH_Zoey_3.1.99.zip
```
The list shows all available dog databases in `snpEff` at the time of writing this. It shows the database name, genome name and source data (where was the genome reference data obtained from). In this project, if you are using the ENSEMBL dog reference genome `CanFam3.1` (https://m.ensembl.org/Canis_lupus_familiaris/Info/Annotation). You should use the lastest version of annotation `CanFam3.1.99`, in which `99` is the version number of annotation release. 

## Transcriptome/RNAseq data set
This data is from [a gene expression study](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA503894). We will supply you with fastq files that you can use to generate assembled RNA transcripts as in the transcriptomic tutorial and practical. **You will need to request the dataset from Dave**. You are expected to carry out QC on the raw fastq files in order to generate clean data for transcript assembly. You will denovo assemble the transcripts with [trinity](https://github.com/trinityrnaseq/trinityrnaseq/wiki) and genome-guided assemble the transcripts using [STAR](https://github.com/alexdobin/STAR) and [StringTie](https://github.com/gpertea/stringtie) then determine how many complete protein coding sequences you have assembled by comparison to *Swissprot* and the quality of your transcript assemblies using [BUSCO](https://pubmed.ncbi.nlm.nih.gov/26059717/). You will compare the assemblies you will generate using different amounts of input data to answer your research question *What is the effect of increasing read numbers on the quality of Trinity assembly?*. 

Additional information regarding determing how many complete protein coding sequences you have assembled. You can follow this documentation (https://github.com/trinityrnaseq/trinityrnaseq/wiki/Counting-Full-Length-Trinity-Transcripts) to count the full length transcripts. SwissProt protein database can be downloaded from this [link](https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz). In your VM, if `Trinity` has been installed in `3000` conda environment, the `analyze_blastPlus_topHit_coverage.pl` script can be directly invoked from shell under your `3000` conda environment.

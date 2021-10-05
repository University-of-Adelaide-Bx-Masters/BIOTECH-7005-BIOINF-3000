# Assignment 5 - [*23 marks total*]

**Due before Tuesday 12th October**
## Correct submission [*3 Marks*]
Your answers to all questions should be submitted to myUni as a `.zip` file containing your two bash scripts and your html `knitted` R markdown file for part 4 that contains the answers to the additional questions (can be in any readable format, but my order of preference is; html, pdf, Rmd). [*1 mark for submitting as a single zip file*]

The `.zip` filename must start with your student number [*1 mark*] and your bash script must be able to run without errors.
Meaningful comments are strongly advised.

For all scripts, please use the directory `~/Assignment5` as the parent directory for all downloads and analysis. Use the following directory structure [*1 mark*]: 
```
.
├── Assignment5
   ├── 0_rawData
   │   ├── fastq
   │   └── FastQC
   ├── 1_trimmedData
   │   ├── fastq
   │   ├── FastQC
   │   └── log
   ├── 2_alignedData
   │   ├── bam
   │   ├── FastQC
   │   ├── log
   │   ├── sorted_bam
   │   └── vcf
   ├── bash
   ├── R
   ├── Ref

```
**It is important you use this directory structure and construct your scripts to use the relative paths as invoked from `~/Assignment5`.  Scripts are marked based on the fact that they run. The expectation is that they will run within the relative path `~/Assignment5`. You will be expected to hard code this into your script.**


## Practical questions [*13 marks*]


This assignment will use _C elegans_ sequence data from [BioProject_598355]. These are 100nt paired-end reads. 

0. Install the `SRA Toolkit` and the latest version of `freebayes` using conda `conda install -c bioconda sra-tools freebayes=1.3.2` this will give you the `fastq-dump` tool that you will use to get the data and an updated `freebayes` that will generate a VCFv4.2 file. (no marks for this)
    
1. Write a __script__ with informative and useful comments [*1 mark*] __to download and trim and clean__ (This should be easy as you can re-purpose a previous script - __Note that I give you examples of commands with paths, but you should always ensure that your script uses the correct path - I do not guarantee that my commands as listed will all have the correct path or be typo free__): 
    + download the fastq read data using sra-tools. To download the dataset from the Short Read Archive (SRA) at NCBI use the following command: `fastq-dump --gzip --split-files SRR3452285 -O ~/Assignment5/0_rawData/fastq/` This will generate two fastq.gz files, one for each end [*1 mark*]. Note that in your script you may use variables to specify directories for input and output (see below)
    + Because of the constraints imposed by your VMs, we will only download the reference sequence for ChrI and limit our analysis to to this chromosome. Because of limitations with downloading from Box with `wget/curl` we have put these files into `/home/student/data/Assignment_5` on your VMs.  [*1 mark*]. 
    + run fastqc and trim the reads as for previous practicals/assignment. The command for trimming will be: `cutadapt -m 35 -q 30 -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT -o ./2_trimmedData/fastq/SRR3452285_1.fastq.gz -p ./2_trimmedData/fastq/SRR3452285_2.fastq.gz ./0_rawData/fastq//SRR3452285_1.fastq.gz ./0_rawData/fastq/SRR3452285_2.fastq.gz > ./2_trimmed_Data/log/cutadapt.log` [*1 mark*]
2. Write a __script__ with informative and useful comments [*1 mark*] __to align the cleaned reads to ChrI and call variants__:
    + Build the index for `bwa` for the chr1.fna sequence [*1 mark*]
    + Align the reads `bwa mem -t 2 ./Ref/WBcel235_chrI 1_trimmedData/fastq/SRR3452285_1.fastq.gz 1_trimmedData/fastq/SRR3452285_2.fastq.gz | samtools view -bhS -F4 - > 2_alignedData/bam/SRR3452285_chrI.bam` [*1 mark*]
    + Sort the `bam` file `samtools sort ./2_alignedData/bam/SRR3452285_chrI.bam -o ./2_alignedData/sorted_bam/SRR3452285_chrI.bam` [*1 mark*]
    + Index the `bam` file `samtools index ./2_alignedData/sorted_bam/SRR3452285_chrI.bam ` [*1 mark*]
    + Deduplicate with picard `picard MarkDuplicates I=./2_alignedData/sorted_bam/SRR3452285_chrI.bam O=./2_alignedData/sorted_bam/SRR3452285_chrI_rmdup.bam M=dups.metrics.txt REMOVE_DUPLICATES=true` [*1 mark*]
    + Get the samtools stats for the aligned, de-duplicated data `samtools stats ./2_alignedData/sorted_bam/SRR3452285_chrI_rmdup.bam > ./2_alignedData/log/SRR3452285_chrI_rmdup.stats` [*1 mark*]
    + Index the `bam` file `SRR3452285_chrI_rmdup.bam` [*1 mark*]
    + Use freebayes to call variants `freebayes -f ./Ref/WBcel235_chrI.fa ./2_alignedData/sorted_bam/SRR3452285_chrI_rmdup.bam > ./2_alignedData/vcf/SRR3452285_chrI.vcf` [*1 mark*]
3. Use RStudio to analyse your variants. __Create a new R project called `Assignment5.Rproj`__ in `~/Assignment5/R` and __create a new R markdown document called `Assignment5.Rmd`__.  You can then copy the code below 

~~~text
# Assessment Cover Sheet
<!--

Replace the comments - this is a comment - in the table with the relevant information.
Do not alter the first two lines of the table.

-->
| \  | \  |
|:--- |:--- |
| Student Name     | <!-- Your name --> |
| Student ID       | <!-- Your student ID -->  |
| Assessment Title | <!-- Assessment title -->  |	 
| Course/Program   | Bioinformatics and Systems Biology |
| Lecturer/Tutor   | <!-- Section instructor --> |
| Date Submitted   | <!-- Date --> |

## KEEP A COPY

Please be sure to make a copy of your work. If you have submitted assessment work electronically make sure you have a backup copy.

## PLAGIARISM AND COLLUSION

Plagiarism: using another person’s ideas, designs, words or works without appropriate acknowledgement.

Collusion: another person assisting in the production of an assessment submission without the express requirement, or consent or knowledge of the assessor.

## CONSEQUENCES OF PLAGIARISM AND COLLUSION

The penalties associated with plagiarism and collusion are designed to impose sanctions on offenders that reflect the seriousness of the University’s commitment to academic integrity. Penalties may include: the requirement to revise and resubmit assessment work, receiving a result of zero for the assessment work, failing the course, expulsion and/or receiving a financial penalty.

I declare that all material in this assessment is my own work except where there is clear acknowledgement and reference to the work of others. I have read the University Policy Statement on Academic Honesty & Assessment Obligations (http://www.adelaide.edu.au/policies/230).

I give permission for my assessment work to be reproduced and submitted to other academic staff for the purposes of assessment and to be copied, submitted to and retained by the University's plagiarism detection software provider for the purposes of electronic checking of plagiarism.
 
 
<!-- Replace a1234567 with your student ID and add the date. -->
| Signed   | Date       |
| -------- | ---------- |
| a1234567 | yyyy/mm/dd |

~~~


into your `Assignment5.Rmd` source window just below:
```
---
title: "Assignment5"
author: "Your Name"
date: "The date"
output: html_document
---

```
## So that it looks like this:
~~~text
---
title: "Assignment5"
author: "Dave Adelson"
date: "30/09/2020"
output: html_document
---

# Assessment Cover Sheet
<!--

Replace the comments - this is a comment - in the table with the relevant information.
Do not alter the first two lines of the table.

-->
| \  | \  |
|:--- |:--- |
| Student Name     | <!-- Your name --> |
| Student ID       | <!-- Your student ID -->  |
| Assessment Title | <!-- Assessment title -->  |	 
| Course/Program   | Bioinformatics and Systems Biology |
| Lecturer/Tutor   | <!-- Section instructor --> |
| Date Submitted   | <!-- Date --> |

## KEEP A COPY

Please be sure to make a copy of your work. If you have submitted assessment work electronically make sure you have a backup copy.

## PLAGIARISM AND COLLUSION

Plagiarism: using another person’s ideas, designs, words or works without appropriate acknowledgement.

Collusion: another person assisting in the production of an assessment submission without the express requirement, or consent or knowledge of the assessor.

## CONSEQUENCES OF PLAGIARISM AND COLLUSION

The penalties associated with plagiarism and collusion are designed to impose sanctions on offenders that reflect the seriousness of the University’s commitment to academic integrity. Penalties may include: the requirement to revise and resubmit assessment work, receiving a result of zero for the assessment work, failing the course, expulsion and/or receiving a financial penalty.

I declare that all material in this assessment is my own work except where there is clear acknowledgement and reference to the work of others. I have read the University Policy Statement on Academic Honesty & Assessment Obligations (http://www.adelaide.edu.au/policies/230).

I give permission for my assessment work to be reproduced and submitted to other academic staff for the purposes of assessment and to be copied, submitted to and retained by the University's plagiarism detection software provider for the purposes of electronic checking of plagiarism.
 
 
<!-- Replace a1234567 with your student ID and add the date. -->
| Signed   | Date       |
| -------- | ---------- |
| a1234567 | yyyy/mm/dd |

---
~~~

## Paste this text just below the setup code chunk.

~~~text

# Assignment 5 - [*23 marks*]

**Due before Tuesday 20th October**
## Correct submission [*3 Marks*]
Your answers to all questions should be submitted to myUni as a `.zip` file containing your two bash scripts, this html `knitted` R markdown file that contains the answers to the theoretical questions (can be in any readable format, but my order of preference is; html, pdf, plain text). [*1 mark for submitting as a single zip file*]

The `.zip` filename must start with your student number [*1 mark*] and your bash script must be able to run without errors.
Meaningful comments are strongly advised.

For all scripts, please use the directory `~/Assignment5` as the parent directory for all downloads and analysis. Use the following directory structure [*1 mark*]: 
```
.
├── Assignment5
   ├── 0_rawData
   │   ├── fastq
   │   └── FastQC
   ├── 1_trimmedData
   │   ├── fastq
   │   ├── FastQC
   │   └── log
   ├── 2_alignedData
   │   ├── bam
   │   ├── FastQC
   │   ├── log
   │   ├── sorted_bam
   │   └── vcf
   ├── bash
   ├── R
   ├── Ref

```
**It is important you use this directory structure and construct your scripts to use the relative paths as invoked from `~/Assignment5`.  Scripts are marked based on the fact that they run. The expectation is that they will run within the relative path `~/Assignment5`. You will be expected to hard code this into your script.**


## Practical questions [*13 marks*]

This assignment will use _C elegans_ sequence data from [BioProject_598355]. These are 100nt paired-end reads. 

0. Install the `SRA Toolkit` and the latest version of `freebayes` using conda `conda install -c bioconda sra-tools freebayes=1.3.2` this will give you the `fastq-dump` tool that you will use to get the data and an updated `freebayes` that will generate a VCFv4.2 file. (no marks for this)
    
1. Write a __script__ with informative and useful comments [*1 mark*] __to download and trim and clean__ (This should be easy as you can re-purpose a previous script - __Note that I give you examples of commands with paths, but you should always ensure that your script uses the correct path - I do not guarantee that my commands as listed will all have the correct path or be typo free__): 
    + download the fastq read data using sra-tools. To download the dataset from the Short Read Archive (SRA) at NCBI use the following command: `fastq-dump --gzip --split-files SRR3452285 -O ~/Assignment5/0_rawData/fastq/` This will generate two fastq.gz files, one for each end [*1 mark*]. Note that in your script you may use variables to specify directories for input and output (see below)
    + Because of the constraints imposed by your VMs, we will only download the reference sequence for ChrI and limit our analysis to to this chromosome. download the genome reference and gff3 files for _C elegans_ chromosome I to ~/Assignment5/Ref/ (see `https://universityofadelaide.box.com/shared/static/ct8alnh7a7t0z6x6m6g17t6n23raysbo.gz` and `https://universityofadelaide.box.com/shared/static/ssvrdjx5imz6ktrw5d0owuumx1ry158h.gz` ) [*1 mark*]. 
    + run fastqc and trim the reads as for previous practicals/assignment. The command for trimming will be: `cutadapt -m 35 -q 30 -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT -o ./2_trimmedData/fastq/SRR3452285_1.fastq.gz -p ./2_trimmedData/fastq/SRR3452285_2.fastq.gz ./0_rawData/fastq//SRR3452285_1.fastq.gz ./0_rawData/fastq/SRR3452285_2.fastq.gz > ./2_trimmed_Data/log/cutadapt.log` [*1 mark*]
2. Write a __script__ with informative and useful comments [*1 mark*] __to align the cleaned reads to ChrI and call variants__:
    + Build the index for `bwa` for the chr1.fna sequence [*1 mark*]
    + Align the reads `bwa mem -t 2 ./Ref/WBcel235_chrI 1_trimmedData/fastq/SRR3452285_1.fastq.gz 1_trimmedData/fastq/SRR3452285_2.fastq.gz | samtools view -bhS -F4 - > 2_alignedData/bam/SRR3452285_chrI.bam` [*1 mark*]
    + Sort the `bam` file `samtools sort ./2_alignedData/bam/SRR3452285_chrI.bam -o ./2_alignedData/sorted_bam/SRR3452285_chrI.bam` [*1 mark*]
    + Index the `bam` file `samtools index ./2_alignedData/sorted_bam/SRR3452285_chrI.bam ` [*1 mark*]
    + Deduplicate with picard `picard MarkDuplicates I=./2_alignedData/sorted_bam/SRR3452285_chrI.bam O=./2_alignedData/sorted_bam/SRR3452285_chrI_rmdup.bam M=dups.metrics.txt REMOVE_DUPLICATES=true` [*1 mark*]
    + Get the samtools stats for the aligned, de-duplicated data `samtools stats ./2_alignedData/sorted_bam/SRR3452285_chrI_rmdup.bam > ./2_alignedData/log/SRR3452285_chrI_rmdup.stats` [*1 mark*]
    + Index the `bam` file `SRR3452285_chrI_rmdup.bam` [*1 mark*]
    + Use freebayes to call variants `freebayes -f ./Ref/WBcel235_chrI.fa ./2_alignedData/sorted_bam/SRR3452285_chrI_rmdup.bam > ./2_alignedData/vcf/SRR3452285_chrI.vcf` [*1 mark*]
3. Use RStudio to analyse your variants.

[BioProject_598355]:https://www.ncbi.nlm.nih.gov/bioproject/598355
[here]:https://university-of-adelaide-bx-masters.github.io/BIOTECH-7005/COVERSHEET.html


## Getting RStudio set up
##############
# Packages   #
##############
```{r}
library(magrittr)
library(vcfR)
# library(ape)
library(tidyverse)
```
###############
# Functions   #
###############
```{r}
'%nin%' <- Negate('%in%')

ulength <- function(x) {x %>% unique() %>% length()}
fsummary <- function(x) {x %>% as.factor() %>% summary()}
```
############
# Set Up   #
############
```{r}
DirIn <- "/home/student/Assignment5/2_alignedData/vcf"
DirOut <- "/home/student/Assignment5/R/Out"
if(! dir.exists(DirOut)) {dir.create(DirOut)}

DirPlot <- file.path(DirOut, "Plots")
if(! dir.exists(DirPlot)) {dir.create(DirPlot)}
```
###############
# Load VCF    #
###############
## Load file - This loads the vcf file created from the de-duplicated, sorted `.bam` file. 
```{r}
FileIn <- "SRR3452285_chrI.vcf" 
VCF <- read.vcfR(file.path(DirIn, FileIn), verbose = TRUE)

VCF@meta %>% as_tibble() # metadata
VCF@fix %>% as_tibble() # co-ordinates etc and INFO ie data about this position in the genome
VCF@gt %>% as_tibble() # FORMAT ie genotypes

INFO2df(VCF[1:20,]) %>% as_tibble() # converts INFO key value pairs to df
```
#######################
# Analyse Genotypes   #
#######################
## Genotype info stored in separate columns for each sample - only one sample reported here
## Keys are stored in FORMAT col

## Convert FORMAT to DF, add co-ords etc
```{r}
VCF@gt[,1] %>% unique() # all keys are present in all records, usually not the case

ColNames <- VCF@gt[1,1] %>% strsplit(":") %>% unlist() # get col names
GenomeData <- VCF@gt %>% as_tibble %>% dplyr::select(-FORMAT) %>% separate(unknown, ColNames, sep = ":") # convert strings to cols
GenomeData$DP %<>% as.integer()
GenomeData$RO %<>% as.integer()
GenomeData$QR %<>% as.integer()
GenomeData$AO %<>% as.integer()
GenomeData$QA %<>% as.integer()
GenomeData %<>% select(GT, DP, AD, RO, AO, QR, QA, GL)

X <- VCF@fix %>% as_tibble() %>% dplyr::select(CHROM, POS, REF, ALT, QUAL) # get required INFO & other cols
colnames(X) <- c("Chr", "Pos", "Ref", "Alt", "Qual") # rename

GenomeData <- bind_cols(X, GenomeData) # merge INFO and FORMAT data
GenomeData$Pos %<>% as.integer()
GenomeData$Qual %<>% as.numeric()

rm(X)

## Extract P(Best) vs P(NextBest)
GenomeData %<>% rowwise() %>%
    mutate(GLB = GL %>% strsplit(",") %>% unlist()  %>% as.numeric() %>% sort() %>% .[2]) %>%
    ungroup() # split, sort, take the 2nd value (first always 0)

## View col definitions and data to see what each col means
VCF@meta %>% as_tibble() %>% filter(grepl("FORMAT", value))
GenomeData

```

## Plots

```{r}
# Positions with an alt allele, only 2 alleles, qual < 500 (so plot is not squashed)
Title <- sprintf("C elegans Subset Alt Read Fraction vs Likelihood of Alt Allele")
GenomeData %>% filter(GT %in% c("0/1", "1/1"), str_count(AD, ",") == 1, Qual < 500) %>% mutate(AF = AO/(RO+AO)) %>%
    ggplot(aes(x = AF, y = Qual, group = GT, colour = GT)) +
    geom_point(size = 0.3) +
    scale_colour_manual(values = c("pink","dark red")) +
    labs(title = Title, x = "Alt Allele Read Fraction", y = "Likelihood of Alt Allele") +
    theme(axis.text.x = element_text(angle = 0, hjust = 0.5),
          plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5),
          legend.position = "right", legend.title=element_blank())
FileOut <- sprintf("%s%s", Title, ".jpeg")
ggsave(filename = FileOut, path = DirPlot, device = "jpeg")
```

```{r}
Title <- sprintf("C elegans Reference Score vs Genotype Call Quality Score")
GenomeData %>% filter(GT %in% c("0/0", "0/1", "1/1")) %>%
    ggplot(aes(x = GLB, y = Qual, group = GT, colour = GT)) +
    geom_point(size = 0.3) +
    scale_colour_manual(values = c("blue","pink","dark red")) +
    labs(title = Title, x = "GLB", y = "Qual") +
    theme(axis.text.x = element_text(angle = 0, hjust = 0.5),
          plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5),
          legend.position = "right")
FileOut <- sprintf("%s%s", Title, ".jpeg")
ggsave(filename = FileOut, path = DirPlot, device = "jpeg")
```

```{r}
Title <- sprintf("C elegans Reference Score vs Genotype Call Quality Score (Trunc)")
GenomeData %>% filter(GT %in% c("0/0", "0/1", "1/1")) %>% filter(GLB > -100, Qual < 1000) %>%
    ggplot(aes(x = GLB, y = Qual, group = GT, colour = GT)) +
    geom_point(size = 0.3) +
    scale_colour_manual(values = c("blue","pink","dark red")) +
    labs(title = Title, x = "GLB", y = "Qual") +
    theme(axis.text.x = element_text(angle = 0, hjust = 0.5),
          plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5),
          legend.position = "right", legend.title=element_blank())
FileOut <- sprintf("%s%s", Title, ".jpeg")
ggsave(filename = FileOut, path = DirPlot, device = "jpeg")
```
## Additional questions [*7 marks*] 

4. For this line in the R code `GenomeData %<>% select(GT, DP, AD, RO, AO, QR, QA, GL)`, what are the VCF descriptions found in the `vcf` file header for these fields? [*3 Marks*]

5. For the plot titled `C elegans Subset Alt Read Fraction vs Likelihood of Alt Allele` is the distribution what you expect?  Why? [*3 Marks*]

6. How many properly paired reads were there in the sorted, deduplicated `.bam` file used to generate the `.vcf`? [*1 Mark*]
~~~
## This will give you a .Rmd file that will allow you to generate the relevant output by running the code and it includes the all of the questions you need to answer. 
## These additional questions worth 7 marks are reproduced below for clarity.

4. For this line in the R code `GenomeData %<>% select(GT, DP, AD, RO, AO, QR, QA, GL)`, what are the VCF descriptions found in the `vcf` file header for these fields? 

5. For the plot titled `C elegans Subset Alt Read Fraction vs Likelihood of Alt Allele` is the distribution what you expect?  Why? 

6. How many properly paired reads were there in the sorted, deduplicated `.bam` file used to generate the `.vcf`? 


[BioProject_598355]:https://www.ncbi.nlm.nih.gov/bioproject/598355
[here]:https://university-of-adelaide-bx-masters.github.io/BIOTECH-7005/COVERSHEET.html

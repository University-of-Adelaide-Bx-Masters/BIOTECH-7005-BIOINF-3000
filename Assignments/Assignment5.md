# Assignment 5 - [*25 marks*]

**Due before Tuesday 20th October**

Your answers to all questions should be submitted to myUni as a `.zip` file containing your bash scripts, your html rendered R markdown file for part 4, and a single file containing the answers to the theoretical questions (can be in any readable format, but my order of preference is; html, pdf, plain text). [*1 mark for submitting as a single zip file*]

The `.zip` filename must start with your student number [*1 mark*] and your bash script must be able to run without errors.
Meaningful comments are strongly advised [*1 mark*]

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


## Practical questions [*12 marks*]

This assignment will use _C elegans_ sequence data from [BioProject_598355]. These are 100nt paired-end reads. 

1. Install the SRA Toolkit using conda `conda install -c bioconda sra-tools` this will give you the `fastq-dump` tool that you will use to get the data. (no marks for this)
    
2. Write a script with informative and useful comments [*1 mark*] download and trim and clean : 
    + download the fastq read data using sra-tools. To download the dataset from the Short Read Archive (SRA) at NCBI use the following command: `fastq-dump --gzip --split-files SRR3452285 -O ~/Assignment5/0_rawData/fastq/` This will generate two fastq.gz files, one for each end [*1 mark*]
    + Because of the constraints imposed by your VMs, we will only download the reference sequence for ChrI and limit our analysis to to this chromosome. download the genome reference and gff3 files for _C elegans_ chromosome I to ~/Assignment5/Ref/ (see `https://universityofadelaide.box.com/shared/static/ct8alnh7a7t0z6x6m6g17t6n23raysbo.gz` and `https://universityofadelaide.box.com/shared/static/ssvrdjx5imz6ktrw5d0owuumx1ry158h.gz` ) [*1 mark*]. 
    + run fastqc and trim the reads as for previous practicals/assignment. The command for trimming will be (note I am giving you a hint for your script with this command): `cutadapt -m 35 -q 30 -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT -o ${TRIMFQ}/SRR3452285_1.fastq.gz -p ${TRIMFQ}/SRR3452285_2.fastq.gz ${RAWFQ}/SRR3452285_1.fastq.gz ${RAWFQ}/SRR3452285_2.fastq.gz > ${TRIMLOG}/cutadapt.log` [*1 mark*]
3. Write a script with informative and useful comments [*1 mark*] to align the cleaned reads to ChrI and call variants:
    + Build the index for `bwa` for the chr1.fna sequence [*1 mark*]
    + Align the reads `bwa mem -t 2 ./Ref/WBcel235_chrI 1_trimmedData/fastq/SRR3452285_1.fastq.gz 1_trimmedData/fastq/SRR3452285_2.fastq.gz | samtools view -bhS -F4 - > 2_alignedData/bam/SRR3452285_chrI.bam` [*1 mark*]
    + Get the samtools stats for the aligned data `samtools stats 2_alignedData/bam/SRR2003569_chI.bam > 2_alignedData/log/SRR2003569_chI.stats` [*1 mark*]
    + Sort the `bam` file `samtools sort ./2_alignedData/bam/SRR3452285_chrI.bam -o ./2_alignedData/sorted_bam/SRR3452285_chrI.bam` [*1 mark*]
    + Index the `bam` file `samtools index ./2_alignedData/sorted_bam/SRR3452285_chrI.bam ` [*1 mark*]
    + Deduplicate with picard `java -jar /path/to/picard/tools/picard.jar MarkDuplicates I=./2_alignedData/sorted_bam/SRR3452285_chrI.bam O=./2_alignedData/sorted_bam/SRR3452285_chrI_rmdup.bam M=./2alignedData/log/dups.metrics.txt REMOVE_DUPLICATES=true` [*1 mark*]
    + Index the `bam` file `SRR3452285_chrI_rmdup.bam` [*1 mark*]
    + Use freebayes to call variants `freebayes -f ./Ref/WBcel235_chrI.fa ./2_alignedData/sorted_bam/SRR3452285_chrI_rmdup.bam > ./2_alignedData/vcf/SRR3452285_chrI.vcf` [*1 mark*]
4. Use RStudio to analyse your variants. Create a new R project called `Assignment5` in `~/Assignment5/R` and create a new R markdown document called `Assignment5`.  You can then copy the code below into your `Assignment5.Rmd` source window just below the setup code chunk.

```

## Getting RStudio set up
##############
## Packages ##
##############
```{r}
library(magrittr)
library(vcfR)
# library(ape)
library(tidyverse)
```
###############
## Functions ##
###############
```{r}
'%nin%' <- Negate('%in%')

ulength <- function(x) {x %>% unique() %>% length()}
fsummary <- function(x) {x %>% as.factor() %>% summary()}
```
############
## Set Up ##
############
```{r}
DirIn <- "/home/a1071750/student/Assignment5/2_alignedData/vcf"
DirOut <- "/home/a1071750/student/Assignment5/R/Out"
if(! dir.exists(DirOut)) {dir.create(DirOut)}

DirPlot <- file.path(DirOut, "Plots")
if(! dir.exists(DirPlot)) {dir.create(DirPlot)}
```
###############
## Load VCF ##
###############
# Load file - This loads the vcf file created from the de-duplicated, sorted `.bam` file. 
```{r}
FileIn <- "SRR3452285_chrI.vcf" 
VCF <- read.vcfR(file.path(DirIn, FileIn), verbose = TRUE)

VCF@meta %>% as_tibble() # metadata
VCF@fix %>% as_tibble() # co-ordinates etc and INFO ie data about this position in the genome
VCF@gt %>% as_tibble() # FORMAT ie genotypes

INFO2df(VCF[1:20,]) %>% as_tibble() # converts INFO key value pairs to df
```
#######################
## Analyse Genotypes ##
#######################
# Genotype info stored in separate columns for each sample - only one sample reported here
# Keys are stored in FORMAT col

# Convert FORMAT to DF, add co-ords etc
```{r}
VCF@gt[,1] %>% unique() # all keys are present in all records, usually not the case

ColNames <- VCF@gt[1,1] %>% strsplit(":") %>% unlist() # get col names
GenomeData <- VCF@gt %>% as_tibble %>% dplyr::select(-FORMAT) %>% separate(unknown, ColNames, sep = ":") # convert strings to cols
GenomeData$DP %<>% as.integer()
GenomeData$RO %<>% as.integer()
GenomeData$QR %<>% as.integer()
GenomeData$AO %<>% as.integer()
GenomeData$QA %<>% as.integer()

X <- VCF@fix %>% as_tibble() %>% dplyr::select(CHROM, POS, REF, ALT, QUAL) # get required INFO & other cols
colnames(X) <- c("Chr", "Pos", "Ref", "Alt", "Qual") # rename

GenomeData <- bind_cols(X, GenomeData) # merge INFO and FORMAT data
GenomeData$Pos %<>% as.integer()
GenomeData$Qual %<>% as.numeric()

rm(X)

# Extract P(Best) vs P(NextBest)
GenomeData %<>% rowwise() %>%
    mutate(GLB = GL %>% strsplit(",") %>% unlist()  %>% as.numeric() %>% sort() %>% .[2]) %>%
    ungroup() # split, sort, take the 2nd value (first always 0)

# View col definitions and data to see what each col means
VCF@meta %>% as_tibble() %>% filter(grepl("FORMAT", value))
GenomeData

```

## Plots
```{r}
Title <- sprintf("C elegans Reference Score vs Genotype Call Quality Score")
GenomeData %>% filter(GT %in% c("0/0", "0/1", "1/1")) %>%
    ggplot(aes(x = GLB, y = Qual, group = GT, colour = GT)) +
    geom_point(size = 0.3) +
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
    labs(title = Title, x = "GLB", y = "Qual") +
    theme(axis.text.x = element_text(angle = 0, hjust = 0.5),
          plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5),
          legend.position = "right", legend.title=element_blank())
FileOut <- sprintf("%s%s", Title, ".jpeg")
ggsave(filename = FileOut, path = DirPlot, device = "jpeg")
```

```
## Theoretical questions [*12 marks*]


[BioProject_598355]:https://www.ncbi.nlm.nih.gov/bioproject/598355

### Genes For Question 3
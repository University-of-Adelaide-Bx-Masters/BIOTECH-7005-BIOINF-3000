# Week 8 Practical - VCF Analysis - Dr Rick Tearle


# VCF Analysis practical

As with all previous practicals, start by creating a `~/Practical_8/` subdirectory on your VM. You will need to move the data and R script and Rproj files (see below) from your `~/data/variant_calling/` directory to the `~/Practical_8/` directory prior to working on it. 

## About this practical

This practical is designed to give you experience with VCF files and their analysis and will use the genome and NGS data from Assignment 4. It includes an Rscript and BAM and VCF files based on that data.  

The files you will need are in your `~/student/data/variant_calling` directory. These are the files you should have:

```
-rw-r--r-- 1 student student 9.1M Sep 12 06:07 Arabidopsis_thaliana.TAIR10.48.gff3.gz
-rw-r--r-- 1 student student  35M Sep 12 06:07 Arabidopsis_thaliana.TAIR10.dna.toplevel.fa.gz
-rw-rw-r-- 1 student student 1.1G Sep 12 06:27 SRR5882792_Athaliana_TAIR10.bam
-rw-rw-r-- 1 student student 351K Sep 14 01:17 SRR5882792_Athaliana_TAIR10.bam.bai
-rw-rw-r-- 1 student student  15M Sep 12 06:28 SRR5882792_Athaliana_TAIR10_picard.vcf
-rwxr-xr-x 1 student student 9.3K Sep 14 05:47 VCF_Analysis.R
-rwxr-xr-x 1 student student  234 Sep 14 05:47 VCF_Analysis.Rproj

```

Move these files to your `~/Practical_8/` directory. You will then be able to open the Rproj and Rscript files to start working. 

The information for reads containing variants is stored in the BAM files and variant information is stored in the VCF files. We will visualise both to see the evidence that supports the variant calls. You will use an R script to summarise the data in the VCF file and to identify interesting variants to interrogate further. 

For the visualisation of variant data you will use IGV [Integrated Genome Viewer][1] and you will need to install this on your own computer as it will not run on a browser from your VM. There are versions of IGV available for download for any operating system you might be using: Windowns, MacOS or Linux. Visualisation of your results will require you to download the R output of the analysis to your computer. 

Specific instructions will be provided during the practical to allow you to carry out the steps required. In previous practicals you have been able to access background material, along with a detailed description of practical steps and the reasons behind them. For this practical, all of that information will be provided verbally during the scheduled practical session, and if time permits it will be incorporated into a revised version of this document. In contrast to previous practical sessions, we will record this practical and upload it to MyUni to provide an accessible record of the instructions and procedures. 

## Set up

First load the packages we need.

```{r, results='hide', eval = FALSE}
library(magrittr)
library(vcfR)
library(tidyverse)
```
Then load the functions we need.

```{r, results='hide', eval = FALSE}
'%nin%' <- Negate('%in%')

ulength <- function(x) {x %>% unique() %>% length()}
fsummary <- function(x) {x %>% as.factor() %>% summary()}
```

Now set up the paths and environment for input and output. 

```{r results='hide', eval=FALSE} 
DirIn <- "/home/student/Practical_8" 
DirOut <- "/home/student/Practical_8/Out" 
if(! dir.exists(DirOut)) {dir.create(DirOut)} 

DirPlot <- file.path(DirOut, "Plots") 
if(! dir.exists(DirPlot)) {dir.create(DirPlot)} 
```


## Some notes
The R package vcfR contains code to load and manipulate VCF files. We will only use a small proportion of its functionality but if you find yourself interrogating VCF files regularly, it will be worth your while reading up on the functions in this package.

I have provided code to load reference sequences and annotation (gff) files in R as well, see sections **Load GFF File** and **Load Reference File**. We will not use them here: they are provided for the sake of completeness.

After loading packages, we will start with the section **IGV**. This just contains a lists of co-ordinates for you to look at in IGV. We will inspect through these locations and I will have some comments on each.

Because IGV cannot deal with `*.gz` files, you will have to decompress your genomic data file

```bash
pigz -d ~/student/Practical_8/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa.gz
```
You will also need to create an index for that file

```bash
samtools faidx ~/Project_8/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa
```

Download the following files to your local computer using RStudio's File browser.
Simply select 1 file at a time by checking the checkbox and click "More" >> "Export...".
Click the "Download" button and save it somewhere obvious.

* `~/student/Practical_8/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa`
* `~/student/Practical_8/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa.fai`
* `~/student/Practical_8/SRR5882792_Athaliana_TAIR10.bam`
* `~/student/Practical_8/SRR5882792_Athaliana_TAIR10.bam.bai`

Visit, [IGV-web](https://igv.org/app/) and load the genome from a `Local File ...` by selecting both the `Arabidopsis_thaliana.TAIR10.dna.toplevel.fa` and `Arabidopsis_thaliana.TAIR10.dna.toplevel.fa.fai` files.
Once the reference genome is loaded, load a "Track" from a `Local File ...` by selecting both the `SRR5882792_Athaliana_TAIR10.bam` and `SRR5882792_Athaliana_TAIR10.bam.bai` files.

![IGV-web SRR11140748 Illumina](images/SRR11140748-illumina.png)

Next, go to **Load VCF**. After loading the VCF file we will look at how to access some of the data.

The next section is **Analyse Genotypes**. We will extract a subset of the columns from the VCF file as a data frame and reformat some of them. Data frames are a convenient way to store and interrogate this information.

We will dig in to the data to see what evidence is being used to call variants, in **Plots**. We will look at the relationship between the genotypes called, the number of reads and the strength of evidence. When we find an interesting relationship we will take examples and look at them in IGV, where we can get a sense of the reads that support a call.

Finally, I would like you to take the data and the code, and generate your own plots etc, to see if you can find a threshold to balance True and False Positives. The columns Qual and GLB will be useful here.

[1]: https://software.broadinstitute.org/software/igv/

# Week 8 Practical - VCF Analysis - Dr Rick Tearle
{:.no_toc}

* TOC
{:toc}

# VCF Analysis practical

As with all previous practicals, start by creating a `~/Practical_8/` subdirectory on your VM. You will need to move the data and R script and Rproj files (see below) from your `~/data/variant_calling/` directory to the `~/Practical_8/` directory prior to working on it. 

## About this practical

This practical is designed to give you experience with VCF files and their analysis and will use the genome and NGS data from Assignment 4. It includes an Rscript and BAM and VCF files based on that data.  

The files you will need are in your `~/sudent/data/variant_calling` directory. These are the files you should have:

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

For the visualisation of variant data you will use IGV [Integreated Genome Viewer][1] and you will need to install this on your own computer as it will not run on a browser from your VM. There are versions of IGV available for download for any operating system you might be using: Windowns, MacOS or Linux. Visualisation of your results will require you to download the R output of the analysis to your computer. 

Specific instructions will be provided during the practical to allow you to carry out the steps required. In previous practicals you have been able to access background material, along with a detailed description of practical steps and the reasons behind them. For this practical, all of that information will be provided verbally during the scheduled practical session, and if time permits it will be incorporated into a revised version of this document. In contrast to previous practical sessions, we will record this practical and upload it to MyUni to provide an accessible record of the instructions and procedures. 

[1]: https://software.broadinstitute.org/software/igv/
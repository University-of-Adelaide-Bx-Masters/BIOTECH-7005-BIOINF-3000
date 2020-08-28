#!/bin/bash

INIT=$(pwd)

# update conda & make sure r isn't there
chown -R student:student /home/student/miniconda3/*
echo -e "Updating conda\n"
CONDA=/home/student/miniconda3/bin/conda
${CONDA} update -y -n base -c defaults conda

# add bioconda channels
${CONDA} config --add channels conda-forge 
${CONDA} config --add channels bioconda 
${CONDA} config --add channels maxibor  # Required for adapterremoval
${CONDA} config --add channels biocore  # Required for mafft

echo -e "Installing conda pkgs\n"
${CONDA} install --yes samtools bowtie2 picard jalview freebayes cutadapt sabre star fastqc bwa sambamba hisat2 bedtools subread bcftools


# Now the apt ones
#echo -e "Installing apt pkgs\n"
#apt update
#apt -y upgrade
#apt install -y bwa picard sambamba bedtools hisat2 subread bcftools default-jre default-jdk

# This was from 2019 and the zip archive was not in /opt
# Fix up the sabre installation
#echo -e "Correcting sabre installation\n"
#cd /opt 
#unzip master.zip
#mv sabre-master sabre
#cd sabre
#make
#cd /usr/local/bin
#rm sabre
#ln -s /opt/sabre/sabre ./sabre
#rm /opt/master.zip

# Install STAR
#echo -e "Installing STAR\n"
#cd /opt
#git clone https://github.com/alexdobin/STAR.git
#cd STAR/source
#make STAR
#cd /usr/local/bin
#ln -s /opt/STAR/bin/Linux_x86_64/STAR ./STAR

# Install FastQC
#cd /opt
#wget https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.8.zip
#unzip fastqc_v0.11.8.zip
#rm fastqc_v0.11.8.zip
#chmod 755 FastQC
#cd /usr/local/bin
#ln -s /opt/FastQC/fastqc ./fastqc

# Cleanup
echo -e "Cleaning up installations\n"
apt autoremove -y
${CONDA} clean --all --yes

echo -e "Returning to initial directory\n"
cd ${INIT}

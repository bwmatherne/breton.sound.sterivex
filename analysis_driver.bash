#!/bin/bash

############################################################################

# Analysis Driver for Sterivex filter data using Riffomonas.org instructions
# Settings are inabled for running on the LSU HPC

############################################################################

# 1. Download mothur with miniconda

conda install -c bioconda mothur
conda install -c bioconda/label/cf201901 mothur

#2. Create make file
make.file(inputdir=data/raw, type=gz, prefix=stability, numcols=4)
make.contigs(file=stability,  oligos=barcode.oligos, checkorient=t, inputdir=data/raw, outputdir=data/mothur/, processors=20)


Undetermined_S0_L001_I1_001.fastq.gz
metadata.tsv                       Undetermined_S0_L001_R1_001.fastq.gz
mothur.1660174282.logfile          Undetermined_S0_L001_R2_001.fastq.gz

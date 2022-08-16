#!/bin/bash

############################################################################

# QIIME2 Analysis Driver for breton.sound.sterivex filter data using Riffomonas.org instructions
# Settings are inabled for running on the LSU HPC

############################################################################

#Demultiplex and denoise each Illumina run
# Note - R settings for dada2 are found in the seq1 code script

qsub code/seq1.pbs
qsub code/seq2.pbs
qsub code/seq3.pbs
qsub code/seq4.pbs

# Combine output files from each run and make a single table, rept seqs, and tree files
qsub code/combine_table.pbs

# Taxonomic assignment based on classifier computed in Spring 2022
qsub code/trainset.pbs

# Run script to remove mitochondria, chloroplasts, and sequences with no phyla assignement
# Script will also produce new phylogenetic trees
qsub code/filter-taxa.pbs

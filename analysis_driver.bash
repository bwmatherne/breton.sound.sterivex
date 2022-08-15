#!/bin/bash

############################################################################

# Analysis Driver for breton.sound.sterivex filter data using Riffomonas.org instructions
# Settings are inabled for running on the LSU HPC

############################################################################

#Demultiplex and denoise each Illumina run
# Note - R settings for dad2 are foudn in the seq1 code script

qsub code/seq1.pbs
qsub code/seq2.pbs
qsub code/seq3.pbs
qsub code/seq4.pbs

# Combine Metadata files and filter out diversion information


#Load packages
library(dplyr)
library(ggplot2)
library(tidyr)
library(lubridate)
library(stringr)
library(ggfortify)
library(prediction)
library(tidyverse)

# clear R's brain
rm(list = ls())

seq3 <- read_tsv('data/references/metadata.seq123.tsv')
seq4 <- read_tsv('data/references/metadata.seq4.tsv')

join <- full_join(seq3, seq4)

# Remove NA values from table for qiime
df <- as.data.frame(join)
df <- sapply(df, as.character)
df[is.na(df)] <- " "
df <- as.data.frame(df)

#Select only Diversion Event Samples
diversion <- df %>% filter(Event == 'Diversion' | Event == 'categorical')


write_tsv(diversion,'data/references/diversion.metadata.tsv')

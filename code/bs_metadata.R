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

full <- read_tsv('data/raw/combined.metadata.tsv')

# Remove NA values from table for qiime
df <- as.data.frame(full)
df <- sapply(df, as.character)
df[is.na(df)] <- " "
df <- as.data.frame(df)

#Select only Diversion Event Samples
bs <- df %>% filter(Site == 'BS' | Site == 'categorical')


write_tsv(bs,'data/process/bs.metadata.tsv')

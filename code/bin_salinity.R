# R Script to Bin Salinity Data

#Load dependencies

if (!requireNamespace("devtools", quietly = TRUE)){install.packages("devtools")}
devtools::install_github("jbisanz/qiime2R")


library(qiime2R)
library(tidyverse)

#Importing metadata
bin <- read_q2metadata("data/process/bs.metadata.tsv")

hist(bin$Salinity)

bin %>% group_by(Location) %>%
  mutate(avg = mean(Salinity))

# subset by location
inshore <- bin %>% filter(Location == "Inshore")
hist(inshore$Salinity)
summary(inshore) # Highest salinity is 6.21 ppt



coast<- bin %>% filter(Location == "Coastal") 
hist(coast$Salinity)
filter(coast, Salinity <= 6.21) %>% count() #22 samples are <= 6.21

#Create bin column for three groups based on salinity findings

coast.low <-  filter(coast, Salinity <= 6.21) %>% add_column(bin_sal = "coast.low")
coast.hi <- filter(coast, Salinity > 6.21) %>% add_column(bin_sal = "coast.hi")
inshore <- add_column(inshore, bin_sal= "inshore") 

stack <- bind_rows(coast.hi,coast.low,inshore)

write_tsv(stack, file = "data/process/bin.metadata.tsv")

# Combine with bs.metadata file for use in qiime2

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

BS <- read_tsv('data/process/bs.metadata.tsv') 
SAL <- read_tsv('data/process/bin.metadata.tsv')
SAL <- select(SAL, 
              SampleID,
              bin_sal)

join <- left_join(BS, SAL, by = c("#SampleID" = "SampleID"))

join$bin_sal <- replace_na(join$bin_sal,"categorical")

# Remove NA values from table for qiime
df <- as.data.frame(join)
df <- sapply(df, as.character)
df[is.na(df)] <- " "
df <- as.data.frame(df)



write_tsv(df,'data/process/sal.metadata.tsv')






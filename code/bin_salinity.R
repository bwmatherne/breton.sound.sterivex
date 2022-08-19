# R Script to Bin Salinity Data

#Load dependencies

if (!requireNamespace("devtools", quietly = TRUE)){install.packages("devtools")}
devtools::install_github("jbisanz/qiime2R")


library(qiime2R)
library(tidyverse)

#Importing metadata
bin <- read_q2metadata("data/process/bs.metadata.tsv")


#perform binning with specific number of bins
bin1 <- bin %>% mutate(bin_sal = ntile(Salinity, n=3))
head(bin1) # Note, evenly distributes samples but not specific by salinty value

#perform binning with specific bin ranges

bin %>% mutate(bin_sal = cut(Salinity, breaks=c(0,5,15))) %>% group_by(bin_sal) %>% count()
bin %>% mutate(bin_sal = cut(Salinity, breaks=3)) %>% head()

glimpse(bin)
write_tsv(bin, file = "data/process/bin.metadata.tsv")


sal <- pull(bin,Salinity)
hist(bin$Salinity)

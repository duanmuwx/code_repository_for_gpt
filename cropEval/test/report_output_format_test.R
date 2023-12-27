library(readxl)
library(writexl)
library(biosim)
library(tidyverse)
library(reshape2)
library(stringi)
library(ggplot2)
library(ggpmisc)
data('popmean')
data('adv')
data('hybobs')
data('blupebv')
data('Aov')
data('pheno')
region = 'SpML'
traits <- c('YLD','MST','TWT','PHT','EHT','PMD','LDG','NCLB','GLS')
ctrlnsample = 'ZMN00545'
ctrlssample = 'ZMN00080'
outebv <- del_result(popmean, region,
                     traits, ctrlnsample,
                     ctrlssample, adv,
                     hybobs, blupebv,
                     Aov)

print(names(outebv))
print(head(outebv[['SS_GCA']]))

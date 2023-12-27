
#This script produces a demo dataset to be used in cropEval package. The values of
#traits have been modified and cannot be used for real analyses.

rm(list = ls())
pdata <- readRDS('F:/projects/geneticAnalysis/datafiles/phenotypes_P_all.rds')
names(pdata)
ml <- pdata$SpML
su <- pdata$SU
su <- su[,c("Location","Field","Range" ,"Pass","Set","Rep","Type","IsCK","X22plot" ,
              "label","Name","GenoID","FGenoID","MGenoID","weight","MST","YLD","PHT",
              "EHT","PMD.","LDG.")]
names(su)[names(su)=='PMD.'] <- 'PMD'
names(su)[names(su)=='LDG.'] <- 'LDG'
names(su)[names(su)=='X22plot.'] <- 'X22Plot'
names(su)[names(su)=='label'] <- 'Label'

head(su)
ml <- ml[,c("Location","Field","Range" ,"Pass","Set","Rep","Type","IsCK","X22Plot" ,
              "Label","Name","GenoID","FGenoID","MGenoID","weight","MST","YLD","PHT",
              "EHT","PMD","LDG")]
names(su)==names(ml)

traits <- c("weight","MST","PHT","EHT","PMD","LDG")
for (i in 1:length(traits)){
    rand <- runif(nrow(su))-0.5
    xbar <- mean(su[,traits[i]],na.rm = T)*0.25*rand
    su[,traits[i]] <- ifelse(!is.na(su[,traits[i]]),su[,traits[i]] + xbar,NA)
}
su$YLD <- ifelse(!is.na(su$YLD),su$weight*(100-su$MST)*1.292636,NA)
for (i in 1:length(traits)){
    rand <- (runif(nrow(ml))-0.5)
    xbar <- mean(ml[,traits[i]],na.rm = T)*0.25*rand
    ml[,traits[i]] <-ifelse(!is.na(ml[,traits[i]]),ml[,traits[i]] + xbar,NA)
}
ml$YLD <- ifelse(!is.na(ml$YLD),ml$weight*(100-ml$MST)*1.292636,NA)

#plot(YLD2~YLD,ml)
hist(ml$YLD)
phdemo <- list(SpML=pdata$SpML,SU=pdata$SU)
rm( "i", "ml","pdata" ,"rand" ,"su" , "traits","xbar")

save.image('data/phenotypes_P_alltrails.rda')
#
# library(readxl)
# library(writexl)
# library(biosim)
# library(tidyverse)
# library(reshape2)
# library(stringi)
# library(asreml)


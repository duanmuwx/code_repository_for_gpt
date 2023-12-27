


rm(list=ls())
missingAdding <- function(w){
    m <- table(w$Range,w$Pass)
    row.names(m)
    names(m)
    m <- as.matrix(m)
    colnames(m)
    
    nr=nrow(m)
    nc=ncol(m)
    ndiff <- nc*nr-sum(m)
    
    if (ndiff>0){
        missing <- c()
        for (i in 1:nr){
            for (j in 1:nc){
                if (m[i,j]==0){
                    rn <- rownames(m)[i]
                    cn <- colnames(m)[j]
                    missing <- rbind(missing,                                                      c(as.numeric(as.character(rn)),
                                                                                                     as.numeric(as.character(cn))))
                }
            }
        }
        
        zeros <- data.frame(matrix(NA,nrow(missing),ncol(w)))
        names(zeros) <- names(w)
        for (i in 1:nrow(missing)){
            zeros$Pass[i] <- missing[i,2]
            zeros$Range[i] <- missing[i,1]
            zeros$Name[i] <- paste0('DUMMY',i)
            zeros$Field[i] <- w$Field[1]
            zeros$Location[i] <- w$Location[1]
        }
        
        ww <- rbind(w,zeros)
        table(ww$Range,ww$Pass)
    } else {
        ww <- w
    }
    
    
    return (ww)
    
}

replicationCheck <- function(df){
    q <- data.frame(table(df$Name))
    return (nrow(q[q$Freq>1,])/nrow(q))
}

blueEstimation <- function(datdf,repct,trait){
    if (all(c('Name','Range','Pass','Location','Field','Set',
              'Type','trait') %in% names(datdf))){
        datdf$Range <- as.factor(datdf$Range)
        datdf$Pass <- as.factor(datdf$Pass)
        datdf$Name <- as.factor(datdf$Name)
        datdf$Set <- as.factor(datdf$Set)
        datdf$Location <- as.factor(datdf$Location)
        
        if (nrow(datdf[!is.na(datdf$trait),])<130){
          cat('asr <- asreml(trait ~ 1,random = ~idv(Name),
                              residual = ~idv(units), data = datdf)\n')
          asr <- asreml(trait ~ Range + Pass,random = ~idv(Name),
                        residual = ~idv(units), data = datdf)
          
        } else if (nrow(datdf[!is.na(datdf$trait),])<150){
          cat('asr <- asreml(trait ~ 1,random = ~idv(Name),
                              residual = ~ar1(Range):ar1(Pass), data = datdf)\n')
          asr <- asreml(trait ~ 1,random = ~idv(Name),
                        residual = ~ar1(Range):ar1(Pass), data = datdf)
        
            
        } else if(trait=='TWT'){# & datdf$Field[1]==2 & datdf$Type[1]=='PLate'){
            cat('asr <- asreml(trait ~ Pass,random = ~idv(Name),
                          residual = ~ar1(Range):ar1(Pass), data = datdf)\n')
            asr <- asreml(trait ~ Pass,random = ~idv(Name),# + idv(units),
                          residual = ~ar1(Range):ar1(Pass), data = datdf)
            
        } else {
            cat('asr <- asreml(trait ~ Range + Pass,random = ~idv(Name),
                          residual = ~ar1(Range):ar1(Pass), data = datdf)\n')
            asr <- asreml(trait ~ Range + Pass,random = ~idv(Name),# + idv(units),
                          residual = ~ar1(Range):ar1(Pass), data = datdf)
            
        }
        asr <- update(asr)
        if(!asr$converge)asr <- update(asr)
        if(!asr$converge)asr <- update(asr)
        if(!asr$converge)asr <- update(asr)
        print(wald(asr))
        summary(asr)$varcomp
        H2 <- vpredict(asr,H2~V1/(V1+V2))
        H2$Location <- unique(datdf$Location)[1]
        H2$Field <- unique(datdf$Field)[1]
        H2$Type <- unique(datdf$Type)[1]
        H2$N <- nrow(datdf)
        H2
        
        # if (repct>0.1){ #BLUE will be estimated 
        #                 #if there are replications
        #     asr <- asreml(YLD ~ Name,
        #                   residual = ~ar1(Range):ar1(Pass), 
        #                   data = datdf)
        #     if(!asr$converge)asr <- update(asr)
        #     if(!asr$converge)asr <- update(asr)
        #     if(!asr$converge)asr <- update(asr)
        #     
        # } 
        
        blue <- predict(asr, classify = "Name",pworkspace="800mb",
                        data=datdf)$pvals
       
        if(nrow(blue[!is.na(blue$predicted.value),]!=0)){
          blue <- blue[blue$Name!='0',]
          blue$Location <- unique(datdf$Location)[1]
          blue$Field <- unique(datdf$Field)[1]
          blue$Type <- unique(datdf$Type)[1]
          blue$N <- nrow(datdf)
          blue$converge <- asr$converge
        } else {
          blue <- data.frame(asr$coefficients$random)
          names(blue)[1] <- 'predicted.value'
          newcol <- stringsplit(row.names(blue),'_')
          blue$Name <- newcol$X2
          blue <- blue[,c(2,1)]
          blue$std.error <- 'NA'
          blue$status <- 'NA'
          blue$Location <- unique(datdf$Location)[1]
          blue$Field <- unique(datdf$Field)[1]
          blue$N <- nrow(datdf)
          blue$converge <- asr$converge
          
        }
        head(blue)
        
       
        
         blue <- blue[blue$Name!='0',]
        blue$Location <- unique(datdf$Location)[1]
        blue$Field <- unique(datdf$Field)[1]
        blue$Type <- unique(datdf$Type)[1]
        blue$N <- nrow(datdf)
        blue$converge <- asr$converge
        head(blue)
        
        title <- paste(trait,'BLUE:  Location =',
                       unique(datdf$Location)[1],
                       ' Field =',unique(datdf$Field)[1],
                       ' H2 =',round(H2[1,1],2),'N =',nrow(datdf),
                       '\n')
        qqnorm(blue$predicted.value,main=title)
        qqline(blue$predicted.value)
        hist(blue$predicted.value,main = title)
    } else {
        cat('Not all necessary columns are included in data frame\n')
        head(datdf)
    }
    outlist <- list(H2=H2,blue=blue,asr=asr)
    
    return (outlist)
    
}


library(readxl)
library(writexl)
library(biosim)
library(tidyverse)
library(reshape2)
library(stringi) 
library(asreml)

# source('scripts/functions.R')
# source('scripts/blueEstimation.R')

#filename <- '22YTSPML-C.xlsx'
#sheets <- c('P-SpML')
outfile <- 'Output/blue_ML.csv'

# 22P-ML.xlsx -------------------------------------------------------

#vartype <- c('numeric','text',rep('numeric',5),rep('text',8),
            # rep('numeric',12),'text',rep('numeric',2),'text')
#d <- data.frame(read_excel(filename,sheet = sheets[1],col_type=vartype))
#names(d)[names(d)=='P.'] <- 'Phash'
#names(d)[names(d)=='E.'] <- 'Ehash'
#names(d)[names(d)=='PMD.'] <- 'PMDhash'
#names(d)[names(d)=='LDG.'] <- 'LDGhash'
#library("readxl")
#d<-read_excel('22YTSPML-C.xlsx',sheet = 1)

d<-read.csv("E:/P1项目/pheno/blue_estimation/22YTSPML-C.csv",header=T,sep=",",check.names = F)
d <- d[order(d$Location,d$Type,d$Field,d$Range,d$Pass),]
table(d$Location) 
table(d$Field,d$Type)
ped <- unique(d[,c('Name','FGenoID','MGenoID')])

dim(ped)  #1902
#d$Type[d$Type %in% c('PME','Pshort')] <- 'PME-Pshort'
unique(d$Type)

traits <- c('MST','TWT','YLD','PHT','EHT','GLS'#'Phash','Ehash','PMDhash',
           )# 'PMD',"LDG" "NCLB""GLS"
trait <- 'GLS'

dat <- d

dat <- dat[order(dat$Location,dat$Type,dat$Field,
                 dat$Range,dat$Pass),]

locations <- unique(dat$Location)
table(dat$Field,dat$Type)

repctage <- c()
blue <- c()
h2 <- list()
blue <- list()
for (trait in traits){
    blue[[trait]] <- c()
    h2[[trait]] <- c()
}

trait <- 'YLD'
loc <- 1
i=1
j=1

cat('Outliers checking\n')
for (trait in traits){
  w <- dat
  names(w)[names(w)==trait] <- 'trait'
  
  std=sd(w$trait,na.rm = T)
  mn=mean(w$trait,na.rm = T)
  w$abnormal <- ifelse(w$trait<mn-std*3 | w$trait>mn+std*3,
                       2,1)
  if (trait=="MST") w$abnomal[w$trait>=43] <- 2
  if (trait=="YLD") w$abnomal[w$trait<200] <- 2
  if (trait=="YLD") w$abnomal[w$MST>=43] <- 2#特殊异常值
  if (trait=="PMD") w$abnomal<-1
  if (trait=="GLS") w$abnomal<-1
  if (trait=="NCLB") w$abnomal<-1
  title <- paste(trait,'Phenotype ','N =',nrow(w[!is.na(w$trait),]),'\n')
  
  qqnorm(w$trait,main=title,col=w$abnormal)
  qqline(w$trait)
  #hist(w$trait,main=title)
  
  
  #if (trait=='MST')w$trait[w$trait>=43] <- NA
  #if (trait=='YLD')w$trait[w$trait<200] <- NA
  #if (trait=='YLD')w$trait[w$MST>=43] <- NA
  #if (trait=='Phash')w$trait[w$trait>=43] <- NA
  
 # w$trait[w$trait<mn-std*3 | w$trait>mn+std*3] <- NA
  w$trait[w$abnomal==2]<-NA
  title <- paste(trait,'Phenotype ','N =',nrow(w[!is.na(w$trait),]),
                 ' Outliers removed\n')
  
  qqnorm(w$trait,main=title)
  qqline(w$trait)
  #hist(w$trait,main=title)
}


cat('\n\nBLUE calculated for each field and outliers checked separately\n')
for (trait in traits){
    for (loc in 1:length(locations)){

        df <- dat[dat$Location==locations[loc],]
        fields <- unique(df$Field)
        fields <- fields[order(fields)]

        for (i in 1:length(fields)){
            dff <- df[df$Field==fields[i],]
            types <- unique(dff$Type)

            for (j in 1:length(types)){
                cat('loc=',loc,'Field=',i,'Type=',j,'\n')
                w <- dff[dff$Type==types[j],]
                
                repct <- replicationCheck(w)
                repctage <- rbind(repctage,c(locations[loc],
                                             fields[i],types[j],repct))
                cat('\nLocation:',locations[loc],'Field:',
                    fields[i],'Type',types[j],'N:',nrow(w),
                    'Replications:',repct,'\n')
                if (locations[loc]=='6JLDF' & fields[i]==2 &
                    types[j]=='P1') {
                    w$Range[w$Set==25] <- 50
                    w$Range[w$Set==26] <- 51
                }

                #print(table(w$Range,w$Pass))
                ww <- missingAdding(w)
                ww <- ww[order(ww$Range,ww$Pass),]
                #print(table(ww$Range,ww$Pass))

                for (trait in traits){
                if (trait=='YLD'){
                    datw <- ww[,c('Name','Range','Pass','Location',
                              'Field','Type','Set',trait,'MST')]
                } else {
                    datw <- ww[,c('Name','Range','Pass','Location',
                                  'Field','Type','Set',trait)]
                }
                names(datw)[8] <- 'trait'
                
                std=sd(datw$trait,na.rm = T)
                mn=mean(datw$trait,na.rm = T)
                datw$abnormal <- ifelse(datw$trait<mn-std*3 | datw$trait>mn+std*3,
                                        2,1)
                
                datw$trait[datw$trait<mn-std*3] <- NA
                datw$trait[datw$trait>mn+std*3] <- NA
                if (trait=='MST')datw$trait[datw$trait>=43] <- NA
                if (trait=='YLD')datw$trait[datw$trait<200] <- NA
                if (trait=='YLD')datw$trait[datw$MST>=43] <- NA
                if (trait=='P#')datw$trait[datw$trait>=43] <- NA
                

                if (nrow(datw[!is.na(datw$trait),])>20){
                    title <- paste(trait,'Phenotype:  Location =',
                                   unique(datw$Location)[1],'Type =',types[j],
                                   ' Field =',datw$Field[1],
                                   'N =',nrow(datw),
                                   '\n')

                    qqnorm(datw$trait,main=title,col=datw$abnormal)
                    qqline(datw$trait)
                    #hist(datw$trait)

                    cat(types[j],trait,nrow(datw[!is.na(datw$trait),]),
                        '\n')
                    head(datw)
                    
                    blueout <- blueEstimation(datdf=datw,
                                              trait=trait,repct = repct)
                    blue[[trait]] <- rbind(blue[[trait]],blueout$blue)
                    h2[[trait]] <- rbind(h2[[trait]],blueout$H2)
                }
                }
            }
        }
    }
}

repctage <- as.data.frame(repctage)
names(repctage) <- c('Location','Field','Type','Repct')
repctage$Repct <- as.numeric(repctage$Repct)
bluefile <- 'Output/blue_ML.xlsx'
h2file <- 'Output/h2_ML.xlsx'
rdsfile <- 'Output/blue_ML.rds'
writexl::write_xlsx(blue,path=bluefile)
writexl::write_xlsx(h2,path=h2file)
saveRDS(blue,rdsfile)

ped <- unique(dat[,c('Name','FGenoID','MGenoID')])
dim(ped)

pedfile <- 'Output/blue_ML_ped.rds'
saveRDS(ped,pedfile)

write.csv(ped,'Output/pedigree_ML.csv',quote = F,row.names = F)





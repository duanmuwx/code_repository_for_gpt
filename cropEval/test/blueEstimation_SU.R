rm(list = ls())
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
        datdf <- datdf[!is.na(datdf$Name),]
        
        if (nrow(datdf[!is.na(datdf$trait),])<0){
            cat('BLUE not estimable due to only one Range','\n')
            H2 <- data.frame(Estimate=NA,SE=NA,Location=unique(datdf$Location)[1],
                             Field=unique(datdf$Field)[1],Type=datdf$Type[1],N=nrow(datdf))
            blue <- data.frame(Name=datdf$Name,predicted.value=datdf$trait, std.error=NA,
                               status='unestimable',Location=unique(datdf$Location)[1],
                               Field=unique(datdf$Field)[1],Type=datdf$Type[1],
                               N=nrow(datdf),converge='NA')
            print(head(blue))
            asr <- NULL
            title <- paste(trait,'BLUE is not estimable:  Location =',
                           unique(datdf$Location)[1],
                           ' Field =',unique(datdf$Field)[1],
                           'N =',nrow(datdf),
                           '\n')
            
        } else {
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
                asr <- asreml(trait ~ Pass,random = ~idv(Name),
                              residual = ~ar1(Range):ar1(Pass), data = datdf)
                
            } else {
                cat('asr <- asreml(trait ~ Range + Pass,random = ~idv(Name),
                              residual = ~ar1(Range):ar1(Pass), data = datdf)\n')
                asr <- asreml(trait ~ Range + Pass,random = ~idv(Name),
                              residual = ~ar1(Range):ar1(Pass), data = datdf)
                
            }
    
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
        }
        
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

filename <- 'datafiles/22YTSu-filte.xlsx'
sheets <- c('P')
outfile <- 'Output/blue_ML.csv'

# 22P-SU.xlsx -------------------------------------------------------

vartype <- c('numeric','text',rep('numeric',5),rep('text',8),
             rep('numeric',12),'text',rep('text',2))
d <- data.frame(read_excel(filename,sheet = sheets[1],col_type=vartype))
names(d)[names(d)=='name'] <- 'Name'
names(d)[names(d)=='P.'] <- 'Phash'
names(d)[names(d)=='E.'] <- 'Ehash'
names(d)[names(d)=='PMD.'] <- 'PMD%'
names(d)[names(d)=='LDG.'] <- 'LDG%'


d <- d[order(d$Location,d$Type,d$Field,d$Range,d$Pass),]
table(d$Location) 
table(d$Field,d$Type)
ped <- unique(d[,c('Name','FGenoID','MGenoID')])

dim(ped)  #1999
d$Type[d$Type %in% c('PME','Pshort')] <- 'PME-Pshort'
unique(d$Type)

d <- d[!is.na(d$Range) & !is.na(d$Pass),]
dim(d) #26630

#品种名含有HNXZLL的数据排除
d <- d[substr(d$Name,1,6)!='HNXZLL',]

traits <- c('MST','TWT','YLD','PHT','EHT',#'Phash','Ehash','PMDhash','LDGhash',
            'PMD%','LDG%')
trait <- 'YLD'

summ <- descriptivePlus(dat = d,groupedby = c('Location'),
                        traitName = traits)

# Trait Location    N        Mean         SD        Min       Max
# 1    MST    7HBHD 2354  27.6690102   3.846965   4.400000   39.9100
# 2    MST   7HBSJZ   95  32.5646316   2.177354  25.230000   35.2500
# 3    MST    7HBXL 3402  30.0960817   2.975760   3.510000   35.4500
# 4    MST    7SDCQ 2100  23.4753857   4.190668   0.000000   33.6600
# 5    MST    7SDQD 2642  26.4373883   2.956409  11.110000   35.4000
# 6    MST    8AHSZ 1977  18.5983359   4.096569   4.540000   35.6000
# 7    MST    8HNAY 2870  25.5903937   3.679044   4.420000   37.7100
# 8    MST    8HNSQ 2690  16.7167695   3.281781   6.910000   25.1600
# 9    MST    8HNXM 2391  22.0746215   3.868924   0.000000   46.0000
# 10   MST    8HNXX 2306  25.3581526   2.659301  15.700000   37.6000
# 11   MST   8HNZMD 3155  19.4012456   2.702959   7.050000   29.6100
# 12   TWT    8HNXM 2291  72.5189873   3.488657  61.400000   80.2000
# 13   TWT    8HNXX 2305  71.2487852   2.123267  56.300000   89.6000
# 14   YLD    7HBHD 2354 530.3174444 185.351044  30.893992 1053.3047
# 15   YLD   7HBSJZ   95 717.7922027  90.032677 465.336105  943.6547
# 16   YLD    7HBXL 3402 710.0052593 169.175533 126.903317 1194.2176
# 17   YLD    7SDCQ 2100 568.6677863 157.525846   5.092584 1048.6543
# 18   YLD    7SDQD 2642 619.7884565 124.972948  78.077042 1109.8900
# 19   YLD    8AHSZ 1977 492.5749952 109.870955  54.293800  842.8142
# 20   YLD    8HNAY 2870 648.2418254 110.756098 219.798475 1062.5310
# 21   YLD    8HNSQ 2690 637.5112949 129.148289 156.788303 1122.2748
# 22   YLD    8HNXM 2391 571.9216659 116.507143  89.540097  956.4056
# 23   YLD    8HNXX 2306 542.8356869  91.243061   0.000000  849.5118
# 24   YLD   8HNZMD 3155 645.2000754 116.057343 110.538446  973.3753
# 25   PHT    7HBHD 2228 287.7890485  18.803415 198.000000  385.0000
# 26   PHT   7HBSJZ   96 292.9270833  18.498503 246.000000  338.0000
# 27   PHT    7HBXL 3402 277.3665491  22.500050 198.000000  388.0000
# 28   PHT    7SDCQ 1964 299.6664969  20.767085 220.000000  362.0000
# 29   PHT    7SDQD 2648 304.9856495  21.723991 215.000000  365.0000
# 30   PHT    8AHSZ 1786 282.7469205  22.046368 200.000000  345.0000
# 31   PHT    8HNXM 2282 291.3716039  18.474278 215.000000  346.0000
# 32   PHT    8HNXX 2739 277.5735670  24.247159 172.000000  894.0000
# 33   PHT   8HNZMD 3036 286.8152174  23.062292 207.000000  377.0000
# 34   EHT    7HBHD 2228 123.5570018  14.110917  72.000000  316.0000
# 35   EHT   7HBSJZ   96 115.2604167  13.550098  86.000000  149.0000
# 36   EHT    7HBXL 3399 110.1959400  14.908446  50.000000  177.0000
# 37   EHT    7SDCQ 1964 111.1155804  14.412891  64.000000  161.0000
# 38   EHT    7SDQD 2648 110.8191088  13.577526  66.000000  158.0000
# 39   EHT    8AHSZ 1784  99.3587444  13.197461  56.000000  177.0000
# 40   EHT    8HNXM 2282 105.9907975  13.769264  59.000000  170.0000
# 41   EHT    8HNXX 2738 110.2954711  15.256879  19.000000  165.0000
# 42   EHT   8HNZMD 3036 108.7193676  12.866196  64.000000  170.0000
# 43   PMD    8HNSQ 2743   4.0962450   7.488396   0.000000   50.0000
# 44   PMD    8HNXM  926   6.8131749   7.116990   1.000000   48.0000
# 45   PMD    8HNXX 2308  12.6083189  12.130556   0.000000   47.0000
# 46   LDG    7HBHD 2417  24.0579230  15.740548   0.000000   53.0000
# 47   LDG   7HBSJZ   96   6.5208333   8.974794   2.000000   52.0000
# 48   LDG    7HBXL 3402   0.9188713   4.503245   0.000000   50.0000
# 49   LDG    7SDCQ 2099   6.1796093  12.066966   0.000000   48.0000
# 50   LDG    7SDQD 2648   8.4418429  12.443982   0.000000   50.0000
# 51   LDG    8AHSZ 1997   4.4536805  12.661837   0.000000   52.0000
# 52   LDG    8HNAY 2902   0.6078567   1.669086   0.000000   28.0000
# 53   LDG    8HNSQ 2743  20.4261757  16.923507   0.000000   53.0000
# 54   LDG    8HNXM 2291   3.4980358   6.302448   0.000000   52.0000
# 55   LDG    8HNXX  432   0.0000000   0.000000   0.000000    0.0000
# 56   LDG   8HNZMD 3155   6.2614897   9.151672   0.000000   47.0000

summ2 <- dcast(summ, formula= Location ~ Trait, 
               value.var = 'N')
#summ2[is.na(summ2)] <- ''
summ2


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
loc <- 3
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
    if (trait=='MST')w$abnormal[w$trait>=43] <- 2
    if (trait=='YLD')w$abnormal[w$trait<200] <- 2
    if (trait=='YLD')w$abnormal[w$MST>=43] <- 2
    if (trait=='Phash')w$abnormal[w$trait>=43] <- 2
    title <- paste(trait,'Phenotype ','N =',nrow(w[!is.na(w$trait),]),'\n')
    
    qqnorm(w$trait,main=title,col=w$abnormal)
    qqline(w$trait)
    #hist(w$trait,main=title)
    
    
    if (trait=='MST')w$trait[w$trait>=43] <- NA
    if (trait=='YLD')w$trait[w$trait<200] <- NA
    if (trait=='YLD')w$trait[w$MST>=43] <- NA
    if (trait=='Phash')w$trait[w$trait>=43] <- NA
    
    w$trait[w$trait<mn-std*3 | w$trait>mn+std*3] <- NA
                         
    title <- paste(trait,'Phenotype ','N =',nrow(w[!is.na(w$trait),]),
                   ' Outliers removed\n')

    qqnorm(w$trait,main=title)
    qqline(w$trait)
    #hist(w$trait,main=title)
}

trait <- 'LDG'
cat('\n\nBLUE calculated for each field and outliers checked separately\n')
for (trait in traits){
    dd <-  dat[!is.na(dat[[trait]]),]
    locations <- unique(dd$Location)
    locations <- locations[order(locations)]
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
                
                if (nrow(w[!is.na(w[[trait]]),])>0){
                    #print(table(w$Range,w$Pass))
                    ww <- missingAdding(w)
                    ww <- ww[order(ww$Range,ww$Pass),]
                    #print(table(ww$Range,ww$Pass))
    
                    #for (trait in traits){
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
                    if (trait=='Phash')datw$trait[datw$trait>=43] <- NA
                    

                #if (nrow(datw[!is.na(datw$trait),])>40){
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
                    
                    blue <- blue[substr(blue$Name,1,5)!='DUMMY',]
                    
                    blueout <- blueEstimation(datdf=datw,
                                              trait=trait,repct = repct)
                    blue[[trait]] <- rbind(blue[[trait]],blueout$blue)
                    h2[[trait]] <- rbind(h2[[trait]],blueout$H2)
                #}
                } else {
                    cat('No observations available.\n')
                }
            }
        }
    }
}

repctage <- as.data.frame(repctage)
names(repctage) <- c('Location','Field','Type','Repct')
repctage$Repct <- as.numeric(repctage$Repct)
bluefile <- 'Output/blue_SU.xlsx'
h2file <- 'Output/h2_SU.xlsx'
rdsfile <- 'Output/blue_SU.rds'
writexl::write_xlsx(blue,path=bluefile)
writexl::write_xlsx(h2,path=h2file)
saveRDS(blue,rdsfile)

ped <- unique(dat[,c('Name','FGenoID','MGenoID')])
dim(ped)

pedfile <- 'Output/blue_SU_ped.rds'
saveRDS(ped,pedfile)

writexl::write_xlsx(ped,path='Output/pedigree_SU.xlsx')





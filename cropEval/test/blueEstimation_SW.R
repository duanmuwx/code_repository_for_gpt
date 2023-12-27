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
            
            if(nrow(blue[!is.na(blue$predicted.value),])!=0){
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
                blue$Type <- unique(datdf$Type)[1]
                blue$N <- nrow(datdf)
                blue$converge <- asr$converge
                
            }
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
    head(blue)
    head(H2)
    outlist <- list(H2=H2,blue=blue,asr=asr)
    cat('ASREML estimation is successful...\n')
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

filename <- 'datafiles/22YTSW-C.xlsx'
sheets <- c('22YTSW')
outfile <- 'Output/blue_SW.csv'

# 22P-SU.xlsx -------------------------------------------------------

vartype <- c('numeric','text',rep('numeric',5),rep('text',8),
             rep('numeric',16))
d <- data.frame(read_excel(filename,sheet = sheets[1],col_type=vartype))
names(d)[names(d)=='P.'] <- 'Phash'
names(d)[names(d)=='E.'] <- 'Ehash'
names(d)[names(d)=='PMD.'] <- 'PMDhash'
names(d)[names(d)=='ER.'] <- 'ERhash'


d <- d[order(d$Location,d$Type,d$Field,d$Range,d$Pass),]
table(d$Location)
# 9GZLPS  9SCST  9YNKM  9YNQJ 
# 323    312    312    312
table(d$Field,d$Type)
ped <- unique(d[,c('Name','FGenoID','MGenoID')])

dim(ped)  #1999
# d$Type[d$Type %in% c('PME','Pshort')] <- 'PME-Pshort'
unique(d$Type)

#d <- d[!is.na(d$Range) & !is.na(d$Pass),]
dim(d) #26630

#品种名含有HNXZLL的数据排除

traits <- c('MST','YLD','PHT','EHT','NCLB','Rust','SCLB','ER') #'PMD','FEW','ER',
trait <- 'YLD'

summ <- descriptivePlus(dat = d,groupedby = c('Location'),
                        traitName = traits)
summ
#    Trait Location   N       Mean          SD        Min       Max
# 1    MST   9GZLPS 323  31.533437   3.4252718  22.620000   40.8200
# 2    MST    9SCST 311   9.405402   1.8962718   4.050000   17.7100
# 3    MST    9YNKM 308  28.747078   2.6761508  17.900000   38.6000
# 4    MST    9YNQJ 312  19.256731   2.5801280  14.100000   28.3000
# 5    YLD   9GZLPS 323 672.457897 107.6204231 415.844252  983.1005
# 6    YLD    9SCST 311 584.666790  95.8275361 306.542730  838.3475
# 7    YLD    9YNKM 308 649.000347 147.8349850 318.848015 1187.5111
# 8    YLD    9YNQJ 312 718.317358 151.8537217 326.354310 1101.1265
# 9    PHT   9GZLPS 323 220.526316  18.0028057 165.000000  270.0000
# 10   PHT    9SCST 312 299.375000  23.3457946 230.000000  375.0000
# 11   PHT    9YNKM 312 255.157051  26.1324418   9.000000  325.0000
# 12   PHT    9YNQJ 312 285.858974  20.5294702 230.000000  350.0000
# 13   EHT   9GZLPS 323  74.752322  12.2576109  45.000000  110.0000
# 14   EHT    9SCST 312 123.733974  16.3291033  55.000000  170.0000
# 15   EHT    9YNKM 312  97.852564  18.1509584  60.000000  230.0000
# 16   EHT    9YNQJ 312 125.756410  16.7129623  80.000000  190.0000
# 17   PMD   9GZLPS 323   8.745787   0.6753760   2.263158    9.0000
# 18   PMD    9SCST 312   8.232146   1.6617533   1.000000    9.0000
# 19   FEW    9YNKM 308   8.644026   2.0302003   4.190000   16.2100
# 20  NCLB    9YNKM 312   7.810897   1.2878894   3.000000    9.0000
# 21    ER    9YNKM 308   8.636023   0.6556621   4.459459    9.0000
# 22    ER    9YNQJ 312   7.820513   1.1812424   3.000000    9.0000
# 23  Rust    9YNQJ 312   3.608974   1.7051981   1.000000    7.0000
# 24  SCLB    9YNQJ 312   5.907051   1.8706607   3.000000    9.0000

summ2 <- dcast(summ, formula= Location ~ Trait, 
               value.var = 'N')
#summ2[is.na(summ2)] <- ''
summ2

#   Location EHT  ER MST NCLB PHT Rust SCLB YLD
# 1   9GZLPS 323  NA 323   NA 323   NA   NA 323
# 2    9SCST 312  NA 311   NA 312   NA   NA 311
# 3    9YNKM 312 308 308  312 312   NA   NA 308
# 4    9YNQJ 312 312 312   NA 312  312  312 312
dat <- d

dat <- dat[order(dat$Location,dat$Type,dat$Field,
                 dat$Range,dat$Pass),]

locations <- unique(dat$Location)
table(dat$Field,dat$Type)



cat('Outliers checking\n')
for (trait in traits){
    w <- dat
    names(w)[names(w)==trait] <- 'trait'
    
    std=sd(w$trait,na.rm = T)
    mn=mean(w$trait,na.rm = T)
    w$abnormal <- ifelse(w$trait<mn-std*3 | w$trait>mn+std*3,
                            2,1)
    if (trait=='MST')w$abnormal[w$trait>=43 | w$trait<10] <- 2
    if (trait=='YLD')w$abnormal[w$trait<200] <- 2
    if (trait=='YLD')w$abnormal[w$MST>=43] <- 2
    if (trait=='Phash')w$abnormal[w$trait>=43] <- 2
    title <- paste(trait,'Phenotype ','N =',nrow(w[!is.na(w$trait),]),'\n')
    
    qqnorm(w$trait,main=title,col=w$abnormal)
    qqline(w$trait)
    #hist(w$trait,main=title)
    
    
    if (trait=='MST')w$trait[w$trait>=43| w$trait<10] <- NA
    if (trait=='YLD')w$trait[w$trait<200] <- NA
    if (trait=='YLD')w$trait[w$MST>=43] <- NA
    if (trait=='Phash')w$trait[w$trait>=43] <- NA
    
    if (!trait %in% c('NCLB')) w$trait[w$trait<mn-std*3 | w$trait>mn+std*3] <- NA
                         
    title <- paste(trait,'Phenotype ','N =',nrow(w[!is.na(w$trait),]),
                   ' Outliers removed\n')

    qqnorm(w$trait,main=title)
    qqline(w$trait)
    #hist(w$trait,main=title)
}


trait <- 'MST'
loc <- 2
i=1
j=1

repctage <- c()
h2 <- list()
blueall <- list()
for (trait in traits){
    blueall[[trait]] <- c()
    h2[[trait]] <- c()
}

trait <- 'LDG'
cat('\n\nBLUE calculated for each field and outliers checked separately\n')
for (trait in traits){
    dd <-  dat[!is.na(dat[[trait]]),]
    locations <- unique(dd$Location)
    locations <- locations[order(locations)]
    locations
    for (loc in 1:length(locations)){

        df <- dat[dat$Location==locations[loc],]
        fields <- unique(df$Field)
        fields <- fields[order(fields)]
        fields
        
        for (i in 1:length(fields)){
            dff <- df[df$Field==fields[i],]
            types <- unique(dff$Type)
            types
            
            for (j in 1:length(types)){
                cat('############ loc=',loc,'Field=',i,'Type=',j,'###########\n')
                w <- dff[dff$Type==types[j],]
                
                repct <- replicationCheck(w)
                repctage <- rbind(repctage,c(locations[loc],
                                             fields[i],types[j],repct))
                cat('\nLocation:',locations[loc],'Field:',
                    fields[i],'Type',types[j],'N:',nrow(w),
                    'Replications:',repct,'\n')
                # if (locations[loc]=='6JLDF' & fields[i]==2 &
                #     types[j]=='P1') {
                #     w$Range[w$Set==25] <- 50
                #     w$Range[w$Set==26] <- 51
                # }
                
                if (nrow(w[!is.na(w[[trait]]),])>0){
                    #print(table(w$Range,w$Pass))
                    ww <- missingAdding(w)
                    ww <- ww[order(ww$Range,ww$Pass),]
                    print(table(ww$Range,ww$Pass))
    
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
                    if (!trait %in% c('NCLB')) datw$abnormal <- 
                        ifelse(datw$trait<mn-std*3 | datw$trait>mn+std*3,2,1)
                    
                    datw$trait[datw$trait<mn-std*3] <- NA
                    datw$trait[datw$trait>mn+std*3] <- NA
                    if (trait=='MST')datw$trait[datw$trait>=43 | datw$trait<10] <- NA
                    if (trait=='YLD')datw$trait[datw$trait<200] <- NA
                    if (trait=='YLD')datw$trait[datw$MST>=43] <- NA
                    if (trait=='Phash')datw$trait[datw$trait>=43] <- NA
                    
                    if (!trait %in% c('NCLB')) {
                        datw$trait[datw$trait<mn-std*3 | datw$trait>mn+std*3] <- NA
                    }
                    
                #if (nrow(datw[!is.na(datw$trait),])>40){
                    title <- paste(trait,'Phenotype:  Location =',
                                   unique(datw$Location)[1],'Type =',types[j],
                                   ' Field =',datw$Field[1],
                                   'N =',nrow(datw),
                                   '\n')

                    qqnorm(datw$trait,main=title)
                    qqline(datw$trait)
                    #hist(datw$trait)

                    cat(types[j],trait,nrow(datw[!is.na(datw$trait),]),
                        '\n')
                    head(datw)
                    
                    blueout <- blueEstimation(datdf=datw,
                                              trait=trait,repct = repct)
                    head(blueout$blue)
                    head(blueout$H2)
                    blueall[[trait]] <- rbind(blueall[[trait]],blueout$blue)
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
writexl::write_xlsx(blueall,path=bluefile)
writexl::write_xlsx(h2,path=h2file)
saveRDS(blueall,rdsfile)

ped <- unique(dat[,c('Name','FGenoID','MGenoID')])
dim(ped)

pedfile <- 'Output/blue_SU_ped.rds'
saveRDS(ped,pedfile)

writexl::write_xlsx(ped,path='Output/pedigree_SU.xlsx')





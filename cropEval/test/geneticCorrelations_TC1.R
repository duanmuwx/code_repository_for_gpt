extracting <- function(asr,var,individual,term,console=FALSE) {
    ebv <- summary(asr,coef=T)$coef.random
    ebv <- as.data.frame(ebv)
    bvalue <- as.data.frame(ebv[grepl(term,row.names(ebv)),])
    
    head(bvalue)
    dim(bvalue)
    #print(dim(bvalue))
    bvalue$id <- 0
    bvalue$id <- 0
    for (i in 1:dim(bvalue)[1]){
        #print(i)
        s <- strsplit(row.names(bvalue)[i],"[_]")[[1]]
        bvalue$id[i] <- s[2]
    }
    v <- (1-bvalue$'std.error'^2/var)
    bvalue$rel <- ifelse(v<=0,0,v)
    bvalue$acc <- sqrt(bvalue$rel)
    
    
    bvalue <- bvalue[bvalue$id %in% individual,]
    
    bvalue <- bvalue[,c(4,1,2,5,6)]
    names(bvalue) <- c("ID","EBV","se","rel","acc")
    if (console) {
        print(head(bvalue))
        cat('No of individuals:',dim(bvalue),'\n')
    }
    return(bvalue)
    
}

library(readxl)
library(writexl)
library(biosim)
library(tidyverse)
library(reshape2)
library(stringi) 
library(asreml)
library(openxlsx)
library(ggplot2)

library(biosim)
library(matrixcalc)
library(asreml)
library(MCMCglmm)
library(ggpmisc)


pedfile <- 'F:/projects/geneticAnalysis/datafiles/TC1_ped.xlsx'
ped <- data.frame(read_excel(pedfile,sheet = 1))
names(ped) <- c('ID','Dam','Sire')
ped <- ped[!(ped$ID=='22TC1XM13626' & ped$Dam==0),]  #Duplicated, but pedigree info are not exactly the same. Remove one and set Dam=0
ped$Dam[ped$ID=='22TC1XM13626'] <- '0'
ped[ped$ID=='22TC1XM13626',]
ped$ID <- toupper(ped$ID)
ped$Sire <- toupper(ped$Sire)
ped$Dam <- toupper(ped$Dam)
ped$ID[ped$ID=='C7'] <- 'C7-2'
ped$Dam[ped$Dam=='C7'] <- 'C7-2'

# ainv <- ainverse(ped)
# saveRDS(ainv,   'F:/projects/geneticAnalysis/Output/ainverse_TC1.rds')
ainv <- readRDS('F:/projects/geneticAnalysis/Output/ainverse_TC1.rds')

parents <- data.frame(read_excel('F:\\projects\\pedReconst\\output\\TC1_parents.xlsx',sheet = 1))
control <- parents[is.na(parents$FGenoID),]
newcol <- stringsplit(control$GenoID,'/')
control$FGenoID <- newcol$X1
control$MGenoID <- newcol$X2
control
parents$FGenoID[parents$name=='DH605'] <- control$FGenoID[control$name=='DH605']
parents$MGenoID[parents$name=='DH605'] <- control$MGenoID[control$name=='DH605']
parents$FGenoID[parents$name=='DK159'] <- control$FGenoID[control$name=='DK159']
parents$MGenoID[parents$name=='DK159'] <- control$MGenoID[control$name=='DK159']
parents$FGenoID[parents$name=='XY335'] <- control$FGenoID[control$name=='XY335']
parents$MGenoID[parents$name=='XY335'] <- control$MGenoID[control$name=='XY335']
parents$FGenoID[parents$name=='ZD958'] <- control$FGenoID[control$name=='ZD958']
parents$MGenoID[parents$name=='ZD958'] <- control$MGenoID[control$name=='ZD958']
parents <- parents[parents$name!='Filler',]
dim(parents) #6047

traits <- c('MST','TWT','YLD','PHT','EHT','NCLB','GLS') #,'PMD','PMD%','LDG','LDG%'

trait <- 'YLD'

bluelist <- readRDS('F:/projects/geneticAnalysis/Output/blue_TC1.rds')
# for (trait in traits){
blue <- bluelist[[trait]]
blue$Name <- toupper(blue$Name)
# head(blue)
blue$Field <- as.factor(blue$Field)
blue$LF <- paste0(blue$Location,'_',blue$Field)
blue$LF <- as.factor(blue$LF)
blue$Name <- toupper(blue$Name)
blue$Name <- as.factor(blue$Name)
blue$site <- ifelse(blue$Location %in% c("6JLGZL","7LNHC","7LNTL"),
                    'Spring','Summer')
blue$site <- factor(blue$site)
blue <- blue[order(blue$site),]

blue2 <- merge(blue,ped,by.x='Name',by.y='ID',all.x=T)
blue2$Sire <- factor(blue2$Sire)
blue2$Dam <- factor(blue2$Dam)

# filename <- 'F:/projects/geneticAnalysis/datafiles/22TC1.xlsx'
# sheets <- c('TC1-yld-HC')
# 
# vartype <- c('numeric','text',rep('numeric',5),rep('text',8),
#              rep('numeric',8))
# yld <- data.frame(read_excel(filename,sheet = sheets[1],col_type=vartype))
# yld <- yld[,1:21]
# yld$Location <- factor(yld$Location)
# yld$Range <- factor(yld$Range)
# yld$Pass <- factor(yld$Pass)
# yld$FGenoID <- factor(toupper(yld$FGenoID))
# yld$MGenoID <- factor(toupper(yld$MGenoID))
# yld$LF <- paste0(yld$Location,'_',yld$Field)
# yld$LF <- factor(yld$LF)
# yld$name <- toupper(yld$name)


traits <- c('MST','TWT','YLD','PHT','EHT','NCLB','GLS')
trait <- 'YLD'

idlist <- c()
for (i in 1:(length(traits))){
    blue <- bluelist[[traits[i]]] 
    blue$ID <- paste0(blue$Name,'_',blue$Location,'_',blue$Field,
                       '_',blue$Type)
    idlist <- unique(c(idlist,blue$ID))
}
length(idlist)
blues <- stringsplit(idlist,'_')
names(blues) <- c('Name','Location','Field')

blues <- merge(ped,blues,by.x = 'ID',by.y='Name',all.y=T)
blues <- blues[blues$ID!='',]
names(blues)[1] <- 'Name'
# blues$Fam <- paste0(blues$SIRE,'_',blues$DAM)
# fam <- unique(blues[,c('Fam','SIRE','DAM')])
# fam$FamNum <- 1:nrow(fam)
# blues <- merge(blues,fam[,c('Fam','FamNum')],by='Fam')
blues$ID <- paste0(blues$Name,'_',blues$Location,'_',blues$Field)
head(blues)

for (i in 1:(length(traits))){
    print(traits[i])
    blue <- bluelist[[traits[i]]] 
    names(blue)[2] <- traits[i]
    blue$ID <- paste0(blue$Name,'_',blue$Location,'_',blue$Field)
    
    blues <- merge(blues,blue[,c('ID',traits[i])],by='ID',all.x=T)
}

head(blues)


trait <- 'YLD'

genetCorr <- c()
corr <- as.data.frame(matrix(1,length(traits),length(traits)))
names(corr) <- traits
row.names(corr) <- traits
for (i in 1:(length(traits)-1)){
    for (j in (i+1):length(traits)){
        trait1 <- traits[i]
        trait2 <- traits[j]
        cat(trait1,trait2,'\n')
        
        blue <- blues[,c('Name','Location','Field',trait1,trait2)]
        blue2 <- bluelist[[traits[j]]]
        names(blue)[4] <- 'trait1'
        names(blue)[5] <- 'trait2'

        blue$Name <- as.factor(blue$Name)
        blue$Field <- as.factor(blue$Field)
        blue$Location <- as.factor(blue$Location)
        
        mod = asreml(cbind(trait1,trait2) ~ trait + trait:Location +
                         trait:Field, 
                     random = ~ us(trait):vm(Name,ainv),
                     residual = ~ units:us(trait), 
                     data=blue)
        mod <- update(mod)
        summary(mod)$varcomp
        rg <- vpredict(mod,rg~V2/sqrt(V1*V3))
        rg$trait1 <- trait1
        rg$trait2 <- trait2
        rg$converge <- mod$converge
        genetCorr <- rbind(genetCorr,rg)
        corr[i,j] <- rg$Estimate[1]
        corr[j,i] <- rg$Estimate[1]
    }
    
}
corr
write.csv(corr,'Output/geneticCorrelation_TC1_matrix.csv',quote = F,
          row.names = F)
write.csv(genetCorr,'Output/geneticCorrelation_TC1.csv',quote = F,
          row.names = F)
saveRDS(corr,'Output/geneticCorrelation_TC1_matrix.rds')

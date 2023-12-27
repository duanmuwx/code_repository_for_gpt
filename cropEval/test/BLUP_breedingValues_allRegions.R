
rm(list=ls())

extracting <- function(asr,var,individual,term) {
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
    #print(head(bvalue))
    
    print(dim(bvalue))
    return(bvalue)
    
}

GCA_SCA_Estimation <- function(bluehb,ainv,ped){
    
    bluehb <- merge(bluehb,ped,by='Name',all.x=T)
    
    bluehb$SIRE <- as.factor(bluehb$SIRE)
    bluehb$DAM <- as.factor(bluehb$DAM)  
    bluehb$SIRE[is.na(bluehb$SIRE)] <- '0'
    bluehb$DAM[is.na(bluehb$DAM)] <- '0'

    
    bluehb$Fam <- paste0(bluehb$SIRE,'_',bluehb$DAM)
    
    bluehb$Fam[bluehb$Fam=='NA_NA'] <- '0_0'
    
    fam <- unique(bluehb[,c('Fam','SIRE','DAM')])
    fam$FamNum <- 1:nrow(fam)
    bluehb <- merge(bluehb,fam[,c('Fam','FamNum')],by='Fam')
    head(bluehb)
    #writecsv(bluehb,'temp.csv')
    bluehb$FamNum <- as.factor(bluehb$FamNum)
    bluehb$Location  <- as.factor(bluehb$Location )
    bluehb$Field  <- as.factor(bluehb$Field)
    bluehb$LF <- paste0(bluehb$Location,'_',bluehb$Field)
    bluehb$LF  <- as.factor(bluehb$LF)
    
    asr2 <- asreml(predicted.value ~ LF,
                   random = ~ vm(SIRE,ainv) + vm(DAM,ainv) + idv(FamNum), 
                   residual = ~idv(units),
                   data = bluehb)
    
    summary(asr2)$varcomp
    vs <- summary(asr2)$varcomp$component[2]
    vd <- summary(asr2)$varcomp$component[3]
    vf <- summary(asr2)$varcomp$component[1]
    
    h2 <- vpredict(asr2,h2~(V2+V3)/(V1+V2+V3+V4))
    print(h2)
    
    d2 <- vpredict(asr2,h2~V1/(V1+V2+V3+V4))
    print(d2)
    loglik2 <- asr2$loglik
    asr3 <- asreml(predicted.value ~ LF,
                   random = ~ vm(SIRE,ainv) + vm(DAM,ainv),# + idv(FamNum), 
                   residual = ~idv(units),
                   data = bluehb)
    loglik3 <- asr3$loglik
    
    d2$LRTvalue <- (loglik2-loglik3)*2
    
    bvalue <- summary(asr2,coef=T)$coef.random
    head(bvalue)
    #write.csv(bvalue,'temp.csv')
    ind <- unique(c(bluehb$SIRE,bluehb$DAM))
    ebv <- extracting(asr2,var=vs,individual = ind,term='vm(Name, ainv)*')
    head(ebv)
    febv <- extracting(asr2,var=vf,individual = fam$FamNum,term='FamNum*')
    febv <- merge(fam,febv,by.x='FamNum',by.y='ID')
    febv <- febv[order(-febv$EBV),]
    head(febv)
    #write.csv(febv,'SCA.csv')
    
    ind <- unique(c(bluehb$SIRE))
    sebv <- extracting(asr2,var=vs,individual = ind,term='vm(SIRE, ainv)*')
    sebv <- sebv[sebv$EBV!=0,]
    sebv$par='SIRE'
    #sebv <- sebv[1:315,]
    ind <- unique(c(bluehb$DAM))
    debv <- extracting(asr=asr2,var=vd,individual = ind,term='vm(DAM, ainv)*')
    debv <- debv[debv$EBV!=0,]
    debv$par='DAM'
    
    sebv <- sebv[!substr(row.names(sebv),1,6)=='vm(DAM',]
    debv <- debv[!substr(row.names(debv),1,7)=='vm(SIRE',]
    
    outlist <- list(h2=h2,d2=d2,varcomp=summary(asr2)$varcomp,Mebv=sebv,Febv=debv,
                    SCAebv=febv)
    return(outlist)
}

programRun <- function(bluelist,ainv,ped){
    h2 <- c()
    d2 <- c()
    varc <- c()
    SCAebv <- c()
    
    Mebv <- c()
    Febv <- c()
    
    for (i in 1:length(traits)){
        trait <- traits[i]
        blue <- bluelist[[traits[i]]]
        asrout <- GCA_SCA_Estimation(bluehb=blue,ainv=ainv,ped=ped)
        
        asrout$h2
        asrout$d2
        asrout$varcomp
        if (i==1){
            SCAebv <- asrout$SCAebv[,c(1:6,8)]
            Mebv <- asrout$Mebv[,c(1:3,5)]
            Febv <- asrout$Febv[,c(1:3,5)]
            names(SCAebv)[names(SCAebv)=='EBV'] <- trait
            names(Mebv)[names(Mebv)=='EBV'] <- trait
            names(Febv)[names(Febv)=='EBV'] <- trait
            names(SCAebv)[names(SCAebv)=='se'] <- paste0(trait,'se')
            names(Mebv)[names(Mebv)=='se'] <- paste0(trait,'se')
            names(Febv)[names(Febv)=='se'] <- paste0(trait,'se')
            
            names(SCAebv)[names(SCAebv)=='acc'] <- paste0(trait,'acc')
            names(Mebv)[names(Mebv)=='acc'] <- paste0(trait,'acc')
            names(Febv)[names(Febv)=='acc'] <- paste0(trait,'acc')
            
        } else {
            names(asrout$SCAebv)[names(asrout$SCAebv)=='EBV'] <- trait
            names(asrout$Mebv)[names(asrout$Mebv)=='EBV'] <- trait
            names(asrout$Febv)[names(asrout$Febv)=='EBV'] <- trait
            
            names(asrout$SCAebv)[names(asrout$SCAebv)=='se'] <- paste0(trait,'se')
            names(asrout$Mebv)[names(asrout$Mebv)=='se'] <- paste0(trait,'se')
            names(asrout$Febv)[names(asrout$Febv)=='se'] <- paste0(trait,'se')
            
            names(asrout$SCAebv)[names(asrout$SCAebv)=='acc'] <- paste0(trait,'acc')
            names(asrout$Mebv)[names(asrout$Mebv)=='acc'] <- paste0(trait,'acc')
            names(asrout$Febv)[names(asrout$Febv)=='acc'] <- paste0(trait,'acc')
            
            SCAebv <- merge(SCAebv,asrout$SCAebv[,c(1,5,6,8)],by='FamNum')
            
            Mebv<-merge(Mebv,asrout$Mebv[,c(1:3,5)],by='ID')
            Febv<-merge(Febv,asrout$Febv[,c(1:3,5)],by='ID')
        }
        
        asrout$h2$trait <- trait
        asrout$d2$trait <- trait
        asrout$varcomp$trait <- trait
        #filename=paste0("Febv_",trait,".csv")
        #write.csv(asrout$Febv,filename)
        #filename=paste0("Mebv_",trait,".csv")
        #write.csv(asrout$Mebv,filename)
        
        h2 <- rbind(h2,asrout$h2)
        d2 <- rbind(d2,asrout$d2)
        varc <- rbind(varc,asrout$varcomp)
    }
    
    #dim(ebv)
    h2
    d2
    #head(febv)
    varc$item <- row.names(varc)
    varc <- varc[,c(6,7,1:5)]
    
    varc
    exout <- list(h2=h2,d2=d2,varcomp=varc,Mebv=Mebv,Febv=Febv,SCAebv=SCAebv)
    
    return(exout)
}

library(biosim)
library(matrixcalc)
library(asreml)
library(MCMCglmm)
library(openxlsx)
library(readxl)

#setwd(choose.dir())

ainv <- readRDS('Output/pedigree_P1_inbredlines_ainverse.rds')
ped <- readRDS('Output/pedigree_P1_inbredlines.rds')
names(ped)[names(ped)=='FGenoID'] <- 'SIRE'
names(ped)[names(ped)=='MGenoID'] <- 'DAM'
head(ped)

bluedat <- readRDS('Output/bluelist_Pexp.rds')
regions <- names(bluedat)
regions
#[1] "EM" "ML" "SU" "SW"

ebv <- list()

for (i in 1:3){
    reg <- regions[i]  
    if (reg=='EM') traits <- c('MST','YLD','PHT','EHT',"LDG","NCLB")
    if (reg=='ML') traits <- c('MST','TWT','YLD','PHT','EHT',"GLS")
    if (reg=='SU') traits <- c('MST','TWT','YLD','PHT','EHT',"PMD%")
    if (reg=='SW') traits <- c('MST','YLD','PHT','EHT')#"NCLB","Rust","SCLB","Rust","SCLB","
    print(traits)
    ebv[[reg]] <- programRun(bluelist=bluedat[[reg]],ainv,ped)
    
    writexl::write_xlsx(ebv[[reg]],path=paste0('Output/blup_EBV_',reg,'.xlsx'))
  
}

saveRDS(ebv,'Output/BLUP_breedingValues_allRegions.rds')




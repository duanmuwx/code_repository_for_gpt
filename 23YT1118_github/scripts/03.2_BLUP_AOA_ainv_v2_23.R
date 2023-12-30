rm(list = ls())
GCA_SCA_Estimation <- function(bluehb, ainv, ped){
  
  bluehb <- merge(bluehb,ped,by='Name',all.x=T) # Merge Blue files and genealogy files by using left join
  #openxlsx::write.xlsx(bluehb, file = "Output/YLD_MMSP.xlsx")
 
  bluehb$SIRE <- as.character(bluehb$SIRE)
  bluehb$DAM <- as.character(bluehb$DAM)
  bluehb$SIRE[is.na(bluehb$SIRE)] <- '0' # Missing values are represented by 0
  bluehb$DAM[is.na(bluehb$DAM)] <- '0'   # Missing values are represented by 0
  bluehb$SIRE <- as.factor(bluehb$SIRE)
  bluehb$DAM <- as.factor(bluehb$DAM)
  
  
  bluehb$Fam <- paste0(bluehb$DAM, '/', bluehb$SIRE) # Connecting parent sample and parent sample.
  bluehb$Fam[bluehb$Fam=='NA/NA'] <- '0/0'           # Missing values are represented by 0_0.
  
  # Group according to paternal group and maternal group, and mark with numbers.
  fam <- unique(bluehb[,c("Fam",'DAM','SIRE')])
  fam$FamNum <- 1:nrow(fam)
  bluehb <- merge(bluehb, fam[,c('Fam','FamNum')], by='Fam')
  head(bluehb)
  #writecsv(bluehb,'temp.csv')
  
  bluehb$FamNum <- as.factor(bluehb$FamNum)
  bluehb$Location  <- as.factor(bluehb$Location )
  bluehb$Field  <- as.factor(bluehb$Field)
  bluehb$LF <- paste0(bluehb$Location,'_',bluehb$Field) # Combine Location and Field as one factor to get LF.
  bluehb$LF  <- as.factor(bluehb$LF)
  bluehb$Year  <- as.factor(bluehb$Year)
  # LF as a fixed factor, maternal and paternal and maternal_paternal as a random factor,
  # establish a linear model
  asr2 <- asreml(predicted.value ~ LF ,#+ Year
                 random = ~ vm(SIRE,ainv) + vm(DAM,ainv) + idv(FamNum),
                 residual = ~idv(units),
                 data = bluehb)
  iter <- 1
  while(!asr2$converge & (iter < 10)){
    message("更新 mgxe 迭代第 ", iter, " 次")
    asr2  <- update(asr2)
    iter <- iter + 1
  }
  # Extracting variance components from the model.
  fixed<-asr2$coefficients$fixed
  popmean <- fixed[row.names(fixed)=='(Intercept)',]
  
  # The generalized heritability(H2) is calculated according to the variance component of each factor.
  summary(asr2)$varcomp
  vs <- summary(asr2)$varcomp$component[2]
  vd <- summary(asr2)$varcomp$component[3]
  vf <- summary(asr2)$varcomp$component[1]
  h2 <- vpredict(asr2,h2~(V2+V3)/(V1+V2+V3+V4))
  print(h2)
  d2 <- vpredict(asr2,h2~V1/(V1+V2+V3+V4))
  print(d2)
  
  # The maximum likelihood value is obtained by establishing the model,
  # and whether the factor(FamNum) is significant is tested.
  loglik2 <- asr2$loglik
  asr3 <- asreml(predicted.value ~ LF ,#+ Year
                 random = ~ vm(SIRE,ainv) + vm(DAM,ainv),# + idv(FamNum),
                 residual = ~idv(units),
                 data = bluehb)
  loglik3 <- asr3$loglik
  d2$LRTvalue <- abs(loglik2-loglik3)*2
  
  
  #bvalue <- summary(asr2,coef=T)$coef.random
  #head(bvalue)
  #write.csv(bvalue,'temp.csv')

  
  # Extracting breeding value from the model, that is, non-additive effect
  febv <- extracting(asr2,var=vf,individual = fam$FamNum,term='FamNum*')
  febv <- merge(fam,febv,by.x='FamNum',by.y='ID')
  febv = full_join(febv,bluehb[,c(1,2)])
  febv=unique(febv)
  febv <- febv[order(-febv$EBV),]
  head(febv)
  #write.csv(febv,'SCA.csv')
  
  # The breeding value is extracted from the model, that is, the additive effect of the maternal line.
  ind <- unique(c(bluehb$SIRE))
  sebv <- extracting(asr=asr2,var=vs,individual = ind,term='SIRE*')
  sebv <- sebv[sebv$EBV!=0,]
  sebv$par='SIRE'
  #sebv <- sebv[1:315,]
  
  # The breeding value is extracted from the model, that is, the additive effect of the paternal line.
  ind <- unique(c(bluehb$DAM))
  debv <- extracting(asr=asr2,var=vd,individual = ind,term='DAM*')
  debv <- debv[debv$EBV!=0,]
  debv$par='DAM'
  
  #sebv <- sebv[!substr(row.names(sebv),1,6)=='vm(DAM',]
  #debv <- debv[!substr(row.names(debv),1,7)=='vm(SIRE',]
  
  outlist <- list(h2=h2,d2=d2,varcomp=summary(asr2)$varcomp,Mebv=sebv,Febv=debv,
                  SCAebv=febv,popmean=popmean)
  return(outlist)
}

extracting <- function(asr, var, individual, term, gxe=FALSE) {
  # Extract some important statistics from the model,
  # such as variance component, standard error, etc.
  ebv <- summary(asr, coef=T)$coef.random
  ebv <- as.data.frame(ebv)
  bvalue <- as.data.frame(ebv[grepl(term,row.names(ebv)),])
  
  head(bvalue)
  dim(bvalue)

  bvalue$id <- 0
  for (i in 1:dim(bvalue)[1]){
    s <- strsplit(row.names(bvalue)[i],"[_]")[[1]]
    if (gxe){
      bvalue$id[i] <- s[3]
    } else {
      bvalue$id[i] <- s[2]
    }
  }
  v <- (1-bvalue$'std.error'^2/var)
  bvalue$rel <- ifelse(v<=0,0,v)
  bvalue$acc <- sqrt(bvalue$rel)
  
  bvalue <- bvalue[bvalue$id %in% individual,]
  bvalue <- bvalue[,c(4,1,2,5,6)]
  row.names(bvalue) <- seq(1:nrow(bvalue))
  names(bvalue) <- c("ID","EBV","se","rel","acc")
  print(head(bvalue))
  
  return(bvalue)
  
}


    library(asreml)
    library(tidyverse)

    

    pedigree_P_inbredlines <- readRDS("datafiles/pedigree_19_23P_1031.rds")
    
    colnames(pedigree_P_inbredlines) <- c("Name", "DAM", "SIRE")
    ainv <- ainverse(pedigree_P_inbredlines)

    # all.blue <- readRDS("../data/19-22BLUE_AOA.rds")
    bluedat <- readRDS("Output/23P_datalist_bluelist.rds")
    all.blue = bluedat$bluelist$blue
    
    all_trait_name <- c("YLD14",	"MST","PHT","EHT")#,	"TWT"


    h2.all.trait <- NULL
    d2.all.trait <- NULL
    popMean.all.trait <- NULL
    popMean.tmp <- data.frame(popmean=0)
    varcomp.all.trait<-NULL

    
    for (i in 1:length(all_trait_name)){
        trait=all_trait_name[i]
        blue.value <- all.blue[[which(names(all.blue) == trait)]]
        #colnames(blue.value)[10] <- "Year"
        # blue.value$AOA[blue.value$Station == "CC"] <- "SP"
        # blue.value$AOA[blue.value$Station == "NW"] <- "SP"
        # blue.value$AOA[blue.value$Station == "JN"] <- "SU"
        # blue.value$AOA[blue.value$Station == "XX"] <- "SU"

        AOA.vec <- names(table(blue.value$AOA))
        #AOA.vec=c("EMSP","MMSP","LMSP")
        SCAebv.all.trait <- NULL
        Mebv.all.trait <- NULL
        Febv.all.trait <- NULL

        for (aoa in AOA.vec){
            cau.data <- blue.value %>% dplyr::filter(AOA %in% aoa)
            asrout <- GCA_SCA_Estimation(bluehb = cau.data, ainv = ainv, ped = pedigree_P_inbredlines)
            
            popMean.tmp$popmean <- asrout$popmean   
            popMean.tmp$AOA <- aoa
            popMean.tmp$Trait <- trait
            popMean.all.trait <- rbind(popMean.all.trait,popMean.tmp)
            
            varcomp.tmp<-asrout$varcomp
            varcomp.tmp$AOA <- aoa
            varcomp.tmp$Trait <- trait
            varcomp.all.trait <- rbind(varcomp.all.trait, varcomp.tmp)

            h2.tmp <-  asrout$h2
            h2.tmp$AOA <- aoa
            h2.tmp$Trait <- trait
            h2.all.trait <- rbind(h2.all.trait, h2.tmp)

            d2.tmp <-  asrout$d2
            d2.tmp$AOA <- aoa
            d2.tmp$Trait <- trait
            d2.all.trait <- rbind(d2.all.trait, d2.tmp)

            # popMean$popmean[i] <- asrout$popmean

            SCAebv <- asrout$SCAebv[,c(9,1:6,8)]
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

            SCAebv$AOA <- aoa
            Mebv$AOA <- aoa
            Febv$AOA <- aoa

            SCAebv.all.trait <- rbind(SCAebv.all.trait, SCAebv)
            Mebv.all.trait <- rbind(Mebv.all.trait, Mebv)
            Febv.all.trait <- rbind(Febv.all.trait, Febv)
       
        }

        outdata <- list(SCAebv.tmp = SCAebv.all.trait,
                        Mebv.tmp = Mebv.all.trait,
                        Febv.tmp = Febv.all.trait)

        #res.all.2 <- c(res.all.2, list(outdata))
        
        if (i==1){
          SCAebv.all <-outdata$SCAebv.tmp
          Mebv.all<-outdata$Mebv.tmp
          Febv.all<-outdata$Febv.tmp
        } else{
        SCAebv.all <- merge(SCAebv.all,outdata$SCAebv.tmp,by=c("Name","Fam","AOA","DAM","SIRE"),all=TRUE)
        SCAebv.all <- SCAebv.all[, -grep("^FamNum", names(SCAebv.all))]
        Mebv.all<-merge(Mebv.all,outdata$Mebv.tmp,by=c("ID","AOA"),all=TRUE)
        Febv.all<-merge(Febv.all,outdata$Febv.tmp,by=c("ID","AOA"),all=TRUE)
        }

        #names(outdata) <- paste0(trait, "-", names(outdata))

        #res.all <- c(res.all, outdata)

    }


    varcomp.all.trait$item<-row.names(varcomp.all.trait)
    
    idx<-list(h2 = cbind(Type = rep("h2", dim(h2.all.trait)[1]),h2.all.trait),
              d2 = cbind(Type = rep("d2", dim(d2.all.trait)[1]),d2.all.trait),
              popMean = cbind(Type = rep("popmean", dim(popMean.all.trait)[1]),popMean.all.trait),
              varcomp = varcomp.all.trait,
              SCAebv=SCAebv.all,
              Mebv=Mebv.all,
              Febv=Febv.all
              )


    #openxlsx::write.xlsx(idx, file = "Output/1031GCA_SCA_by_23.xlsx")

    saveRDS(idx, "Output/1118GCA_SCA_by_23.rds")






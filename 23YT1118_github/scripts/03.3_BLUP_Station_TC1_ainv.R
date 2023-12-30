rm(list=ls())
GCA_SCA_Estimation_TC <- function(bluehb, ainv, ped){
  
  bluehb <- merge(bluehb,ped,by='Name',all.x=T) # Merge Blue files and genealogy files by using left join
  #openxlsx::write.xlsx(bluehb, file = "Output/YLD_MMSP.xlsx")
  
  bluehb$SIRE <- as.character(bluehb$SIRE)
  bluehb$DAM <- as.character(bluehb$DAM)
  bluehb$SIRE[is.na(bluehb$SIRE)] <- '0' # Missing values are represented by 0
  bluehb$DAM[is.na(bluehb$DAM)] <- '0'   # Missing values are represented by 0
  bluehb$SIRE <- as.factor(bluehb$SIRE)
  bluehb$DAM <- as.factor(bluehb$DAM)
  
  
  bluehb$Location  <- as.factor(bluehb$Location )
  bluehb$Field  <- as.factor(bluehb$Field)
  bluehb$LF <- paste0(bluehb$Location,'_',bluehb$Field) # Combine Location and Field as one factor to get LF.
  bluehb$LF  <- as.factor(bluehb$LF)
  bluehb$Name <- as.factor(bluehb$Name)
  # LF as a fixed factor, maternal and paternal and maternal_paternal as a random factor,
  # establish a linear model
  asr2 <- asreml(predicted.value ~ LF+SIRE ,
                 random = ~ vm(DAM,ainv),
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
  va <- summary(asr2)$varcomp$component[1]

  h2 <- vpredict(asr2,h2~(V1)/(V1+V2))
  print(h2)

 
  # The breeding value is extracted from the model, that is, the additive effect of the maternal line.
  # ind <- unique(c(bluehb$SIRE))
  # sebv <- extracting(asr=asr2,var=va,individual = ind,term='Name*')
  # sebv <- sebv[sebv$EBV!=0,]
  # sebv$par='SIRE'
  
  ind <- unique(c(bluehb$DAM))
  debv <- extracting(asr=asr2,var=va,individual = ind,term='DAM*')
  debv <- debv[debv$EBV!=0,]
  debv$par='DAM'

  outlist <- list(h2=h2,varcomp=summary(asr2)$varcomp,ebv=debv,
                  popmean=popmean)
  return(outlist)
}

extracting <- function(asr, var, individual, term, gxe=FALSE) {
  ebv <- summary(asr, coef=T)$coef.random
  ebv <- as.data.frame(ebv)
  dim(ebv)
  bvalue <- as.data.frame(ebv[grepl(term, row.names(ebv)), ])

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
  bvalue$rel <- ifelse(v<=0, 0, v)
  bvalue$acc <- sqrt(bvalue$rel)
  
  bvalue <- bvalue[bvalue$id %in% individual,]
  
  bvalue <- bvalue[,c(4,1,2,5,6)]
  row.names(bvalue) <- seq(1:nrow(bvalue))
  names(bvalue) <- c("ID", "EBV", "se", "rel", "acc")
  print(head(bvalue))
  
  return(bvalue)
  
}


library(tidyverse)

library(asreml)
library(qgenet)

# ped <- readxl::read_excel('datafiles/Filter_Pedigree_Split_TC1_2023.10.28.xlsx')
# colnames(ped) <- c("ID", "DAM", "SIRE")
# idlist=ped[substr(ped$ID,1,7)=="23TC1XM",]$ID
# pedtc1=pedigreeTrimming(ped,idlist)
# 
# ped3$ID[ped3$DAM=='ZMN00557' & ped3$SIRE=='ZMN00005/ZMN00381//!-017' ]  <- '23TC2XM4158E'
# ped3$ID[ped3$DAM=='ZMN00811' & ped3$SIRE=='ZMN00005/ZMN00381//!-018' ]  <- '23TC2XM4159E'
# ped3$ID[ped3$DAM=='ZMN00451' & ped3$SIRE=='ZMN00806	']  <- '23TC2XM4163'
# ped3$ID[ped3$DAM=='ZMN00461' & ped3$SIRE=='ZMN00806	' ]  <- '23TC2XM4167'
# ped3$ID[ped3$DAM=='ZMN00557' & ped3$SIRE=='ZMN00005/ZMN00381//!-012' ]  <- '23TC2XM4152E'
# ped3$ID[ped3$DAM=='ZMN00811' & ped3$SIRE=='ZMN00005/ZMN00381//!-011' ]  <- '23TC2XM4163E'
# ped3$ID[ped3$DAM=='ZMN00811' & ped3$SIRE=='ZMN00005/ZMN00381//!-015' ]  <- '23TC2XM4167E'
# ped3$ID[ped3$DAM=='ZMN00824' & ped3$SIRE=='ZMN01541	' ]  <- '23TC2XM4163E'
# 
# openxlsx::write.xlsx(pedtc1,file ="datafiles/23TC1ped_2023.10.30.xlsx") 

# matches <- data.frame(read.table(textConnection(
#   'Name DAM SIRE
#     DH605M  0 0
#     DH605F  0 0
#     DH605  DH605F DH605M
#     Chang72 0 0
#     Zheng58 0 0
#     ZD958  Zheng58 Chang72
#     DMY1M 0 0
#     DMY1F 0 0
#     DMY1  DMY1F DMY1M
#     XY1219M 0 0
#     XY1219F 0 0
#     XY1219  XY1219F XY1219M
#     XY1483 PH2GAA PH26JA
#     PH2GAA 0 0
#     PH26JA 0 0
#     XY1225 PHHJC PH1CRW
#     PHHHJC 0 0
#     WG737 WG4582 WG603
#     WG4582 0 0
#     WG603 PH4CV WG5603
#     WC009 0 0
#     V76-1 0 0
#     HY187	V76-1	WC009
#     SCML0849  0 0
#     ZNC442  0 0
#     CD99  ZNC442  SCML0849
#     ZMN00545  0 0
#     ZMN00080  0 0
#     DK159 ZMN00080  ZMN00545
#     ZMN00333  0 0
#     ZMN00154  0 0
#     XY335 ZMN00333  ZMN00154
#     C1563M  0 0
#     C1563F  0 0
#     C1563 C1563F  C1563M
#     '
# ), header = T, sep = ''))
# 
# ped <- readxl::read_excel('datafiles/23TC1ped_2023.10.30.xlsx')
# colnames(ped) <- c("Name", "DAM", "SIRE")
# ped=unique(rbind(matches,ped))
# saveRDS(ped,'datafiles/23TC1ped_2023.10.31.rds')

#ped=readRDS('datafiles/23TC1ped_2023.10.31.rds')
ped <- readxl::read_excel('datafiles/23TC1ped_2023.10.30.xlsx',sheet = 1)
colnames(ped) <- c("Name", "DAM", "SIRE")
ped=unique(ped)
ainv <- ainverse(ped)


data <- readRDS("Output/23TC1_datalist_bluelist1118.rds") 

all.blue <- data$blue

all_trait_name <- c("YLD14","MST",  "PHT", "EHT")#,"TWT"


h2.all.trait <- NULL
d2.all.trait <- NULL
popMean.all.trait <- NULL
popMean.tmp <- data.frame(popmean=0)
varcomp.all.trait<-NULL


for (i in 1:length(all_trait_name)){
  trait=all_trait_name[i]
  blue.value <- all.blue[[which(names(all.blue) == trait)]]
  #colnames(blue.value)[10] <- "Year"
  blue.value=blue.value[!substr(blue.value$Name,6,8)=="Liu",]
  blue.value=blue.value[!substr(blue.value$Name,1,3)=="PRF",]
  Station.vec <- names(table(blue.value$Station))
  
  
  # Mebv.all.trait <- NULL
  ebv.all.trait <- NULL
  
  for (station in Station.vec){
    cau.data <- blue.value %>% dplyr::filter(Station %in% station)
    asrout <- GCA_SCA_Estimation_TC(bluehb = cau.data, ainv = ainv, ped = ped)
    
    popMean.tmp$popmean <- asrout$popmean   
    popMean.tmp$Station <- station
    popMean.tmp$Trait <- trait
    popMean.all.trait <- rbind(popMean.all.trait,popMean.tmp)
    
    varcomp.tmp<-asrout$varcomp
    varcomp.tmp$Station <- station
    varcomp.tmp$Trait <- trait
    varcomp.all.trait <- rbind(varcomp.all.trait, varcomp.tmp)
    
    h2.tmp <-  asrout$h2
    h2.tmp$Station <- station
    h2.tmp$Trait <- trait
    h2.all.trait <- rbind(h2.all.trait, h2.tmp)
    
    
    
    # popMean$popmean[i] <- asrout$popmean
    
    
    # Mebv <- asrout$Mebv[,c(1:3,5)]
    ebv <- asrout$ebv[,c(1:3,5)]
    
    
    # names(Mebv)[names(Mebv)=='EBV'] <- trait
    names(ebv)[names(ebv)=='EBV'] <- trait
    
    # names(Mebv)[names(Mebv)=='se'] <- paste0(trait,'se')
    names(ebv)[names(ebv)=='se'] <- paste0(trait,'se')
    
    
    # names(ebv)[names(ebv)=='acc'] <- paste0(trait,'acc')
    names(ebv)[names(ebv)=='acc'] <- paste0(trait,'acc')
    
    
    # Mebv$Station <- station
    ebv$Station <- station
    
    
    # Mebv.all.trait <- rbind(Mebv.all.trait, Mebv)
    ebv.all.trait <- rbind(ebv.all.trait, ebv)
    
  }
  
  outdata <- list(
    # matchesMebv.tmp = Mebv.all.trait,
    ebv.tmp = ebv.all.trait)
  
  #res.all.2 <- c(res.all.2, list(outdata))
  
  if (i==1){
    
    # Mebv.all<-outdata$Mebv.tmp
    ebv.all<-outdata$ebv.tmp
  } else{
    
    # Mebv.all<-merge(Mebv.all,outdata$Mebv.tmp,by=c("ID","Station"),all=TRUE)
    ebv.all<-merge(ebv.all,outdata$ebv.tmp,by=c("ID","Station"),all=TRUE)
  }
  
}


varcomp.all.trait$item<-row.names(varcomp.all.trait)

idx<-list(h2 = cbind(Type = rep("h2", dim(h2.all.trait)[1]),h2.all.trait),
          
          popMean = cbind(Type = rep("popmean", dim(popMean.all.trait)[1]),popMean.all.trait),
          varcomp = varcomp.all.trait,
          # Mebv=Mebv.all,
          ebv=ebv.all
)


#openxlsx::write.xlsx(idx, file = "Output/1118GCA_TC1_by_station.xlsx")

saveRDS(idx, "Output/1118GCA_23TC1_by_station.rds")

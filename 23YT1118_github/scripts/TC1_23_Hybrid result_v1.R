#####diease mean#####
#~Name#
rm(list = ls())
library(asreml)
library(tidyverse)
library(reshape2)
library(dplyr)


# diease.tmp <- readRDS('datafiles/23phenos_TC1_1103.rds')
# 
#  xx_zq=readRDS("datafiles/seedAllocation_Sites_TC1_XX_20230418_ZQreallocated(1).rds")
#  jn_xx=readxl::read_xlsx("datafiles/seedAllocation_Sites_TC1_JN_20230418_XXreallocated(1).xlsx",sheet = "新乡")
#  name=xx_zq$章丘$Name
#  namejn=jn_xx$Name
# 
#  zq=diease.tmp[diease.tmp$Location %in% c("SDZQ1","SDZQ2"),]
#  tmp=zq[!(zq$Name %in% name),]
#  zq$Station="XX"
# 
#  
# 
#  xx=diease.tmp[diease.tmp$Location %in% c("HNXX1","HNXX3"),]
#  tmp.xx=xx[!(xx$Name %in% namejn),]
#  xx$Station="JN"
# 
#  all=diease.tmp[!(diease.tmp$Location %in% c("HNXX1","HNXX3","SDZQ1","SDZQ2")),]
# 
#  diease.tmp.o=rbind(zq,tmp,xx,tmp.xx,all)
#  diease.tmp.o=unique(diease.tmp.o)
#  saveRDS(diease.tmp.o,'datafiles/23phenos_TC1_1107.rds')
#  

diease.tmp <- readRDS('datafiles/23phenos_TC1_1118reallocated.rds')
table(diease.tmp$Station,diease.tmp$Location)

diease.tmp=diease.tmp[!diease.tmp$Type=="Filler",]
#diease.tmp=diease.tmp[substr(diease.tmp$Type,1,2)=="P1",]



diease.tmp.1=diease.tmp[diease.tmp$Station %in% c("CC","NW"),]
diease.tmp.1$Station="SP"

diease.tmp.2=diease.tmp[diease.tmp$Station %in% c("JN","XX"),]
diease.tmp.2$Station="SU"

diease=rbind(diease.tmp,diease.tmp.1,diease.tmp.2)
diease=unique(diease)


result.all.trait=NULL
diease_filtered <- diease[!(diease$PlotDiscarded %in% c("Yes")), ]

all_trait_name <- c("RTLPCT", "TDPPCT", "GSPPCT","STKLPCT","ERTLPCT","LRTLPCT", "STKRPCT",
                    "PMDPCT",'NCLB','SCLB','GLS',"RSTSOU",
                    "SHBLSC","PMDSC")#"YLD14","MST","PHT","EHT",

for (i in 1:length(all_trait_name)){
  
  trait=all_trait_name[i]
  blue.value <- diease_filtered[,c("Year","Location","Field","Range","Pass","Name","Station",trait)]
  
  names(blue.value)[names(blue.value)==trait]<-"trait"
  
  aoa.vec <- names(table(blue.value$Station))
  #aoa.vec=c("SU","SP")
  result.all=NULL
  
  for (aoa in aoa.vec){
    cau.data <- blue.value %>% dplyr::filter(Station %in% aoa)
    cau.data=cau.data[!is.na(cau.data$trait),]
    # 检查 trait 列是否全为空
    if (all(is.na(cau.data$trait) )) {
      # 如果全为空，则跳过计算
      cat(paste0(trait,aoa,"is empty,skipping\n"))
      
      
    } else {
      result <- aggregate(trait ~ Name, cau.data, FUN = function(x) c(length(x),mean=mean(x), max(x), min(x)))
      result$Station=aoa
      
      new_df <- cbind(result$Name, result$Station,result$trait)
      
      colnames(new_df) <- c("Name","Station",paste0("Obs_",trait),paste0("mean_",trait),paste0("max_",trait),paste0("min_",trait))
      
      
      result.all <- rbind(result.all,new_df)
      
    }
    
  }
  if (!is.null(result.all)) {
    if (is.null(result.all.trait)){
      result.all.trait <-result.all
      
    } else{
      if (!is.null(result.all.trait) && !is.null(result.all)) {
        result.all.trait <- merge(result.all.trait,result.all,by=c("Name","Station"),all=TRUE)
        # 将第3列到最后一列设置为数值型
        result.all.trait <- result.all.trait %>%
          mutate_at(vars(3:ncol(.)), as.numeric)
      }
    }
  }
}



saveRDS(result.all.trait, "Output/1118diease_station_by_23TC1.rds")

#~Name+Location#
rm(list = ls())
library(asreml)
library(tidyverse)
library(reshape2)
library(dplyr)

diease.tmp <- readRDS('datafiles/23phenos_TC1_1118reallocated.rds')
#diease <- readRDS('datafiles/23phenos_TC2_1031.rds')

diease.tmp=diease.tmp[!diease.tmp$Type=="Filler",]
#diease.tmp=diease.tmp[substr(diease.tmp$Type,1,2)=="P1",]



diease.tmp.1=diease.tmp[diease.tmp$Station %in% c("CC","NW"),]
diease.tmp.1$Station="SP"

diease.tmp.2=diease.tmp[diease.tmp$Station %in% c("JN","XX"),]
diease.tmp.2$Station="SU"

diease=rbind(diease.tmp,diease.tmp.1,diease.tmp.2)
diease=unique(diease)

result.all.trait=NULL
# diease$AOA[diease$AOA %in% c("MMSP","LMSP")]="SP"
# diease$AOA[diease$AOA %in% c("NCSU","MCSU","SCSU")]="SU"
diease_filtered <- diease[!(diease$PlotDiscarded %in% c("Yes")), ]

all_trait_name <- c("RTLPCT", "TDPPCT", "GSPPCT","STKLPCT","ERTLPCT","LRTLPCT", "STKRPCT",
                    "PMDPCT",'NCLB','SCLB','GLS',"RSTSOU",
                    "SHBLSC","PMDSC")#"YLD14","MST","PHT","EHT",

for (i in 1:length(all_trait_name)){
  
  trait=all_trait_name[i]
  blue.value <- diease_filtered[,c("Year","Location","Field","Range","Pass","Name","Station",trait)]
  
  names(blue.value)[names(blue.value)==trait]<-"trait"
  
  AOA.vec <- names(table(blue.value$Station))
  
  result.all=NULL
  
  for (aoa in AOA.vec){
    cau.data <- blue.value %>% dplyr::filter(Station %in% aoa)
    cau.data=cau.data[!is.na(cau.data$trait),]
    # 检查 trait 列是否全为空
    if (all(is.na(cau.data$trait) )) {
      # 如果全为空，则跳过计算
      cat(paste0(trait,aoa,"is empty,skipping\n"))
      
      
    } else {
      result <- aggregate(trait ~ Name+Location, cau.data, FUN = function(x) c(mean=mean(x), max(x), min(x),length(x)))
      result$Station=aoa
      new_df <- cbind(result$Name, result$Station,result$Year,result$Location,result$trait)
      
      colnames(new_df) <- c("Name","Station","Location",paste0("mean_",trait),paste0("max_",trait),
                            paste0("min_",trait),paste0("Obs_",trait))
      
      
      result.all <- rbind(result.all,new_df)
      
    }
    
  }
  if (!is.null(result.all)) {
    if (is.null(result.all.trait)){
      result.all.trait <-result.all
      
    } else{
      if (!is.null(result.all.trait) && !is.null(result.all)) {
        result.all.trait <- merge(result.all.trait,result.all,by=c("Name","Station","Location"),all=TRUE)
        # 将第3列到最后一列设置为数值型
        result.all.trait <- result.all.trait %>%
          mutate_at(vars(5:ncol(.)), as.numeric)
      }
    }
  }
}



saveRDS(result.all.trait, "Output/1118diease_loc_by_23TC1.rds")


tmp.1=readRDS("Output/1118diease_station_by_23TC1.rds")
dim(tmp.1)

tmp.1$Location="All"
names(tmp.1)

dim(result.all.trait)

tmp.1=tmp.1[,names(result.all.trait)]



result=rbind(tmp.1,result.all.trait)
result=unique(result)

#openxlsx::write.xlsx(result, file = "Output/1103diease_by_23P1result.xlsx")
saveRDS(result, "Output/1118diease_23TC1result.rds")


#####datclean mean#####
#~ Name
rm(list = ls())
library(asreml)
library(tidyverse)
library(reshape2)
library(dplyr)
data <- readRDS("Output/23TC1_datalist_bluelist1118.rds") 

all.blue <- data$datclean
# ckblue <- readRDS("Output/23P_datalist_bluelist.rds") 
# ckblue<- ckblue$bluelist$blue

result.all.trait=NULL
all_trait_name <- c( "YLD14","MST", "PHT", "EHT")#,"TWT"

for (i in 1:length(all_trait_name)){
  trait=all_trait_name[i]
  blue.value.tmp <- all.blue[[which(names(all.blue) == trait)]]
  #colnames(blue.value)[10] <- "Year"
  blue.value.tmp=blue.value.tmp[!substr(blue.value.tmp$Name,6,8)=="Liu",]
  blue.value.tmp=blue.value.tmp[!substr(blue.value.tmp$Name,1,3)=="PRF",]
  
  
  blue.value.tmp.1=blue.value.tmp[blue.value.tmp$Station %in% c("CC","NW"),]
  blue.value.tmp.1$Station="SP"
  
  blue.value.tmp.2=blue.value.tmp[blue.value.tmp$Station %in% c("JN","XX"),]
  blue.value.tmp.2$Station="SU"
  
  blue.value=rbind(blue.value.tmp,blue.value.tmp.1,blue.value.tmp.2)
  blue.value=unique(blue.value)
  
  # blue.value$AOA[blue.value$AOA %in% c("MMSP","LMSP")]="SP"
  # blue.value$AOA[blue.value$AOA %in% c("NCSU","MCSU","SCSU")]="SU"
  # ck=ckblue[[which(names(ckblue) == trait)]]
  # ck=ck[ck$AOA=="EMSP",]
  # ck=ck[ck$Name %in% c("DMY1","DMY3","C1563"),]
  # blue.value=rbind(blue.value,ck)
  # blue.value$Station[blue.value$Location=="HLJHULCC"]="CC"
  
  names(blue.value)[names(blue.value)==trait]<-"trait"
  
  aoa.vec <- names(table(blue.value$Station))
  
  
  
  # Mebv.all.trait <- NULL
  result.all=NULL
  
  for (aoa in aoa.vec){
    
    cau.data <- blue.value %>% dplyr::filter(Station %in% aoa)
    cau.data$Name=as.character(cau.data$Name)
    cau.data=cau.data[!is.na(cau.data$trait),]
    if (all(is.na(cau.data$trait) )) {
      # 如果全为空，则跳过计算
      cat(paste0(trait,station,"is empty,skipping\n"))
      
      
    } else {
      result <- aggregate(trait ~ Name, cau.data, FUN = function(x) c(length(x),mean=mean(x)))
      result$Station=aoa
      
      new_df <- cbind(result$Name, result$Station,result$trait)
      
      colnames(new_df) <- c("Name","Station",paste0("Obs_",trait),paste0("Pheno_",trait))
      
      
      result.all <- rbind(result.all,new_df)
      
    }
    
  }
  if (!is.null(result.all)) {
    if (is.null(result.all.trait)){
      result.all.trait <-result.all
      
    } else{
      if (!is.null(result.all.trait) && !is.null(result.all)) {
        result.all.trait <- merge(result.all.trait,result.all,by=c("Name","Station"),all=TRUE)
        # 将第3列到最后一列设置为数值型
        result.all.trait <- result.all.trait %>%
          mutate_at(vars(3:ncol(.)), as.numeric)
      }
    }
  }
}
saveRDS(result.all.trait, "Output/1118datclean_station_by_23TC1.rds")


#~ Name+Location#
rm(list = ls())
library(asreml)
library(tidyverse)
library(reshape2)
library(dplyr)
data <- readRDS("Output/23TC1_datalist_bluelist1118.rds") 

all.blue <- data$datclean
table(all.blue$YLD14$Location,all.blue$YLD14$Station)
# ckblue <- readRDS("Output/23P_datalist_bluelist.rds") 
# ckblue<- ckblue$bluelist$blue


result.all.trait=NULL
all_trait_name <- c( "YLD14","MST", "PHT", "EHT")#,"TWT"

for (i in 1:length(all_trait_name)){
  trait=all_trait_name[i]
  blue.value.tmp <- all.blue[[which(names(all.blue) == trait)]]
  #colnames(blue.value)[10] <- "Year"
  blue.value.tmp=blue.value.tmp[!substr(blue.value.tmp$Name,6,8)=="Liu",]
  blue.value.tmp=blue.value.tmp[!substr(blue.value.tmp$Name,1,3)=="PRF",]
  
  blue.value.tmp.1=blue.value.tmp[blue.value.tmp$Station %in% c("CC","NW"),]
  blue.value.tmp.1$Station="SP"
  
  blue.value.tmp.2=blue.value.tmp[blue.value.tmp$Station %in% c("JN","XX"),]
  blue.value.tmp.2$Station="SU"
  
  blue.value=rbind(blue.value.tmp,blue.value.tmp.1,blue.value.tmp.2)
  
  # ck=ckblue[[which(names(ckblue) == trait)]]
  # ck=ck[ck$AOA=="EMSP",]
  # ck=ck[ck$Name %in% c("DMY1","DMY3","C1563"),]
  # blue.value=rbind(blue.value,ck)
  # blue.value$Station[blue.value$Location=="HLJHULCC"]="CC"
  
  names(blue.value)[names(blue.value)==trait]<-"trait"
  
  aoa.vec <- names(table(blue.value$Station))
  
  
  # Mebv.all.trait <- NULL
  result.all=NULL
  
  for (aoa in aoa.vec){
    cau.data <- blue.value %>% dplyr::filter(Station %in% aoa)
    cau.data$Name=as.character(cau.data$Name)
    cau.data=cau.data[!is.na(cau.data$trait),]
    
    if (all(is.na(cau.data$trait) )) {
      # 如果全为空，则跳过计算
      cat(paste0(trait,aoa,"is empty,skipping\n"))
      
      
    } else {
      result <- aggregate(trait ~ Name+Location, cau.data, FUN = function(x) c(length(x),mean=mean(x)))
      result$Station=aoa
      new_df <- cbind(result$Name,result$Station,result$Location,result$trait)
      
      colnames(new_df) <- c("Name","Station","Location",paste0("Obs_",trait),paste0("Pheno_",trait))
      
      
      result.all <- rbind(result.all,new_df)
      
    }
    
  }
  if (!is.null(result.all)) {
    if (is.null(result.all.trait)){
      result.all.trait <-result.all
      
    } else{
      if (!is.null(result.all.trait) && !is.null(result.all)) {
        result.all.trait <- merge(result.all.trait,result.all,by=c("Name","Station","Location"),all=TRUE)
        # 将第3列到最后一列设置为数值型
        result.all.trait <- result.all.trait %>%
          mutate_at(vars(4:ncol(.)), as.numeric)
      }
    }
  }
}

saveRDS(result.all.trait, "Output/1118dataclean_loc_by_23TC1.rds")

p=readRDS("Output/1118datclean_station_by_23TC1.rds")
p$Location="All"


dim(p)
names(p)
#p=p[,c(1:2,7,3:6)]

dim(result.all.trait)
names(result.all.trait)
p=p[,names(result.all.trait)]

p_yld=rbind(result.all.trait,p)

p_yld=unique(p_yld)

p_yld <- p_yld %>% arrange(Name, Station)


saveRDS(p_yld, "Output/1118dataclean_by_23TC1result.rds")

#####blue mean#####
#~ Name#
rm(list = ls())
library(asreml)
library(tidyverse)
library(reshape2)
library(dplyr)
data <- readRDS("Output/23TC1_datalist_bluelist1118.rds") 

all.blue <- data$blue
# ckblue <- readRDS("Output/23P_datalist_bluelist.rds") 
# ckblue<- ckblue$bluelist$blue

result.all.trait=NULL
all_trait_name <- c( "YLD14","MST", "PHT", "EHT")#,"TWT"

for (i in 1:length(all_trait_name)){
  trait=all_trait_name[i]
  blue.value.tmp <- all.blue[[which(names(all.blue) == trait)]]
  #colnames(blue.value)[10] <- "Year"
  blue.value.tmp=blue.value.tmp[!substr(blue.value.tmp$Name,6,8)=="Liu",]
  blue.value.tmp=blue.value.tmp[!substr(blue.value.tmp$Name,1,3)=="PRF",]
  
  blue.value.tmp.1=blue.value.tmp[blue.value.tmp$Station %in% c("CC","NW"),]
  blue.value.tmp.1$Station="SP"
  
  blue.value.tmp.2=blue.value.tmp[blue.value.tmp$Station %in% c("JN","XX"),]
  blue.value.tmp.2$Station="SU"
  
  blue.value=rbind(blue.value.tmp,blue.value.tmp.1,blue.value.tmp.2)
  blue.value=unique(blue.value)
  
  aoa.vec <- names(table(blue.value$Station))
  
  
  # Mebv.all.trait <- NULL
  result.all=NULL
  
  for (aoa in aoa.vec){
    
    cau.data <- blue.value %>% dplyr::filter(Station %in% aoa)
    cau.data$Name=as.character(cau.data$Name)
    if (all(is.na(cau.data$predicted.value) )) {
      # 如果全为空，则跳过计算
      cat(paste0(trait,station,"is empty,skipping\n"))
      
      
    } else {
      result <- aggregate(predicted.value ~ Name, cau.data, FUN = function(x) c(mean=mean(x)))
      result$Station=aoa
      
      new_df <- cbind(result$Name, result$Station,result$predicted.value)
      
      colnames(new_df) <- c("Name","Station",paste0("PredVal_",trait))
      
      
      result.all <- rbind(result.all,new_df)
      
    }
    
  }
  if (!is.null(result.all)) {
    if (is.null(result.all.trait)){
      result.all.trait <-result.all
      
    } else{
      if (!is.null(result.all.trait) && !is.null(result.all)) {
        result.all.trait <- merge(result.all.trait,result.all,by=c("Name","Station"),all=TRUE)
        # 将第3列到最后一列设置为数值型
        result.all.trait <- result.all.trait %>%
          mutate_at(vars(3:ncol(.)), as.numeric)
      }
    }
  }
}


saveRDS(result.all.trait, "Output/1118blue_station_by_23TC1.rds")

#~ Name+Location#
rm(list = ls())
library(asreml)
library(tidyverse)
library(reshape2)
library(dplyr)
data <- readRDS("Output/23TC1_datalist_bluelist1118.rds") 

all.blue <- data$blue
# ckblue <- readRDS("Output/23P_datalist_bluelist.rds") 
# ckblue<- ckblue$bluelist$blue


result.all.trait=NULL
all_trait_name <- c( "YLD14","MST", "PHT", "EHT")#,"TWT"

for (i in 1:length(all_trait_name)){
  trait=all_trait_name[i]
  blue.value.tmp <- all.blue[[which(names(all.blue) == trait)]]
  #colnames(blue.value)[10] <- "Year"
  blue.value.tmp=blue.value.tmp[!substr(blue.value.tmp$Name,6,8)=="Liu",]
  blue.value.tmp=blue.value.tmp[!substr(blue.value.tmp$Name,1,3)=="PRF",]
  
  blue.value.tmp.1=blue.value.tmp[blue.value.tmp$Station %in% c("CC","NW"),]
  blue.value.tmp.1$Station="SP"
  
  blue.value.tmp.2=blue.value.tmp[blue.value.tmp$Station %in% c("JN","XX"),]
  blue.value.tmp.2$Station="SU"
  
  blue.value=rbind(blue.value.tmp,blue.value.tmp.1,blue.value.tmp.2)
  blue.value=unique(blue.value)
  
  aoa.vec <- names(table(blue.value$Station))
  
  
  # Mebv.all.trait <- NULL
  result.all=NULL
  
  for (aoa in aoa.vec){
    cau.data <- blue.value %>% dplyr::filter(Station %in% aoa)
    cau.data$Name=as.character(cau.data$Name)
    
    if (all(is.na(cau.data$predicted.value) )) {
      # 如果全为空，则跳过计算
      cat(paste0(trait,aoa,"is empty,skipping\n"))
      
      
    } else {
      result <- aggregate(predicted.value ~ Name+Location, cau.data, FUN = function(x) c(mean=mean(x)))
      result$Station=aoa
      new_df <- cbind(result$Name,result$Station,result$Location,result$predicted.value)
      
      colnames(new_df) <- c("Name","Station","Location",paste0("PredVal_",trait))
      
      
      result.all <- rbind(result.all,new_df)
      
    }
    
  }
  if (!is.null(result.all)) {
    if (is.null(result.all.trait)){
      result.all.trait <-result.all
      
    } else{
      if (!is.null(result.all.trait) && !is.null(result.all)) {
        result.all.trait <- merge(result.all.trait,result.all,by=c("Name","Station","Location"),all=TRUE)
        # 将第3列到最后一列设置为数值型
        result.all.trait <- result.all.trait %>%
          mutate_at(vars(4:ncol(.)), as.numeric)
      }
    }
  }
}

saveRDS(result.all.trait, "Output/1118blue_loc_by_23TC1.rds")




p=readRDS("Output/1118blue_station_by_23TC1.rds")
p$Location="All"

dim(p)
names(p)

dim(result.all.trait)
names(result.all.trait)

p=p[,names(result.all.trait)]

p1blue=rbind(p,result.all.trait)

p1blue=unique(p1blue)

p1blue <- p1blue %>% arrange(Name, Station)


saveRDS(p1blue, "Output/1118blue_by_23TC1result.rds")


#####Compare CK#####
rm(list = ls())
library(stringr)
library(dplyr)
del_result <- function(region, traits, ctrl,data,idx) {
  tc1=data[data[[idx]] == region,]
  tc.tmp=NULL
  
  locations=unique(tc1$Location)
  for (loc in locations) {
    tc2.tmp=tc1[tc1$Location==loc,]
    
    ctrlf <- tc2.tmp[tc2.tmp$Name == ctrl, c("Name", paste0( "PredVal_", traits))]
    if(nrow(ctrlf)==0){
      tc3.tmp=tc1[tc1$Location=="All",]
      
      ctrlf <- tc3.tmp[tc3.tmp$Name == ctrl, c("Name", paste0( "PredVal_", traits))]
    }
    
    for (i in 1:length(traits)) {
      
      tc2.tmp[[paste0( "CKpct_", traits[i])]] <- 
        round((tc2.tmp[[paste0( "PredVal_", traits[i])]] )/(ctrlf[[paste0( "PredVal_", traits[i])]]) * 100, 1)
      
    }
    tc2.tmp <- tc2.tmp %>% arrange(Name,-CKpct_YLD14)
    tc.tmp=rbind(tc.tmp,tc2.tmp)   
  }
  result = tc.tmp
  return (result)
}





tcblue=readRDS("Output/1118blue_by_23TC1result.rds")

tcpheno=readRDS("Output/1118dataclean_by_23TC1result.rds")

diease= readRDS("Output/1118diease_23TC1result.rds")

tc.1=full_join(tcblue,tcpheno)
tc=full_join(tc.1,diease)

tc=unique(tc)

traits=c( "YLD14","MST", "PHT", "EHT")

idx="Station"
region_all <-  unique(tc$Station)#
outebv_SP.tmp=NULL
outebv_SU.tmp=NULL


for (region in region_all) {
  if(region=="HRB"){
    
    
    ctrl = 'DMY3'  # female
    outebv_EMSP = del_result(region, traits, ctrl ,
                             data=tc, idx="Station")
    
  }
  
  if (region %in% c('CC',"NW","SP")){
    
    
    ctrl = 'XY1483'  # female
    outebv_SP=del_result(region, traits,ctrl, 
                         data=tc, idx="Station")
    
    outebv_SP.tmp=rbind(outebv_SP.tmp,outebv_SP)
    
    
  }
  if (region %in% c("JN" ,"XX"  , "SU")){
    
    ctrl= 'DH605'  # female
    outebv_SU=del_result(region, traits, ctrl,
                         data=tc, idx="Station")
    
    outebv_SU.tmp=rbind(outebv_SU.tmp,outebv_SU)
    
  }
  if (region == 'SW'){
    
    ctrl= 'CD99'  # female
    outebv_SW=del_result(region, traits, ctrl, 
                         data=tc, idx="Station")
    
  }
  
}  
result=rbind(outebv_EMSP,outebv_SP.tmp,outebv_SU.tmp,outebv_SW)
result <-  result[, colSums(is.na(result)) != nrow(result)]
result <- result %>% mutate(across(where(is.numeric), ~round(.x, 1)))
names(result)
#result=result[,c(1:4,80)]
result <- result %>% arrange(Name,Station,Location,-CKpct_YLD14)

traitdiease=c("RTLPCT", "TDPPCT", "GSPPCT","STKLPCT","ERTLPCT","LRTLPCT", "STKRPCT",
              'NCLB','SCLB','GLS',"RSTSOU",
              "SHBLSC","PMDSC")#"PMDPCT",
cn <- c("mean_", "max_", "min_","Obs_")
colname <- c()
for (i in 1:length(traitdiease)) {
  colname <- c(colname, paste0(cn, traitdiease[i]))
  
}

traits <- c('YLD14', 'MST', 'PHT', 'EHT')

# for (trait in traits) {
#   result=result[,-grep(paste0("max_",trait),names(result),value = FALSE)]
#   result=result[,-grep(paste0("min_",trait),names(result),value = FALSE)]
#   names(result)[names(result)==paste0("mean_",trait)]=paste0("Pheno_",trait)
# }


traits_1 = unlist(lapply(traits, function(x){c(str_c('Pheno_',x),str_c('PredVal_', x), 
                                               str_c('CKpct_',x),str_c('Obs_',x)
)}))

result=result[,c("Name","Station","Location" ,traits_1,colname)]

pheno <- readRDS('datafiles/23phenos_TC1_1118reallocated.rds')
pheno=pheno[,c("Name","GenoID","Pedigree")]
pheno=pheno[!is.na(pheno$Name),]
pheno=pheno[!is.na(pheno$GenoID),]
pheno=pheno[!is.na(pheno$Pedigree),]
pheno=pheno[!pheno$Pedigree=="FL",]
pheno=pheno[!pheno$Pedigree=="INF",]
pheno=unique(pheno)
#pheno <- pheno[!duplicated(pheno$Name),]
result.tmp=merge(result,pheno,by="Name",all.x = TRUE)#添加系谱
dim(result.tmp)
names(result.tmp)

result.tmp=result.tmp[,c(1,73,2:71)]
#openxlsx::write.xlsx(result,file ="Output/1107result_by_23TC1.xlsx") 
saveRDS(result,"Output/1118result_by_23TC1.rds")



#####Sort#####
library(dplyr)

rm(list = ls())
#setwd("D:/YangZhenyu/Project/00/Breeding_Anova")
# P3P4.Anova.res <- openxlsx::read.xlsx(xlsxFile = "./datafiles/1101result_by_23TC1.xlsx")
#P3P4.Anova.res <- readRDS("Output/1104result_by_23P1.rds")
# pheno.data     <- openxlsx::read.xlsx(xlsxFile = "./datafiles/Yield Trial Master Catalog-2023.10.31xlsx.xlsx")
P3P4.Anova.res=readRDS("Output/1118result_by_23TC1.rds")
P3P4.Anova.res <- unique(P3P4.Anova.res)

AOA.vec <- sort(unique(P3P4.Anova.res$Station))

length(unique(paste0(P3P4.Anova.res$Location, "_", P3P4.Anova.res$Name, "_", P3P4.Anova.res$AOA)))

#AOA.vec <- c('SP','MMSP','LMSP','SU','NCSU','MCSU','SCSU','EMSP','SWCN')
AOA.vec <- c('JN','XX','SU','CC','NW','SP','HRB','SW')

res.all   <- NULL
ped.error <- NULL

P1.Anova.res <- P3P4.Anova.res
res2 <- NULL
for (aoa in AOA.vec) {
  P3.Anova.Aoa.res <- P1.Anova.res %>% dplyr::filter(Station == aoa)
  if(dim(P3.Anova.Aoa.res)[1] != 0){
    
    P3.Anova.Aoa.all.res    <- P3.Anova.Aoa.res %>% dplyr::filter(Location == "All")
    P3.Anova.Aoa.no.all.res <- P3.Anova.Aoa.res %>% dplyr::filter(Location != "All")
    
    order.res <- P3.Anova.Aoa.all.res[order(P3.Anova.Aoa.all.res$PredVal_YLD14, decreasing = T), ]
    res1 <- NULL
    for (j in 1:dim(order.res)[1]) {
      data.2 <- P3.Anova.Aoa.no.all.res %>% dplyr::filter(Name == order.res[j,1])
      data.3 <- data.2[order(data.2$PredVal_YLD14, decreasing = T), ]
      data.4 <- rbind(order.res[j,], data.3)
      
      res1.2 <- data.4
      
      res1.3 <- res1.2[1, ]
      res1.2 <- res1.2[-1,]
      # res1.2 <- res1.2[order(res1.2$BlueAvg_YLD14, decreasing = T), ]
      res1.2 <- res1.2[order(res1.2$Location), ]
      res1.3 <- rbind(res1.2, res1.3)
      res1   <- rbind(res1, res1.3)
    }
    res2 <- rbind(res2, res1)
  }
}
res.all <- rbind(res.all, res2)

openxlsx::write.xlsx(res.all, "Output/1118result_by_23TC1-Sort.xlsx")


# write.xlsx(ped.error, "./output/Pedigree重名的杂交种.xlsx")


rm(list = ls())
library(asreml)
library(tidyverse)
library(reshape2)
library(dplyr)

diease <- readRDS('datafiles/19_23_P_phenos_1031.rds')
#diease <- readRDS('datafiles/23phenos_1031.rds')

#diease <- readRDS('datafiles/23phenos_TC1_1031.rds')

#diease <- readRDS('datafiles/23phenos_TC2_1031.rds')

diease=diease[!diease$Type=="Filler",]
# diease=diease[substr(diease$Type,1,2)=="P1",]
diease=diease[substr(diease$Type,1,1)=="P",]
#diease=diease[diease$Year %in%c("2022","2023"),]
result.all.trait=NULL
diease$AOA[diease$AOA %in% c("MMSP","LMSP")]="SP"
diease$AOA[diease$AOA %in% c("NCSU","MCSU","SCSU")]="SU"
diease_filtered <- diease[!(diease$PlotDiscarded %in% c("Yes")), ]

all_trait_name <- c("RTLPCT", "TDPPCT", "GSPPCT","STKLPCT","ERTLPCT","LRTLPCT", "STKRPCT",
                    "PMDPCT",'NCLB','SCLB','GLS',"RSTSOU",
                    "SHBLSC","PMDSC","EARTSC")

for (i in 1:length(all_trait_name)){
  
  trait=all_trait_name[i]
  blue.value <- diease_filtered[,c("Year","Location","Field","Range","Pass","Name","AOA",trait)]
  
  names(blue.value)[names(blue.value)==trait]<-"trait"
  
  AOA.vec <- names(table(blue.value$AOA))
  
  result.all=NULL
  
  for (aoa in AOA.vec){
    cau.data <- blue.value %>% dplyr::filter(AOA %in% aoa)
    # 检查 trait 列是否全为空
    if (all(is.na(cau.data$trait) )) {
      # 如果全为空，则跳过计算
      cat(paste0(trait,aoa,"is empty,skipping\n"))
      
      
    } else {
      result <- aggregate(trait ~ Name+Year+Location, cau.data, FUN = function(x) c(mean=mean(x), max(x), min(x),length(x)))
      result$AOA=aoa
      new_df <- cbind(result$Name, result$AOA,result$Year,result$Location,result$trait)
      
      colnames(new_df) <- c("Name","AOA","Year","Location",paste0("mean_",trait),paste0("max_",trait),
                            paste0("min_",trait),paste0("Obs_",trait))
      
      
      result.all <- rbind(result.all,new_df)
      
    }
    
  }
  if (!is.null(result.all)) {
    if (is.null(result.all.trait)){
      result.all.trait <-result.all
      
    } else{
      if (!is.null(result.all.trait) && !is.null(result.all)) {
        result.all.trait <- merge(result.all.trait,result.all,by=c("Name","AOA","Year","Location"),all=TRUE)
        # 将第3列到最后一列设置为数值型
        result.all.trait <- result.all.trait %>%
          mutate_at(vars(5:ncol(.)), as.numeric)
      }
    }
  }
}



# openxlsx::write.xlsx(result.all.trait, file = "Output/1031diease_by_19_23.xlsx")
# saveRDS(result.all.trait, "Output/1031diease_by_19_23.rds")

#openxlsx::write.xlsx(result.all.trait, file = "Output/1102P_diease_by_22_23.xlsx")
saveRDS(result.all.trait, "Output/1102P_year_loc_diease_by_19_23.rds")


yearloc=readRDS("Output/1102P_year_loc_diease_by_19_23.rds")
tmp.2=rbind(result.all.trait,yearloc)
dim(tmp.2)
# openxlsx::write.xlsx(result.all.trait, file = "Output/1031diease_by_23TC1.xlsx")
# saveRDS(result.all.trait, "Output/1031diease_by_23TC1.rds")


# openxlsx::write.xlsx(result.all.trait, file = "Output/1031diease_by_23TC2.xlsx")     
# saveRDS(result.all.trait, "Output/1031diease_by_23TC2.rds") 

tmp.1=readRDS("Output/1102P_diease_by_19_23_.tmp.rds")
dim(tmp.1)
tmp.1$Year="All"
tmp.1$Location="All"
names(tmp.1)
names(tmp.2)
tmp.1=tmp.1[,names(tmp.2)]



result=rbind(tmp.1,tmp.2)
openxlsx::write.xlsx(result, file = "Output/1102P34_diease_by_19_23result.xlsx")
saveRDS(result, "Output/1102P34_diease_by_19_23result.rds")

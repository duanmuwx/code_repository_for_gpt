rm(list = ls())
library(asreml)
library(tidyverse)
library(reshape2)
library(dplyr)

# diease <- readRDS('datafiles/19_23_P_phenos_1031.rds')
diease <- readRDS('datafiles/23phenos_1118.rds')

#diease <- readRDS('datafiles/23phenos_TC1_1031.rds')

#diease <- readRDS('datafiles/23phenos_TC2_1031.rds')

diease=diease[!diease$Type=="Filler",]
#diease=diease[substr(diease$Type,1,1)=="P",]
result.all.trait=NULL
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
      result <- aggregate(trait ~ Name, cau.data, FUN = function(x) c(mean=mean(x), max(x), min(x)))
      result$AOA=aoa
      new_df <- cbind(result$Name,result$trait, result$AOA)
      
      colnames(new_df) <- c("Name",paste0("mean_",trait),paste0("max_",trait),paste0("min_",trait),"AOA")
      
     
      result.all <- rbind(result.all,new_df)
      
      }

  }
  if (!is.null(result.all)) {
    if (is.null(result.all.trait)){
      result.all.trait <-result.all
      
    } else{
      if (!is.null(result.all.trait) && !is.null(result.all)) {
        result.all.trait <- merge(result.all.trait,result.all,by=c("Name","AOA"),all=TRUE)
        # 将第3列到最后一列设置为数值型
        result.all.trait <- result.all.trait %>%
          mutate_at(vars(3:ncol(.)), as.numeric)
      }
    }
  }
}



# openxlsx::write.xlsx(result.all.trait, file = "Output/1031diease_by_19_23.xlsx")
# saveRDS(result.all.trait, "Output/1031diease_by_19_23.rds")

#openxlsx::write.xlsx(result.all.trait, file = "Output/1118diease_by_23.xlsx")
saveRDS(result.all.trait, "Output/1118diease_by_23.rds")

# openxlsx::write.xlsx(result.all.trait, file = "Output/1031diease_by_23TC1.xlsx")
# saveRDS(result.all.trait, "Output/1031diease_by_23TC1.rds")


# openxlsx::write.xlsx(result.all.trait, file = "Output/1031diease_by_23TC2.xlsx")     
# saveRDS(result.all.trait, "Output/1031diease_by_23TC2.rds") 


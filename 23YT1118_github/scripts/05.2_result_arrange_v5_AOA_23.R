rm(list = ls())

#' Parent pairing frequency
#'
#' @param parent
#'
#' @return
#' @export
#'
#' @examples
#'
parent_Obs <- function(pheno,parent,traits){
  
  
  Geno = unique(pheno[[parent]])
  Geno_count = NULL
  hybrid_Ob = NULL
  
  
  result.tmp = NULL
  
  for (j in 1:length(Geno)) {
    result = NULL
    name=Geno[j]
    ZMN=pheno[pheno[[parent]]==name,]
    #ZMN=ZMN[!is.na(ZMN$Name),]
    
    hybrid=length(unique(ZMN$Name))
    
    temp = data.frame(GenoID = name, HybridCount = hybrid)
    
    Geno_count = rbind(Geno_count, temp)
    
    for (i in 1:length(traits)) {
      
      ZMN.tmp=ZMN
      trait = traits[i]
      
      names(ZMN.tmp)[names(ZMN.tmp) == trait] <- 'trait'
      ZMN.tmp=ZMN.tmp[!is.na(ZMN.tmp$trait),]
      
      if (all(is.na(ZMN.tmp$trait))) {
        # 如果全为空，则跳过计算
        cat(paste0(name," ",trait, "is empty, skipping\n"))
        result$GenoID=name
        result[paste0("Obs_",trait)] <- NA
        result[paste0("Loc_",trait)] <- NA
        result[paste0("mean_",trait)] <- NA
        result[paste0("max_",trait)] <- NA
        result[paste0("min_",trait)] <- NA
      } else {
        
        hybrid_Ob = data.frame(table(ZMN.tmp[[parent]]))
        hybrid_Ob$N_loc=length(unique(ZMN.tmp$Location))
        # if (trait %in% c("STKRPCT","LDGPCT",'GSPPCT','RTLPCT',"TDPPCT",'NCLB','GLS',"PMDPCT")){
        hybrid_Ob$mean=mean(ZMN.tmp$trait)
        hybrid_Ob$max=max(ZMN.tmp$trait)
        hybrid_Ob$min=min(ZMN.tmp$trait)
      
        colnames(hybrid_Ob) = c("GenoID",paste0("Obs_",trait),paste0("Loc_",trait),
                                paste0("mean_",trait),paste0("max_",trait),paste0("min_",trait))
       # } 
        
        
        
        if (is.null(result)) {
          result = hybrid_Ob
        } else {
          result = merge(result, hybrid_Ob, by = 'GenoID', all = TRUE)
        }
      }
    }
    result.tmp=rbind(result.tmp,result)
  }
  
  result.tmp.tmp=merge(Geno_count,result.tmp, by = 'GenoID', all = TRUE)
  return(result.tmp.tmp)
}



#' Statistical observation of hybrid phenotype
#'
#' @param pheno a data.frame
#' @param traits Multiple traits can be selected
#'
#' @return
#' @export
#'
#' @examples
hybrid_Obs <- function(pheno, traits) {
  
  
  result <- NULL
  
  for (j in 1:length(traits)) {
    sub_list = pheno
    trait = traits[j]
    
    names(sub_list)[names(sub_list) == trait] <- 'trait'
    sub_list = sub_list[!is.na(sub_list$trait),]
    
    if (all(is.na(sub_list$trait))) {
      # 如果全为空，则跳过计算
      cat(paste0(trait, " is empty, skipping\n"))
      
      result[[paste0("Loc_",trait)]] <- NA
      result[[paste0("Obs_",trait)]] <- NA
    } else {
      
      hybrid_Ob = data.frame(table(sub_list$Name,sub_list$Location))
      
      
      # 统计同一个 Name 测定的表型数和 Location 数
      library(dplyr)
      
      result.tmp <- hybrid_Ob %>%
        group_by(Var1) %>%
        summarise(
          N_Loc = n_distinct(Var2),
          N_Ob = sum(Freq)
        )
      
      names(result.tmp) = c("Hybrid",paste0("Loc_",trait),paste0("Obs_",trait))
     
      
      if (is.null(result)) {
        result = result.tmp
      } else {
        result = merge(result, result.tmp, by = 'Hybrid', all = TRUE)
      }
    }
    
  }
  
  return(result)  
}


#' Collation of genetic evaluation results
#'
#' @param popmean Population average
#' @param region  Maize varieties
#' @param traits  phenotypic traits
#' @param ctrlmsample  Reference male parent
#' @param ctrlfsample Reference female parent
#' @param adv
#' @param hybobs
#' @param blupebv breeding value
#' @param Aov Variance test table
#'
#' @return
#'
#' @examples

del_result <- function(popmean, region, traits, ctrlmsample, ctrlfsample, blupebv,idx) {

  traits_1 = unlist(lapply(traits, function(x){c(x,  str_c('CKpct_', x), 
                                                 str_c('Loc_',x),str_c('Obs_',x)
                                                 
                                                 )}))
  
  # 处理常规性状  gca
  Fobs=parent_Fall[parent_Fall[[idx]] == region,]
  Mobs=parent_Mall[parent_Mall[[idx]] == region,]
  
  MGCA <- blupebv$Mebv[blupebv$Mebv[[idx]]==region,]
  FGCA <- blupebv$Febv[blupebv$Febv[[idx]]==region,]
  
  names(MGCA)[1] <- "GenoID"
  names(FGCA)[1] <- "GenoID"
  
  
  MGCA <- MGCA[, -grep("se", names(MGCA), value = FALSE)]
  FGCA <- FGCA[, -grep("se", names(FGCA), value = FALSE)]
  MGCA <- MGCA[, -grep("acc", names(MGCA), value = FALSE)]
  FGCA <- FGCA[, -grep("acc", names(FGCA), value = FALSE)]
  
  MGCA <- merge(MGCA, Mobs, by = c("GenoID",idx), all.x = T)
  
  names(MGCA)[names(MGCA)=='HybridCount'] <- "nHybrids"
  
  FGCA <- merge(FGCA, Fobs, by = c("GenoID",idx), all.x = T)
  
  names(FGCA)[names(FGCA)=='HybridCount'] <- "nHybrids"
  
  
  # 计算基准值
  
  ctrlm <- MGCA[MGCA$GenoID == ctrlmsample, c("GenoID", traits)]
  ctrlf <- FGCA[FGCA$GenoID == ctrlfsample, c("GenoID", traits)]
  

  
  for (i in 1:length(traits)) {
    mu <- popmean$mean[popmean[[idx]]==region & popmean$Trait == traits[i]]
    mu <- as.numeric(mu)
    if(nrow(ctrlm)==1){
      MGCA[[paste0( "CKpct_", traits[i])]] <- 
      round((MGCA[[traits[i]]] + mu)/(ctrlm[[traits[i]]] + mu) * 100, 1)
    }else{
      MGCA[[paste0( "CKpct_", traits[i])]] <- NA
    }
    
    if(nrow(ctrlf)==1){
      FGCA[[paste0( "CKpct_", traits[i])]] <- 
        round((FGCA[[traits[i]]] + mu)/(ctrlf[[traits[i]]] + mu) * 100, 1)
    }else{
      FGCA[[paste0( "CKpct_", traits[i])]] <- NA
    }
    
  
  }
  
  # 处理常规性状  sca
  sca <- blupebv$SCAebv[blupebv$SCAebv[[idx]]==region,]
  sca <- sca[, -grep("se", names(sca), value = FALSE)]
  sca <- sca[, -grep("acc", names(sca), value = FALSE)]
  
 
  names(sca)[names(sca) == "SIRE"] <- "MGenoID" # 父本
  names(sca)[names(sca) == "DAM"] <- "FGenoID"  # 母本
  names(sca)[names(sca) == "Fam"] <- c("GenoID")

  sca <- sca[sca$GenoID != "0/0", ]   
  
  # 计算基准值

  ctrl <- sca$GenoID[sca$GenoID == paste0(ctrlfsample, '/', ctrlmsample)]
  
  
  for (i in 1:length(traits)) {
    names(sca)[names(sca) == traits[i]] <- paste0("SCA_", traits[i])
  }
  
  m_gca <- MGCA[, c("GenoID", traits)]
  
  
  names(m_gca)[2:(length(traits) + 1)] <- 
    paste0("MGenoID_GCA_", names(m_gca)[2: (length(traits) + 1)])
  
  
  f_gca <- FGCA[, c("GenoID", traits)]
  
  names(f_gca)[2: (length(traits) + 1)] <- 
    paste0("FGenoID_GCA_", names(f_gca)[2: (length(traits) + 1)])
  
  sca <- merge(sca, f_gca, by.x = "FGenoID", by.y = "GenoID", all.x = TRUE)
  sca <- merge(sca, m_gca, by.x = "MGenoID", by.y = "GenoID", all.x = TRUE)

  
  for(i in 1:length(traits)){
    mu <- popmean$mean[popmean[[idx]]==region & popmean$Trait == traits[i]]
    mu <- as.numeric(mu)
    
    MGCA[[traits[i]]] =  MGCA[[traits[i]]]+ mu
    FGCA[[traits[i]]] =  FGCA[[traits[i]]]+ mu
    
    sca[paste0("Perf_", traits[i])] <- 
      sca[paste0("MGenoID_GCA_",traits[i])] + 
      sca[paste0("FGenoID_GCA_", traits[i])] + 
      sca[paste0("SCA_", traits[i])]+mu
    
    sca[paste0("CKpct_", traits[i])] <- 
      round(sca[paste0("Perf_", traits[i])]/
              (sca[sca$GenoID==ctrl,paste0("Perf_", traits[i])]) * 100, 1)
 
  }

  MGCA <- MGCA[order(-MGCA$CKpct_YLD), ]
  MGCA$Rank_YLD=1:nrow(MGCA)
  FGCA <- FGCA[order(-FGCA$CKpct_YLD), ]
  FGCA$Rank_YLD=1:nrow(FGCA)

  # 排序
  colname <- c("GenoID",idx, "nHybrids","Rank_YLD", traits_1)

  MGCA <- MGCA[, colname]
  FGCA <- FGCA[, colname]
  #MGCA <-  MGCA[, colSums(is.na(MGCA)) != nrow(MGCA)]
  #FGCA <-  FGCA[, colSums(is.na(FGCA)) != nrow(FGCA)]
  #保留三位小数
  MGCA <- MGCA %>% mutate(across(where(is.numeric), ~round(.x, 1)))
  FGCA <- FGCA %>% mutate(across(where(is.numeric), ~round(.x, 1)))
  
  # 合并疾病性状和常规性状
  

  sca <- sca %>% mutate(across(where(is.numeric), ~round(.x, 1)))
  
  sca.tmp=sca
  
  hybobs=hyb_obs[hyb_obs[[idx]]==region, ] 
  
  sca.tmp = merge(sca.tmp,hybobs, by.x = c("Name",idx),by.y = c("Hybrid",idx), all.x = T)
  sca.tmp <- sca.tmp[order(-sca.tmp$CKpct_YLD), ]
  sca.tmp$Rank_YLD=1:nrow(sca.tmp)
  
  
  colname <- c()
  
  cn <- c("FGenoID_GCA_", "MGenoID_GCA_", "SCA_", "Perf_", "CKpct_","Loc_","Obs_")
  
  for (i in 1:length(traits)) {
    colname <- c(colname, paste0(cn, traits[i]))
    
  }
  
  
  
  sca.tmp <- sca.tmp[, c("Name",idx,"GenoID","FGenoID","MGenoID", "Rank_YLD",colname)]
  

  outebv <- list(FP_GCA = FGCA, 
                 MP_GCA = MGCA, 
                 SCA_GCA = sca.tmp
                 )

  return (outebv)

}




#' popmean
#'
#' @param pheno_datclean a list,Phenotypic data after data quality control
#' @param traits 
#'
#' @return
#' @export
#'
#' @examples
popmean<-function(pheno_datclean,traits,idx,idx.vec){
  mean.all=c()
  for (trait in traits) {
    pheno_dat <- pheno_datclean[[which(names(pheno_datclean) == trait)]]
    #idx.vec <- names(table(pheno_dat[[idx]]))
    for (tmp in idx.vec){
      cau.data <- pheno_dat[pheno_dat[[idx]]==tmp,]
      names(cau.data)[names(cau.data) == trait]="trait"
      cau.data=cau.data[!is.na(cau.data$trait),]
      mean=data.frame(mean(cau.data$trait))
      names(mean)="mean"
      mean$Trait=trait
      mean[[idx]]=tmp
      mean.all=rbind(mean.all,mean)
    }
  }
  return(mean.all)
}





#' modify excel files
#'
#' @param data a list
#' @param sheets name of list
#'
#' @return
#' @export
#'
#' @examples
modify=function(data,sheets){
  library(openxlsx)  
  wb <- createWorkbook()
  for(sheet in sheets){
    addWorksheet(wb, sheet)
    writeData(wb, 
              sheet, 
              data[[sheet]], 
              startCol = 1, 
              startRow = 1)
    #rowNames = TRUE)
    modifyBaseFont(wb, fontSize = 12, fontName = "Arial")
    setColWidths(wb, sheet = sheet, cols = 1:(ncol(data[[sheet]])), widths = "4" 
    )
    
    
    # style <- createStyle(fgFill = "#008B8B",
    #                      fontColour = "#FFFFFF",wrapText = TRUE)
    # addStyle(wb,
    #          sheet = sheet,
    #          style=style,
    #          rows = 1,
    #          cols = 1:(ncol(data[[sheet]])),
    #          gridExpand = TRUE)
    style2 <-  createStyle(
      borderColour = "#4F81BD",
      border='TopBottomLeftRight')
    addStyle(wb,
             sheet = sheet,
             style=style2,
             rows = 2:(nrow(data[[sheet]])+1),
             cols = 1:(ncol(data[[sheet]])),
             gridExpand = TRUE)
    # style3 <-  createStyle(
    #   borderColour = "#4F81BD",
    #   border='bottom',
    #   borderStyle='double')
    # addStyle(wb,
    #          sheet = sheet,
    #          style=style3,
    #          rows = nrow(data[[sheet]])+1,
    #          cols = 1:(ncol(data[[sheet]])),
    #          gridExpand = TRUE,
    #          stack=TRUE)
    # style4 <- createStyle(fgFill = "#F0FFFF")
    # addStyle(wb,
    #         sheet = sheet,
    #         style=style4,
    #         rows = 2:(nrow(data[[sheet]])+1),
    #         cols = 1:(ncol(data[[sheet]])),
    #     gridExpand = TRUE,
    #     stack=TRUE)
    # showGridLines(wb,
    #              sheet,
    #             showGridLines = FALSE)
    
    style <- createStyle(
      fgFill = "#008B8B",
      fontColour = "#FFFFFF",
      textDecoration='bold',
      halign='center')
    addStyle(wb,
             sheet = sheet,
             style=style,
             rows = 1,
             cols = 1:(ncol(data[[sheet]])),
             gridExpand = TRUE)
    # style2 <-  createStyle(
    #   borderColour = "#4F81BD",
    #   border='TopBottomLeftRight',
    #   halign='center')
    # addStyle(wb,
    #          sheet = sheet,
    #          style=style2,
    #          rows = 2:(nrow(data[[sheet]])+1),
    #          cols = 1:(ncol(data[[sheet]])),
    #          gridExpand = TRUE)
    freezePane(
      wb,
      sheet,
      #firstActiveRow = NULL,
      #firstActiveCol = NULL,
      firstRow = TRUE,
      firstCol = TRUE)
    
  }
  
  saveWorkbook(wb, 
               file=file_path,
               overwrite = T)
  # saveRDS(wb, file = file_path1)
}


#-------------------------------------------------------------------------------#

library(tidyverse)
library(reshape2)
library(stringi)
library(writexl)
library(openxlsx)


# raw_data (a data frame) for Hybrid_Observations and parent_Observations
# traits <- c('YLD14', 'MST', 'PHT', 'EHT',"TWT","LDGPCT",'GSPPCT','RTLPCT',"TDPPCT",'NCLB','SCLB','GLS',"RSTSOU","SHBLSC","STKRPCT","PMDPCT","PMDSC","EARTSC")
traits <- c('YLD14', 'MST', 'PHT', 'EHT',"RTLPCT", 
            "TDPPCT", "GSPPCT","STKLPCT","ERTLPCT","LRTLPCT", "STKRPCT",
            "PMDPCT",'NCLB','SCLB','GLS',"RSTSOU",
            "SHBLSC","PMDSC")#,"TWT"

#pheno <- readRDS("datafiles/19_23_P_phenos_1028.rds")   ## raw_data (a data frame) for Hybrid_Observations and parent_Observations

pheno <- readRDS("datafiles/23phenos_1118.rds")

pheno=pheno[substr(pheno$Type,1,1)=="P",]
pheno <- pheno[!(pheno$PlotDiscarded %in% c("Yes")), ]
# pheno$AOA[pheno$Station =="CC"] <- "SP"
# pheno$AOA[pheno$Station =="NW"] <- "SP"
# 
# pheno$AOA[pheno$Station == "JN"] <- "SU"
# pheno$AOA[pheno$Station == "XX"] <- "SU"



pedigree_P_inbredlines <- readRDS("datafiles/pedigree_19_23P_1031.rds")
colnames(pedigree_P_inbredlines) <- c("Name", "DAM", "SIRE")
pheno1<-merge(pheno,pedigree_P_inbredlines,by="Name",all.x = T)
pheno1=pheno1[!is.na(pheno1$Name),]
pheno1=pheno1[!is.na(pheno1$DAM),]
pheno1=pheno1[!is.na(pheno1$SIRE),]
pheno1=pheno1[!pheno1$SIRE==0,]
pheno1=pheno1[!pheno1$DAM==0,]
# pheno1$SIRE[is.na(pheno1$SIRE)] <- '0' # Missing values are represented by 0
# pheno1$DAM[is.na(pheno1$DAM)] <- '0'   # Missing values are represented by 0
# 

parent_Fall=c()
parent_Mall=c()
hyb_obs=c()
popmeand=data.frame()
idx="AOA" #? can modify eg. Station
region_all <-  c('EMSP',"MMSP","LMSP","SCSU","MCSU","SWCN","NCSU")#
#region_all <-  c('HEB', 'JN', 'CC', 'XX', 'SW')
#AOA.vec <- names(table(pheno1[[idx]]))

for (aoa in region_all) {
  cau.data = pheno1[pheno1[[idx]]==aoa,]
  
  parent_F=parent_Obs(pheno=cau.data,parent = "DAM",traits)
  parent_F[[idx]]=aoa
  parent_M=parent_Obs(pheno=cau.data,parent = "SIRE",traits)
  parent_M[[idx]]=aoa
  hybs=hybrid_Obs(pheno=cau.data,traits)
  hybs[[idx]]=aoa
  parent_Fall=rbind(parent_Fall,parent_F)
  parent_Mall=rbind(parent_Mall,parent_M)
  hyb_obs=rbind(hyb_obs,hybs)
  
}
#out=list(parent_F=parent_Fall,parent_M=parent_Mall,hybrid=hyb_obs)
#openxlsx::write.xlsx(out, file = "Output/parent_hybrid_obs23.xlsx")
#saveRDS(out, file = "Output/parent_hybrid_obs23.rds")

#calculate population mean
pheno_datclean=readRDS("Output/23P_datalist_bluelist.rds") #for calculate population mean
pheno_datclean=pheno_datclean$datalist$datclean

traits <- c('YLD14', 'MST', 'PHT', 'EHT')#,"TWT"

popmean=popmean(pheno_datclean,traits,idx.vec=region_all,idx="AOA")#? idx can modify eg. Station



#breeding value result 
blupebv = readRDS("Output/1118GCA_SCA_by_23.rds") #MST,YLD,PHT,EHT,TWT
diease <- readRDS("Output/1118diease_by_23.rds")#diease traits

# blupebv = readRDS("Output/1028GCA_SCA_by_23.rds") #MST,YLD,PHT,EHT,TWT
# diease <- readRDS("Output/1028diease_by_23.rds")#diease traits

#diease <- diease %>% mutate(across(where(is.numeric), ~round(.x, 3)))
#traits <- c('YLD14', 'MST', 'PHT', 'EHT',"TWT","LDGPCT",'GSPPCT','RTLPCT',"TDPPCT",'NCLB','SCLB','GLS',"RSTSOU","SHBLSC","STKRPCT","PMDPCT")

region_all <-  c('EMSP',"MMSP","LMSP","SWCN","SCSU","MCSU","NCSU")#
outebv_SP.tmp=NULL
outebv_SU.tmp=NULL

for (region in region_all) {
  if(region=="EMSP"){
  
    ctrlmsample = 'C1563M'  # male
    ctrlfsample = 'C1563F'  # female
    outebv_EMSP = del_result(popmean, region, traits, ctrlmsample, ctrlfsample, 
                         blupebv,idx="AOA")

  }
  
  if (region %in% c('MMSP',"LMSP")){
    
    ctrlmsample = 'ZMN00545'  # male
    ctrlfsample = 'ZMN00080'  # female
    outebv_SP=del_result(popmean, region, traits,ctrlmsample, ctrlfsample, 
                         blupebv, idx="AOA")
    
    outebv_SP.tmp$FP_GCA=rbind(outebv_SP.tmp$FP_GCA,outebv_SP$FP_GCA)
    outebv_SP.tmp$MP_GCA=rbind(outebv_SP.tmp$MP_GCA,outebv_SP$MP_GCA)
    outebv_SP.tmp$SCA_GCA=rbind(outebv_SP.tmp$SCA_GCA,outebv_SP$SCA_GCA)
    
   }
  if (region %in% c( 'SCSU',"MCSU","NCSU")){
    ctrlmsample = 'Chang72'  # male
    ctrlfsample = 'Zheng58'  # female
    outebv_SU=del_result(popmean, region, traits, ctrlmsample, ctrlfsample, 
                         blupebv, idx="AOA")
    
    outebv_SU.tmp$FP_GCA=rbind(outebv_SU.tmp$FP_GCA,outebv_SU$FP_GCA)
    outebv_SU.tmp$MP_GCA=rbind(outebv_SU.tmp$MP_GCA,outebv_SU$MP_GCA)
    outebv_SU.tmp$SCA_GCA=rbind(outebv_SU.tmp$SCA_GCA,outebv_SU$SCA_GCA)
  }
  if (region == 'SWCN'){
    ctrlmsample = 'SCML0849'  # male
    ctrlfsample = 'ZNC442'  # female
    outebv_SW=del_result(popmean, region, traits, ctrlmsample, ctrlfsample, 
                         blupebv, idx="AOA")
    
  }
   
}  

h2=blupebv$h2
h2$LRTvalue = NA
h2=h2[,c(5,1:4,6)]
h2=h2 %>% mutate(across(where(is.numeric), ~round(.x,3)))

d2=blupebv$d2
d2=d2[,c(6,1:3,5,4)]
d2=d2 %>% mutate(across(where(is.numeric), ~round(.x, 3)))
h2_d2 = rbind(h2,d2)

va=blupebv$varcomp
va=va[,c(7,8,1,2,6)]
va=va %>% mutate(across(where(is.numeric), ~round(.x, 3)))
va=va[!(substr(va$item,1,7)=="units!R" ),]
va$item=rep(c("SCA方差","父本GCA方差","母本GCA方差","剩余方差"),nrow(va)/4)

popMean=popmean
popMean =popMean[!is.na(popMean$mean),]
popMean = popMean[,c(2,1,3)]
popMean=popMean %>% mutate(across(where(is.numeric), ~round(.x, 1)))

result=list()

traitdiease=c("RTLPCT", "TDPPCT", "GSPPCT","STKLPCT","ERTLPCT","LRTLPCT", "STKRPCT",
              "PMDPCT",'NCLB','SCLB','GLS',"RSTSOU",
              "SHBLSC","PMDSC")
traits_1 = unlist(lapply(traitdiease, function(x){c(str_c('mean_', x), 
                                                    str_c('max_',x),str_c('min_',x),
                                                    str_c('Loc_',x),str_c('Obs_',x)
)}))

result$FP_GCA=rbind(outebv_EMSP$FP_GCA,outebv_SP.tmp$FP_GCA,outebv_SU.tmp$FP_GCA,outebv_SW$FP_GCA)#
result$FP_GCA=merge(result$FP_GCA,parent_Fall[,c("GenoID",traits_1,idx)],by=c("GenoID",idx),all.x = TRUE)


result$FP_GCA <- result$FP_GCA %>% arrange(AOA, Rank_YLD)
# result$FP_GCA[is.na(result$FP_GCA)] <- ""
# result$FP_GCA <- result$FP_GCA[, colSums(result$FP_GCA != "") > 0]
#
columns_to_convert <- names(result$FP_GCA )[3:ncol(result$FP_GCA )]

# 循环遍历需要转换的列，并将其转换为数值型
for (col in columns_to_convert) {
  result$FP_GCA[[col]] <- as.numeric(result$FP_GCA[[col]])
}
result$FP_GCA <- result$FP_GCA %>% mutate(across(where(is.numeric), ~round(.x, 1)))
result$FP_GCA <-  result$FP_GCA[, colSums(is.na(result$FP_GCA)) != nrow(result$FP_GCA)]


result$MP_GCA=rbind(outebv_EMSP$MP_GCA,outebv_SP.tmp$MP_GCA,outebv_SU.tmp$MP_GCA,outebv_SW$MP_GCA)#
result$MP_GCA=merge(result$MP_GCA,parent_Mall[,c("GenoID",traits_1,idx)],by=c("GenoID",idx),all.x = TRUE)

# result$MP_GCA <- result$MP_GCA %>% mutate(across(where(is.numeric), ~round(.x, 1)))
result$MP_GCA <- result$MP_GCA %>% arrange(AOA, Rank_YLD)
# result$MP_GCA[is.na(result$MP_GCA)] <- ""
# result$MP_GCA <- result$MP_GCA[, colSums(result$MP_GCA != "") > 0]
columns_to_convert <- names(result$MP_GCA )[3:ncol(result$MP_GCA )]

# 循环遍历需要转换的列，并将其转换为数值型
for (col in columns_to_convert) {
  result$MP_GCA[[col]] <- as.numeric(result$MP_GCA[[col]])
}
result$MP_GCA <- result$MP_GCA %>% mutate(across(where(is.numeric), ~round(.x, 1)))
result$MP_GCA <-  result$MP_GCA[, colSums(is.na(result$MP_GCA)) != nrow(result$MP_GCA)]



result$SCA_GCA=rbind(outebv_EMSP$SCA_GCA,outebv_SP.tmp$SCA_GCA,outebv_SU.tmp$SCA_GCA,outebv_SW$SCA_GCA)#

diease.tmp=merge(diease,hyb_obs,by.x = c("Name",idx),by.y = c("Hybrid",idx),all.x = TRUE)
result$SCA_GCA=merge(result$SCA_GCA,diease.tmp[,c("Name",traits_1,idx)],by=c("Name",idx),all.x = TRUE)

#result$SCA_GCA=merge(result$SCA_GCA,hyb_obs[,c("Hybrid",traits_2,idx)],by.x = c("Name",idx),by.y = c("Hybrid",idx),all.x = TRUE)

#result$SCA_GCA <- result$SCA_GCA[order(-result$SCA_GCA$CKpct_YLD), ]

# columns_to_convert <- names(result$MP_GCA )[3:ncol(result$MP_GCA )]
# 
# # 循环遍历需要转换的列，并将其转换为数值型
# for (col in columns_to_convert) {
#   result$MP_GCA[[col]] <- as.numeric(result$MP_GCA[[col]])
# }

result$SCA_GCA <-  result$SCA_GCA[, colSums(is.na(result$SCA_GCA)) != nrow(result$SCA_GCA)]
result$SCA_GCA <- result$SCA_GCA %>% mutate(across(where(is.numeric), ~round(.x, 1)))

###添加类型
emsptype <- readxl::read_excel('datafiles/23YT-TrialType.xlsx', sheet = 1)
EMSP=result$SCA_GCA[substr(result$SCA_GCA$AOA,1,4)=="EMSP",]
for (i in 1:nrow(emsptype)) {
  EMSP$Type[EMSP$GenoID==emsptype$GenoID[i]] = emsptype$Type[i]
}

dim(EMSP)
EMSP=EMSP[,c(1,2,105,3:104)]

sptype <- readxl::read_excel('datafiles/23YT-TrialType.xlsx', sheet = 2)
SP=result$SCA_GCA[substr(result$SCA_GCA$AOA,1,4) %in% c("LMSP","MMSP"),]
for (i in 1:nrow(sptype)) {
  SP$Type[SP$GenoID==sptype$GenoID[i]] = sptype$Type[i]
}

dim(SP)
SP=SP[,c(1,2,105,3:104)]



sutype <- readxl::read_excel('datafiles/23YT-TrialType.xlsx', sheet = 3)
SU=result$SCA_GCA[substr(result$SCA_GCA$AOA,3,4)=="SU",]
for (i in 1:nrow(sutype)) {
  SU$Type[SU$GenoID==sutype$GenoID[i]] = sutype$Type[i]
}

dim(SU)
SU=SU[,c(1,2,105,3:104)]


swtype <- readxl::read_excel('datafiles/23YT-TrialType.xlsx', sheet =4)
SW=result$SCA_GCA[substr(result$SCA_GCA$AOA,1,4)=="SWCN",]
for (i in 1:nrow(swtype)) {
  SW$Type[SW$GenoID==swtype$GenoID[i]] = swtype$Type[i]
}

dim(SW)
SW=SW[,c(1,2,105,3:104)]

result$SCA_GCA=rbind(EMSP,SP,SU,SW)
result$SCA_GCA <- result$SCA_GCA %>% arrange(AOA, Rank_YLD)




result$h2_d2 = h2_d2
result$varcomp = va
result$popMean=popMean

file_path <- "Output/23Presult_GCA_SCA_1118.xlsx"
#modify(data=result,sheets=names(result))
openxlsx::write.xlsx(result,file =  "Output/23Presult_GCA_SCA_1118.xlsx") 

# file_path <- "Output/23P_GCA_SCA_1028.xlsx"


# file_path1 <- "Output/P_allTraits1028_SP.1.rds"
# 

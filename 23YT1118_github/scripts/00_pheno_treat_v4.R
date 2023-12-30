######从postgres数据库直接读取23年数据#####
# rm(list=ls())
# install.packages("RPostgres")
# install.packages("RPostgreSQL")
# library(DBI)
# 
# # Connect to a specific postgres database 
# con <- dbConnect(RPostgres::Postgres(),dbname = 'postgres', 
#                  host = '10.100.0.189', 
#                  port = 1235, # or any other port specified by your DBA
#                  user = 'postgres_readonly',
#                  password = 'qyng@2009')
# res <- dbSendQuery(con, "SELECT * FROM dbt.\"Yield_Trial_Master_Catalog-O-23-10-23\"")
# res <- dbFetch(res)
# res
# res[50332, ]
# res[74337,]


#########合并2023年表型数据和19-22年数据#######
#19-22年数据#
rm(list = ls())
library(dplyr)
p4 <- data.frame(read.table(textConnection(
  ' Name               Ped    Hybrid
QTY702 ZMN01045/ZMN00739   D208930  
QTY704 ZMN00662/ZMN00735   D208240
QTY705 ZMN00662/ZMN01212  D2181668
QTY706 ZMN00662/ZMN01187  D2181543
QTY707 ZMN00081/ZMN01236  D2181787
QTY805 ZMN01027/ZMN00154   D208897
QTY806 ZMN00081/ZMN01260   D208575
QTY808 ZMN01019/ZMN00739   D208902
QTY809 ZMN01045/ZMN01260   D218760
QTY811 ZMN00661/ZMN00392   D188186
QTY814 ZMN01128/ZMN01260  D2181251
QTY815 ZMN01021/ZMN00741   D217865
QTY902 ZMN00661/ZMN00547  D1980035
QTY801  ZMN00080/ZMN00762 D1980837
QTY703 ZMN00662/ZMN00739  D208254
QTY701  ZMN00659/ZMN00748 D208277
QTY901  ZMN00661/ZMN00762 D208309
QTY804  ZMN00430/ZMN00762 D1980833
QTY803  ZMN00430/ZMN01084 D208681
QTY807  ZMN00932/ZMN01263 D208564
QTY802  ZMN00932/ZMN01261 D208574
QTY810  ZMN01027/ZMN01307 D218958
QTY813  ZMN00662/ZMN01307 D218957
QTY812  ZMN00625/ZMN00762 D1980839'
),header = T,sep=''))



pheno1=readRDS('datafiles/19_22_P_AOA_Station.rds')
names(pheno1)
pheno1 <- pheno1 %>% select(-c('Set',"Rep","FEW","PMD","LDG"))
names(pheno1)[names(pheno1)=="Plot"]= "PlotID" 
names(pheno1)[names(pheno1)=="weight"]= "PLOTWT"
names(pheno1)[names(pheno1)=="YLD"]= "YLD14"
names(pheno1)[names(pheno1)=="FNSC"]= "FNSCNT"
names(pheno1)[names(pheno1)=="HARVEAR" ]= "EARCNT" 
names(pheno1)[names(pheno1)=="PMDP"  ]= "PMDCNT" 
names(pheno1)[names(pheno1)=="LDGP"  ]= "TDPCNT" 
names(pheno1)[names(pheno1)=="ERP"  ]= "EARTCNT"
names(pheno1)[names(pheno1)=="ER"  ]= "EARTSC"
names(pheno1)[names(pheno1)=="Discard"  ]= "PlotDiscarded"
names(pheno1)[names(pheno1)== "Rust"   ]= "RSTSOU"
# pheno1$AOA[pheno1$Station=="HRB"]="EMSP"
# pheno1$AOA[pheno1$Station=="CC"]="MLSP"
# pheno1$AOA[pheno1$Station=="NW"]="MLSP"
# pheno1$AOA[pheno1$Station=="JN"]="SU"
# pheno1$AOA[pheno1$Station=="XX"]="SU"
table(pheno1$AOA,pheno1$Station)



# pheno23$RTLPCT=round(pheno23$RTLCNT/pheno23$ERSCNT*100,3)
pheno1$TDPPCT=round(pheno1$TDPCNT/pheno1$FNSCNT*100,3)
pheno1$PMDPCT=(pheno1$PMDCNT/pheno1$FNSCNT)*100


traits <- c("TDPPCT", "PMDPCT")

# phenos_filtered.1 <- phenos_filtered %>% 
#   filter(if_any(all_of(traits), ~ . < 0))
# 
# phenos_new <- anti_join(phenos, phenos_filtered)

pheno1 <- pheno1 %>% 
  mutate(across(all_of(traits), ~ if_else(. > 100, 100, .)))

#2023年数据#

col_types <- c(rep("text",2),rep("numeric",2),rep("text",1),rep("numeric",3),rep("text",7),
               rep("numeric",6),rep("text",2),rep("numeric",39),rep("text",2))

# col_types <- c(rep("text",6),rep("numeric",5),rep("text",6),
#                rep("numeric",2),rep("text",1),
#                rep("numeric",6),rep("text",2),rep("numeric",58),rep("text",1))
pheno23 <- readxl::read_excel('datafiles/2023年产量测试数据-2023.11.18.xlsx', sheet = 1,col_types = col_types)

pheno23 <- pheno23 %>% select(-c( "PedID", "RecID", "Serp"))

pheno23$YLD14[pheno23$BookName=="SCST"]=pheno23$YLD14_BE[pheno23$BookName=="SCST"]
pheno23 <- pheno23 %>% select(-c('YLD14_BE'))

names(pheno23)[names(pheno23)=="pass"]="Pass"
names(pheno23)[names(pheno23)=="BookPrj"]="Station"
names(pheno23)[names(pheno23)=="TierNo"]="Field"
names(pheno23)[names(pheno23)=="TrialType"]="Type"
names(pheno23)[names(pheno23)=="Location"]="Location_C"
names(pheno23)[names(pheno23)=="BookName"]="Location"
pheno23$Location[pheno23$Location_C=="呼兰CC"]="HLJHULCC"
pheno23$Station[pheno23$Location_C=="呼兰CC"]="CC"

pheno23$Year=2023
pheno23$AOA=NA
pheno23$AOA <- ifelse(pheno23$Location %in% c("HLJJD", "HLJJMS", "HLJMS", "HLJWK", "HLJHUL", "HLJFY", "HLJHL"),
                      "EMSP",
                      pheno23$AOA)
pheno23$AOA <- ifelse(pheno23$Location %in% c("JLCLTC2", "JLDHTC2", "JLGZL", "JLYS", "JLNA", "JLSY","NMNM",
                                              "HLJHULCC","SXCZ","SXGP","SXXZA","SXXZB"),
                      "MMSP",pheno23$AOA)
pheno23$AOA <- ifelse(pheno23$Location %in% c("LNCT2", "LNCT1", "NMKEQ", "SXTG", "SXYC","SXQX"),
                      "LMSP",
                      pheno23$AOA)


pheno23$AOA <- ifelse(pheno23$Location %in% c("HNQF", "HNQB1", "HNQB2", "HNXX1", "HNXX2",
                                              "HNXX3","SDCX","HNYJ","HBCX","SDCQ1","SDCQ2","SDST",
                                              "SDZQ1","SDZQ2","HBRZ","SDGX1","SDGX2","SDWS1","SDWS2",
                                              "HBNJ2","HBNJ1","SDSG"),
                      "MCSU",
                      pheno23$AOA)

pheno23$AOA <- ifelse(pheno23$Location %in% c("HNZC", "HNXM1", "HNXM2", "HNXM3", "HNXM4",
                                              "HNXM5","HNXP","HNSP","HNXH","AHSX2"),
                      "SCSU",
                      pheno23$AOA)

pheno23$AOA <- ifelse(pheno23$Location %in% c("SCCD","SCST","GZAS","YNFM","YNZY"),
                      "SWCN",
                      pheno23$AOA)
pheno23$AOA <- ifelse(pheno23$Location %in% c("HBWQ1","HBWQ2","HBZX","HBJZ"),
                      "NCSU",
                      pheno23$AOA)





#pheno23$Name[pheno23$GenoID=='ZMN00824/ZMN00547']  <- 'D178070A'
pheno23$Name[pheno23$GenoID=='ZMN00786/ZMN00805']  <- 'D219005SW'
pheno23$Name[pheno23$GenoID=='ZMN01116/ZMN00435']  <- 'D219110SW'
pheno23$Name[pheno23$GenoID=='ZMN00805/ZMN00436']  <- 'D219118SW'
pheno23$Name[pheno23$GenoID=='ZMN00805/ZMN00789']  <- 'D219165SW'
pheno23$Name[pheno23$GenoID=='ZMN00375/ZMN00789']  <- 'D219166SW'
pheno23$Name[pheno23$GenoID=='ZMN01116/ZMN00789']  <- 'D219169SW'
pheno23$Name[pheno23$GenoID=='ZMN00081/ZMN00789']  <- 'D219171SW'
pheno23$Name[pheno23$GenoID=='ZMN00662/ZMN00789']  <- 'D219172SW'
pheno23$Name[pheno23$GenoID=='ZMN01027/ZMN00789']  <- 'D219175SW'
pheno23$Name[pheno23$GenoID=='ZMN00727/ZMN00789']  <- 'D219176SW'

# pheno23$AOA[pheno23$Station=="HRB"]="EMSP"
# pheno23$AOA[pheno23$Station=="CC"]="MLSP"
# pheno23$AOA[pheno23$Station=="NW"]="MLSP"
# pheno23$AOA[pheno23$Station=="JN"]="SU"
# pheno23$AOA[pheno23$Station=="XX"]="SU"
# pheno23$AOA[pheno23$Station=="SW"]="SWCN"

table(pheno23$AOA,pheno23$Station)
pheno23$FGenoID=NA
pheno23$MGenoID=NA
pheno23$EARTCNT=NA
pheno23$EARTSC=NA

pheno23$FGenoID <- as.character(pheno23$FGenoID)
pheno23$MGenoID <- as.character(pheno23$MGenoID)
pheno23$EARTCNT <- as.numeric(pheno23$EARTCNT)
pheno23$EARTSC <- as.numeric(pheno23$EARTSC)
#pheno23$RTLCNT[pheno23$Location=="SXYC"]=NA
#pheno23$TDPCNT=pheno23$LRTLCNT+pheno23$STKLCNT

pheno23$TDPCNT <- ifelse(is.na(pheno23$LRTLCNT), pheno23$STKLCNT, 
                         ifelse(is.na(pheno23$STKLCNT), pheno23$LRTLCNT, 
                                pheno23$LRTLCNT + pheno23$STKLCNT))

pheno23$ERSCNT=ifelse(is.na(pheno23$ERSCNT),pheno23$FNSCNT,pheno23$ERSCNT)
# pheno23$RTLPCT=round(pheno23$RTLCNT/pheno23$ERSCNT*100,3)
pheno23$TDPPCT=round(pheno23$TDPCNT/pheno23$ERSCNT*100,3)


pheno23$NCLB[pheno23$NCLB==0]=NA
pheno23$GLS[pheno23$GLS==0]=NA


traits <- c("RTLPCT", "TDPPCT", "GSPPCT","STKLPCT","ERTLPCT","LRTLPCT", "STKRPCT", "PMDPCT")


phenos_new <- pheno23 %>% 
  mutate(across(all_of(traits), ~ if_else(. > 100, 100, .)))



for (i in 1:nrow(p4)){
  phenos_new$Name[phenos_new$Name==p4$Hybrid[i]] <- p4$Name[i]
  
}

phenos_new=phenos_new[!is.na(phenos_new$Name),]

saveRDS(phenos_new,"datafiles/23phenos_1118.rds")
writexl::write_xlsx(phenos_new, path='datafiles/23phenos_1118.xlsx')


pheno23.TC1=phenos_new[phenos_new$Type=="TC1",]
pheno23.TC1=pheno23.TC1[!substr(pheno23.TC1$Name,6,8)=="Liu",]
pheno23.TC1=pheno23.TC1[!is.na(pheno23.TC1$Name),]
saveRDS(pheno23.TC1,"datafiles/23phenos_TC1_1118.rds")
writexl::write_xlsx(pheno23.TC1, path='datafiles/23phenos_TC1_1118.xlsx')


pheno23.TC2=phenos_new[phenos_new$Type=="TC2",]
pheno23.TC2=pheno23.TC2[!substr(pheno23.TC2$Name,6,8)=="Liu",]
TC2CK=phenos_new[phenos_new$IsCK=="CK"&phenos_new$AOA=="EMSP"&phenos_new$Type=="P1",]
TC2CK=TC2CK[!is.na(TC2CK$Name),]
pheno23.TC2=rbind(pheno23.TC2,TC2CK)
pheno23.TC2=pheno23.TC2[!is.na(pheno23.TC2$Name),]
pheno23.TC2=unique(pheno23.TC2)
saveRDS(pheno23.TC2,"datafiles/23phenos_TC2_1118.rds")
writexl::write_xlsx(pheno23.TC2, path='datafiles/23phenos_TC2_1118.xlsx')


##合并##
setdiff(names(phenos_new), names(pheno1))#在 pheno23 中但不在 pheno1 中的列名

for(i in c("Location_C" ,"Pedigree", "Sub","PLTDAT","EMGDAT" ,
           "SDVSC","EVG" ,"ERSCNT" ,"苗期综合评价","SLKDAT",
           "SLKDAY" , "POLDAT" ,"POLDAY", "GSPCNT","GSPPCT" ,
           "STKLCNT" ,"STKLPCT","ERTLCNT" ,"ERTLPCT","LRTLCNT" ,
           "LRTLPCT","RTLCNT" , "RTLPCT", "SHBLSC" , "STKRCNT","STKRPCT" , "PMDSC", "DISESC"        
           ,"MLKDAT" , "MLKDAY", "BREDSC",  "HarvNote"      
           ,"收获前综合评价")){
  pheno1[[i]] <- NA
}
#pheno1=pheno1[,names(phenos_new)]
phenos_new=phenos_new[,names(pheno1)]

names(phenos_new)
names(pheno1)
phenos=rbind(pheno1,phenos_new)


for (i in 1:nrow(p4)){
  phenos$Name[phenos$Name==p4$Hybrid[i]] <- p4$Name[i]
  
}

phenos=phenos[!is.na(phenos$Name),]
saveRDS(phenos,"datafiles/19_23_P_phenos_1118.rds")
writexl::write_xlsx(phenos, path='datafiles/19_23_P_phenos_1118.xlsx')






#########TC1互种问题#######
rm(list = ls())
diease.tmp <- readRDS('datafiles/23phenos_TC1_1118.rds')

 xx_zq=readRDS("datafiles/seedAllocation_Sites_TC1_XX_20230418_ZQreallocated(1).rds")
 jn_xx=readxl::read_xlsx("datafiles/seedAllocation_Sites_TC1_JN_20230418_XXreallocated(1).xlsx",sheet = "新乡")
 name=xx_zq$章丘$Name
 namejn=jn_xx$Name

 zq=diease.tmp[diease.tmp$Location %in% c("SDZQ1","SDZQ2"),]
 tmp=zq[!(zq$Name %in% name),]
 zq$Station="XX"



 xx=diease.tmp[diease.tmp$Location %in% c("HNXX1","HNXX3"),]
 tmp.xx=xx[!(xx$Name %in% namejn),]
 xx$Station="JN"

 all=diease.tmp[!(diease.tmp$Location %in% c("HNXX1","HNXX3","SDZQ1","SDZQ2")),]

 diease.tmp.o=rbind(zq,tmp,xx,tmp.xx,all)
 diease.tmp.o=unique(diease.tmp.o)
 saveRDS(diease.tmp.o,'datafiles/23phenos_TC1_1118reallocated.rds')


 
#########TC2互种问题########
diease.tmp <- readRDS('datafiles/23phenos_TC2_1118.rds')

xx_zq=readxl::read_xlsx("datafiles/seedAllocation_Sites_TC2_XX_20230404_ZQreallocated.xlsx",sheet = "章丘")
jn_xx=readxl::read_xlsx("datafiles/seedAllocation_Sites_TC2_JN_20230402_XXreallocated(2).xlsx",sheet = "新乡")

name=xx_zq$Name
namejn=jn_xx$Name

zq=diease.tmp[diease.tmp$Location %in% c("SDZQ1","SDZQ2"),]
tmp=zq[!(zq$Name %in% name),]
zq$Station="XX"

xx=diease.tmp[diease.tmp$Location %in% c("HNXX1","HNXX2","HNXX3"),]
tmp.xx=xx[!(xx$Name %in% namejn),]
xx$Station="JN"

all=diease.tmp[!(diease.tmp$Location %in% c("HNXX1","HNXX2","HNXX3","SDZQ1","SDZQ2")),]

diease.tmp.o=rbind(zq,tmp,xx,tmp.xx,all)
diease.tmp.o=unique(diease.tmp.o)

table(diease.tmp$Station,diease.tmp$Location)
saveRDS(diease.tmp.o,'datafiles/23phenos_TC2_1118reallocated.rds')

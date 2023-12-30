rm(list=ls())
setwd("D:/YangZhenyu/Project/00/Breeding_Anova")
source("./src/anova.sup.v5.R")
# Reading data ------------------------------------------------------------

p4 <- data.frame(read.table(textConnection(
    'Name               Ped   Hybrid
 QTY702 ZMN01045/ZMN00739  D208930
 QTY704 ZMN00662/ZMN00735  D208240
 QTY705 ZMN00662/ZMN01212 D2181668
 QTY706 ZMN00662/ZMN01187 D2181543
 QTY707 ZMN00081/ZMN01236 D2181787
 QTY801 ZMN00080/ZMN00762 D1980837
 QTY804 ZMN00430/ZMN00762 D1980833
 QTY805 ZMN01027/ZMN00154  D208897
 QTY806 ZMN00081/ZMN01260  D208575
 QTY808 ZMN01019/ZMN00739  D208902
 QTY809 ZMN01045/ZMN01260  D218760
 QTY810 ZMN01027/ZMN01307	  D218958
 QTY811 ZMN00661/ZMN00392  D188186
 QTY812 ZMN00625/ZMN00762 D1980839
 QTY813 ZMN00662/ZMN01307	  D218957
 QTY814 ZMN01128/ZMN01260 D2181251
 QTY815 ZMN01021/ZMN00741  D217865
 QTY902 ZMN00661/ZMN00547 D1980035
 DK159  ZMN00545/ZMN00080    DK159
 XY335  ZMN00333/ZMN00154    XY335
 XY1483 XXX                 XY1483
 ZD958  XXX                  ZD958
 CD99   XXX                  CD99
 ZY335  XXX                  ZY335
 DH605  XXX                  DH605
 XY1446 XXX                 XY1446
 DMY1   XXX                   DMY1
 DMY3   XXX                   DMY3
 C1563  XXX                  C1563'
),header = T,sep=''))

p4

bluelist.all <- readRDS('datafiles/2023_11_03_19-23P_BLUE1103/23P_datalist_bluelist.rds')  #7个熟期的数据都包括了
bluelist <- bluelist.all$bluelist$blue
table(bluelist$YLD14$AOA)

# EMSP  LMSP  MCSU  MMSP  NCSU  SCSU  SWCN
# 2906 17888 50690 36977  8427 30562   805

phenlist <- readRDS('datafiles/2023_11_03_19-23P_BLUE1103/19-23P_datclean1103.rds')  #7个熟期的数据都包括了
aoas <- unique(bluelist$YLD14$AOA);aoas <- aoas[order(aoas)]
table(phenlist$YLD14$AOA)

# EMSP  LMSP  MCSU  MMSP  NCSU  SCSU  SWCN
# 6317 26974 84050 57012 15637 46080  1848

p3list <- read.csv('datafiles/2023_11_03_19-23P_BLUE1103/23YT DHB-HHH-last_20230311_P3P4(1).csv')

p3list <- p3list[p3list$Category %in% c('P3', 'TD'),]
# p3 <- c(p3,'DK159','XY1483','XY335','ZD958','CD99','ZY335','DH605','XY1446','DMY1','DMY3','C1563')
dim(p3list) #212

# Comparision with BLUP EBV
blupEBV <- data.frame(read_excel('./output/P_allTraits1026_SP_diease_YLD(1).xlsx',sheet = 3))

traits <- c("YLD14","MST","PHT","EHT")

# aoas <- c('LMSP','MMSP','MCSU','SCSU','SP','SU')  #EMSP,NCSU,SWCN have singularity problem due to limited number of observations
# aoas <- c('EMSP','SWCN')  #EMSP,NCSU,SWCN have singularity problem due to limited number of observations
# aoa <- 'SP'
# aoa <- 'EMSP'

aoas <- c('SP','MMSP','LMSP','SU','NCSU','MCSU','SCSU','EMSP','SWCN')  #EMSP,NCSU,SWCN have singularity problem due to limited number of observations

# controls <- data.frame(AOA=c('SWCN','EMSP','LMSP','MMSP','NCSU','MCSU','SCSU','SP','SU'),
#                        Name=c('CD99','DMY3','DK159','DK159','ZD958','ZD958','ZD958','DK159','ZD958'))

controls <- data.frame(AOA=c('SWCN', 'EMSP','LMSP', 'MMSP', 'NCSU', 'MCSU', 'SCSU',  'SP',   'SU'),
                           Name=c('CD99','DMY3','DK159','DK159','DH605','DH605','DH605','DK159','DH605'))


# aoa <- "SP"
"C1563"
"DMY3"

pdf(file = "./output/YLD_MST.Compare.pdf", width = 8, height = 5)

res.all  <- NULL
res.name <- NULL
for (aoa in aoas){
    if (aoa %in% c('LMSP','MMSP','SP')){
        p3 <- p3list$HY[p3list$SP=='SP']
    }else if(aoa %in% c('MCSU','SCSU','SU')){
        p3 <- unique(p3list$HY[p3list$JN=='JN' | p3list$XX=='XX'])
    } else {
        p3 <- p3list$HY[p3list$HLJ=='HLJ']
    }
    if(length(p3) != 0){
        p3 <- unique(c(p3,'DK159','XY1483','XY335','ZD958','CD99','ZY335','DH605','XY1446','DMY1','DMY3','C1563'))
        p3.tmp <- P3Estimation(p3 = p3,
                     bluelist = bluelist,
                     phenlist = phenlist,
                     traits = traits,
                     aoa = aoa,
                     blupEBV =blupEBV,
                     control = controls$Name[controls$AOA==aoa])
        res.all <- c(res.all, list(p3.tmp))
        res.name <- c(res.name, paste0(aoa, "_P3"))
    }
    if(length(p4) != 0){
        p4.tmp <- P4Estimation(p4 = p4,
                               bluelist = bluelist,
                               phenlist = phenlist,
                               traits   = traits,
                               aoa = aoa,
                               blupEBV = blupEBV,
                               control = controls$Name[controls$AOA==aoa])
        res.all <- c(res.all, list(p4.tmp))
        res.name <- c(res.name, paste0(aoa, "_P4"))
    }
}

dev.off()

names(res.all) <- res.name

openxlsx::write.xlsx(res.all, file = "./output/2023年_P3试验_杂交种_方差分析-20231104-K1-1.xlsx")
#===================================================================================================
#添加系谱
library(DT)
library(dplyr)
library(openxlsx)
rm(list = ls())

bind.res <- NULL
all.sheet.name <- openxlsx::getSheetNames(file = "./output/2023年_P3试验_杂交种_方差分析-20231104-K1-1.xlsx")
for (sheet in 1:length(all.sheet.name)) {
    res.tmp <- openxlsx::read.xlsx(xlsxFile = "./output/2023年_P3试验_杂交种_方差分析-20231104-K1-1.xlsx", sheet = sheet)
    split.res.tmp <- unlist(strsplit(all.sheet.name[sheet], "_"))
    res.tmp$AOA   <- split.res.tmp[1]
    res.tmp$Type  <- split.res.tmp[2]
    bind.res <- rbind(bind.res, res.tmp)
}

P3P4.Anova.all <- cbind(Name = bind.res[,1],
                        bind.res[,20:21],
                        bind.res[,2:19])

P3P4.Anova.all.0 <- unique(P3P4.Anova.all[,1:21])

P3P4.Anova.all.0[which(is.na(P3P4.Anova.all.0$Phen_YLD14)), 9] <- NA
P3P4.Anova.all.0[which(is.na(P3P4.Anova.all.0$Phen_MST)), 13]  <- NA
P3P4.Anova.all.0[which(is.na(P3P4.Anova.all.0$Phen_PHT)), 17]  <- NA
P3P4.Anova.all.0[which(is.na(P3P4.Anova.all.0$Phen_EHT)), 21]  <- NA

P3P4.Anova.all.0$ID <- paste0(P3P4.Anova.all.0$Name, "_", P3P4.Anova.all.0$Type)

pheno.all  <- openxlsx::read.xlsx("./datafiles/2023_11_03_19-23P_BLUE1103/2023年产量测试数据-2023.11.3.xlsx")
hybred_ped <- data.frame(ID = paste0(pheno.all$Name, "_", pheno.all$TrialType),
                         Pedigree = pheno.all$Pedigree)

hybred_ped <- unique(hybred_ped)

P3P4.merge.ped <- merge(P3P4.Anova.all.0, hybred_ped, by = "ID", all.x = T)

P3P4.merge.ped <- P3P4.merge.ped[,-1]

P3P4.merge.ped <- cbind(P3P4.merge.ped[,1:5],
                        Pedigree = P3P4.merge.ped[,22],
                        P3P4.merge.ped[,6:21])

openxlsx::write.xlsx(P3P4.merge.ped, "./output/2023年-P3试验-杂交种方差分析及病害统计-20231104-K2-1.xlsx")


#====================================================================================================
#合并病害数据+排序
library(DT)
library(dplyr)
library(openxlsx)
rm(list = ls())

pheno.all      <- readRDS(file = "./datafiles/2023_11_03_19-23P_BLUE1103/1103P_diease_by_23result.rds")
P3P4.Anova.all.0 <-  openxlsx::read.xlsx(xlsxFile = "./output/2023年-P3试验-杂交种方差分析及病害统计-20231104-K2-1.xlsx")

pheno.all$NAYL        <- paste0(pheno.all$Name, "_", pheno.all$AOA, "_", pheno.all$Year, "_", pheno.all$Location)
pheno.all.0 <- pheno.all[,-c(1:4)]



P3P4.Anova.all.0$NAYL <- paste0(P3P4.Anova.all.0$Name, "_", P3P4.Anova.all.0$AOA, "_", P3P4.Anova.all.0$Year, "_", P3P4.Anova.all.0$Location)

merge.res <- merge(P3P4.Anova.all.0, pheno.all.0, by = "NAYL", all.x = T)
merge.res <- unique(merge.res)

P3P4.Anova.res <- merge.res[,-1]

AOA.vec <- c('SP','MMSP','LMSP','SU','NCSU','MCSU','SCSU','EMSP','SWCN')

res.all   <- NULL
ped.error <- NULL

for (P in c("P3")) {
    P3.Anova.res <- P3P4.Anova.res %>% dplyr::filter(Type == P)
    res2 <- NULL
    for (aoa in AOA.vec) {
        P3.Anova.Aoa.res <- P3.Anova.res %>% dplyr::filter(AOA == aoa)
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
                res1.2 <- res1.2[order(res1.2$Year, res1.2$Location), ]
                res1.3 <- rbind(res1.2, res1.3)
                res1   <- rbind(res1, res1.3)
            }
            res2 <- rbind(res2, res1)
        }
    }
    res.all <- rbind(res.all, res2)
}

openxlsx::write.xlsx(res.all, "./output/2023年-P3试验-杂交种方差分析及病害统计-20231104-K3-1.xlsx")








library(DT)
library(stringr)
library(dplyr)
rm(list = ls())
setwd("D:/Programs/Rproject/Breeding_Analysis")
Rcpp::sourceCpp(file = "./pedigreeCleaning.cpp")

Pedigree.all <- openxlsx::read.xlsx(xlsxFile = "./data/Yield Trial Master Catalog-2023.10.28-1.xlsx")

TC1.Pedigree <- Pedigree.all %>% dplyr::filter(TrialType == "TC1")

TC1.All.num <- dim(TC1.Pedigree)[1]

#
xxx.idx <- which(grepl(pattern = "\\(", x = TC1.Pedigree$Pedigree))

TC1.Pedigree$Pedigree[xxx.idx]

# str_replace_all("ZJSS084(M753M/DMY3M)/ZJNS050(DK517F/51417B)", "(?<=\\().+?(?=\\))", "xxx")

#处理带括号的系谱
#-,=,%代表3种不同类型的错误
unusual.ped.fun <- function(Name, Ped){

    JString <- Ped
    MString.df <- str_match_all(JString, "(?<=\\().+?(?=\\))")[[1]]
    LString.df <- str_locate_all(JString, "(?<=\\().+?(?=\\))")[[1]]

    for (jml in 1:dim(MString.df)[1]) {
        substring(JString, first = LString.df[jml,1], last = LString.df[jml,2]) <- gsub("/", "|", MString.df[jml,1])
        # print(JString)
    }
    ped.split  <- pedigree(name = Name, GenoID = JString)
    ped.len    <- dim(ped.split)[1]
    ped.new    <- ped.split[ped.len, ]
    ped.new$P1 <- gsub("\\|", "/", ped.new$P1)
    ped.new$P2 <- gsub("\\|", "/", ped.new$P2)
    for (i in 1:(ped.len-1)) {
        if(grepl("\\|", ped.split[i, ]$name)){
            Xname <- gsub("\\|", "/", ped.split[i, ]$name)
            XString.df <- str_match_all(Xname, "(?<=\\().+?(?=\\))")[[1]]
            ped.tmp1z  <- pedigree(name = Xname, GenoID = XString.df)
            ped.new    <- rbind(ped.tmp1z, ped.new)
        }else{
            ped.new <- rbind(ped.split[i, ], ped.new)
        }
    }
    return(ped.new)
}

#拆分每个样本的系谱
TC1.Ped.all = NULL
for (i in 1:TC1.All.num) {

    TC1.Pedigree[i,]$Pedigree <- gsub("（", "(", TC1.Pedigree[i,]$Pedigree)
    TC1.Pedigree[i,]$Pedigree <- gsub("）", ")", TC1.Pedigree[i,]$Pedigree)

    TC1.Pedigree[i,]$Pedigree <- toupper(TC1.Pedigree[i,]$Pedigree)

    braces.R.num <- str_count(TC1.Pedigree[i,]$Pedigree, "\\(")
    braces.L.num <- str_count(TC1.Pedigree[i,]$Pedigree, "\\)")

    if(is.na(TC1.Pedigree[i,]$Pedigree)){
        tc1.ped.split <- data.frame(name = TC1.Pedigree[i,]$Name,
                                    P1   = "-",
                                    P2   = "-")
        print(paste0("Error -", TC1.Pedigree[i,]$Name, "   ", TC1.Pedigree[i,]$Pedigree))
    }else{
        if((braces.R.num == braces.L.num) & (braces.R.num >= 1) & (braces.L.num >= 1)){
            tc1.ped.split <- unusual.ped.fun(Name = TC1.Pedigree[i,]$Name,
                                             Ped  = TC1.Pedigree[i,]$Pedigree)
        }else if(braces.R.num != braces.L.num){
            tc1.ped.split <- data.frame(name = TC1.Pedigree[i,]$Name,
                                        P1   = "=",
                                        P2   = "=")
            print(paste0("Error =",TC1.Pedigree[i,]$Name, "   ", TC1.Pedigree[i,]$Pedigree))
        }else{

            tc1.ped.split <- pedigree(name   = TC1.Pedigree[i,]$Name,
                                      GenoID = TC1.Pedigree[i,]$Pedigree)

            xxxx.idx <- dim(tc1.ped.split)[1]
            if(grepl("-|/", tc1.ped.split[xxxx.idx, ]$P1) | grepl("-|/", tc1.ped.split[xxxx.idx, ]$P2) | grepl("liu", tc1.ped.split[xxxx.idx, ]$name) | tc1.ped.split[xxxx.idx, ]$P1 == 0 | tc1.ped.split[xxxx.idx, ]$P2 == 0 | grepl("Liu", tc1.ped.split[xxxx.idx, ]$name) | !grepl("//!",  TC1.Pedigree[i,]$Pedigree)){
                #对追溯中的父母本进行校正
                error.idx1 <- which((tc1.ped.split$P1!=0)&(tc1.ped.split$P2==0))
                if(length(error.idx1)!=0){
                    err1.ped <- tc1.ped.split[error.idx1, ]
                    tc1.ped.split[error.idx1, 2] <- 0
                    tc1.ped.split <- tc1.ped.split %>% dplyr::filter(!name %in% err1.ped$P1)
                }
                error.idx2 <- which((tc1.ped.split$P2!=0)&(tc1.ped.split$P1==0))
                if(length(error.idx2)!=0){
                    err2.ped <- tc1.ped.split[error.idx1, ]
                    tc1.ped.split[error.idx2, 3] <- 0
                    tc1.ped.split <- tc1.ped.split %>% dplyr::filter(!name %in% err2.ped$P2)
                }
            }else{
                tc1.ped.split <- data.frame(name = TC1.Pedigree[i,]$Name,
                                            P1   = "%",
                                            P2   = "%")
                print(paste0("Error %",TC1.Pedigree[i,]$Name, "   ", TC1.Pedigree[i,]$Pedigree))
            }
        }
    }
    TC1.Ped.all <- rbind(TC1.Ped.all, tc1.ped.split)
}
#去除重复
TC1.Ped.all.unique <- unique(TC1.Ped.all)
openxlsx::write.xlsx(TC1.Ped.all.unique, file = "./output/Pedigree_Split_TC1_2023.10.28.xlsx")

#---------------检查是否每个测验种都进行了系谱拆分---------------
which(!TC1.Pedigree$Name %in% TC1.Ped.all.unique$name)

#---------------检查是否有重名但系谱不一致的测验种---------------
duplicate_name <-  TC1.Ped.all.unique %>% dplyr::group_by(name) %>% dplyr::summarise(freq = n()) %>% dplyr::filter(freq > 1) %>% select(name)
duplicate_name

#---------------提取没有系谱信息（-/-）的测验种---------------
Ped.error.D1  <- TC1.Ped.all.unique %>% dplyr::filter(P1 == "-")
Err.D1.output <- Pedigree.all  %>% dplyr::filter(Name %in% Ped.error.D1$name)
openxlsx::write.xlsx(Err.D1.output, file = "./output/2023_TC1_没有系谱信息.xlsx")

#---------------提取系谱中括号不完整（=/=）的测验种---------------
Ped.error.D2  <- TC1.Ped.all.unique %>% dplyr::filter(P1 == "=")
Err.D2.output <- Pedigree.all  %>% dplyr::filter(Name %in% Ped.error.D2$name)
openxlsx::write.xlsx(Err.D2.output, file = "./output/2023_TC1_括号不完整.xlsx")

#---------------提取系谱可能不完整（系谱显示是个DH系）的测验种---------------
Ped.error.D3  <- TC1.Ped.all.unique %>% dplyr::filter(P1 == "%")
Err.D3.output <- Pedigree.all  %>% dplyr::filter(Name %in% Ped.error.D3$name)
openxlsx::write.xlsx(Err.D3.output, file = "./output/2023_TC1_缺少Tester.xlsx")

#---------------检查是否有的父母本没有系谱信息，如果没有则其父母本为0/0---------------
Ped.final <- TC1.Ped.all.unique %>% dplyr::filter(P1 != "-") %>% dplyr::filter(P1 != "=") %>% dplyr::filter(P1 != "%")

rm.zero.ped <- Ped.final %>% dplyr::filter( (P1 != 0)&(P2 != 0) )
rm.zero.ped[which(!rm.zero.ped$P1 %in% TC1.Ped.all.unique$name), ]
rm.zero.ped[which(!rm.zero.ped$P2 %in% TC1.Ped.all.unique$name), ]

openxlsx::write.xlsx(Ped.final, file = "./output/Filter_Pedigree_Split_TC1_2023.10.28.xlsx")



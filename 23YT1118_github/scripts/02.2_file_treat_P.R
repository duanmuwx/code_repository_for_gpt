rm(list = ls())

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
QTY812  ZMN00625/ZMN00762 D1980839
QTY813  ZMN00662/ZMN01307 D218957'
),header = T,sep=''))


###########提取19-22年P试验去除异常值的数据###############
data <- readRDS("datafiles/19_22P_datalist_bluelist.rds")
datclean <- lapply(data, function(x) x$datalist$datclean)

years = names(data)
features = c("YLD14", "MST","PHT", "EHT", "TWT")  # 特征名称

#按性状合并多年份的数据
datclean_merge<-list()
for (year in years) {
  for (feature in features){
    datclean_merge[[feature]] <- rbind(datclean_merge[[feature]],
                                    datclean[[year]][[feature]]
                                   
    )
  }
}

####加入23年P试验数据####
data23=readRDS("Output/23P_datalist_bluelist.rds")

datclean23 <- data23$datalist$datclean
datclean_merge.tmp=list()
for (feature in features){
  datclean_merge.tmp[[feature]] <- rbind(datclean_merge[[feature]],
                                     datclean23[[feature]])
   # datclean_merge.tmp[[feature]]$AOA[datclean_merge.tmp[[feature]]$Station == "CC"] = "MLSP"
   # datclean_merge.tmp[[feature]]$AOA[datclean_merge.tmp[[feature]]$Station == "NW"] = "MLSP"
   # datclean_merge.tmp[[feature]]$AOA[datclean_merge.tmp[[feature]]$Station == "HRB"] = "EMSP"                              
   # datclean_merge.tmp[[feature]]$AOA[datclean_merge.tmp[[feature]]$Station == "JN"] = "SU"
   # datclean_merge.tmp[[feature]]$AOA[datclean_merge.tmp[[feature]]$Station == "XX"] = "SU"
   # datclean_merge.tmp[[feature]]$AOA[datclean_merge.tmp[[feature]]$Station == "SW"] = "SWCN"
}

for (i in 1:nrow(p4)){
  datclean_merge.tmp$YLD14$Name[datclean_merge.tmp$YLD14$Name==p4$Hybrid[i]] <- p4$Name[i]
  datclean_merge.tmp$MST$Name[datclean_merge.tmp$MST$Name==p4$Hybrid[i]] <- p4$Name[i]
  datclean_merge.tmp$PHT$Name[datclean_merge.tmp$PHT$Name==p4$Hybrid[i]] <- p4$Name[i]
  datclean_merge.tmp$EHT$Name[datclean_merge.tmp$EHT$Name==p4$Hybrid[i]] <- p4$Name[i]
  datclean_merge.tmp$TWT$Name[datclean_merge.tmp$TWT$Name==p4$Hybrid[i]] <- p4$Name[i]
}

saveRDS(datclean_merge.tmp,"Output/19-23P_datclean1118.rds") 

#################提取P试验BLUE的数据##############
#data <- readRDS("Output/19_22P_datalist_bluelist.rds")
bluedat <- lapply(data, function(x) x$bluelist$blue)

years = names(data)
features = c("YLD14","MST",  "PHT", "EHT", "TWT")  # 特征名称

#按性状合并多年份的数据
bluedat_merge<-list()
for (year in years) {
  for (feature in features){
    bluedat_merge[[feature]] <- rbind(bluedat_merge[[feature]],
                                      bluedat[[year]][[feature]] )
 
  }
}
####加入23年P试验数据####
bluedat23 <- data23$bluelist$blue
bluedat_merge.tmp=list()
for (feature in features){
  bluedat_merge.tmp[[feature]] <- rbind(bluedat_merge[[feature]],
                                         bluedat23[[feature]] )
  # bluedat_merge.tmp[[feature]]$AOA[bluedat_merge.tmp[[feature]]$Station == "CC"] = "MLSP"
  # bluedat_merge.tmp[[feature]]$AOA[bluedat_merge.tmp[[feature]]$Station == "NW"] = "MLSP"
  # bluedat_merge.tmp[[feature]]$AOA[bluedat_merge.tmp[[feature]]$Station == "HRB"] = "EMSP" 
  # bluedat_merge.tmp[[feature]]$AOA[bluedat_merge.tmp[[feature]]$Station == "JN"] = "SU"
  # bluedat_merge.tmp[[feature]]$AOA[bluedat_merge.tmp[[feature]]$Station == "XX"] = "SU"
  # bluedat_merge.tmp[[feature]]$AOA[bluedat_merge.tmp[[feature]]$Station == "SW"] = "SWCN"
}



for (i in 1:nrow(p4)){
  bluedat_merge.tmp$YLD14$Name[bluedat_merge.tmp$YLD14$Name==p4$Hybrid[i]] <- p4$Name[i]
  bluedat_merge.tmp$MST$Name[bluedat_merge.tmp$MST$Name==p4$Hybrid[i]] <- p4$Name[i]
  bluedat_merge.tmp$PHT$Name[bluedat_merge.tmp$PHT$Name==p4$Hybrid[i]] <- p4$Name[i]
  bluedat_merge.tmp$EHT$Name[bluedat_merge.tmp$EHT$Name==p4$Hybrid[i]] <- p4$Name[i]
  bluedat_merge.tmp$TWT$Name[bluedat_merge.tmp$TWT$Name==p4$Hybrid[i]] <- p4$Name[i]
}

saveRDS(bluedat_merge.tmp,"Output/19-23P_BLUE1118.rds") 





###########提取19-22年P试验异常值标签的数据###############
data <- readRDS("datafiles/19_22P_datalist_bluelist.rds")
datLabelled <- lapply(data, function(x) x$datalist$datLabelled)

years = names(data)
features = c("YLD14", "MST","PHT", "EHT", "TWT","PLOTWT")  # 特征名称

#按性状合并多年份的数据
datLabelled_merge<-list()
for (year in years) {
  for (feature in features){
    datLabelled_merge[[feature]] <- rbind(datLabelled_merge[[feature]],
                                          datLabelled[[year]][[feature]]
                                       
    )
  }
}

####加入23年P试验异常值标签的数据####
data23=readRDS("Output/23P_datalist_bluelist.rds")

datLabelled23 <- data23$datalist$datLabelled
datLabelled_merge.tmp=list()
for (feature in features){
  datLabelled_merge.tmp[[feature]] <- rbind(datLabelled_merge[[feature]],
                                         datLabelled23[[feature]])
  datLabelled_merge.tmp[[feature]] =datLabelled_merge.tmp[[feature]][!is.na(datLabelled_merge.tmp[[feature]]$Name),]
  # datclean_merge.tmp[[feature]]$AOA[datclean_merge.tmp[[feature]]$Station == "CC"] = "MLSP"
  # datclean_merge.tmp[[feature]]$AOA[datclean_merge.tmp[[feature]]$Station == "NW"] = "MLSP"
  # datclean_merge.tmp[[feature]]$AOA[datclean_merge.tmp[[feature]]$Station == "HRB"] = "EMSP"                              
  # datclean_merge.tmp[[feature]]$AOA[datclean_merge.tmp[[feature]]$Station == "JN"] = "SU"
  # datclean_merge.tmp[[feature]]$AOA[datclean_merge.tmp[[feature]]$Station == "XX"] = "SU"
  # datclean_merge.tmp[[feature]]$AOA[datclean_merge.tmp[[feature]]$Station == "SW"] = "SWCN"
}

for (i in 1:nrow(p4)){
  datLabelled_merge.tmp$YLD14$Name[datLabelled_merge.tmp$YLD14$Name==p4$Hybrid[i]] <- p4$Name[i]
  datLabelled_merge.tmp$MST$Name[datLabelled_merge.tmp$MST$Name==p4$Hybrid[i]] <- p4$Name[i]
  datLabelled_merge.tmp$PHT$Name[datLabelled_merge.tmp$PHT$Name==p4$Hybrid[i]] <- p4$Name[i]
  datLabelled_merge.tmp$EHT$Name[datLabelled_merge.tmp$EHT$Name==p4$Hybrid[i]] <- p4$Name[i]
  datLabelled_merge.tmp$TWT$Name[datLabelled_merge.tmp$TWT$Name==p4$Hybrid[i]] <- p4$Name[i]
}

#saveRDS(datLabelled_merge.tmp,"Output/19-23P_datLabelled1118.rds") 
writexl::write_xlsx(datLabelled_merge.tmp, path='datafiles/19_23_P_datLabelled1118.xlsx')

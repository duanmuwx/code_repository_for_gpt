
#########merge pdf#####
library(qpdf)
# 载入pdftk库


# 定义要合并的PDF文件路径
pdf_files <- list.files("Output/", pattern = "_plot_", full.names = TRUE)

# 使用pdftk将PDF文件合并为一个PDF文件
output_file <- "Output/23数据分布图.pdf"
pdf_combine(pdf_files, output_file)





#########merge P blue#####


rm(list=ls())
# 从文件夹中读取所有符合条件的 RDS 文件
file_names <- list.files("Output/", pattern = "_P_datalist_bluelist.rds$", full.names = TRUE)

# 创建一个空列表，用于存储读取的 RDS 文件内容
data_list <- list()

# 读取每个 RDS 文件，并将其内容添加到 data_list 列表中
for (file in file_names) {
  tmp <- readRDS(file)
  data_list[[length(data_list) + 1]] <- tmp
}

# 将 data_list 中所有内容合并到一个数据框中
traits=c("YLD14",	"MST",	"TWT","PHT","EHT")

datalist = list()
bluelist = list()

# 循环遍历合并后的 outblue
for (i in 1:length(data_list)) {
  # 获取当前的 outblue
  outblue = data_list[[i]]
  
  
  # 合并 datclean 和 datLabelled 的行
  datclean =  outblue$datalist$datclean
  datLabelled =outblue$datalist$datLabelled
  
  for (trait in traits) {
    # 将合并后的行添加到 datalist 中
    datalist$datclean[[trait]] = rbind(datalist$datclean[[trait]], datclean[[trait]])
    datalist$datLabelled[[trait]] = rbind(datalist$datLabelled[[trait]], datLabelled[[trait]])
  }
  
  
  # 合并 blue 和 h2 的行
  blue = outblue$bluelist$blue
  h2 =  outblue$bluelist$h2
  for (trait in traits) {
    # 将合并后的行添加到 bluelist 中
    bluelist$blue[[trait]] = rbind(bluelist$blue[[trait]], blue[[trait]])
    bluelist$h2[[trait]] = rbind(bluelist$h2[[trait]], h2[[trait]])
  }
  
}
idx=list(datalist=datalist,bluelist=bluelist)

table(idx$bluelist$blue$YLD14$Location)


saveRDS(idx,file = "Output/23P_datalist_bluelist.rds")




#########merge TC1 blue#######
rm(list=ls())
# 从文件夹中读取所有符合条件的 RDS 文件
file_names <- list.files("Output/", pattern = "_TC1_datalist_bluelist.rds$", full.names = TRUE)

# 创建一个空列表，用于存储读取的 RDS 文件内容
data_list <- list()

# 读取每个 RDS 文件，并将其内容添加到 data_list 列表中
for (file in file_names) {
  tmp <- readRDS(file)
  data_list[[length(data_list) + 1]] <- tmp
}

# 将 data_list 中所有内容合并到一个数据框中
traits=c("YLD14",	"MST",	"TWT","PHT","EHT")

datalist = list()
bluelist = list()

# 循环遍历合并后的 outblue
for (i in 1:length(data_list)) {
  # 获取当前的 outblue
  outblue = data_list[[i]]
  
  
  # 合并 datclean 和 datLabelled 的行
  datclean =  outblue$datalist$datclean
  datLabelled =outblue$datalist$datLabelled
  
  for (trait in traits) {
    # 将合并后的行添加到 datalist 中
    datalist$datclean[[trait]] = rbind(datalist$datclean[[trait]], datclean[[trait]])
    datalist$datLabelled[[trait]] = rbind(datalist$datLabelled[[trait]], datLabelled[[trait]])
    datalist$datclean[[trait]]$Station[datalist$datclean[[trait]]$Location=="HLJHULCC"]="CC"
    datalist$datLabelled[[trait]]$Station[datalist$datLabelled[[trait]]$Location=="HLJHULCC"]="CC"
    
     }
  
  
  # 合并 blue 和 h2 的行
  blue = outblue$bluelist$blue
  h2 =  outblue$bluelist$h2
  for (trait in traits) {
    # 将合并后的行添加到 bluelist 中
    bluelist$blue[[trait]] = rbind(bluelist$blue[[trait]], blue[[trait]])
    bluelist$h2[[trait]] = rbind(bluelist$h2[[trait]], h2[[trait]])
    bluelist$blue[[trait]]$Station[bluelist$blue[[trait]]$Location=="HLJHULCC"]="CC"
    }
  
}
idx=list(datalist=datalist,bluelist=bluelist)

table(idx$bluelist$blue$YLD14$Location)


saveRDS(idx,file = "Output/23TC1_datalist_bluelist1103.rds")



#互种问题#
tc1=readRDS("Output/23TC1_datalist_bluelist1103.rds")

xx_zq=readRDS("datafiles/seedAllocation_Sites_TC1_XX_20230418_ZQreallocated(1).rds")
jn_xx=readxl::read_xlsx("datafiles/seedAllocation_Sites_TC1_JN_20230418_XXreallocated(1).xlsx",sheet = "新乡")

name=xx_zq$章丘$Name
namejn=jn_xx$Name
blue=list()
#traits=names(idx)
traits=c("YLD14","MST","PHT","EHT")

idx=tc1$bluelist$blue
for (trait in traits) {
  zq=idx[[trait]][idx[[trait]]$Location %in% c("SDZQ1","SDZQ2"),]
  tmp=zq[!(zq$Name %in% name),]
  zq$Station="XX"

  
  
  xx=idx[[trait]][idx[[trait]]$Location %in% c("HNXX1","HNXX3"),]
  tmp.xx=xx[!(xx$Name %in% namejn),]
  xx$Station="JN"
  
  all=idx[[trait]][!(idx[[trait]]$Location %in% c("HNXX1","HNXX3","SDZQ1","SDZQ2")),]
  
  
  
  blue[[trait]]=rbind(zq,tmp,xx,tmp.xx,all)
  blue[[trait]]=unique(blue[[trait]])
}

datclean=list()
idx=tc1$datalist$datclean
for (trait in traits) {
  zq=idx[[trait]][idx[[trait]]$Location %in% c("SDZQ1","SDZQ2"),]
  tmp=zq[!(zq$Name %in% name),]
  zq$Station="XX"

  xx=idx[[trait]][idx[[trait]]$Location %in% c("HNXX1","HNXX3"),]
  tmp.xx=xx[!(xx$Name %in% namejn),]
  xx$Station="JN"
  
  all=idx[[trait]][!(idx[[trait]]$Location %in% c("HNXX1","HNXX3","SDZQ1","SDZQ2")),]

  
  datclean[[trait]]=rbind(zq,tmp,xx,tmp.xx,all)
  datclean[[trait]]=unique(datclean[[trait]])
  
  
}

vec=list(datclean=datclean,blue=blue)
saveRDS(vec,"Output/23TC1_datalist_bluelist1118.rds")


#########merge TC2 blue#######
rm(list=ls())
# 从文件夹中读取所有符合条件的 RDS 文件
file_names <- list.files("Output/", pattern = "_TC2_datalist_bluelist.rds$", full.names = TRUE)

# 创建一个空列表，用于存储读取的 RDS 文件内容
data_list <- list()

# 读取每个 RDS 文件，并将其内容添加到 data_list 列表中
for (file in file_names) {
  tmp <- readRDS(file)
  data_list[[length(data_list) + 1]] <- tmp
}

# 将 data_list 中所有内容合并到一个数据框中
traits=c("YLD14",	"MST",	"TWT","PHT","EHT")

datalist = list()
bluelist = list()

# 循环遍历合并后的 outblue
for (i in 1:length(data_list)) {
  # 获取当前的 outblue
  outblue = data_list[[i]]
  
  
  # 合并 datclean 和 datLabelled 的行
  datclean =  outblue$datalist$datclean
  datLabelled =outblue$datalist$datLabelled
  
  for (trait in traits) {
    # 将合并后的行添加到 datalist 中
    datalist$datclean[[trait]] = rbind(datalist$datclean[[trait]], datclean[[trait]])
    datalist$datLabelled[[trait]] = rbind(datalist$datLabelled[[trait]], datLabelled[[trait]])
    datalist$datclean[[trait]]$Station[datalist$datclean[[trait]]$Location=="HLJHULCC"]="CC"
    datalist$datLabelled[[trait]]$Station[datalist$datLabelled[[trait]]$Location=="HLJHULCC"]="CC"
    
    }
  
  
  # 合并 blue 和 h2 的行
  blue = outblue$bluelist$blue
  h2 =  outblue$bluelist$h2
  for (trait in traits) {
    # 将合并后的行添加到 bluelist 中
    bluelist$blue[[trait]] = rbind(bluelist$blue[[trait]], blue[[trait]])
    bluelist$h2[[trait]] = rbind(bluelist$h2[[trait]], h2[[trait]])
    bluelist$blue[[trait]]$Station[bluelist$blue[[trait]]$Location=="HLJHULCC"]="CC"
    
     }
  
}
idx=list(datalist=datalist,bluelist=bluelist)

table(idx$bluelist$blue$YLD14$Location)


saveRDS(idx,file = "Output/23TC2_datalist_bluelist1103.rds")


#互种问题#
tc2=readRDS("Output/23TC2_datalist_bluelist1103.rds")


xx_zq=readxl::read_xlsx("datafiles/seedAllocation_Sites_TC2_XX_20230404_ZQreallocated.xlsx",sheet = "章丘")
jn_xx=readxl::read_xlsx("datafiles/seedAllocation_Sites_TC2_JN_20230402_XXreallocated(2).xlsx",sheet = "新乡")


name=xx_zq$Name
namejn=jn_xx$Name
blue=list()
#traits=names(idx)
traits=c("YLD14","MST","PHT","EHT")

idx=tc2$bluelist$blue
for (trait in traits) {
  zq=idx[[trait]][idx[[trait]]$Location %in% c("SDZQ1","SDZQ2"),]
  tmp=zq[!(zq$Name %in% name),]
  zq$Station="XX"
  
  
  
  xx=idx[[trait]][idx[[trait]]$Location %in% c("HNXX1","HNXX2","HNXX3"),]
  tmp.xx=xx[!(xx$Name %in% namejn),]
  xx$Station="JN"
  
  all=idx[[trait]][!(idx[[trait]]$Location %in% c("HNXX1","HNXX2","HNXX3","SDZQ1","SDZQ2")),]
  
  
  
  blue[[trait]]=rbind(zq,tmp,xx,tmp.xx,all)
  blue[[trait]]=unique(blue[[trait]])
}

datclean=list()
idx=tc2$datalist$datclean
for (trait in traits) {
  zq=idx[[trait]][idx[[trait]]$Location %in% c("SDZQ1","SDZQ2"),]
  tmp=zq[!(zq$Name %in% name),]
  zq$Station="XX"
  
  xx=idx[[trait]][idx[[trait]]$Location %in% c("HNXX1","HNXX2","HNXX3"),]
  tmp.xx=xx[!(xx$Name %in% namejn),]
  xx$Station="JN"
  
  all=idx[[trait]][!(idx[[trait]]$Location %in% c("HNXX1","HNXX2","HNXX3","SDZQ1","SDZQ2")),]
  
  
  datclean[[trait]]=rbind(zq,tmp,xx,tmp.xx,all)
  datclean[[trait]]=unique(datclean[[trait]])
  
  
}

vec=list(datclean=datclean,blue=blue)
saveRDS(vec,"Output/23TC2_datalist_bluelist1118.rds")


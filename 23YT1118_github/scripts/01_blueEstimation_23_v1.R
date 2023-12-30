rm(list=ls())
outliersRemoval <- function(datw, trait, model) {

  
  names(datw)[names(datw)==trait] <- 'trait'
  #datw=datw[!is.na(datw$trait),]
  # 3SD方法处理异常值
  if (model == '3SD') {
    std <- sd(datw$trait, na.rm = TRUE)
    mn <- mean(datw$trait, na.rm = TRUE)
   
    datw$abnormal <- ifelse(datw$trait < mn - std * 3 | datw$trait > mn + std * 3, 2, 1)

  }
  
  #boxplot方法处理异常值
  if (model == 'boxplot') {
    bx <- boxplot(datw$trait,ylab=trait)
    min <- bx$stats[1, 1]
    max <- bx$stats[5, 1]
    c(min, max)
    datw$abnormal <- ifelse(datw$trait < min | datw$trait > max, 2, 1)
  }
  
  #其他常识性设定范围处理异常值
  
  if (trait == 'MST') datw$abnormal[datw$trait < 10] <- 3
  if (trait == 'MST') datw$abnormal[datw$trait >= 40] <- 3
  if (trait == 'YLD14') datw$abnormal[datw$trait < 200] <- 3
  if (trait == 'YLD14') datw$abnormal[datw$MST >= 40] <- 3
  if (trait == 'YLD14') datw$abnormal[datw$MST < 10] <- 3
  
  if (trait == 'PLOTWT') datw$abnormal[datw$trait  < 2] <- 3
  if (trait == 'YLD14') datw$abnormal[datw$PLOTWT  < 2] <- 3
  if (trait == 'MST') datw$abnormal[datw$PLOTWT  < 2] <- 3
  
  
  title <- paste0(trait,' Loc=',
                 unique(datw$Location)[1],' Type=',type,
                
                 ' N=',nrow(datw[!is.na(datw$trait),])
                )
  
  
  qqnorm(datw$trait,ylab = trait,col=datw$abnormal,main=title)
  qqline(datw$trait)
  #频数分布直方图
  hh<-hist(datw$trait,breaks = 25,plot=F)
  plot(hh,
       freq=T,col=rgb(0.7,0.2,0.2,alpha=0.5),border = "white",
       xlab = trait,main=title,
       ylab="Frequency",labels=T)
  
  # qqnorm(datw$trait,main=title,col=datw$abnormal)
  # qqline(datw$trait)
  
  datLabelled <- datw
  #将标记为异常值的部分处理为NA
  datw$trait[datw$abnormal %in% c(2,3)] <- NA
  #cat(nrow(outliers),'records were identified as outliers for trait',trait,'.\n')
  
  
  # title <- paste(year,trait,'Location =',
  #                unique(datw$Location)[1],'Type =',datw$Type[1],
  #                'Field =',datw$Field[1],
  #                'N =',nrow(datw[!is.na(datw$trait),]),
  #                 'Ol rm\n')
  # 
  # qqnorm(datw$trait,main=title)
  # qqline(datw$trait)
  

  names(datw)[names(datw)=="trait"] <- trait  

  names(datLabelled)[names(datLabelled)=="trait"] <- trait 
  names(datLabelled)[names(datLabelled) == "abnormal"] <- paste0(trait, "_abnormal")
  
  #将处理后数据的被处理列名称改回原名称
  
  datw <- subset(datw, select = -abnormal)
  # 将结果存储到列表中
  outlist <- list(datclean = datw, datLabelled = datLabelled)
  #names(datLabelled)
  return(outlist)
  
}
missingAdding <- function(w){
  m <- table(w$Range,w$Pass)
  row.names(m)
  names(m)
  m <- as.matrix(m)
  colnames(m)
  
  nr=nrow(m)
  nc=ncol(m)
  ndiff <- nc*nr-sum(m)
  
  if (ndiff>0){
    missing <- c()
    for (i in 1:nr){
      for (j in 1:nc){
        if (m[i,j]==0){
          rn <- rownames(m)[i]
          cn <- colnames(m)[j]
          missing <- rbind(missing, c(as.numeric(as.character(rn)),
                                      as.numeric(as.character(cn))))
        }
      }
    }
    
    zeros <- data.frame(matrix(NA,nrow(missing),ncol(w)))
    names(zeros) <- names(w)
    for (i in 1:nrow(missing)){
      zeros$Pass[i] <- missing[i,2]
      zeros$Range[i] <- missing[i,1]
      zeros$Name[i] <- paste0('DUMMY',i)
      zeros$Field[i] <- w$Field[1]
      zeros$Location[i] <- w$Location[1]
      zeros$Type[i] <- w$Type[1]
      zeros$Year[i] <- w$Year[1]
      zeros$AOA[i] <- w$AOA[1]
      zeros$Station[i] <- w$Station[1]
    }
    
    ww <- rbind(w,zeros)
    table(ww$Range,ww$Pass)
  } else {
    ww <- w
  }
  
  
  return (ww)
  
}

replicationCheck <- function(df){
  q <- data.frame(table(df$Name))
  return (nrow(q[q$Freq>1,])/nrow(q))
}

blueEstimation <- function(datdf,repct,trait){
  if (length(setdiff(c("Name", "Field", "Range", "Pass", 
                       "Location", "Type", "trait", "Year", "AOA", "Station"), 
                     names(datdf))) == 0){
    datdf$Range <- as.factor(datdf$Range)
    datdf$Pass <- as.factor(datdf$Pass)
    datdf$Name <- as.factor(datdf$Name)
    #datdf$Set <- as.factor(datdf$Set)
    datdf <- datdf[!is.na(datdf$Name),]
    
    if (nrow(datdf[!is.na(datdf$trait),])<100){
      cat('BLUE not estimable due to only one Range','\n')
      result <- aggregate(trait ~ Name, datdf, mean)
      H2 <- data.frame(Estimate=NA,SE=NA,Location=unique(datdf$Location)[1],
                       Field=unique(datdf$Field)[1],Type=datdf$Type[1],N=nrow(datdf),
                       Year = unique(datdf$Year)[1],
                       AOA = unique(datdf$AOA)[1],
                       Station = unique(datdf$Station)[1],
                       Trait=trait)
      blue <- data.frame(Name=result$Name,
                         predicted.value=result$trait, 
                         std.error=NA,
                         status='unestimable',
                         Location=unique(datdf$Location)[1],
                         Field=unique(datdf$Field)[1],
                         Type=datdf$Type[1],
                         N=nrow(datdf),
                         converge='NA',
                         Year = unique(datdf$Year)[1],
                         AOA = unique(datdf$AOA)[1],
                         Station = unique(datdf$Station)[1],
                         Trait=trait)
      print(head(blue))
      asr <- NULL
      title <- paste(trait,'BLUE is not estimable: Loc =',
                     unique(datdf$Location)[1],
                     ' Type = ',type,
                     'Field =',unique(datdf$Field)[1],
                     'N =',nrow(datdf[!is.na(datdf$trait),])
                    )
      
    } else {
      if (nrow(datdf[!is.na(datdf$trait),])<130){
        cat('asr <- asreml(trait ~ Range + Pass,random = ~idv(Name),
                              residual = ~idv(units), data = datdf)\n')
        asr <- asreml(trait ~ Range + Pass,random = ~idv(Name),
                      residual = ~idv(units), data = datdf)
        
      } else if (nrow(datdf[!is.na(datdf$trait),])<150){
        cat('asr <- asreml(trait ~ 1,random = ~idv(Name),
                              residual = ~ar1(Range):ar1(Pass), data = datdf)\n')
        asr <- asreml(trait ~ 1,random = ~idv(Name),
                      residual = ~ar1(Range):ar1(Pass), data = datdf)
        
      } else if(trait=='TWT'){# & datdf$Field[1]==2 & datdf$Type[1]=='PLate'){
        cat('asr <- asreml(trait ~ Pass,random = ~idv(Name),
                              residual = ~ar1(Range):ar1(Pass), data = datdf)\n')
        asr <- asreml(trait ~ Pass,random = ~idv(Name),
                      residual = ~ar1(Range):ar1(Pass), data = datdf)
        
      } else {
        cat('asr <- asreml(trait ~ Range + Pass,random = ~idv(Name),
                              residual = ~ar1(Range):ar1(Pass), data = datdf)\n')
        asr <- asreml(trait ~ Range + Pass,random = ~idv(Name),
                      residual = ~ar1(Range):ar1(Pass), data = datdf)
        
      }
      asr <- update(asr)
      if(!asr$converge)asr <- update(asr)
      if(!asr$converge)asr <- update(asr)
      if(!asr$converge)asr <- update(asr)
      print(wald(asr))
      summary(asr)$varcomp

      #plot(varioGram(asr))
      
      H2 <- vpredict(asr,H2~V1/(V1+V2))
      H2$Location <- unique(datdf$Location)[1]
      H2$Field <- unique(datdf$Field)[1]
      H2$Type <- unique(datdf$Type)[1]
      H2$N <- nrow(datdf)
      H2$Year = unique(datdf$Year)[1]
      H2$AOA = unique(datdf$AOA)[1]
      H2$Station = unique(datdf$Station)[1]
      H2$Trait=trait
      H2
      
      # if (repct>0.1){ #BLUE will be estimated 
      #                 #if there are replications
      #     asr <- asreml(YLD ~ Name,
      #                   residual = ~ar1(Range):ar1(Pass), 
      #                   data = datdf)
      #     if(!asr$converge)asr <- update(asr)
      #     if(!asr$converge)asr <- update(asr)
      #     if(!asr$converge)asr <- update(asr)
      #     
      # } 
      
      blue <- predict(asr, classify = "Name",pworkspace="800mb",
                      data=datdf)$pvals
      if(nrow(blue[!is.na(blue$predicted.value),])!=0){
        blue <- blue[blue$Name!='0',]
        blue$Location <- unique(datdf$Location)[1]
        blue$Field <- unique(datdf$Field)[1]
        blue$Type <- unique(datdf$Type)[1]
        blue$N <- nrow(datdf)
        blue$converge <- asr$converge
        blue$Year = 2023
        blue$AOA <- unique(datdf$AOA)[1]
        blue$Station <- unique(datdf$Station)[1]
        blue$Trait=trait
      } else {
        blue <- data.frame(asr$coefficients$random)
        fixed=asr$coefficients$fixed
        mu<-fixed[row.names(fixed)=="(Intercept)",]
        
        names(blue)[1] <- 'predicted.value'
        blue$predicted.value <- blue$predicted.value + mu
        newcol <- stringsplit(row.names(blue),'_')
        blue$Name <- newcol$X2
        blue <- blue[blue$Name!='0',]
        blue <- blue[,c(2,1)]
        blue$std.error <- 'NA'
        blue$status <- 'NA'
        blue$Location <- unique(datdf$Location)[1]
        blue$Field <- unique(datdf$Field)[1]
        blue$Type <- unique(datdf$Type)[1]
        blue$N <- nrow(datdf)
        blue$converge <- asr$converge
        blue$Year = 2023
        blue$AOA <- unique(datdf$AOA)[1]
        blue$Station <- unique(datdf$Station)[1]
        blue$Trait=trait
      }
      head(blue)
      
      title <- paste0(trait,' BLUE:Loc=',
                      unique(datdf$Location)[1],
                      ' Type=',type,
                      ' H2=',round(H2[1,1],2),' N=',nrow(datdf[!is.na(datdf$trait),])
                     )
    }

    #blue=blue[!is.na(blue$predicted.value),]
    
    boxplot(blue$predicted.value,ylab=paste0(trait,"_BLUE"),main=title)
    qqnorm(blue$predicted.value,main=title,ylab=paste0(trait,"_BLUE"))
    qqline(blue$predicted.value)
    #频数分布直方图
    hh<-hist(blue$predicted.value,breaks = 25,plot=F)
    plot(hh,
         freq=T,col=rgb(0.7,0.2,0.2,alpha=0.5),border = "white",
         xlab = paste0(trait,"_BLUE"),main=title,
         ylab="Frequency",labels=T)
    
  } else {
    cat('Not all necessary columns are included in data frame\n')
    head(datdf)
  }
  
  blue<-subset(blue,substr(blue$Name,1,5)!='DUMMY')
  
  
  
  outlist <- list(H2=H2,blue=blue,asr=asr)
  cat('ASREML estimation is successful...\n')
  
  return (outlist)
  
}


programRun<-function(dat,traits){
  
  dat <- dat[order(dat$Location,dat$Field,dat$Range,dat$Pass),]#dat$Type,
  
  dat <- dat[!is.na(dat$Range) & !is.na(dat$Pass),]
  dat <- dat[!is.na(dat$Name),]
  
  
  
  repctage <- c()
  
  h2 <- list()
  blue <- list()
  
  datclean<-list()
  
  datLabelled<-list()
  
  
  for (trait in traits){
    blue[[trait]] <- c()
    h2[[trait]] <- c()
    
    datLabelled[[trait]]<-c()
    datclean[[trait]]<-c()
  }

  for (trait in traits) {
     dat=dat[!is.na(dat[[trait]]),]
     repct <- replicationCheck(dat)
     repctage <- rbind(repctage,c(repct))
      # cat("\nLocation:",locations[loc],'Field:',
      #         fields[i],'Type',types[j],'N:',nrow(w),
      #         'Replications:',repct,'\n')
   
      
          if (nrow(dat[!is.na(dat[[trait]]),])>0){
            #table(ww$Range,ww$Pass)
            #处理缺失值与排序
            ww <- missingAdding(dat)
            ww <- ww[order(ww$Range, ww$Pass),]
            
            if (trait=='YLD14'){
              dataw <- ww[,c('Name','Field','Range','Pass','Location',
                             'Type',trait,'MST',"PLOTWT","Year","AOA","Station")]
            }else if(trait=='MST') {
              dataw <- ww[,c('Name','Field','Range','Pass','Location',
                             'Type',trait,"PLOTWT","Year","AOA","Station")]
            } else {
              dataw <- ww[,c('Name','Field','Range','Pass','Location',
                             'Type',trait,"Year","AOA","Station")]
            }

            datatreat <-outliersRemoval(datw=dataw, trait, model="boxplot")
            
            datLabelled[[trait]]<-rbind(datLabelled[[trait]],datatreat$datLabelled)
            datclean[[trait]]<-rbind(datclean[[trait]],datatreat$datclean)
            
            
            names(datatreat$datclean)[names(datatreat$datclean)==trait] <- 'trait'
            blueout <- blueEstimation(datdf=datatreat$datclean,
                                      trait=trait,repct = repct)
            
            blue[[trait]] <- rbind(blue[[trait]],blueout$blue)
            h2[[trait]] <- rbind(h2[[trait]],blueout$H2)
            
            
            #}
          } else {
            cat('No observations available.\n')
          }

  }
  datalist <- list(datclean = datclean,datLabelled=datLabelled)
  bluelist<- list(blue=blue, h2=h2)
  
  return(list(datalist=datalist, bluelist=bluelist))
  
}

merge.fun=function(combined_outblue,traits){
  datalist = list()
  bluelist = list()
  
  # 循环遍历合并后的 outblue
  for (i in 1:length(combined_outblue)) {
    # 获取当前的 outblue
      outblue = combined_outblue[[i]]
    
  
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
  return(idx)
}

drawplot=function(idx,traits){
  for (trait in traits) {
    
    trait.tmp=idx$datalist$datLabelled[[trait]]
    title <- paste0(trait,' Loc=',
                    unique(trait.tmp$Location)[1],' Type=',type,
                    
                    ' N=',nrow(trait.tmp[!is.na(trait.tmp[[trait]]),]))
    
    col=paste0(trait,"_abnormal")
    boxplot(trait.tmp[[trait]],ylab=trait)
    qqnorm(trait.tmp[[trait]],ylab = trait,col=trait.tmp[[col]],main=title)
    qqline(trait.tmp[[trait]])
    #频数分布直方图
    hh<-hist(trait.tmp[[trait]],breaks = 25,plot=F)
    plot(hh,
         freq=T,col=rgb(0.7,0.2,0.2,alpha=0.5),border = "white",
         xlab = trait,main="",
         ylab="Frequency",labels=F)
    text(hh$mids, hh$counts, labels = hh$counts, cex = 0.8)
    
    
    
    title <- paste0(trait,' BLUE:Loc=',
                    unique(idx$bluelist$blue[[trait]]$Location)[1],
                    ' Type=',type,
                    ' N=',nrow(idx$bluelist$blue[[trait]][!is.na(idx$bluelist$blue[[trait]]$predicted.value),])
    )
    
    
    #blue=blue[!is.na(blue$predicted.value),]
    
    boxplot(idx$bluelist$blue[[trait]]$predicted.value,ylab=paste0(trait,"_BLUE"))
    qqnorm(idx$bluelist$blue[[trait]]$predicted.value,main=title,ylab=paste0(trait,"_BLUE"))
    qqline(idx$bluelist$blue[[trait]]$predicted.value)
    #频数分布直方图
    hh<-hist(idx$bluelist$blue[[trait]]$predicted.value,breaks = 25,plot=F)
    plot(hh,
         freq=T,col=rgb(0.7,0.2,0.2,alpha=0.5),border = "white",
         xlab = paste0(trait,"_BLUE"),main="",
         ylab="Frequency",labels=F)
    text(hh$mids, hh$counts, labels = hh$counts, cex = 0.8)
  }
}

#统计表型值、异常值、BLUE值数量、及表型平均值、最大值、最小值
analysis=function(idx,traits){
  
count=data.frame(Trait = character(0), N_pheno = numeric(0),Pheno_mean= numeric(0),Pheno_max= numeric(0),
                 Pheno_min= numeric(0),
                 N_outlier = numeric(0),Pct_outlier=numeric(0),
                 N_BLUE= numeric(0),BLUE_mean= numeric(0),BLUE_max= numeric(0),
                 BLUE_min= numeric(0),
                 Type=character(0),Location=character(0),
                 AOA=character(0),Station=character(0))

 for (trait in traits) {
  trait.tmp=idx$datalist$datLabelled[[trait]]
  colname=paste0(trait,"_abnormal")
  trait.tmp=trait.tmp[!is.na(trait.tmp[[colname]]),]
  
  trait.tmp.out =trait.tmp[trait.tmp[[colname]] %in% c(2,3),]
  
  blue.tmp=idx$bluelist$blue[[trait]]
  
  blue.tmp=blue.tmp[!is.na(blue.tmp$predicted.value),]
  
  count=rbind(count,data.frame(Trait=trait,
                               N_pheno=nrow(trait.tmp),
                               #N_discard=N_discard,
                               #Pct_discard=round(N_discard/nrow(trait.tmp)*100,2),
                               Pheno_mean=round(mean(trait.tmp[[trait]]),3),
                               Pheno_max= round(max(trait.tmp[[trait]]),3),
                               Pheno_min= round(min(trait.tmp[[trait]]),3),
                               N_outlier=nrow(trait.tmp.out),
                               Pct_outlier=round(nrow(trait.tmp.out)/nrow(trait.tmp)*100,3),
                               N_BLUE=nrow(blue.tmp),
                               BLUE_mean=round(mean(blue.tmp$predicted.value),3),
                               BLUE_max= round(max(blue.tmp$predicted.value),3),
                               BLUE_min= round(min(blue.tmp$predicted.value),3),
                               Type=unique(blue.tmp$Type)[1],
                               Location=unique(blue.tmp$Location)[1],
                               AOA=unique(blue.tmp$AOA)[1],
                               Station=unique(blue.tmp$Station)[1]))
  
 }
return(count)
}

#library(readxl)
library(writexl)
library(biosim)
library(tidyverse)
library(reshape2)
library(stringi) 
library(asreml)

#asreml.license.activate()

# phenotype -------------------------------------------------------
#getwd()
#setwd("E:/P1 project/")
#phenos <- readRDS('datafiles/19_22_P_AOA_Station.rds')
# col_types <- c(rep("numeric",2),rep("text",3),rep("numeric",3),rep("text",7),
#                rep("numeric",6),rep("text",2),rep("numeric",29),rep("text",2))

# col_types <- c(rep("text",2),rep("numeric",2),rep("text",1),rep("numeric",3),rep("text",7),
#                rep("numeric",6),rep("text",2),rep("numeric",39),rep("text",2))


# col_types <- c(rep("text",6),rep("numeric",5),rep("text",6),
#                rep("numeric",2),rep("text",1),
#                rep("numeric",6),rep("text",2),rep("numeric",58),rep("text",1))

#phenos1 <- readxl::read_excel('datafiles/Yield Trial Master Catalog-10.30.xlsx', sheet = 1,col_types = col_types)
phenos1=readRDS("datafiles/23phenos_1118.rds")

# names(phenos1)[names(phenos1)=="pass"]="Pass"
# names(phenos1)[names(phenos1)=="BookPrj"]="Station"
# names(phenos1)[names(phenos1)=="TierNo"]="Field"
# names(phenos1)[names(phenos1)=="TrialType"]="Type"
# names(phenos1)[names(phenos1)=="Location"]="Location_C"
# names(phenos1)[names(phenos1)=="BookName"]="Location"

phenos.tmp=phenos1[phenos1$Name %in% c("INF","PRF","INF-4","discard"),]
phenos.tmp$Name=paste0("PRF",1:nrow(phenos.tmp))
phenos2=phenos1[!(phenos1$Name %in% c("INF","PRF","INF-4","discard")),]

phenos=rbind(phenos2,phenos.tmp)
phenos=unique(phenos)
par(mfrow=c(2,3),mar = c(5, 4, 4, 2) + 0.1,cex.main = 0.8) 


##############################CC and NW Station################################
# JLCLTC2 -----------------------------------------------------------------

cau.data <- phenos %>% dplyr::filter(Location %in% "JLCLTC2")
cau.data$AOA="MMSP"
cau.data$Year=2023


target_cols <- c( "HarvNote")#"苗期综合评价",,"收获前综合评价"

# 搜索包含字符的行，并返回逻辑向量
has_char <- apply(cau.data[target_cols], 1, function(x) any(grepl("[[:alnum:][:punct:]]+", x)))


discard=cau.data[has_char,]
discard[,c("HarvNote")]

ph_filtered=cau.data[!has_char,]
N_discard=nrow(ph_filtered[ph_filtered$PlotDiscarded %in% c("Yes"), ])

ph_filtered <- ph_filtered[!(ph_filtered$PlotDiscarded %in% c("Yes")), ]
names(ph_filtered)
ph=ph_filtered

traits=c("YLD14",	"MST","PLOTWT","TWT")

combined_outblue = list()

# 循环遍历不同的 type
for (i in 1:2) {
  # 根据type选择特定的数据子集
  if (i == 1) {
    
    type="Field1_2_TC2"
    subset_ph = ph[ph$Pass %in% 1:48 & ph$Range %in% 2:43, ]
  } else if (i == 2) {
    type="Field3_TC2"
    subset_ph = ph[ph$Pass %in% 49:60 & ph$Range %in% 2:37, ]#5号地不进行空间分析，数据加入到总的中。导致总数据BLUE值不符合正态分布
  } 
  # 运行程序并生成 outblue
  outblue.tmp <-  programRun(dat = subset_ph, traits)
  
  # 将当前的 outblue 添加到 combined_outblue 列表中
  combined_outblue = c(combined_outblue, list(outblue.tmp))
  
}

idx=merge.fun(combined_outblue,traits)

count=analysis(idx,traits)
#count$Type="P"
write.table(count,"Output/count.csv",sep=",",row.names=FALSE,col.names=TRUE,append=TRUE)

#####绘图

pdf("Output/JLCL_plot_TC2.pdf")

par(mfrow=c(2,3),mar = c(5, 4, 4, 2) + 0.1,cex.main = 0.8) 
type="TC2"
drawplot(idx,traits)
dev.off()

outblue_name <- paste0("Output/23JLCL_TC2_datalist_bluelist.rds")

saveRDS(idx,file = outblue_name)





# JLDH_TC2--------------------------------------------------------------------

cau.data <- phenos %>% dplyr::filter(Location %in% "JLDHTC2")
cau.data$AOA="MMSP"
cau.data$Year=2023


target_cols <- c( "HarvNote")#"苗期综合评价",,"收获前综合评价"

# 搜索包含字符的行，并返回逻辑向量
has_char <- apply(cau.data[target_cols], 1, function(x) any(grepl("[[:alnum:][:punct:]]+", x)))


discard=cau.data[has_char,]
discard[,c("HarvNote")]

ph_filtered=cau.data[!has_char,]
N_discard=nrow(ph_filtered[ph_filtered$PlotDiscarded %in% c("Yes"), ])

ph_filtered <- ph_filtered[!(ph_filtered$PlotDiscarded %in% c("Yes")), ]
names(ph_filtered)
ph=ph_filtered

traits=c("YLD14",	"MST","PLOTWT","TWT")

combined_outblue = list()

# 循环遍历不同的 type
for (i in 1:2) {
  # 根据type选择特定的数据子集
  if (i == 1) {
    
    type="Field1_2_TC2"
    subset_ph = ph[ph$Pass %in% 1:48 & ph$Range %in% 2:39, ]
  } else if (i == 2) {
    type="Field3_TC2"
    subset_ph = ph[ph$Pass %in% 49:72 & ph$Range %in% 2:28, ]
  } 
  # 运行程序并生成 outblue
  outblue.tmp <-  programRun(dat = subset_ph, traits)
  
  # 将当前的 outblue 添加到 combined_outblue 列表中
  combined_outblue = c(combined_outblue, list(outblue.tmp))
  
}

idx=merge.fun(combined_outblue,traits)

count=analysis(idx,traits)
count$Type="TC2"
write.table(count,"Output/count.csv",sep=",",row.names=FALSE,col.names=FALSE,append=TRUE)

#####绘图

pdf("Output/JLDH_plot_TC2.pdf")

par(mfrow=c(2,3),mar = c(5, 4, 4, 2) + 0.1,cex.main = 0.8) 
type="TC2"
drawplot(idx,traits)
dev.off()

outblue_name <- paste0("Output/23JLDH_TC2_datalist_bluelist.rds")

saveRDS(idx,file = outblue_name)







# JLGZL--------------------------------------------------------------------
cau.data <- phenos %>% dplyr::filter(Location %in% "JLGZL")
cau.data$AOA="MMSP"
cau.data$Year=2023


target_cols <- c("HarvNote")#"收获前综合评价", ,"苗期综合评价"

# 搜索包含字符的行，并返回逻辑向量
has_char <- apply(cau.data[target_cols], 1, function(x) any(grepl("[[:alnum:][:punct:]]+", x)))


discard=cau.data[has_char,]
discard[,c("HarvNote")]

ph_filtered=cau.data[!has_char,]
N_discard=nrow(ph_filtered[ph_filtered$PlotDiscarded %in% c("Yes"), ])
N_discard
ph_filtered <- ph_filtered[!(ph_filtered$PlotDiscarded %in% c("Yes")), ]
names(ph_filtered)
ph=ph_filtered

####JLGZL_P####
traits=c("YLD14",	"MST","PLOTWT","PHT","EHT","TWT")#	

combined_outblue <- list()

# 循环遍历不同的 type
for (i in 1:6) {
  # 根据type选择特定的数据子集
  if (i == 1) {
    
    type="Field1_P3_P4"
    subset_ph = ph[ph$Pass %in% 1:23 & ph$Range %in% 2:14, ]
    
  } else if (i == 2){
    type="Field1_P1"
    subset_ph = ph[ph$Pass %in% 1:24 & ph$Range %in% 15:66, ]
    
  } else if (i == 3){
    type="Field2_P1"
    subset_ph = ph[ph$Pass %in% 25:48 & ph$Range %in% 2:66, ]
    
  } else if (i == 4){
    type="Field3_P1"
    subset_ph = ph[ph$Pass %in% 49:72 & ph$Range %in% 2:34, ]
    
  } else if (i == 5){
    type="Field6_P1"
    subset_ph = ph[ph$Pass %in% 121:144 & ph$Range %in% 58:66, ]
    
  } else if (i == 6){
    type="Field7_8_P1"
    subset_ph = ph[ph$Pass %in% 145:180 & ph$Range %in% 2:66, ]
    
  } 
  # 运行程序并生成 outblue
  outblue.tmp = programRun(dat = subset_ph, traits)
  
  # 将当前的 outblue 添加到 combined_outblue 列表中
  combined_outblue = append(combined_outblue, list(outblue.tmp))
  
}

idx=merge.fun(combined_outblue,traits)

count=analysis(idx,traits)
count$Type="P"
write.table(count,"Output/count.csv",sep=",",row.names=FALSE,col.names=FALSE,append=TRUE)
#####绘图

pdf("Output/JLGZL_plot_P.pdf")

par(mfrow=c(2,3),mar = c(5, 4, 4, 2) + 0.1,cex.main = 0.8) 
type="P"
drawplot(idx,traits)
dev.off()

outblue_name <- paste0("Output/23JLGZL_P_datalist_bluelist.rds")

saveRDS(idx,file = outblue_name )

####JLGZL_TC1####
traits=c("YLD14",	"MST","PLOTWT","PHT","EHT","TWT")#	

combined_outblue <- list()

# 循环遍历不同的 type
for (i in 1:3) {
  # 根据type选择特定的数据子集
  if (i == 1) {
    
    type="Field3_TC1"
    subset_ph = ph[ph$Pass %in% 49:72 & ph$Range %in% 35:66, ]
    #subset_ph =subset_ph[substr(subset_ph$Type,1,1)=="P",]
  } else if (i==2){
    type="Field4_TC1"
    subset_ph = ph[ph$Pass %in% 73:96 & ph$Range %in% 2:66, ]
  } else if (i==3){
    type="Field5_TC1"
    subset_ph = ph[ph$Pass %in% 97:120 & ph$Range %in% 2:49, ]
  } 
  # 运行程序并生成 outblue
  outblue.tmp = programRun(dat = subset_ph, traits)
  
  # 将当前的 outblue 添加到 combined_outblue 列表中
  combined_outblue = append(combined_outblue, list(outblue.tmp))
  
}

idx=merge.fun(combined_outblue,traits)

count=analysis(idx,traits)
#count$Type="P"
write.table(count,"Output/count.csv",sep=",",row.names=FALSE,col.names=FALSE,append=TRUE)
#####绘图

pdf("Output/JLGZL_plot_TC1.pdf")

par(mfrow=c(2,3),mar = c(5, 4, 4, 2) + 0.1,cex.main = 0.8) 
type="TC1"
drawplot(idx,traits)
dev.off()

outblue_name <- paste0("Output/23JLGZL_TC1_datalist_bluelist.rds")

saveRDS(idx,file = outblue_name )


####JLGZL_TC2####
traits=c("YLD14",	"MST","PLOTWT","PHT","EHT","TWT")#	

combined_outblue <- list()

# 循环遍历不同的 type
for (i in 1:2) {
  # 根据type选择特定的数据子集
  if (i == 1) {
    
    type="Field5_TC2"
    subset_ph = ph[ph$Pass %in% 97:120 & ph$Range %in% 50:66, ]
    #subset_ph =subset_ph[substr(subset_ph$Type,1,1)=="P",]
  } else if (i==2){
    type="Field6_TC2"
    subset_ph = ph[ph$Pass %in% 121:144 & ph$Range %in% 2:57, ]
  } 
  # 运行程序并生成 outblue
  outblue.tmp = programRun(dat = subset_ph, traits)
  
  # 将当前的 outblue 添加到 combined_outblue 列表中
  combined_outblue = append(combined_outblue, list(outblue.tmp))
  
}

idx=merge.fun(combined_outblue,traits)

count=analysis(idx,traits)
#count$Type="P"
write.table(count,"Output/count.csv",sep=",",row.names=FALSE,col.names=FALSE,append=TRUE)
#####绘图

pdf("Output/JLGZL_plot_TC2.pdf")

par(mfrow=c(2,3),mar = c(5, 4, 4, 2) + 0.1,cex.main = 0.8) 
type="TC2"
drawplot(idx,traits)
dev.off()

outblue_name <- paste0("Output/23JLGZL_TC2_datalist_bluelist.rds")

saveRDS(idx,file = outblue_name )



# JLYS_P--------------------------------------------------------------------
cau.data <- phenos %>% dplyr::filter(Location %in% "JLYS")
cau.data$AOA="MMSP"
cau.data$Year=2023


target_cols <- c("HarvNote")#"收获前综合评价", ,"苗期综合评价"

# 搜索包含字符的行，并返回逻辑向量
has_char <- apply(cau.data[target_cols], 1, function(x) any(grepl("[[:alnum:][:punct:]]+", x)))


discard=cau.data[has_char,]
discard[,c("HarvNote")]

ph_filtered=cau.data[!has_char,]
N_discard=nrow(ph_filtered[ph_filtered$PlotDiscarded %in% c("Yes"), ])
N_discard
ph_filtered <- ph_filtered[!(ph_filtered$PlotDiscarded %in% c("Yes")), ]
names(ph_filtered)
ph=ph_filtered


traits=c("YLD14",	"MST","PLOTWT","PHT","EHT","TWT")#	

combined_outblue <- list()

# 循环遍历不同的 type
for (i in 1:3) {
  # 根据type选择特定的数据子集
  if (i == 1) {
    
    type="Field1_P3_P4"
    subset_ph = ph[ph$Pass %in% 1:23 & ph$Range %in% 2:16, ]
    
  } else if (i == 2){
    type="Field1_P1"
    subset_ph = ph[ph$Pass %in% 1:24 & ph$Range %in% 17:50, ]
    
  } else if (i == 3){
    type="Field2_3_4_5_6_P1"
    subset_ph = ph[ph$Pass %in% 25:80 & ph$Range %in% 2:51, ]
    
  } 
  # 运行程序并生成 outblue
  outblue.tmp = programRun(dat = subset_ph, traits)
  
  # 将当前的 outblue 添加到 combined_outblue 列表中
  combined_outblue = append(combined_outblue, list(outblue.tmp))
  
}

idx=merge.fun(combined_outblue,traits)

count=analysis(idx,traits)
count$Type="P"
write.table(count,"Output/count.csv",sep=",",row.names=FALSE,col.names=FALSE,append=TRUE)
#####绘图

pdf("Output/JLYS_plot_P.pdf")

par(mfrow=c(2,3),mar = c(5, 4, 4, 2) + 0.1,cex.main = 0.8) 
type="P"
drawplot(idx,traits)
dev.off()

outblue_name <- paste0("Output/23JLYS_P_datalist_bluelist.rds")

saveRDS(idx,file = outblue_name )






# JLNA --------------------------------------------------------------------


cau.data <- phenos %>% dplyr::filter(Location %in% "JLNA")
cau.data$AOA="MMSP"
cau.data$Year=2023


target_cols <- c("HarvNote")#"收获前综合评价", ,"苗期综合评价"

# 搜索包含字符的行，并返回逻辑向量
has_char <- apply(cau.data[target_cols], 1, function(x) any(grepl("[[:alnum:][:punct:]]+", x)))


discard=cau.data[has_char,]
discard[,c("HarvNote")]

ph_filtered=cau.data[!has_char,]
N_discard=nrow(ph_filtered[ph_filtered$PlotDiscarded %in% c("Yes"), ])
N_discard
ph_filtered <- ph_filtered[!(ph_filtered$PlotDiscarded %in% c("Yes")), ]
names(ph_filtered)
ph=ph_filtered


####JLNA_P####
traits=c("YLD14",	"MST",	"TWT","PLOTWT")#,"PHT","EHT"

combined_outblue <- list()

# 循环遍历不同的 type
for (i in 1:6) {
  # 根据type选择特定的数据子集
  if (i == 1) {
    
    type="Field1_P3_P4"
    subset_ph = ph[ph$Pass %in% 1:23 & ph$Range %in% 5:17, ]
  } else if (i == 2) {
    type="Field1_P1"
    subset_ph = ph[ph$Pass %in% 1:24 & ph$Range %in% 18:66, ]
  } else if (i == 3) {
    type="Field2_P1"
    subset_ph = ph[ph$Pass %in% 25:48 & ph$Range %in% 5:66, ]
  } else if (i == 4) {
    type="Field3_P1"
    subset_ph = ph[ph$Pass %in% 49:72 & ph$Range %in% 5:43, ]
  } else if (i==5){
    type="Field7_P1" 
    subset_ph <- ph[ph$Pass %in% 145:168 & ph$Range %in% 41:63, ]
  } else if (i==6){
    type="Field8_9_P1" 
    subset_ph <- ph[ph$Pass %in% 169:200 & ph$Range %in% 2:64, ]
  }
  
  # 运行程序并生成 outblue
  outblue.tmp = programRun(dat = subset_ph, traits)
  
  # 将当前的 outblue 添加到 combined_outblue 列表中
  combined_outblue = append(combined_outblue, list(outblue.tmp))
  
}

idx=merge.fun(combined_outblue,traits)

count=analysis(idx,traits)
count$Type="P"
write.table(count,"Output/count.csv",sep=",",row.names=FALSE,col.names=FALSE,append=TRUE)
#####绘图

pdf("Output/JLNA_plot_P.pdf")

par(mfrow=c(2,3),mar = c(5, 4, 4, 2) + 0.1,cex.main = 0.8) 
type="P"
drawplot(idx,traits)
dev.off()

outblue_name <- paste0("Output/23JLNA_P_datalist_bluelist.rds")

saveRDS(idx,file = outblue_name )

####JLNA_TC2####
traits=c("YLD14",	"MST",	"TWT","PLOTWT","PHT","EHT")#

combined_outblue <- list()

# 循环遍历不同的 type
for (i in 1:3) {
  # 根据type选择特定的数据子集
  if (i == 1) {
    
    type="Field3_TC2"
    subset_ph = ph[ph$Pass %in% 49:72 & ph$Range %in% 44:66, ]
  } else if (i == 2) {
    type="Field4_TC2"
    subset_ph = ph[ph$Pass %in% 73:96 & ph$Range %in% 5:66, ]
  } else if (i == 3) {
    type="Field5_TC2"
    subset_ph = ph[ph$Pass %in% 97:120 & ph$Range %in% 5:21, ]
  } 
  
  # 运行程序并生成 outblue
  outblue.tmp = programRun(dat = subset_ph, traits)
  
  # 将当前的 outblue 添加到 combined_outblue 列表中
  combined_outblue = append(combined_outblue, list(outblue.tmp))
  
}

idx=merge.fun(combined_outblue,traits)

count=analysis(idx,traits)
#count$Type="P"
write.table(count,"Output/count.csv",sep=",",row.names=FALSE,col.names=FALSE,append=TRUE)
#####绘图

pdf("Output/JLNA_plot_TC2.pdf")

par(mfrow=c(2,3),mar = c(5, 4, 4, 2) + 0.1,cex.main = 0.8) 
type="TC2"
drawplot(idx,traits)
dev.off()

outblue_name <- paste0("Output/23JLNA_TC2_datalist_bluelist.rds")

saveRDS(idx,file = outblue_name )


####JLNA_TC1####
traits=c("YLD14",	"MST",	"TWT","PLOTWT")#,"PHT","EHT"

combined_outblue <- list()

# 循环遍历不同的 type
for (i in 1:3) {
  # 根据type选择特定的数据子集
  if (i == 1) {
    
    type="Field5_TC1"
    subset_ph = ph[ph$Pass %in% 97:120 & ph$Range %in% 22:66, ]
  } else if (i == 2) {
    type="Field6_TC1"
    subset_ph = ph[ph$Pass %in% 121:144 & ph$Range %in% 2:63, ]
  } else if (i == 3) {
    type="Field7_TC1"
    subset_ph = ph[ph$Pass %in% 145:168 & ph$Range %in% 2:40, ]
  } 
  
  # 运行程序并生成 outblue
  outblue.tmp = programRun(dat = subset_ph, traits)
  
  # 将当前的 outblue 添加到 combined_outblue 列表中
  combined_outblue = append(combined_outblue, list(outblue.tmp))
  
}

idx=merge.fun(combined_outblue,traits)

count=analysis(idx,traits)
#count$Type="P"
write.table(count,"Output/count.csv",sep=",",row.names=FALSE,col.names=FALSE,append=TRUE)
#####绘图

pdf("Output/JLNA_plot_TC1.pdf")

par(mfrow=c(2,3),mar = c(5, 4, 4, 2) + 0.1,cex.main = 0.8) 
type="TC1"
drawplot(idx,traits)
dev.off()

outblue_name <- paste0("Output/23JLNA_TC1_datalist_bluelist.rds")

saveRDS(idx,file = outblue_name )







# JLSY_P --------------------------------------------------------------------

cau.data <- phenos %>% dplyr::filter(Location %in% "JLSY")
cau.data$AOA="MMSP"
cau.data$Year=2023


target_cols <- c( "HarvNote")#"苗期综合评价",,"收获前综合评价"

# 搜索包含字符的行，并返回逻辑向量
has_char <- apply(cau.data[target_cols], 1, function(x) any(grepl("[[:alnum:][:punct:]]+", x)))


discard=cau.data[has_char,]
discard[,c("HarvNote")]
#discard[,c(27,28,87)]
ph_filtered=cau.data[!has_char,]
N_discard=nrow(ph_filtered[ph_filtered$PlotDiscarded %in% c("Yes"), ])

ph_filtered <- ph_filtered[!(ph_filtered$PlotDiscarded %in% c("Yes")), ]
names(ph_filtered)
ph=ph_filtered


traits=c("YLD14",	"MST","PLOTWT","TWT","PHT","EHT")#

combined_outblue = list()

# 循环遍历不同的 type
for (i in 1:4) {
  # 根据type选择特定的数据子集
  if (i == 1) {
    
    type="Field1_P3_P4"
    subset_ph = ph[ph$Pass %in% 1:23 & ph$Range %in% 2:16, ]
  } else if(i==2){
    type="Field1_P1"
    subset_ph = ph[ph$Pass %in% 1:24 & ph$Range %in% 17:59, ]
  } else if(i==3){
    type="Field2_P1"
    subset_ph = ph[ph$Pass %in% 25:48 & ph$Range %in% 2:59, ]
  } else if(i==4){
    type="Field3_P1"
    subset_ph = ph[ph$Pass %in% 49:72 & ph$Range %in% 2:50, ]
  }
  # 运行程序并生成 outblue
  outblue.tmp <-  programRun(dat = subset_ph, traits)
  
  # 将当前的 outblue 添加到 combined_outblue 列表中
  combined_outblue = c(combined_outblue, list(outblue.tmp))
  
}

idx=merge.fun(combined_outblue,traits)

count=analysis(idx,traits)
count$Type="P"
write.table(count,"Output/count.csv",sep=",",row.names=FALSE,col.names=FALSE,append=TRUE)

#####绘图

pdf("Output/JLSY_plot_P.pdf")

par(mfrow=c(2,3),mar = c(5, 4, 4, 2) + 0.1,cex.main = 0.8) 
type="P"
drawplot(idx,traits)
dev.off()

outblue_name <- paste0("Output/23JLSY_P_datalist_bluelist.rds")

saveRDS(idx,file = outblue_name)









# LNCT2-------------------------------------------------------------------
cau.data <- phenos %>% dplyr::filter(Location %in% "LNCT2")
cau.data$AOA="LMSP"
cau.data$Year=2023


target_cols <- c( "HarvNote")#"收获前综合评价""苗期综合评价",

# 搜索包含字符的行，并返回逻辑向量
has_char <- apply(cau.data[target_cols], 1, function(x) any(grepl("[[:alnum:][:punct:]]+", x)))


discard=cau.data[has_char,]
discard[,c("HarvNote")]

ph_filtered=cau.data[!has_char,]
N_discard=nrow(ph_filtered[ph_filtered$PlotDiscarded %in% c("Yes"), ])

ph_filtered = ph_filtered[!(ph_filtered$PlotDiscarded %in% c("Yes")), ]
names(ph_filtered)
ph=ph_filtered


#### LNCT2 _P####

traits=c("YLD14",	"MST",	"TWT","PLOTWT")

combined_outblue = list()

# 循环遍历不同的 type
for (i in 1:2) {
  # 根据type选择特定的数据子集
  if (i == 1) {
    
    type="3_P1"
    subset_ph = ph[ph$Field == 3 & ph$Range %in% 2:74, ]
  } else if (i == 2) {
    type="4_P1"
    subset_ph = ph[ph$Field %in% 4 & ph$Range %in% 2:67, ]
  } 
  
  # 运行程序并生成 outblue
  outblue.tmp <-  programRun(dat = subset_ph, traits)
  
  # 将当前的 outblue 添加到 combined_outblue 列表中
  combined_outblue = c(combined_outblue, list(outblue.tmp))
  
}

idx=merge.fun(combined_outblue,traits)

count=analysis(idx,traits)
count$Type="P"
write.table(count,"Output/count.csv",sep=",",row.names=FALSE,col.names=FALSE,append=TRUE)
#####绘图

pdf("Output/LNCT2_plot_P.pdf")

par(mfrow=c(2,3),mar = c(5, 4, 4, 2) + 0.1,cex.main = 0.8) 
type="P"
drawplot(idx,traits)
dev.off()

outblue_name <- paste0("Output/23LNCT2_P_datalist_bluelist.rds")

saveRDS(idx,file = outblue_name)



#### LNCT2_TC1 ####


traits=c("YLD14",	"MST",	"TWT","PLOTWT","PHT","EHT")

combined_outblue <- list()

# 循环遍历不同的 type
for (i in 1) {
  # 根据type选择特定的数据子集
  if (i == 1) {
    
    type="TC1"
    subset_ph = ph[ph$Field %in% c(1,2) & ph$Range %in% 2:74, ]
  } 
  
  # 运行程序并生成 outblue
  outblue.tmp = programRun(dat = subset_ph, traits)
  
  # 将当前的 outblue 添加到 combined_outblue 列表中
  combined_outblue <- append(combined_outblue, list(outblue.tmp))
  
}

idx=merge.fun(combined_outblue,traits)

count=analysis(idx,traits)
write.table(count,"Output/count.csv",sep=",",row.names=FALSE,col.names=FALSE,append=TRUE)
#####绘图

pdf("Output/LNCT2_plot_TC1.pdf")

par(mfrow=c(2,3),mar = c(5, 4, 4, 2) + 0.1,cex.main = 0.8) 
type="TC1"
drawplot(idx,traits)
dev.off()

outblue_name <- paste0("Output/23LNCT2_TC1_datalist_bluelist.rds")

saveRDS(idx,file = outblue_name)





# LNCT1 -------------------------------------------------------------------
cau.data <- phenos %>% dplyr::filter(Location %in% "LNCT1")
cau.data$AOA="LMSP"
cau.data$Year=2023


target_cols <- c( "HarvNote")#"收获前综合评价""苗期综合评价",

# 搜索包含字符的行，并返回逻辑向量
has_char <- apply(cau.data[target_cols], 1, function(x) any(grepl("[[:alnum:][:punct:]]+", x)))


discard=cau.data[has_char,]
discard[,c("HarvNote")]

ph_filtered=cau.data[!has_char,]
N_discard=nrow(ph_filtered[ph_filtered$PlotDiscarded %in% c("Yes"), ])

ph_filtered = ph_filtered[!(ph_filtered$PlotDiscarded %in% c("Yes")), ]
names(ph_filtered)
ph=ph_filtered

####LNCT1_p####

traits=c("YLD14",	"MST",	"TWT","PLOTWT")


combined_outblue = list()

# 循环遍历不同的 type
for (i in 1:4) {
  # 根据type选择特定的数据子集
  if (i == 1) {
    
    type="1_P3_4"
    subset_ph = ph[ph$Field == 1 & ph$Range %in% 2:14, ]
  } else if (i == 2) {
    type="1_P1"
    subset_ph = ph[ph$Field == 1 & ph$Range %in% 15:74, ]
  } else if (i==3){
    type="2_P1"
    subset_ph = ph[ph$Field == 2 & ph$Range %in% 2:74, ]
    
  } else if (i==4){
    type="3_P1"
    subset_ph = ph[ph$Field == 3 & ph$Range %in% 2:18, ]
    
  }
  
  # 运行程序并生成 outblue
  outblue.tmp <-  programRun(dat = subset_ph, traits)
  
  # 将当前的 outblue 添加到 combined_outblue 列表中
  combined_outblue = c(combined_outblue, list(outblue.tmp))
  
}

idx=merge.fun(combined_outblue,traits)

count=analysis(idx,traits)
count$Type="P"
write.table(count,"Output/count.csv",sep=",",row.names=FALSE,col.names=FALSE,append=TRUE)
#####绘图

pdf("Output/LNCT1_plot_P.pdf")

par(mfrow=c(2,3),mar = c(5, 4, 4, 2) + 0.1,cex.main = 0.8) 
type="P"
drawplot(idx,traits)
dev.off()

outblue_name <- paste0("Output/23LNCT1_P_datalist_bluelist.rds")

saveRDS(idx,file = outblue_name)
#write.csv(idx$datalist$datLabelled$TWT,"Output/23P_LNCT1_datLabelled_TWT.csv")

####LNCT1_TC2####
traits=c("YLD14",	"MST",	"TWT","PLOTWT","PHT","EHT")

combined_outblue = list()

# 循环遍历不同的 type
for (i in 1:3) {
  # 根据type选择特定的数据子集
  if (i == 1) {
    
    type="3_TC2"
    subset_ph = ph[ph$Field == 3 & ph$Range %in% 19:74, ]
  } else if (i == 2) {
    type="4_TC2"
    subset_ph = ph[ph$Field == 4 & ph$Range %in% 2:74, ]
  } else if (i==3){
    type="5_TC2"
    subset_ph = ph[ph$Field == 5 & ph$Range %in% 2:63, ]
    
  } 
  
  # 运行程序并生成 outblue
  outblue.tmp <-  programRun(dat = subset_ph, traits)
  
  # 将当前的 outblue 添加到 combined_outblue 列表中
  combined_outblue = c(combined_outblue, list(outblue.tmp))
  
}

idx=merge.fun(combined_outblue,traits)

count=analysis(idx,traits)
count$Type="TC2"
write.table(count,"Output/count.csv",sep=",",row.names=FALSE,col.names=FALSE,append=TRUE)
#####绘图

pdf("Output/LNCT1_plot_TC2.pdf")

par(mfrow=c(2,3),mar = c(5, 4, 4, 2) + 0.1,cex.main = 0.8) 
type="TC2"
drawplot(idx,traits)
dev.off()

outblue_name <- paste0("Output/23LNCT1_TC2_datalist_bluelist.rds")

saveRDS(idx,file = outblue_name)







# NMKEQ -------------------------------------------------------------------

cau.data <- phenos %>% dplyr::filter(Location %in% "NMKEQ")
cau.data$AOA="LMSP"
cau.data$Year=2023


target_cols <- c("苗期综合评价", "HarvNote")#"收获前综合评价"

# 搜索包含字符的行，并返回逻辑向量
has_char <- apply(cau.data[target_cols], 1, function(x) any(grepl("[[:alnum:][:punct:]]+", x)))


discard=cau.data[has_char,]
discard[,c("HarvNote")]

ph_filtered=cau.data[!has_char,]
N_discard=nrow(ph_filtered[ph_filtered$PlotDiscarded %in% c("Yes"), ])+nrow(discard)
ph_filtered <- ph_filtered[!(ph_filtered$PlotDiscarded %in% c("Yes")), ]
names(ph_filtered)

traits=c("YLD14",	"MST",	"TWT","PLOTWT","PHT","EHT")


ph=ph_filtered

combined_outblue = list()

# 循环遍历不同的 type
for (i in 1:4) {
  # 根据type选择特定的数据子集
  if (i == 1) {
    
    type="1_P3_4"
    subset_ph = ph[ph$Pass %in% 1:23 & ph$Range %in% 2:15, ]
  } else if (i == 2) {
    type="1_P1"
    subset_ph = ph[ph$Pass %in% 1:24 & ph$Range %in% 16:34, ]
  } else if (i==3){
    type="P1"
    subset_ph = ph[ph$Pass %in% 37:120 & ph$Range %in% 2:34, ]
    
  } else if (i==4){
    type="1_P1"
    subset_ph = ph[ph$Pass %in% 25:36 & ph$Range %in% 2:34, ]
    
  } 
  
  # 运行程序并生成 outblue
  outblue.tmp <-  programRun(dat = subset_ph, traits)
  
  # 将当前的 outblue 添加到 combined_outblue 列表中
  combined_outblue = c(combined_outblue, list(outblue.tmp))
  
}

idx=merge.fun(combined_outblue,traits)

count=analysis(idx,traits)
count$Type="P"
write.table(count,"Output/count.csv",sep=",",row.names=FALSE,col.names=FALSE,append=TRUE)
#####绘图

pdf("Output/NMKEQ_plot_P.pdf")

par(mfrow=c(2,3),mar = c(5, 4, 4, 2) + 0.1,cex.main = 0.8) 
type="P"
drawplot(idx,traits)
dev.off()

outblue_name <- paste0("Output/23NMKEQ_P_datalist_bluelist.rds")

saveRDS(idx,file = outblue_name)



# NMNM-----------------------------------------------------------------

cau.data <- phenos %>% dplyr::filter(Location %in% "NMNM")
cau.data$AOA="MMSP"
cau.data$Year=2023


target_cols <- c( "HarvNote")#"苗期综合评价",,"收获前综合评价"

# 搜索包含字符的行，并返回逻辑向量
has_char <- apply(cau.data[target_cols], 1, function(x) any(grepl("[[:alnum:][:punct:]]+", x)))


discard=cau.data[has_char,]
discard[,c("HarvNote")]
#discard[,c(27,28,87)]
ph_filtered=cau.data[!has_char,]
N_discard=nrow(ph_filtered[ph_filtered$PlotDiscarded %in% c("Yes"), ])

ph_filtered <- ph_filtered[!(ph_filtered$PlotDiscarded %in% c("Yes")), ]
names(ph_filtered)
ph=ph_filtered

traits=c("YLD14",	"MST","PLOTWT","TWT")#,"PHT","EHT"

combined_outblue = list()

# 循环遍历不同的 type
for (i in 1:4) {
  # 根据type选择特定的数据子集
  if (i == 1) {
    
    type="Field1_P3_P4"
    subset_ph = ph[ph$Pass %in% 1:23 & ph$Range %in% 2:15, ]
  } else if (i == 2) {
    type="Field1_P1"
    subset_ph = ph[ph$Pass %in% 1:24 & ph$Range %in% 16:39, ]#
  } else if (i == 3) {
    type="Field2_3_4_P1"
    subset_ph = ph[ph$Pass %in% 25:96 & ph$Range %in% 2:39, ]#
  } else if (i == 4) {
    type="Field5_P1"
    subset_ph = ph[ph$Pass %in% 97:120 & ph$Range %in% 2:14, ]#
  } 
  # 运行程序并生成 outblue
  outblue.tmp <-  programRun(dat = subset_ph, traits)
  
  # 将当前的 outblue 添加到 combined_outblue 列表中
  combined_outblue = c(combined_outblue, list(outblue.tmp))
  
}

idx=merge.fun(combined_outblue,traits)

count=analysis(idx,traits)
count$Type="P"
write.table(count,"Output/count.csv",sep=",",row.names=FALSE,col.names=FALSE,append=TRUE)

#####绘图

pdf("Output/NMNM_plot_P.pdf")

par(mfrow=c(2,3),mar = c(5, 4, 4, 2) + 0.1,cex.main = 0.8) 
type="P"
drawplot(idx,traits)
dev.off()

outblue_name <- paste0("Output/23NMNM_P_datalist_bluelist.rds")

saveRDS(idx,file = outblue_name)





# SXQX_TC2--------------------------------------------------------------------

cau.data <- phenos %>% dplyr::filter(Location %in% "SXQX")
cau.data$AOA="LMSP"
cau.data$Year=2023


target_cols <- c( "HarvNote")#"苗期综合评价",,"收获前综合评价"

# 搜索包含字符的行，并返回逻辑向量
has_char <- apply(cau.data[target_cols], 1, function(x) any(grepl("[[:alnum:][:punct:]]+", x)))


discard=cau.data[has_char,]
discard[,c("HarvNote")]

ph_filtered=cau.data[!has_char,]
N_discard=nrow(ph_filtered[ph_filtered$PlotDiscarded %in% c("Yes"), ])

ph_filtered <- ph_filtered[!(ph_filtered$PlotDiscarded %in% c("Yes")), ]
names(ph_filtered)
ph=ph_filtered

traits=c("YLD14",	"MST","PLOTWT","PHT","EHT")

combined_outblue = list()

type="Field1_2_3_TC2"
subset_ph = ph[ph$Pass %in% 1:72 & ph$Range %in% 2:35, ]

# 运行程序并生成 outblue
outblue.tmp <-  programRun(dat = subset_ph, traits)

# 将当前的 outblue 添加到 combined_outblue 列表中
combined_outblue = c(combined_outblue, list(outblue.tmp))



idx=merge.fun(combined_outblue,traits)

count=analysis(idx,traits)
write.table(count,"Output/count.csv",sep=",",row.names=FALSE,col.names=FALSE,append=TRUE)

#####绘图

pdf("Output/SXQX_plot_TC2.pdf")

par(mfrow=c(2,3),mar = c(5, 4, 4, 2) + 0.1,cex.main = 0.8) 
type="TC2"
drawplot(idx,traits)
dev.off()

outblue_name <- paste0("Output/23SXQX_TC2_datalist_bluelist.rds")

saveRDS(idx,file = outblue_name)


# SXCZ--------------------------------------------------------------------

cau.data <- phenos %>% dplyr::filter(Location %in% "SXCZ")
cau.data$AOA="MMSP"
cau.data$Year=2023


target_cols <- c( "HarvNote")#"苗期综合评价",,"收获前综合评价"

# 搜索包含字符的行，并返回逻辑向量
has_char <- apply(cau.data[target_cols], 1, function(x) any(grepl("[[:alnum:][:punct:]]+", x)))


discard=cau.data[has_char,]
discard[,c("HarvNote")]

ph_filtered=cau.data[!has_char,]
N_discard=nrow(ph_filtered[ph_filtered$PlotDiscarded %in% c("Yes"), ])

ph_filtered <- ph_filtered[!(ph_filtered$PlotDiscarded %in% c("Yes")), ]
names(ph_filtered)
ph=ph_filtered

traits=c("YLD14",	"MST","PLOTWT","TWT","PHT","EHT")

combined_outblue = list()

# 循环遍历不同的 type
for (i in 1:4) {
  # 根据type选择特定的数据子集
  if (i == 1) {
    
    type="Field1_P3_4"
    subset_ph = ph[ph$Pass %in% 1:23 & ph$Range %in% 2:16, ]
  } else if (i == 2) {
    type="Field1_P1"
    subset_ph = ph[ph$Pass %in% 1:24 & ph$Range %in% 17:35, ]
  } else if (i == 3) {
    type="Field2_3_4_P1"
    subset_ph = ph[ph$Pass %in% 25:96 & ph$Range %in% 2:35, ]
  } else if (i == 4) {
    type="Field5_P1"
    subset_ph = ph[ph$Pass %in% 97:120 & ph$Range %in% 2:33, ]
  } 
  # 运行程序并生成 outblue
  outblue.tmp <-  programRun(dat = subset_ph, traits)
  
  # 将当前的 outblue 添加到 combined_outblue 列表中
  combined_outblue = c(combined_outblue, list(outblue.tmp))
  
}

idx=merge.fun(combined_outblue,traits)

count=analysis(idx,traits)
count$Type="P"
write.table(count,"Output/count.csv",sep=",",row.names=FALSE,col.names=FALSE,append=TRUE)

#####绘图

pdf("Output/SXCZ_plot_P.pdf")

par(mfrow=c(2,3),mar = c(5, 4, 4, 2) + 0.1,cex.main = 0.8) 
type="P"
drawplot(idx,traits)
dev.off()

outblue_name <- paste0("Output/23SXCZ_P_datalist_bluelist.rds")

saveRDS(idx,file = outblue_name)




# SXTG_TC2 --------------------------------------------------------------------

cau.data <- phenos %>% dplyr::filter(Location %in% "SXTG")
cau.data$AOA="LMSP"
cau.data$Year=2023


target_cols <- c( "HarvNote")#"苗期综合评价",,"收获前综合评价"

# 搜索包含字符的行，并返回逻辑向量
has_char <- apply(cau.data[target_cols], 1, function(x) any(grepl("[[:alnum:][:punct:]]+", x)))


discard=cau.data[has_char,]
discard[,c("HarvNote")]

ph_filtered=cau.data[!has_char,]
N_discard=nrow(ph_filtered[ph_filtered$PlotDiscarded %in% c("Yes"), ])

ph_filtered <- ph_filtered[!(ph_filtered$PlotDiscarded %in% c("Yes")), ]
names(ph_filtered)
ph=ph_filtered

traits=c("YLD14",	"MST","PLOTWT","PHT","EHT")

combined_outblue = list()

# 循环遍历不同的 type
for (i in 1:2) {
  # 根据type选择特定的数据子集
  if (i == 1) {
    
    type="Field1_2_3_4_TC2"
    subset_ph = ph[ph$Pass %in% 1:80 & ph$Range %in% 2:31, ]
  } else if (i == 2) {
    type="Field5_TC2"
    subset_ph = ph[ph$Pass %in% 81:88 & ph$Range %in% 2:10, ]#5号地不进行空间分析，数据加入到总的中。导致总数据BLUE值不符合正态分布
  } 
  # 运行程序并生成 outblue
  outblue.tmp <-  programRun(dat = subset_ph, traits)
  
  # 将当前的 outblue 添加到 combined_outblue 列表中
  combined_outblue = c(combined_outblue, list(outblue.tmp))
  
}

idx=merge.fun(combined_outblue,traits)

count=analysis(idx,traits)
#count$Type="P"
write.table(count,"Output/count.csv",sep=",",row.names=FALSE,col.names=FALSE,append=TRUE)

#####绘图

pdf("Output/SXTG_plot_TC2.pdf")

par(mfrow=c(2,3),mar = c(5, 4, 4, 2) + 0.1,cex.main = 0.8) 
type="TC2"
drawplot(idx,traits)
dev.off()

outblue_name <- paste0("Output/23SXTG_TC2_datalist_bluelist.rds")

saveRDS(idx,file = outblue_name)










# SXGP --------------------------------------------------------------------

cau.data <- phenos %>% dplyr::filter(Location %in% "SXGP")
cau.data$AOA="MMSP"
cau.data$Year=2023


target_cols <- c()#"收获前综合评价""苗期综合评价", "HarvNote"

# 搜索包含字符的行，并返回逻辑向量
has_char <- apply(cau.data[target_cols], 1, function(x) any(grepl("[[:alnum:][:punct:]]+", x)))


discard=cau.data[has_char,]
discard[,c("HarvNote")]

ph_filtered=cau.data[!has_char,]
N_discard=nrow(ph_filtered[ph_filtered$PlotDiscarded %in% c("Yes"), ])

ph_filtered <- ph_filtered[!(ph_filtered$PlotDiscarded %in% c("Yes")), ]
names(ph_filtered)
ph=ph_filtered

####SXGP_P####
traits=c("YLD14",	"MST","PLOTWT","PHT","EHT")

combined_outblue = list()

# 循环遍历不同的 type
for (i in 1:4) {
  # 根据type选择特定的数据子集
  if (i == 1) {
    
    type="1_P3_4"
    subset_ph = ph[ph$Field %in% 1 & ph$Range %in% 2:14, ]
  } else if (i == 2) {
    type="1_P1"
    subset_ph = ph[ph$Field %in% 1 & ph$Range %in% 15:53, ]
  } else if (i==3){
    type="2_3_P1"
    subset_ph = ph[ph$Field %in% 2:3 & ph$Range %in% 2:53, ]
    
  } else if (i==4){
    type="4_P1"
    subset_ph = ph[ph$Field %in% 4 & ph$Range %in% 2:9, ]
    
  } 
  
  # 运行程序并生成 outblue
  outblue.tmp <-  programRun(dat = subset_ph, traits)
  
  # 将当前的 outblue 添加到 combined_outblue 列表中
  combined_outblue = c(combined_outblue, list(outblue.tmp))
  
}

idx=merge.fun(combined_outblue,traits)

count=analysis(idx,traits)
count$Type="P"
write.table(count,"Output/count.csv",sep=",",row.names=FALSE,col.names=FALSE,append=TRUE)
#####绘图

pdf("Output/SXGP_plot_P.pdf")

par(mfrow=c(2,3),mar = c(5, 4, 4, 2) + 0.1,cex.main = 0.8) 
type="P"
drawplot(idx,traits)
dev.off()

outblue_name <- paste0("Output/23SXGP_P_datalist_bluelist.rds")

saveRDS(idx,file = outblue_name)







####SXGP_TC2####
traits=c("YLD14",	"MST","PLOTWT","PHT","EHT")

combined_outblue = list()

# 循环遍历不同的 type
for (i in 1:3) {
  # 根据type选择特定的数据子集
  if (i == 1) {
    
    type="Field4_TC2"
    subset_ph = ph[ph$Pass %in% 73:96 & ph$Range %in% 10:53, ]
  } else if (i == 2) {
    type="Field5_TC2"
    subset_ph = ph[ph$Pass %in% 97:120 & ph$Range %in% 2:53, ]
  } else if (i==3){
    type="Field6_TC2"
    subset_ph = ph[ph$Pass %in% 121:144 & ph$Range %in% 2:7, ]
    
  } 
  
  # 运行程序并生成 outblue
  outblue.tmp <-  programRun(dat = subset_ph, traits)
  
  # 将当前的 outblue 添加到 combined_outblue 列表中
  combined_outblue = c(combined_outblue, list(outblue.tmp))
  
}

idx=merge.fun(combined_outblue,traits)

count=analysis(idx,traits)
write.table(count,"Output/count.csv",sep=",",row.names=FALSE,col.names=FALSE,append=TRUE)
#####绘图

pdf("Output/SXGP_plot_TC2.pdf")

par(mfrow=c(2,3),mar = c(5, 4, 4, 2) + 0.1,cex.main = 0.8) 
type="TC2"
drawplot(idx,traits)
dev.off()

outblue_name <- paste0("Output/23SXGP_TC2_datalist_bluelist.rds")

saveRDS(idx,file = outblue_name)


####SXGP_TC1####
traits=c("YLD14",	"MST","PLOTWT","PHT","EHT")

combined_outblue = list()

# 循环遍历不同的 type
for (i in 1:3) {
  # 根据type选择特定的数据子集
  if (i == 1) {
    
    type="Field6_TC1"
    subset_ph = ph[ph$Pass %in% 121:144 & ph$Range %in% 8:53, ]
  } else if (i == 2) {
    type="Field7_8_TC1"
    subset_ph = ph[ph$Pass %in% 145:180 & ph$Range %in% 2:53, ]
  } else if (i==3){
    type="Field9_TC1"
    subset_ph = ph[ph$Pass %in% 181:188 & ph$Range %in% 2:49, ]
    
  } 
  
  # 运行程序并生成 outblue
  outblue.tmp <-  programRun(dat = subset_ph, traits)
  
  # 将当前的 outblue 添加到 combined_outblue 列表中
  combined_outblue = c(combined_outblue, list(outblue.tmp))
  
}

idx=merge.fun(combined_outblue,traits)

count=analysis(idx,traits)
write.table(count,"Output/count.csv",sep=",",row.names=FALSE,col.names=FALSE,append=TRUE)
#####绘图

pdf("Output/SXGP_plot_TC1.pdf")

par(mfrow=c(2,3),mar = c(5, 4, 4, 2) + 0.1,cex.main = 0.8) 
type="TC1"
drawplot(idx,traits)
dev.off()

outblue_name <- paste0("Output/23SXGP_TC1_datalist_bluelist.rds")

saveRDS(idx,file = outblue_name)





# SXXZA-----------------------------------------------------------------

cau.data <- phenos %>% dplyr::filter(Location %in% "SXXZA")
cau.data$AOA="MMSP"
cau.data$Year=2023


target_cols <- c( "HarvNote")#"苗期综合评价",,"收获前综合评价"

# 搜索包含字符的行，并返回逻辑向量
has_char <- apply(cau.data[target_cols], 1, function(x) any(grepl("[[:alnum:][:punct:]]+", x)))


discard=cau.data[has_char,]
discard[,c("HarvNote")]
#discard[,c(27,28,87)]
ph_filtered=cau.data[!has_char,]
N_discard=nrow(ph_filtered[ph_filtered$PlotDiscarded %in% c("Yes"), ])

ph_filtered <- ph_filtered[!(ph_filtered$PlotDiscarded %in% c("Yes")), ]
names(ph_filtered)
ph=ph_filtered

traits=c("YLD14",	"MST","PLOTWT","PHT","EHT")#"TWT",

combined_outblue = list()

# 循环遍历不同的 type
for (i in 1:2) {
  # 根据type选择特定的数据子集
  if (i == 1) {
    
    type="Field1_P3_P4"
    subset_ph = ph[ph$Pass %in% 1:23 & ph$Range %in% 2:15, ]
  } else if (i == 2) {
    type="Field2_3_P1"
    subset_ph = ph[ph$Pass %in% 25:72 & ph$Range %in% 2:51, ]#
  } 
  # 运行程序并生成 outblue
  outblue.tmp <-  programRun(dat = subset_ph, traits)
  
  # 将当前的 outblue 添加到 combined_outblue 列表中
  combined_outblue = c(combined_outblue, list(outblue.tmp))
  
}

idx=merge.fun(combined_outblue,traits)

count=analysis(idx,traits)
count$Type="P"
write.table(count,"Output/count.csv",sep=",",row.names=FALSE,col.names=FALSE,append=TRUE)

#####绘图

pdf("Output/SXXZA_plot_P.pdf")

par(mfrow=c(2,3),mar = c(5, 4, 4, 2) + 0.1,cex.main = 0.8) 
type="P"
drawplot(idx,traits)
dev.off()

outblue_name <- paste0("Output/23SXXZA_P_datalist_bluelist.rds")

saveRDS(idx,file = outblue_name)


# SXXZB-----------------------------------------------------------------

cau.data <- phenos %>% dplyr::filter(Location %in% "SXXZB")
cau.data$AOA="MMSP"
cau.data$Year=2023


target_cols <- c( "HarvNote")#"苗期综合评价",,"收获前综合评价"

# 搜索包含字符的行，并返回逻辑向量
has_char <- apply(cau.data[target_cols], 1, function(x) any(grepl("[[:alnum:][:punct:]]+", x)))


discard=cau.data[has_char,]
discard[,c("HarvNote")]

#discard[,c(27,28,87)]
ph_filtered=cau.data[!has_char,]
N_discard=nrow(ph_filtered[ph_filtered$PlotDiscarded %in% c("Yes"), ])

ph_filtered <- ph_filtered[!(ph_filtered$PlotDiscarded %in% c("Yes")), ]
names(ph_filtered)
ph=ph_filtered

####SXXZB_P####
traits=c("YLD14",	"MST","PLOTWT","PHT","EHT")#"TWT",

combined_outblue = list()

# 循环遍历不同的 type
for (i in 1) {
  # 根据type选择特定的数据子集
  if (i == 1) {
    
    type="Field1_P1"
    subset_ph = ph[ph$Pass %in% 1:24 & ph$Range %in% 2:16, ]
  } 
  # 运行程序并生成 outblue
  outblue.tmp <-  programRun(dat = subset_ph, traits)
  
  # 将当前的 outblue 添加到 combined_outblue 列表中
  combined_outblue = c(combined_outblue, list(outblue.tmp))
  
}

idx=merge.fun(combined_outblue,traits)

count=analysis(idx,traits)
count$Type="P"
write.table(count,"Output/count.csv",sep=",",row.names=FALSE,col.names=FALSE,append=TRUE)

#####绘图

pdf("Output/SXXZB_plot_P.pdf")

par(mfrow=c(2,3),mar = c(5, 4, 4, 2) + 0.1,cex.main = 0.8) 
type="P"
drawplot(idx,traits)
dev.off()

outblue_name <- paste0("Output/23SXXZB_P_datalist_bluelist.rds")

saveRDS(idx,file = outblue_name)


####SXXZB_TC2####
traits=c("YLD14",	"MST","PLOTWT","PHT","EHT")#"TWT",

combined_outblue = list()

# 循环遍历不同的 type
for (i in 1:2) {
  # 根据type选择特定的数据子集
  if (i == 1) {
    
    type="Field1_TC2"
    subset_ph = ph[ph$Pass %in% 1:24 & ph$Range %in% 17:26, ]
  } else if(i==2) {
    type="Field2_3_4_5_TC2"
    subset_ph = ph[ph$Pass %in% 25:108 & ph$Range %in% 1:26, ]
  }
  # 运行程序并生成 outblue
  outblue.tmp <-  programRun(dat = subset_ph, traits)
  
  # 将当前的 outblue 添加到 combined_outblue 列表中
  combined_outblue = c(combined_outblue, list(outblue.tmp))
  
}

idx=merge.fun(combined_outblue,traits)

count=analysis(idx,traits)
#count$Type="TC2"
write.table(count,"Output/count.csv",sep=",",row.names=FALSE,col.names=FALSE,append=TRUE)

#####绘图

pdf("Output/SXXZB_plot_TC2.pdf")

par(mfrow=c(2,3),mar = c(5, 4, 4, 2) + 0.1,cex.main = 0.8) 
type="TC2"
drawplot(idx,traits)
dev.off()

outblue_name <- paste0("Output/23SXXZB_TC2_datalist_bluelist.rds")

saveRDS(idx,file = outblue_name)






# SXYC --------------------------------------------------------------------


cau.data <- phenos %>% dplyr::filter(Location %in% "SXYC")
cau.data$AOA="LMSP"
cau.data$Year=2023


target_cols <- c("HarvNote")#"收获前综合评价", ,"苗期综合评价"

# 搜索包含字符的行，并返回逻辑向量
has_char <- apply(cau.data[target_cols], 1, function(x) any(grepl("[[:alnum:][:punct:]]+", x)))


discard=cau.data[has_char,]
discard[,c("HarvNote")]

ph_filtered=cau.data[!has_char,]
N_discard=nrow(ph_filtered[ph_filtered$PlotDiscarded %in% c("Yes"), ])
N_discard
ph_filtered <- ph_filtered[!(ph_filtered$PlotDiscarded %in% c("Yes")), ]
names(ph_filtered)
ph=ph_filtered


####SXYC_P####
traits=c("YLD14",	"MST","PLOTWT","PHT","EHT")#	"TWT",

combined_outblue <- list()

# 循环遍历不同的 type
for (i in 1:4) {
  # 根据type选择特定的数据子集
  if (i == 1) {
    
    type="Field1_P3_P4"
    subset_ph = ph[ph$Pass %in% 1:23 & ph$Range %in% 2:14, ]
  } else if (i == 2) {
    type="Field1_P1"
    subset_ph = ph[ph$Pass %in% 1:24 & ph$Range %in% 15:77, ]
  } else if (i == 3) {
    type="Field2_P1"
    subset_ph = ph[ph$Pass %in% 25:48 & ph$Range %in% 2:77, ]
  } else if (i == 4) {
    type="Field3_P1"
    subset_ph = ph[ph$Pass %in% 49:72 & ph$Range %in% 2:12, ]
  } 
  
  # 运行程序并生成 outblue
  outblue.tmp = programRun(dat = subset_ph, traits)
  
  # 将当前的 outblue 添加到 combined_outblue 列表中
  combined_outblue = append(combined_outblue, list(outblue.tmp))
  
}

idx=merge.fun(combined_outblue,traits)

count=analysis(idx,traits)
count$Type="P"
write.table(count,"Output/count.csv",sep=",",row.names=FALSE,col.names=FALSE,append=TRUE)
#####绘图

pdf("Output/SXYC_plot_P.pdf")

par(mfrow=c(2,3),mar = c(5, 4, 4, 2) + 0.1,cex.main = 0.8) 
type="P"
drawplot(idx,traits)
dev.off()

outblue_name <- paste0("Output/23SXYC_P_datalist_bluelist.rds")

saveRDS(idx,file = outblue_name )


####SXYC_TC2####
traits=c("YLD14",	"MST","PLOTWT","PHT","EHT")#	"TWT",

combined_outblue <- list()

# 循环遍历不同的 type
for (i in 1) {
  # 根据type选择特定的数据子集
  if (i == 1) {
    
    type="Field5_TC2"
    subset_ph = ph[ph$Pass %in% 97:120 & ph$Range %in% 2:75, ]
  } 
  
  # 运行程序并生成 outblue
  outblue.tmp = programRun(dat = subset_ph, traits)
  
  # 将当前的 outblue 添加到 combined_outblue 列表中
  combined_outblue = append(combined_outblue, list(outblue.tmp))
  
}

idx=merge.fun(combined_outblue,traits)

count=analysis(idx,traits)
#count$Type="P"
write.table(count,"Output/count.csv",sep=",",row.names=FALSE,col.names=FALSE,append=TRUE)
#####绘图

pdf("Output/SXYC_plot_TC2.pdf")

par(mfrow=c(2,3),mar = c(5, 4, 4, 2) + 0.1,cex.main = 0.8) 
type="TC2"
drawplot(idx,traits)
dev.off()

outblue_name <- paste0("Output/23SXYC_TC2_datalist_bluelist.rds")

saveRDS(idx,file = outblue_name )


####SXYC_TC1####
traits=c("YLD14",	"MST","PLOTWT","PHT","EHT")#	"TWT",

combined_outblue <- list()

# 循环遍历不同的 type
for (i in 1:2) {
  # 根据type选择特定的数据子集
  if (i == 1) {
    
    type="Field3_TC1"
    subset_ph = ph[ph$Pass %in% 49:72 & ph$Range %in% 13:77, ]
  } else if (i == 2) {
    type="Field4_TC1"
    subset_ph = ph[ph$Pass %in% 73:96 & ph$Range %in% 2:78, ]
  } 
  
  # 运行程序并生成 outblue
  outblue.tmp = programRun(dat = subset_ph, traits)
  
  # 将当前的 outblue 添加到 combined_outblue 列表中
  combined_outblue = append(combined_outblue, list(outblue.tmp))
  
}

idx=merge.fun(combined_outblue,traits)

count=analysis(idx,traits)
#count$Type="P"
write.table(count,"Output/count.csv",sep=",",row.names=FALSE,col.names=FALSE,append=TRUE)
#####绘图

pdf("Output/SXYC_plot_TC1.pdf")

par(mfrow=c(2,3),mar = c(5, 4, 4, 2) + 0.1,cex.main = 0.8) 
type="TC1"
drawplot(idx,traits)
dev.off()

outblue_name <- paste0("Output/23SXYC_TC1_datalist_bluelist.rds")

saveRDS(idx,file = outblue_name )





#################HRB Station#############
# HLJJD --------------------------------------------------------------------

cau.data <- phenos %>% dplyr::filter(Location %in% "HLJJD")
cau.data$AOA="EMSP"
cau.data$Year=2023


target_cols <- c( "HarvNote")#"苗期综合评价",,"收获前综合评价"

# 搜索包含字符的行，并返回逻辑向量
has_char <- apply(cau.data[target_cols], 1, function(x) any(grepl("[[:alnum:][:punct:]]+", x)))


discard=cau.data[has_char,]
discard[,c("HarvNote")]
#discard[,c(27,28,87)]
ph_filtered=cau.data[!has_char,]
N_discard=nrow(ph_filtered[ph_filtered$PlotDiscarded %in% c("Yes"), ])

ph_filtered <- ph_filtered[!(ph_filtered$PlotDiscarded %in% c("Yes")), ]
names(ph_filtered)
ph=ph_filtered

####HLJJD_P####
traits=c("YLD14",	"MST","PLOTWT","PHT","EHT")#"TWT",

combined_outblue = list()
par(mfrow=c(2,3),mar = c(5, 4, 4, 2) + 0.1,cex.main = 0.8)
# 循环遍历不同的 type
for (i in 1) {
  # 根据type选择特定的数据子集
  if (i == 1) {
    
    type="Field1_P1"
    subset_ph = ph[ph$Pass %in% 1:24 & ph$Range %in% 2:7 &ph$Type %in% c("P1","P3"), ]
  } 
  # 运行程序并生成 outblue
  outblue.tmp <-  programRun(dat = subset_ph, traits)
  
  # 将当前的 outblue 添加到 combined_outblue 列表中
  combined_outblue = c(combined_outblue, list(outblue.tmp))
  
}

idx=merge.fun(combined_outblue,traits)

count=analysis(idx,traits)
count$Type="P"
write.table(count,"Output/count.csv",sep=",",row.names=FALSE,col.names=FALSE,append=TRUE)

#####绘图

pdf("Output/HLJJD_plot_P.pdf")

par(mfrow=c(2,3),mar = c(5, 4, 4, 2) + 0.1,cex.main = 0.8) 
type="P"
drawplot(idx,traits)
dev.off()

outblue_name <- paste0("Output/23HLJJD_P_datalist_bluelist.rds")

saveRDS(idx,file = outblue_name)


####HLJJD_TC2####
traits=c("YLD14",	"MST","PLOTWT","PHT","EHT")#"TWT",

combined_outblue = list()
par(mfrow=c(2,3),mar = c(5, 4, 4, 2) + 0.1,cex.main = 0.8)
# 循环遍历不同的 type
for (i in 1) {
  # 根据type选择特定的数据子集
  if (i == 1) {
    
    type="Field1_TC2"
    subset_ph = ph[ph$Pass %in% 1:24 & ph$Range %in% 2:7 &ph$Type %in% c("TC2"), ]
  } 
  # 运行程序并生成 outblue
  outblue.tmp <-  programRun(dat = subset_ph, traits)
  
  # 将当前的 outblue 添加到 combined_outblue 列表中
  combined_outblue = c(combined_outblue, list(outblue.tmp))
  
}

idx=merge.fun(combined_outblue,traits)

count=analysis(idx,traits)
#count$Type="P"
write.table(count,"Output/count.csv",sep=",",row.names=FALSE,col.names=FALSE,append=TRUE)

#####绘图

pdf("Output/HLJJD_plot_TC2.pdf")

par(mfrow=c(2,3),mar = c(5, 4, 4, 2) + 0.1,cex.main = 0.8) 
type="TC2"
drawplot(idx,traits)
dev.off()

outblue_name <- paste0("Output/23HLJJD_TC2_datalist_bluelist.rds")

saveRDS(idx,file = outblue_name)



# HLJJMS --------------------------------------------------------------------

cau.data <- phenos %>% dplyr::filter(Location %in% "HLJJMS")
cau.data$AOA="EMSP"
cau.data$Year=2023


target_cols <- c( "HarvNote")#"苗期综合评价",,"收获前综合评价"

# 搜索包含字符的行，并返回逻辑向量
has_char <- apply(cau.data[target_cols], 1, function(x) any(grepl("[[:alnum:][:punct:]]+", x)))


discard=cau.data[has_char,]
discard[,c("HarvNote")]
#discard[,c(27,28,87)]
ph_filtered=cau.data[!has_char,]
N_discard=nrow(ph_filtered[ph_filtered$PlotDiscarded %in% c("Yes"), ])

ph_filtered <- ph_filtered[!(ph_filtered$PlotDiscarded %in% c("Yes")), ]
names(ph_filtered)
ph=ph_filtered

####HLJJMS_P####
traits=c("YLD14",	"MST","PLOTWT","PHT","EHT")#"TWT",

combined_outblue = list()
par(mfrow=c(2,3),mar = c(5, 4, 4, 2) + 0.1,cex.main = 0.8)
# 循环遍历不同的 type
for (i in 1) {
  # 根据type选择特定的数据子集
  if (i == 1) {
    
    type="Field1_P1"
    subset_ph = ph[ph$Pass %in% 1:24 & ph$Range %in% 2:7 &ph$Type %in% c("P1","P3"), ]
    subset_ph=subset_ph[!substr(subset_ph$Name,1,5)=="23TC2",]
  } 
  # 运行程序并生成 outblue
  outblue.tmp <-  programRun(dat = subset_ph, traits)
  
  # 将当前的 outblue 添加到 combined_outblue 列表中
  combined_outblue = c(combined_outblue, list(outblue.tmp))
  
}

idx=merge.fun(combined_outblue,traits)

count=analysis(idx,traits)
count$Type="P"
write.table(count,"Output/count.csv",sep=",",row.names=FALSE,col.names=FALSE,append=TRUE)

#####绘图

pdf("Output/HLJJMS_plot_P.pdf")

par(mfrow=c(2,3),mar = c(5, 4, 4, 2) + 0.1,cex.main = 0.8) 
type="P"
drawplot(idx,traits)
dev.off()

outblue_name <- paste0("Output/23HLJJMS_P_datalist_bluelist.rds")

saveRDS(idx,file = outblue_name)


####HLJJD_TC2####
traits=c("YLD14",	"MST","PLOTWT","PHT","EHT")#"TWT",

combined_outblue = list()
par(mfrow=c(2,3),mar = c(5, 4, 4, 2) + 0.1,cex.main = 0.8)
# 循环遍历不同的 type
for (i in 1) {
  # 根据type选择特定的数据子集
  if (i == 1) {
    
    type="Field1_TC2"
    subset_ph = ph[ph$Pass %in% 1:24 & ph$Range %in% 2:7 &ph$Type %in% c("TC2"), ]
  } 
  # 运行程序并生成 outblue
  outblue.tmp <-  programRun(dat = subset_ph, traits)
  
  # 将当前的 outblue 添加到 combined_outblue 列表中
  combined_outblue = c(combined_outblue, list(outblue.tmp))
  
}

idx=merge.fun(combined_outblue,traits)

count=analysis(idx,traits)
#count$Type="P"
write.table(count,"Output/count.csv",sep=",",row.names=FALSE,col.names=FALSE,append=TRUE)

#####绘图

pdf("Output/HLJJMS_plot_TC2.pdf")

par(mfrow=c(2,3),mar = c(5, 4, 4, 2) + 0.1,cex.main = 0.8) 
type="TC2"
drawplot(idx,traits)
dev.off()

outblue_name <- paste0("Output/23HLJJMS_TC2_datalist_bluelist.rds")

saveRDS(idx,file = outblue_name)


# HLJMS --------------------------------------------------------------------

cau.data <- phenos %>% dplyr::filter(Location %in% "HLJMS")
cau.data$AOA="EMSP"
cau.data$Year=2023


target_cols <- c( "HarvNote")#"苗期综合评价",,"收获前综合评价"

# 搜索包含字符的行，并返回逻辑向量
has_char <- apply(cau.data[target_cols], 1, function(x) any(grepl("[[:alnum:][:punct:]]+", x)))


discard=cau.data[has_char,]
discard[,c("HarvNote")]
#discard[,c(27,28,87)]
ph_filtered=cau.data[!has_char,]
N_discard=nrow(ph_filtered[ph_filtered$PlotDiscarded %in% c("Yes"), ])

ph_filtered <- ph_filtered[!(ph_filtered$PlotDiscarded %in% c("Yes")), ]
names(ph_filtered)
ph=ph_filtered
ph$Type[substr(ph$Name,1,5)=="23TC2"]="TC2"

####HLJJMS_P####
traits=c("YLD14",	"MST","PLOTWT","PHT","EHT")#"TWT",

combined_outblue = list()
par(mfrow=c(2,3),mar = c(5, 4, 4, 2) + 0.1,cex.main = 0.8)
# 循环遍历不同的 type
for (i in 1) {
  # 根据type选择特定的数据子集
  if (i == 1) {
    
    type="Field1_P1"
    subset_ph = ph[ph$Pass %in% 1:24 & ph$Range %in% 2:7 &ph$Type %in% c("P1","P3"), ]
    
  } 
  # 运行程序并生成 outblue
  outblue.tmp <-  programRun(dat = subset_ph, traits)
  
  # 将当前的 outblue 添加到 combined_outblue 列表中
  combined_outblue = c(combined_outblue, list(outblue.tmp))
  
}

idx=merge.fun(combined_outblue,traits)

count=analysis(idx,traits)
count$Type="P"
write.table(count,"Output/count.csv",sep=",",row.names=FALSE,col.names=FALSE,append=TRUE)

#####绘图

pdf("Output/HLJMS_plot_P.pdf")

par(mfrow=c(2,3),mar = c(5, 4, 4, 2) + 0.1,cex.main = 0.8) 
type="P"
drawplot(idx,traits)
dev.off()

outblue_name <- paste0("Output/23HLJMS_P_datalist_bluelist.rds")

saveRDS(idx,file = outblue_name)


####HLJJD_TC2####
traits=c("YLD14",	"MST","PLOTWT","PHT","EHT")#"TWT",

combined_outblue = list()
par(mfrow=c(2,3),mar = c(5, 4, 4, 2) + 0.1,cex.main = 0.8)
# 循环遍历不同的 type
for (i in 1) {
  # 根据type选择特定的数据子集
  if (i == 1) {
    
    type="Field1_TC2"
    subset_ph = ph[ph$Pass %in% 1:24 & ph$Range %in% 2:7 &ph$Type %in% c("TC2"), ]
  } 
  # 运行程序并生成 outblue
  outblue.tmp <-  programRun(dat = subset_ph, traits)
  
  # 将当前的 outblue 添加到 combined_outblue 列表中
  combined_outblue = c(combined_outblue, list(outblue.tmp))
  
}

idx=merge.fun(combined_outblue,traits)

count=analysis(idx,traits)
#count$Type="P"
write.table(count,"Output/count.csv",sep=",",row.names=FALSE,col.names=FALSE,append=TRUE)

#####绘图

pdf("Output/HLJMS_plot_TC2.pdf")

par(mfrow=c(2,3),mar = c(5, 4, 4, 2) + 0.1,cex.main = 0.8) 
type="TC2"
drawplot(idx,traits)
dev.off()

outblue_name <- paste0("Output/23HLJMS_TC2_datalist_bluelist.rds")

saveRDS(idx,file = outblue_name)


# HLJWK --------------------------------------------------------------------

cau.data <- phenos %>% dplyr::filter(Location %in% "HLJWK")
cau.data$AOA="EMSP"
cau.data$Year=2023


target_cols <- c( "HarvNote")#"苗期综合评价",,"收获前综合评价"

# 搜索包含字符的行，并返回逻辑向量
has_char <- apply(cau.data[target_cols], 1, function(x) any(grepl("[[:alnum:][:punct:]]+", x)))


discard=cau.data[has_char,]
discard[,c("HarvNote")]
#discard[,c(27,28,87)]
ph_filtered=cau.data[!has_char,]
N_discard=nrow(ph_filtered[ph_filtered$PlotDiscarded %in% c("Yes"), ])

ph_filtered <- ph_filtered[!(ph_filtered$PlotDiscarded %in% c("Yes")), ]
names(ph_filtered)
ph=ph_filtered
ph$Type[substr(ph$Name,1,5)=="23TC2"]="TC2"

####HLJWK_P####
traits=c("YLD14",	"MST","PLOTWT","PHT","EHT")#"TWT",

combined_outblue = list()
par(mfrow=c(2,3),mar = c(5, 4, 4, 2) + 0.1,cex.main = 0.8)
# 循环遍历不同的 type
for (i in 1) {
  # 根据type选择特定的数据子集
  if (i == 1) {
    
    type="Field1_P1"
    subset_ph = ph[ph$Pass %in% 1:24 & ph$Range %in% 2:7 &ph$Type %in% c("P1","P3"), ]
    
  } 
  # 运行程序并生成 outblue
  outblue.tmp <-  programRun(dat = subset_ph, traits)
  
  # 将当前的 outblue 添加到 combined_outblue 列表中
  combined_outblue = c(combined_outblue, list(outblue.tmp))
  
}

idx=merge.fun(combined_outblue,traits)

count=analysis(idx,traits)
count$Type="P"
write.table(count,"Output/count.csv",sep=",",row.names=FALSE,col.names=FALSE,append=TRUE)

#####绘图

pdf("Output/HLJWK_plot_P.pdf")

par(mfrow=c(2,3),mar = c(5, 4, 4, 2) + 0.1,cex.main = 0.8) 
type="P"
drawplot(idx,traits)
dev.off()

outblue_name <- paste0("Output/23HLJWK_P_datalist_bluelist.rds")

saveRDS(idx,file = outblue_name)


####HLJWK_TC2####
traits=c("YLD14",	"MST","PLOTWT","PHT","EHT")#"TWT",

combined_outblue = list()
par(mfrow=c(2,3),mar = c(5, 4, 4, 2) + 0.1,cex.main = 0.8)
# 循环遍历不同的 type
for (i in 1) {
  # 根据type选择特定的数据子集
  if (i == 1) {
    
    type="Field1_TC2"
    subset_ph = ph[ph$Pass %in% 1:24 & ph$Range %in% 2:7 &ph$Type %in% c("TC2"), ]
  } 
  # 运行程序并生成 outblue
  outblue.tmp <-  programRun(dat = subset_ph, traits)
  
  # 将当前的 outblue 添加到 combined_outblue 列表中
  combined_outblue = c(combined_outblue, list(outblue.tmp))
  
}

idx=merge.fun(combined_outblue,traits)

count=analysis(idx,traits)
#count$Type="P"
write.table(count,"Output/count.csv",sep=",",row.names=FALSE,col.names=FALSE,append=TRUE)

#####绘图

pdf("Output/HLJWK_plot_TC2.pdf")

par(mfrow=c(2,3),mar = c(5, 4, 4, 2) + 0.1,cex.main = 0.8) 
type="TC2"
drawplot(idx,traits)
dev.off()

outblue_name <- paste0("Output/23HLJWK_TC2_datalist_bluelist.rds")

saveRDS(idx,file = outblue_name)






# HLJHUL --------------------------------------------------------------------
cau.data <- phenos %>% dplyr::filter(Location_C %in% "呼兰")
cau.data$AOA="EMSP"
cau.data$Year=2023


target_cols <- c("HarvNote")#"收获前综合评价", ,"苗期综合评价"

# 搜索包含字符的行，并返回逻辑向量
has_char <- apply(cau.data[target_cols], 1, function(x) any(grepl("[[:alnum:][:punct:]]+", x)))


discard=cau.data[has_char,]
discard[,c("HarvNote")]

ph_filtered=cau.data[!has_char,]
N_discard=nrow(ph_filtered[ph_filtered$PlotDiscarded %in% c("Yes"), ])
N_discard
ph_filtered <- ph_filtered[!(ph_filtered$PlotDiscarded %in% c("Yes")), ]
names(ph_filtered)
ph=ph_filtered


####HLJHUL_P####
traits=c("YLD14",	"MST","PLOTWT","TWT","PHT","EHT")#	

combined_outblue <- list()

# 循环遍历不同的 type
for (i in 1) {
  # 根据type选择特定的数据子集
  if (i == 1) {
    
    type="Field2_P1"
    subset_ph = ph[ph$Pass %in% 25:48 & ph$Range %in% 45:50, ]
    subset_ph=subset_ph[substr(subset_ph$Type,1,1)=="P",]
  } 
  
  # 运行程序并生成 outblue
  outblue.tmp = programRun(dat = subset_ph, traits)
  
  # 将当前的 outblue 添加到 combined_outblue 列表中
  combined_outblue = append(combined_outblue, list(outblue.tmp))
  
}

idx=merge.fun(combined_outblue,traits)

count=analysis(idx,traits)
count$Type="P"
write.table(count,"Output/count.csv",sep=",",row.names=FALSE,col.names=FALSE,append=TRUE)
#####绘图

pdf("Output/HLJHUL_plot_P.pdf")

par(mfrow=c(2,3),mar = c(5, 4, 4, 2) + 0.1,cex.main = 0.8) 
type="P"
drawplot(idx,traits)
dev.off()

outblue_name <- paste0("Output/23HLJHUL_P_datalist_bluelist.rds")

saveRDS(idx,file = outblue_name )


####HLJHUL_TC2####
traits=c("YLD14",	"MST","PLOTWT","TWT","PHT","EHT")#	

combined_outblue <- list()

# 循环遍历不同的 type
for (i in 1) {
  # 根据type选择特定的数据子集
 if (i == 1) {
    type="Field2_TC2"
    subset_ph = ph[ph$Pass %in% 25:48 & ph$Range %in% 45:50, ]
    subset_ph=subset_ph[substr(subset_ph$Name,1,5)=="23TC2",]
  }
  
  # 运行程序并生成 outblue
  outblue.tmp = programRun(dat = subset_ph, traits)
  
  # 将当前的 outblue 添加到 combined_outblue 列表中
  combined_outblue = append(combined_outblue, list(outblue.tmp))
  
}

idx=merge.fun(combined_outblue,traits)

count=analysis(idx,traits)
#count$Type="P"
write.table(count,"Output/count.csv",sep=",",row.names=FALSE,col.names=FALSE,append=TRUE)
#####绘图

pdf("Output/HLJHUL_plot_TC2.pdf")

par(mfrow=c(2,3),mar = c(5, 4, 4, 2) + 0.1,cex.main = 0.8) 
type="TC2"
drawplot(idx,traits)
dev.off()

outblue_name <- paste0("Output/23HLJHUL_TC2_datalist_bluelist.rds")

saveRDS(idx,file = outblue_name )


####HLJHUL_TC1####
traits=c("YLD14",	"MST","PLOTWT","TWT","PHT","EHT")#	

combined_outblue <- list()

# 循环遍历不同的 type
for (i in 1:2) {
  # 根据type选择特定的数据子集
  if (i == 1) {
    
    type="Field2_TC1"
    subset_ph = ph[ph$Pass %in% 25:48 & ph$Range %in% 51:59, ]
  } else if (i == 2) {
    type="Field3_TC1"
    subset_ph = ph[ph$Pass %in% 49:72 & ph$Range %in% 2:42, ]
  } 
  
  # 运行程序并生成 outblue
  outblue.tmp = programRun(dat = subset_ph, traits)
  
  # 将当前的 outblue 添加到 combined_outblue 列表中
  combined_outblue = append(combined_outblue, list(outblue.tmp))
  
}

idx=merge.fun(combined_outblue,traits)

count=analysis(idx,traits)
#count$Type="P"
write.table(count,"Output/count.csv",sep=",",row.names=FALSE,col.names=FALSE,append=TRUE)
#####绘图

pdf("Output/HLJHUL_plot_TC1.pdf")

par(mfrow=c(2,3),mar = c(5, 4, 4, 2) + 0.1,cex.main = 0.8) 
type="TC1"
drawplot(idx,traits)
dev.off()

outblue_name <- paste0("Output/23HLJHUL_TC1_datalist_bluelist.rds")

saveRDS(idx,file = outblue_name )












# HLJHULCC --------------------------------------------------------------------
cau.data <- phenos %>% dplyr::filter(Location_C %in% "呼兰CC")
cau.data$AOA="MMSP"
cau.data$Year=2023


target_cols <- c("HarvNote")#"收获前综合评价", ,"苗期综合评价"

# 搜索包含字符的行，并返回逻辑向量
has_char <- apply(cau.data[target_cols], 1, function(x) any(grepl("[[:alnum:][:punct:]]+", x)))


discard=cau.data[has_char,]
discard[,c("HarvNote")]

ph_filtered=cau.data[!has_char,]
N_discard=nrow(ph_filtered[ph_filtered$PlotDiscarded %in% c("Yes"), ])
N_discard
ph_filtered <- ph_filtered[!(ph_filtered$PlotDiscarded %in% c("Yes")), ]
names(ph_filtered)
ph=ph_filtered


####HLJHULCC_P####
traits=c("YLD14",	"MST","PLOTWT","TWT","PHT","EHT")#	

combined_outblue <- list()

# 循环遍历不同的 type
for (i in 1:3) {
  # 根据type选择特定的数据子集
   if (i == 1) {
    type="Field3_P1"
    subset_ph = ph[ph$Pass %in% 49:72 & ph$Range %in% 43:59, ]
  } else if (i == 2) {
    type="Field4_P1"
    subset_ph = ph[ph$Pass %in% 73:96 & ph$Range %in% 2:59, ]
  } else if (i == 3) {
    type="Field5_P1"
    subset_ph = ph[ph$Pass %in% 97:120 & ph$Range %in% 2:32, ]
  } 
  
  # 运行程序并生成 outblue
  outblue.tmp = programRun(dat = subset_ph, traits)
  
  # 将当前的 outblue 添加到 combined_outblue 列表中
  combined_outblue = append(combined_outblue, list(outblue.tmp))
  
}

idx=merge.fun(combined_outblue,traits)

count=analysis(idx,traits)
count$Type="P"
write.table(count,"Output/count.csv",sep=",",row.names=FALSE,col.names=FALSE,append=TRUE)
#####绘图

pdf("Output/HLJHULCC_plot_P.pdf")

par(mfrow=c(2,3),mar = c(5, 4, 4, 2) + 0.1,cex.main = 0.8) 
type="P"
drawplot(idx,traits)
dev.off()

outblue_name <- paste0("Output/23HLJHULCC_P_datalist_bluelist.rds")

saveRDS(idx,file = outblue_name )


####HLJHULCC_TC2####
traits=c("YLD14",	"MST","PLOTWT","TWT","PHT","EHT")#	

combined_outblue <- list()

# 循环遍历不同的 type
for (i in 1:2) {
  # 根据type选择特定的数据子集
  if (i == 1) {
    
    type="Field1_TC2"
    subset_ph = ph[ph$Pass %in% 1:24 & ph$Range %in% 2:59, ]
  } else if (i == 2) {
    type="Field2_TC2"
    subset_ph = ph[ph$Pass %in% 25:48 & ph$Range %in% 2:44, ]
  } 
  
  # 运行程序并生成 outblue
  outblue.tmp = programRun(dat = subset_ph, traits)
  
  # 将当前的 outblue 添加到 combined_outblue 列表中
  combined_outblue = append(combined_outblue, list(outblue.tmp))
  
}

idx=merge.fun(combined_outblue,traits)

count=analysis(idx,traits)
#count$Type="P"
write.table(count,"Output/count.csv",sep=",",row.names=FALSE,col.names=FALSE,append=TRUE)
#####绘图

pdf("Output/HLJHULCC_plot_TC2.pdf")

par(mfrow=c(2,3),mar = c(5, 4, 4, 2) + 0.1,cex.main = 0.8) 
type="TC2"
drawplot(idx,traits)
dev.off()

outblue_name <- paste0("Output/23HLJHULCC_TC2_datalist_bluelist.rds")

saveRDS(idx,file = outblue_name )


####HLJHULCC_TC1####
traits=c("YLD14",	"MST","PLOTWT","TWT","PHT","EHT")#	

combined_outblue <- list()

# 循环遍历不同的 type
for (i in 1:2) {
  # 根据type选择特定的数据子集
  if (i == 1) {
    type="Field5_TC1"
    subset_ph = ph[ph$Pass %in% 97:120 & ph$Range %in% 33:59, ]
  } else if (i == 2) {
    type="Field6_TC1"
    subset_ph = ph[ph$Pass %in% 121:124 & ph$Range %in% 2:7, ]
    subset_ph=subset_ph[!substr(subset_ph$Name,1,8)=="23TC2Liu",]
  } 
  
  # 运行程序并生成 outblue
  outblue.tmp = programRun(dat = subset_ph, traits)
  
  # 将当前的 outblue 添加到 combined_outblue 列表中
  combined_outblue = append(combined_outblue, list(outblue.tmp))
  
}

idx=merge.fun(combined_outblue,traits)

count=analysis(idx,traits)
#count$Type="P"
write.table(count,"Output/count.csv",sep=",",row.names=FALSE,col.names=FALSE,append=TRUE)
#####绘图

pdf("Output/HLJHULCC_plot_TC1.pdf")

par(mfrow=c(2,3),mar = c(5, 4, 4, 2) + 0.1,cex.main = 0.8) 
type="TC1"
drawplot(idx,traits)
dev.off()

outblue_name <- paste0("Output/23HLJHULCC_TC1_datalist_bluelist.rds")

saveRDS(idx,file = outblue_name )
















# HLJFY--------------------------------------------------------------------
cau.data <- phenos %>% dplyr::filter(Location %in% "HLJFY")
cau.data$AOA="EMSP"
cau.data$Year=2023


target_cols <- c("HarvNote")#"收获前综合评价", ,"苗期综合评价"

# 搜索包含字符的行，并返回逻辑向量
has_char <- apply(cau.data[target_cols], 1, function(x) any(grepl("[[:alnum:][:punct:]]+", x)))


discard=cau.data[has_char,]
discard[,c("HarvNote")]

ph_filtered=cau.data[!has_char,]
N_discard=nrow(ph_filtered[ph_filtered$PlotDiscarded %in% c("Yes"), ])
N_discard
ph_filtered <- ph_filtered[!(ph_filtered$PlotDiscarded %in% c("Yes")), ]
names(ph_filtered)
ph=ph_filtered

####HLJFY_P####
traits=c("YLD14",	"MST","PLOTWT","PHT","EHT")#	,"TWT"

combined_outblue <- list()

# 循环遍历不同的 type
for (i in 1) {
  # 根据type选择特定的数据子集
  if (i == 1) {
    
    type="Field1_P1"
    subset_ph = ph[ph$Pass %in% 1:24 & ph$Range %in% 2:6, ]
    subset_ph =subset_ph[substr(subset_ph$Type,1,1)=="P",]
  } 
  # 运行程序并生成 outblue
  outblue.tmp = programRun(dat = subset_ph, traits)
  
  # 将当前的 outblue 添加到 combined_outblue 列表中
  combined_outblue = append(combined_outblue, list(outblue.tmp))
  
}

idx=merge.fun(combined_outblue,traits)

count=analysis(idx,traits)
count$Type="P"
write.table(count,"Output/count.csv",sep=",",row.names=FALSE,col.names=FALSE,append=TRUE)
#####绘图

pdf("Output/HLJFY_plot_P.pdf")

par(mfrow=c(2,3),mar = c(5, 4, 4, 2) + 0.1,cex.main = 0.8) 
type="P"
drawplot(idx,traits)
dev.off()

outblue_name <- paste0("Output/23HLJFY_P_datalist_bluelist.rds")

saveRDS(idx,file = outblue_name )

####HLJFY_TC1####
traits=c("YLD14",	"MST","PLOTWT","PHT","EHT")#	,"TWT"

combined_outblue <- list()

# 循环遍历不同的 type
for (i in 1:2) {
  # 根据type选择特定的数据子集
  if (i == 1) {
    
    type="Field1_TC1"
    subset_ph = ph[ph$Pass %in% 1:24 & ph$Range %in% 7:22, ]
    #subset_ph =subset_ph[substr(subset_ph$Type,1,1)=="P",]
  } else if (i==2){
    type="Field2_TC1"
    subset_ph = ph[ph$Pass %in% 25:48 & ph$Range %in% 2:21, ]
  }
  # 运行程序并生成 outblue
  outblue.tmp = programRun(dat = subset_ph, traits)
  
  # 将当前的 outblue 添加到 combined_outblue 列表中
  combined_outblue = append(combined_outblue, list(outblue.tmp))
  
}

idx=merge.fun(combined_outblue,traits)

count=analysis(idx,traits)
#count$Type="P"
write.table(count,"Output/count.csv",sep=",",row.names=FALSE,col.names=FALSE,append=TRUE)
#####绘图

pdf("Output/HLJFY_plot_TC1.pdf")

par(mfrow=c(2,3),mar = c(5, 4, 4, 2) + 0.1,cex.main = 0.8) 
type="TC1"
drawplot(idx,traits)
dev.off()

outblue_name <- paste0("Output/23HLJFY_TC1_datalist_bluelist.rds")

saveRDS(idx,file = outblue_name )






####HLJFY_TC2####
traits=c("YLD14",	"MST","PLOTWT","PHT","EHT")#	,"TWT"

combined_outblue <- list()

# 循环遍历不同的 type
for (i in 1) {
  # 根据type选择特定的数据子集
  if (i == 1) {
    
    type="Field1_TC2"
    subset_ph = ph[ph$Pass %in% 1:24 & ph$Range %in% 2:6, ]
    subset_ph =subset_ph[substr(subset_ph$Type,1,3)=="TC2",]
  } 
  # 运行程序并生成 outblue
  outblue.tmp = programRun(dat = subset_ph, traits)
  
  # 将当前的 outblue 添加到 combined_outblue 列表中
  combined_outblue = append(combined_outblue, list(outblue.tmp))
  
}

idx=merge.fun(combined_outblue,traits)

count=analysis(idx,traits)
#count$Type="P"
write.table(count,"Output/count.csv",sep=",",row.names=FALSE,col.names=FALSE,append=TRUE)
#####绘图

pdf("Output/HLJFY_plot_TC2.pdf")

par(mfrow=c(2,3),mar = c(5, 4, 4, 2) + 0.1,cex.main = 0.8) 
type="TC2"
drawplot(idx,traits)
dev.off()

outblue_name <- paste0("Output/23HLJFY_TC2_datalist_bluelist.rds")

saveRDS(idx,file = outblue_name )
# HLJHL--------------------------------------------------------------------
cau.data <- phenos %>% dplyr::filter(Location %in% "HLJHL")
cau.data$AOA="EMSP"
cau.data$Year=2023


target_cols <- c("HarvNote")#"收获前综合评价", ,"苗期综合评价"

# 搜索包含字符的行，并返回逻辑向量
has_char <- apply(cau.data[target_cols], 1, function(x) any(grepl("[[:alnum:][:punct:]]+", x)))


discard=cau.data[has_char,]
discard[,c("HarvNote")]

ph_filtered=cau.data[!has_char,]
N_discard=nrow(ph_filtered[ph_filtered$PlotDiscarded %in% c("Yes"), ])
N_discard
ph_filtered <- ph_filtered[!(ph_filtered$PlotDiscarded %in% c("Yes")), ]
names(ph_filtered)
ph=ph_filtered

####HLJHL_P####
traits=c("YLD14",	"MST","PLOTWT")#	,"TWT","PHT","EHT"

combined_outblue <- list()

# 循环遍历不同的 type
for (i in 1) {
  # 根据type选择特定的数据子集
  if (i == 1) {
    
    type="Field1_P1"
    subset_ph = ph[ph$Pass %in% 1:24 & ph$Range %in% 2:6, ]
    subset_ph =subset_ph[substr(subset_ph$Type,1,1)=="P",]
  } 
  # 运行程序并生成 outblue
  outblue.tmp = programRun(dat = subset_ph, traits)
  
  # 将当前的 outblue 添加到 combined_outblue 列表中
  combined_outblue = append(combined_outblue, list(outblue.tmp))
  
}

idx=merge.fun(combined_outblue,traits)

count=analysis(idx,traits)
count$Type="P"
write.table(count,"Output/count.csv",sep=",",row.names=FALSE,col.names=FALSE,append=TRUE)
#####绘图

pdf("Output/HLJHL_plot_P.pdf")

par(mfrow=c(2,3),mar = c(5, 4, 4, 2) + 0.1,cex.main = 0.8) 
type="P"
drawplot(idx,traits)
dev.off()

outblue_name <- paste0("Output/23HLJHL_P_datalist_bluelist.rds")

saveRDS(idx,file = outblue_name )

####HLJHL_TC1####
traits=c("YLD14",	"MST","PLOTWT")#	,"TWT","PHT","EHT"

combined_outblue <- list()

# 循环遍历不同的 type
for (i in 1) {
  # 根据type选择特定的数据子集
  if (i == 1) {
    
    type="Field2_TC1"
    subset_ph = ph[ph$Pass %in% 1:24 & ph$Range %in% 7:41, ]
    #subset_ph =subset_ph[substr(subset_ph$Type,1,1)=="P",]
  } 
  # 运行程序并生成 outblue
  outblue.tmp = programRun(dat = subset_ph, traits)
  
  # 将当前的 outblue 添加到 combined_outblue 列表中
  combined_outblue = append(combined_outblue, list(outblue.tmp))
  
}

idx=merge.fun(combined_outblue,traits)

count=analysis(idx,traits)
#count$Type="P"
write.table(count,"Output/count.csv",sep=",",row.names=FALSE,col.names=FALSE,append=TRUE)
#####绘图

pdf("Output/HLJHL_plot_TC1.pdf")

par(mfrow=c(2,3),mar = c(5, 4, 4, 2) + 0.1,cex.main = 0.8) 
type="TC1"
drawplot(idx,traits)
dev.off()

outblue_name <- paste0("Output/23HLJHL_TC1_datalist_bluelist.rds")

saveRDS(idx,file = outblue_name )




####HLJHL_TC2####
traits=c("YLD14",	"MST","PLOTWT")#	,"TWT","PHT","EHT"

combined_outblue <- list()

# 循环遍历不同的 type
for (i in 1) {
  # 根据type选择特定的数据子集
  if (i == 1) {
    
    type="Field1_TC2"
    subset_ph = ph[ph$Pass %in% 1:24 & ph$Range %in% 2:6, ]
    subset_ph =subset_ph[substr(subset_ph$Type,1,2)=="TC",]
  } 
  # 运行程序并生成 outblue
  outblue.tmp = programRun(dat = subset_ph, traits)
  
  # 将当前的 outblue 添加到 combined_outblue 列表中
  combined_outblue = append(combined_outblue, list(outblue.tmp))
  
}

idx=merge.fun(combined_outblue,traits)

count=analysis(idx,traits)
#count$Type="P"
write.table(count,"Output/count.csv",sep=",",row.names=FALSE,col.names=FALSE,append=TRUE)
#####绘图

pdf("Output/HLJHL_plot_TC2.pdf")

par(mfrow=c(2,3),mar = c(5, 4, 4, 2) + 0.1,cex.main = 0.8) 
type="TC2"
drawplot(idx,traits)
dev.off()

outblue_name <- paste0("Output/23HLJHL_TC2_datalist_bluelist.rds")

saveRDS(idx,file = outblue_name )







#########################SU#######################
# HNQF_P --------------------------------------------------------------------
cau.data <- phenos %>% dplyr::filter(Location %in% "HNQF")
cau.data$AOA="MCSU"
cau.data$Year=2023


target_cols <- c("HarvNote")#"收获前综合评价", ,"苗期综合评价"

# 搜索包含字符的行，并返回逻辑向量
has_char <- apply(cau.data[target_cols], 1, function(x) any(grepl("[[:alnum:][:punct:]]+", x)))


discard=cau.data[has_char,]
discard[,c("HarvNote")]

ph_filtered=cau.data[!has_char,]
N_discard=nrow(ph_filtered[ph_filtered$PlotDiscarded %in% c("Yes"), ])

ph_filtered <- ph_filtered[!(ph_filtered$PlotDiscarded %in% c("Yes")), ]
names(ph_filtered)
ph=ph_filtered


traits=c("YLD14",	"MST",	"TWT","PLOTWT")

combined_outblue <- list()

# 循环遍历不同的 type
for (i in 1:5) {
  # 根据type选择特定的数据子集
  if (i == 1) {
    
    type="1_P1"
    subset_ph = ph[ph$Field == 1 & ph$Range %in% 29, ]
  } else if (i == 2) {
    type="P1"
    subset_ph = ph[ph$Field %in% 2:6 & ph$Range %in% 2:29, ]
  } else if (i == 3) {
    type="7_P1"
    subset_ph = ph[ph$Field %in% 7 & ph$Range %in% 2:19, ]
  } else if (i == 4) {
    type="8_P1"
    subset_ph = ph[ph$Field %in% 8, ]
  } else if (i==5){
    type="P3_4" 
    subset_ph <- ph[ph$Field == 1 & ph$Range %in% 2:28, ]
  }
  
  # 运行程序并生成 outblue
  outblue.tmp = programRun(dat = subset_ph, traits)
  
  # 将当前的 outblue 添加到 combined_outblue 列表中
  combined_outblue = append(combined_outblue, list(outblue.tmp))
  
}

idx=merge.fun(combined_outblue,traits)

count=analysis(idx,traits)
count$Type="P"
write.table(count,"Output/count.csv",sep=",",row.names=FALSE,col.names=FALSE,append=TRUE)
#####绘图

pdf("Output/HNQF_plot_P.pdf")

par(mfrow=c(2,3),mar = c(5, 4, 4, 2) + 0.1,cex.main = 0.8) 
type="P"
drawplot(idx,traits)
dev.off()

outblue_name <- paste0("Output/23HNQF_P_datalist_bluelist.rds")

saveRDS(idx,file = outblue_name )

 # a=readRDS("Output/23P_HNQF_datalist_bluelist.rds")
 # write.csv(a$datalist$datLabelled$MST,"Output/23P_HNQF_datLabelled_MST.csv")




# HNQB1_P -------------------------------------------------------------------
cau.data <- phenos %>% dplyr::filter(Location %in% "HNQB1")
cau.data$AOA="MCSU"
cau.data$Year=2023


target_cols <- c( "HarvNote")#"收获前综合评价""苗期综合评价",

# 搜索包含字符的行，并返回逻辑向量
has_char <- apply(cau.data[target_cols], 1, function(x) any(grepl("[[:alnum:][:punct:]]+", x)))


discard=cau.data[has_char,]
discard[,c("HarvNote")]

ph_filtered=cau.data[!has_char,]
N_discard=nrow(ph_filtered[ph_filtered$PlotDiscarded %in% c("Yes"), ])

ph_filtered <- ph_filtered[!(ph_filtered$PlotDiscarded %in% c("Yes")), ]
names(ph_filtered)

traits=c("YLD14",	"MST",	"TWT","PLOTWT")


ph=ph_filtered

combined_outblue = list()

# 循环遍历不同的 type
for (i in 1:3) {
  # 根据type选择特定的数据子集
  if (i == 1) {
    
    type="1_P3_4"
    subset_ph = ph[ph$Field %in% 1 & ph$Range %in% 2:30, ]
  } else if (i == 2) {
    type="1_P1"
    subset_ph = ph[ph$Field %in% 1 & ph$Range %in% 31:42, ]
  } else if (i==3){
    type="P1"
    subset_ph = ph[ph$Field %in% 2:3 & ph$Range %in% 2:42, ]
    
  } 
  
  # 运行程序并生成 outblue
  outblue.tmp <-  programRun(dat = subset_ph, traits)
  
  # 将当前的 outblue 添加到 combined_outblue 列表中
  combined_outblue = c(combined_outblue, list(outblue.tmp))
  
}

idx=merge.fun(combined_outblue,traits)

count=analysis(idx,traits)
count$Type="P"
write.table(count,"Output/count.csv",sep=",",row.names=FALSE,col.names=FALSE,append=TRUE)
#####绘图

pdf("Output/HNQB1_plot_P.pdf")

par(mfrow=c(2,3),mar = c(5, 4, 4, 2) + 0.1,cex.main = 0.8) 
type="P"
drawplot(idx,traits)
dev.off()

outblue_name <- paste0("Output/23HNQB1_P_datalist_bluelist.rds")

saveRDS(idx,file = outblue_name)

# HNQB2 -------------------------------------------------------------------
cau.data <- phenos %>% dplyr::filter(Location %in% "HNQB2")
cau.data$AOA="MCSU"
cau.data$Year=2023


target_cols <- c("苗期综合评价", "HarvNote")#"收获前综合评价"

# 搜索包含字符的行，并返回逻辑向量
has_char <- apply(cau.data[target_cols], 1, function(x) any(grepl("[[:alnum:][:punct:]]+", x)))


discard=cau.data[has_char,]
discard[,c("HarvNote")]

ph_filtered=cau.data[!has_char,]
N_discard=nrow(ph_filtered[ph_filtered$PlotDiscarded %in% c("Yes"), ])

ph_filtered <- ph_filtered[!(ph_filtered$PlotDiscarded %in% c("Yes")), ]
names(ph_filtered)
ph=ph_filtered

####HNQB2_P####
traits=c("YLD14",	"MST",	"TWT","PLOTWT")

combined_outblue = list()

# 循环遍历不同的 type
for (i in 1:3) {
  # 根据type选择特定的数据子集
  if (i == 1) {
    
    type="1_P1"
    subset_ph = ph[ph$Field %in% 1 & ph$Range %in% 2:46, ]
  } else if (i == 2) {
    type="3_P1"
    subset_ph = ph[ph$Field %in% 3 & ph$Range %in% 32:57, ]
  } else if (i==3){
    type="4_P1"
    subset_ph = ph[ph$Field %in% 4 & ph$Range %in% 13:58, ]
    
  } 
  
  # 运行程序并生成 outblue
  outblue.tmp <-  programRun(dat = subset_ph, traits)
  
  # 将当前的 outblue 添加到 combined_outblue 列表中
  combined_outblue = c(combined_outblue, list(outblue.tmp))
  
}

idx=merge.fun(combined_outblue,traits)

count=analysis(idx,traits)
count$Type="P"
write.table(count,"Output/count.csv",sep=",",row.names=FALSE,col.names=FALSE,append=TRUE)
#####绘图

pdf("Output/HNQB2_plot_P.pdf")

par(mfrow=c(2,3),mar = c(5, 4, 4, 2) + 0.1,cex.main = 0.8) 
type="P"
drawplot(idx,traits)
dev.off()

outblue_name <- paste0("Output/23HNQB2_P_datalist_bluelist.rds")

saveRDS(idx,file = outblue_name)





####HNQB2_TC2####
traits=c("YLD14",	"MST",	"TWT","PLOTWT")

combined_outblue = list()

# 循环遍历不同的 type
for (i in 1:3) {
  # 根据type选择特定的数据子集
  if (i == 1) {
    
    type="1_TC2"
    subset_ph = ph[ph$Field %in% 1 & ph$Range %in% 48:57, ]
  } else if (i == 2) {
    type="2_TC2"
    subset_ph = ph[ph$Field %in% 2 & ph$Range %in% 2:57, ]
  } else if (i==3){
    type="3_TC2"
    subset_ph = ph[ph$Field %in% 3 & ph$Range %in% 7:30, ]
    
  } 
  
  # 运行程序并生成 outblue
  outblue.tmp <-  programRun(dat = subset_ph, traits)
  
  # 将当前的 outblue 添加到 combined_outblue 列表中
  combined_outblue = c(combined_outblue, list(outblue.tmp))
  
}

idx=merge.fun(combined_outblue,traits)

count=analysis(idx,traits)

write.table(count,"Output/count.csv",sep=",",row.names=FALSE,col.names=FALSE,append=TRUE)
#####绘图

pdf("Output/HNQB2_plot_TC2.pdf")

par(mfrow=c(2,3),mar = c(5, 4, 4, 2) + 0.1,cex.main = 0.8) 
type="TC2"
drawplot(idx,traits)
dev.off()

outblue_name <- paste0("Output/23HNQB2_TC2_datalist_bluelist.rds")

saveRDS(idx,file = outblue_name)




# HNXX1--------------------------------------------------------------------

cau.data <- phenos %>% dplyr::filter(Location %in% "HNXX1")
cau.data$AOA="MCSU"
cau.data$Year=2023


target_cols <- c( "HarvNote")#"收获前综合评价","苗期综合评价"

# 搜索包含字符的行，并返回逻辑向量
has_char <- apply(cau.data[target_cols], 1, function(x) any(grepl("[[:alnum:][:punct:]]+", x)))


discard=cau.data[has_char,]
discard[,c("HarvNote")]

ph_filtered=cau.data[!has_char,]
N_discard=nrow(ph_filtered[ph_filtered$PlotDiscarded %in% c("Yes"), ])

ph_filtered <- ph_filtered[!(ph_filtered$PlotDiscarded %in% c("Yes")), ]
names(ph_filtered)
ph=ph_filtered

####HNXX1_P####
traits=c("YLD14",	"MST","PLOTWT","PHT","EHT")

combined_outblue = list()

# 循环遍历不同的 type
for (i in 1:4) {
  # 根据type选择特定的数据子集
  if (i == 1) {
    
    type="Field1_P3_4"
    subset_ph = ph[ph$Pass %in% 1:23 & ph$Range %in% 2:27, ]
  } else if (i == 2) {
    type="Field2_P3_4"
    subset_ph = ph[ph$Pass %in% 25:47 & ph$Range %in% 2:9, ]
  } else if (i==3){
    type="Field2_P1"
    subset_ph = ph[ph$Pass %in% 25:48 & ph$Range %in% 10:27, ]
    
  } else if (i==4){
    type="Field3_4_5_P1"
    subset_ph = ph[ph$Pass %in% 49:104 & ph$Range %in% 2:29, ]
    
  } 
  
  # 运行程序并生成 outblue
  outblue.tmp <-  programRun(dat = subset_ph, traits)
  
  # 将当前的 outblue 添加到 combined_outblue 列表中
  combined_outblue = c(combined_outblue, list(outblue.tmp))
  
}

idx=merge.fun(combined_outblue,traits)

count=analysis(idx,traits)
count$Type="P"
write.table(count,"Output/count.csv",sep=",",row.names=FALSE,col.names=FALSE,append=TRUE)
#####绘图

pdf("Output/HNXX1_plot_P.pdf")

par(mfrow=c(2,3),mar = c(5, 4, 4, 2) + 0.1,cex.main = 0.8) 
type="P"
drawplot(idx,traits)
dev.off()

outblue_name <- paste0("Output/23HNXX1_P_datalist_bluelist.rds")

saveRDS(idx,file = outblue_name)

####HNXX1_TC1####
traits=c("YLD14",	"MST","PLOTWT","PHT","EHT")

combined_outblue = list()

  type="Field6_TC1"
  subset_ph = ph[ph$Pass %in% 105:106 & ph$Range %in% 1:30, ]
  
  subset_ph =subset_ph[subset_ph$Type=="TC1",]
  
  # 运行程序并生成 outblue
  outblue.tmp <-  programRun(dat = subset_ph, traits)
  
  # 将当前的 outblue 添加到 combined_outblue 列表中
  combined_outblue = c(combined_outblue, list(outblue.tmp))
  


idx=merge.fun(combined_outblue,traits)

count=analysis(idx,traits)
write.table(count,"Output/count.csv",sep=",",row.names=FALSE,col.names=FALSE,append=TRUE)
#####绘图

pdf("Output/HNXX1_plot_TC1.pdf")

par(mfrow=c(2,3),mar = c(5, 4, 4, 2) + 0.1,cex.main = 0.8) 
type="TC1"
drawplot(idx,traits)
dev.off()

outblue_name <- paste0("Output/23HNXX1_TC1_datalist_bluelist.rds")

saveRDS(idx,file = outblue_name)

####HNXX1_TC2####
traits=c("YLD14",	"MST","PLOTWT","PHT","EHT")

combined_outblue = list()

type="Field6_TC1"
subset_ph = ph[ph$Pass %in% 105:106 & ph$Range %in% 1:30, ]

subset_ph =subset_ph[subset_ph$Type=="TC2",]

# 运行程序并生成 outblue
outblue.tmp <-  programRun(dat = subset_ph, traits)

# 将当前的 outblue 添加到 combined_outblue 列表中
combined_outblue = c(combined_outblue, list(outblue.tmp))



idx=merge.fun(combined_outblue,traits)

count=analysis(idx,traits)
write.table(count,"Output/count.csv",sep=",",row.names=FALSE,col.names=FALSE,append=TRUE)
#####绘图

pdf("Output/HNXX1_plot_TC2.pdf")

par(mfrow=c(2,3),mar = c(5, 4, 4, 2) + 0.1,cex.main = 0.8) 
type="TC2"
drawplot(idx,traits)
dev.off()

outblue_name <- paste0("Output/23HNXX1_TC2_datalist_bluelist.rds")

saveRDS(idx,file = outblue_name)





# HNXX2_P --------------------------------------------------------------------

cau.data <- phenos %>% dplyr::filter(Location %in% "HNXX2")
cau.data$AOA="MCSU"
cau.data$Year=2023


target_cols <- c( "HarvNote")#"苗期综合评价","收获前综合评价"

# 搜索包含字符的行，并返回逻辑向量
has_char <- apply(cau.data[target_cols], 1, function(x) any(grepl("[[:alnum:][:punct:]]+", x)))


discard=cau.data[has_char,]
discard[,c("HarvNote")]

ph_filtered=cau.data[!has_char,]
N_discard=nrow(ph_filtered[ph_filtered$PlotDiscarded %in% c("Yes"), ])

ph_filtered <- ph_filtered[!(ph_filtered$PlotDiscarded %in% c("Yes")), ]
names(ph_filtered)
ph=ph_filtered

traits=c("YLD14",	"MST","PLOTWT","PHT","EHT")

combined_outblue = list()

# 循环遍历不同的 type
for (i in 1:3) {
  # 根据type选择特定的数据子集
  if (i == 1) {
    
    type="Field1_P1"
    subset_ph = ph[ph$Pass %in% 1:24 & ph$Range %in% 2:26, ]
  } else if (i == 2) {
    type="Field2_3_P1"
    subset_ph = ph[ph$Pass %in% 25:72 & ph$Range %in% 2:26, ]
  } else if (i==3){
    type="Field4_P1"
    subset_ph = ph[ph$Pass %in% 73:96 & ph$Range %in% 2:18, ]
    
  } 
  
  # 运行程序并生成 outblue
  outblue.tmp <-  programRun(dat = subset_ph, traits)
  
  # 将当前的 outblue 添加到 combined_outblue 列表中
  combined_outblue = c(combined_outblue, list(outblue.tmp))
  
}

idx=merge.fun(combined_outblue,traits)

count=analysis(idx,traits)
count$Type="P"
write.table(count,"Output/count.csv",sep=",",row.names=FALSE,col.names=FALSE,append=TRUE)

#####绘图

pdf("Output/HNXX2_plot_P.pdf")

par(mfrow=c(2,3),mar = c(5, 4, 4, 2) + 0.1,cex.main = 0.8) 
type="P"
drawplot(idx,traits)
dev.off()

outblue_name <- paste0("Output/23HNXX2_P_datalist_bluelist.rds")

saveRDS(idx,file = outblue_name)



# HNXX3--------------------------------------------------------------------

cau.data <- phenos %>% dplyr::filter(Location %in% "HNXX3")
cau.data$AOA="MCSU"
cau.data$Year=2023


target_cols <- c( "HarvNote")#"苗期综合评价","收获前综合评价"

# 搜索包含字符的行，并返回逻辑向量
has_char <- apply(cau.data[target_cols], 1, function(x) any(grepl("[[:alnum:][:punct:]]+", x)))


discard=cau.data[has_char,]
discard[,c("HarvNote")]

ph_filtered=cau.data[!has_char,]
N_discard=nrow(ph_filtered[ph_filtered$PlotDiscarded %in% c("Yes"), ])

ph_filtered <- ph_filtered[!(ph_filtered$PlotDiscarded %in% c("Yes")), ]
names(ph_filtered)
ph=ph_filtered


####HNXX3_TC2####
traits=c("YLD14",	"MST","PLOTWT","PHT","EHT")

combined_outblue = list()

# 循环遍历不同的 type
for (i in 1:3) {
  # 根据type选择特定的数据子集
  if (i == 1) {
    
    type="Field1_2_TC2"
    subset_ph = ph[ph$Pass %in% 1:48 & ph$Range %in% 6:44, ]
  } else if (i == 2) {
    type="Field3_TC2"
    subset_ph = ph[ph$Pass %in% 49:72 & ph$Range %in% 5:43, ]
  } else if (i==3){
    type="Field4_TC2"
    subset_ph = ph[ph$Pass %in% 73:96 & ph$Range %in% 5:30, ]
    
  } 
  
  # 运行程序并生成 outblue
  outblue.tmp <-  programRun(dat = subset_ph, traits)
  
  # 将当前的 outblue 添加到 combined_outblue 列表中
  combined_outblue = c(combined_outblue, list(outblue.tmp))
  
}

idx=merge.fun(combined_outblue,traits)

count=analysis(idx,traits)

write.table(count,"Output/count.csv",sep=",",row.names=FALSE,col.names=FALSE,append=TRUE)

#####绘图

pdf("Output/HNXX3_plot_TC2.pdf")

par(mfrow=c(2,3),mar = c(5, 4, 4, 2) + 0.1,cex.main = 0.8) 
type="TC2"
drawplot(idx,traits)
dev.off()

outblue_name <- paste0("Output/23HNXX3_TC2_datalist_bluelist.rds")

saveRDS(idx,file = outblue_name)


####HNXX3_TC1####
traits=c("YLD14",	"MST","PLOTWT","PHT","EHT")

combined_outblue = list()

# 循环遍历不同的 type
for (i in 1:4) {
  # 根据type选择特定的数据子集
  if (i == 1) {
    
    type="Field4_TC1"
    subset_ph = ph[ph$Pass %in% 73:96 & ph$Range %in% 31:43, ]
  } else if (i == 2) {
    type="Field5_6_7_TC1"
    subset_ph = ph[ph$Pass %in% 97:168 & ph$Range %in% 4:42, ]
  } else if (i==3){
    type="Field8_9_10_TC1"
    subset_ph = ph[ph$Pass %in% 169:240 & ph$Range %in% 3:42, ]
    
  } else if (i==4){
    type="Field11_12_TC1"
    subset_ph = ph[ph$Pass %in% 241:253 & ph$Range %in% 2:41, ]
    
  } 
  
  # 运行程序并生成 outblue
  outblue.tmp <-  programRun(dat = subset_ph, traits)
  
  # 将当前的 outblue 添加到 combined_outblue 列表中
  combined_outblue = c(combined_outblue, list(outblue.tmp))
  
}

idx=merge.fun(combined_outblue,traits)

count=analysis(idx,traits)
write.table(count,"Output/count.csv",sep=",",row.names=FALSE,col.names=FALSE,append=TRUE)

#####绘图

pdf("Output/HNXX3_plot_TC1.pdf")

par(mfrow=c(2,3),mar = c(5, 4, 4, 2) + 0.1,cex.main = 0.8) 
type="TC1"
drawplot(idx,traits)
dev.off()

outblue_name <- paste0("Output/23HNXX3_TC1_datalist_bluelist.rds")

saveRDS(idx,file = outblue_name)







# HNZC_P --------------------------------------------------------------------

cau.data <- phenos %>% dplyr::filter(Location %in% "HNZC")
cau.data$AOA="SCSU"
cau.data$Year=2023


target_cols <- c( "HarvNote")#"苗期综合评价",,"收获前综合评价"

# 搜索包含字符的行，并返回逻辑向量
has_char <- apply(cau.data[target_cols], 1, function(x) any(grepl("[[:alnum:][:punct:]]+", x)))


discard=cau.data[has_char,]
discard[,c("HarvNote")]
#discard[,c(27,28,87)]
ph_filtered=cau.data[!has_char,]
N_discard=nrow(ph_filtered[ph_filtered$PlotDiscarded %in% c("Yes"), ])

ph_filtered <- ph_filtered[!(ph_filtered$PlotDiscarded %in% c("Yes")), ]
names(ph_filtered)
ph=ph_filtered


traits=c("YLD14",	"MST","PLOTWT","TWT","PHT","EHT")#

combined_outblue = list()

# 循环遍历不同的 type
for (i in 1:4) {
  # 根据type选择特定的数据子集
  if (i == 1) {
    
    type="Field1_P3_P4"
    subset_ph = ph[ph$Pass %in% 1:23 & ph$Range %in% 2:28, ]
  } else if(i==2) {
    type="Field1_P1"
    subset_ph = ph[ph$Pass %in% 1:24 & ph$Range %in% 29:36, ]
  } else if(i==3) {
    type="Field2_3_4_5_P1"
    subset_ph = ph[ph$Pass %in% 25:120 & ph$Range %in% 3:36, ]
  } else if(i==4) {
    type="Field6_P1"
    subset_ph = ph[ph$Pass %in% 121:132 & ph$Range %in% 7:26, ]
  }
  # 运行程序并生成 outblue
  outblue.tmp <-  programRun(dat = subset_ph, traits)
  
  # 将当前的 outblue 添加到 combined_outblue 列表中
  combined_outblue = c(combined_outblue, list(outblue.tmp))
  
}

idx=merge.fun(combined_outblue,traits)

count=analysis(idx,traits)
count$Type="P"
write.table(count,"Output/count.csv",sep=",",row.names=FALSE,col.names=FALSE,append=TRUE)

#####绘图

pdf("Output/HNZC_plot_P.pdf")

par(mfrow=c(2,3),mar = c(5, 4, 4, 2) + 0.1,cex.main = 0.8) 
type="P"
drawplot(idx,traits)
dev.off()

outblue_name <- paste0("Output/23HNZC_P_datalist_bluelist.rds")

saveRDS(idx,file = outblue_name)





# HNXM1 --------------------------------------------------------------------


cau.data <- phenos %>% dplyr::filter(Location %in% "HNXM1")
cau.data$AOA="SCSU"
cau.data$Year=2023


target_cols <- c( "HarvNote")#"苗期综合评价",,"收获前综合评价"

# 搜索包含字符的行，并返回逻辑向量
has_char <- apply(cau.data[target_cols], 1, function(x) any(grepl("[[:alnum:][:punct:]]+", x)))


discard=cau.data[has_char,]
discard[,c("HarvNote")]
#discard[,c(27,28,87)]
ph_filtered=cau.data[!has_char,]
N_discard=nrow(ph_filtered[ph_filtered$PlotDiscarded %in% c("Yes"), ])

ph_filtered <- ph_filtered[!(ph_filtered$PlotDiscarded %in% c("Yes")), ]
names(ph_filtered)
ph=ph_filtered


####HNXM1_P####
traits=c("YLD14",	"MST","PLOTWT","TWT","PHT","EHT")#

combined_outblue = list()

# 循环遍历不同的 type
for (i in 1:8) {
  # 根据type选择特定的数据子集
  if (i == 1) {
    
    type="Field1_P3_P4"
    subset_ph = ph[ph$Pass %in% 1:23 & ph$Range %in% 6:25, ]
  } else if(i==2){
    type="Field2_P3_P4"
    subset_ph = ph[ph$Pass %in% 25:47 & ph$Range %in% 6:11, ]
    
  } else if(i==3){
    type="Field2_P1"
    subset_ph = ph[ph$Pass %in% 25:48 & ph$Range %in% 12:25, ]
    
  } else if(i==4){
    type="Field3_P1"
    subset_ph = ph[ph$Pass %in% 49:72 & ph$Range %in% 1:26, ]
    
  } else if(i==5){
    type="Field4_P1"
    subset_ph = ph[ph$Pass %in% 73:96 & ph$Range %in% 6:26, ]
    
  } else if(i==6){
    type="Field5_6_P1"
    subset_ph = ph[ph$Pass %in% 97:144 & ph$Range %in% 3:26, ]
    
  } else if(i==7){
    type="Field7_P1"
    subset_ph = ph[ph$Pass %in% 145:168 & ph$Range %in% 4:22, ]
    
  } else if(i==8){
    type="Field8_P1"
    subset_ph = ph[ph$Pass %in% 169:192 & ph$Range %in% 5:21, ]
    
  }
  # 运行程序并生成 outblue
  outblue.tmp <-  programRun(dat = subset_ph, traits)
  
  # 将当前的 outblue 添加到 combined_outblue 列表中
  combined_outblue = c(combined_outblue, list(outblue.tmp))
  
}

idx=merge.fun(combined_outblue,traits)

count=analysis(idx,traits)
count$Type="P"
write.table(count,"Output/count.csv",sep=",",row.names=FALSE,col.names=FALSE,append=TRUE)

#####绘图

pdf("Output/HNXM1_plot_P.pdf")

par(mfrow=c(2,3),mar = c(5, 4, 4, 2) + 0.1,cex.main = 0.8) 
type="P"
drawplot(idx,traits)
dev.off()

outblue_name <- paste0("Output/23HNXM1_P_datalist_bluelist.rds")

saveRDS(idx,file = outblue_name)



####HNXM1_TC1####
traits=c("YLD14",	"MST","PLOTWT","TWT","PHT","EHT")#

combined_outblue = list()

# 循环遍历不同的 type
for (i in 1:2) {
  # 根据type选择特定的数据子集
  if (i == 1) {
    
    type="Field9_TC1"
    subset_ph = ph[ph$Pass %in% 193:216 & ph$Range %in% 5:22, ]
  } else if(i==2) {
    type="Field10_TC1"
    subset_ph = ph[ph$Pass %in% 217:240 & ph$Range %in% 6:23, ]
  } 
  # 运行程序并生成 outblue
  outblue.tmp <-  programRun(dat = subset_ph, traits)
  
  # 将当前的 outblue 添加到 combined_outblue 列表中
  combined_outblue = c(combined_outblue, list(outblue.tmp))
  
}

idx=merge.fun(combined_outblue,traits)

count=analysis(idx,traits)
#count$Type="P"
write.table(count,"Output/count.csv",sep=",",row.names=FALSE,col.names=FALSE,append=TRUE)

#####绘图

pdf("Output/HNXM1_plot_TC1.pdf")

par(mfrow=c(2,3),mar = c(5, 4, 4, 2) + 0.1,cex.main = 0.8) 
type="TC1"
drawplot(idx,traits)
dev.off()

outblue_name <- paste0("Output/23HNXM1_TC1_datalist_bluelist.rds")

saveRDS(idx,file = outblue_name)






# HNXM2_TC1 --------------------------------------------------------------------

cau.data <- phenos %>% dplyr::filter(Location %in% "HNXM2")
cau.data$AOA="SCSU"
cau.data$Year=2023


target_cols <- c( "HarvNote")#"苗期综合评价",,"收获前综合评价"

# 搜索包含字符的行，并返回逻辑向量
has_char <- apply(cau.data[target_cols], 1, function(x) any(grepl("[[:alnum:][:punct:]]+", x)))


discard=cau.data[has_char,]
discard[,c("HarvNote")]
#discard[,c(27,28,87)]
ph_filtered=cau.data[!has_char,]
N_discard=nrow(ph_filtered[ph_filtered$PlotDiscarded %in% c("Yes"), ])

ph_filtered <- ph_filtered[!(ph_filtered$PlotDiscarded %in% c("Yes")), ]
names(ph_filtered)
ph=ph_filtered


traits=c("YLD14",	"MST","PLOTWT","TWT","PHT","EHT")#

combined_outblue = list()

# 循环遍历不同的 type
for (i in 1:3) {
  # 根据type选择特定的数据子集
  if (i == 1) {
    
    type="Field1_TC1"
    subset_ph = ph[ph$Pass %in% 1:24 & ph$Range %in% 1:15, ]
  } else if(i==2) {
    type="Field2_TC1"
    subset_ph = ph[ph$Pass %in% 25:48 & ph$Range %in% 1:36, ]
  } else if(i==3) {
    type="Field3_4_5_TC1"
    subset_ph = ph[ph$Pass %in% 49:120 & ph$Range %in% 19:38, ]
  } 
  # 运行程序并生成 outblue
  outblue.tmp <-  programRun(dat = subset_ph, traits)
  
  # 将当前的 outblue 添加到 combined_outblue 列表中
  combined_outblue = c(combined_outblue, list(outblue.tmp))
  
}

idx=merge.fun(combined_outblue,traits)

count=analysis(idx,traits)
#count$Type="P"
write.table(count,"Output/count.csv",sep=",",row.names=FALSE,col.names=FALSE,append=TRUE)

#####绘图

pdf("Output/HNXM2_plot_TC1.pdf")

par(mfrow=c(2,3),mar = c(5, 4, 4, 2) + 0.1,cex.main = 0.8) 
type="TC1"
drawplot(idx,traits)
dev.off()

outblue_name <- paste0("Output/23HNXM2_TC1_datalist_bluelist.rds")

saveRDS(idx,file = outblue_name)




# HNXM3 --------------------------------------------------------------------


cau.data <- phenos %>% dplyr::filter(Location %in% "HNXM3")
cau.data$AOA="SCSU"
cau.data$Year=2023


target_cols <- c( "HarvNote")#"苗期综合评价",,"收获前综合评价"

# 搜索包含字符的行，并返回逻辑向量
has_char <- apply(cau.data[target_cols], 1, function(x) any(grepl("[[:alnum:][:punct:]]+", x)))


discard=cau.data[has_char,]
discard[,c("HarvNote")]
#discard[,c(27,28,87)]
ph_filtered=cau.data[!has_char,]
N_discard=nrow(ph_filtered[ph_filtered$PlotDiscarded %in% c("Yes"), ])

ph_filtered <- ph_filtered[!(ph_filtered$PlotDiscarded %in% c("Yes")), ]
names(ph_filtered)
ph=ph_filtered

####HNXM3_TC1####
traits=c("YLD14",	"MST","PLOTWT","TWT","PHT","EHT")#

combined_outblue = list()

# 循环遍历不同的 type
for (i in 1:2) {
  # 根据type选择特定的数据子集
  if (i == 1) {
    
    type="Field1_2_TC1"
    subset_ph = ph[ph$Pass %in% 1:48 & ph$Range %in% 2:27, ]
  } else if(i==2) {
    type="Field5_TC1"
    subset_ph = ph[ph$Pass %in% 97 & ph$Range %in% 10:22, ]
  } 
  # 运行程序并生成 outblue
  outblue.tmp <-  programRun(dat = subset_ph, traits)
  
  # 将当前的 outblue 添加到 combined_outblue 列表中
  combined_outblue = c(combined_outblue, list(outblue.tmp))
  
}

idx=merge.fun(combined_outblue,traits)

count=analysis(idx,traits)
#count$Type="P"
write.table(count,"Output/count.csv",sep=",",row.names=FALSE,col.names=FALSE,append=TRUE)

#####绘图

pdf("Output/HNXM3_plot_TC1.pdf")

par(mfrow=c(2,3),mar = c(5, 4, 4, 2) + 0.1,cex.main = 0.8) 
type="TC1"
drawplot(idx,traits)
dev.off()

outblue_name <- paste0("Output/23HNXM3_TC1_datalist_bluelist.rds")

saveRDS(idx,file = outblue_name)

####HNXM3_TC2####
traits=c("YLD14",	"MST","PLOTWT","TWT","PHT","EHT")#

combined_outblue = list()

# 循环遍历不同的 type
for (i in 1:3) {
  # 根据type选择特定的数据子集
  if (i == 1) {
    
    type="Field3_TC2"
    subset_ph = ph[ph$Pass %in% 49:72 & ph$Range %in% 2:27, ]
  } else if(i==2) {
    type="Field4_TC2"
    subset_ph = ph[ph$Pass %in% 73:96 & ph$Range %in% 2:22, ]
  } else if(i==3) {
    type="Field5_TC2"
    subset_ph = ph[ph$Pass %in% 97:98 & ph$Range %in% 2:9, ]
  } 
  # 运行程序并生成 outblue
  outblue.tmp <-  programRun(dat = subset_ph, traits)
  
  # 将当前的 outblue 添加到 combined_outblue 列表中
  combined_outblue = c(combined_outblue, list(outblue.tmp))
  
}

idx=merge.fun(combined_outblue,traits)

count=analysis(idx,traits)
#count$Type="P"
write.table(count,"Output/count.csv",sep=",",row.names=FALSE,col.names=FALSE,append=TRUE)

#####绘图

pdf("Output/HNXM3_plot_TC2.pdf")

par(mfrow=c(2,3),mar = c(5, 4, 4, 2) + 0.1,cex.main = 0.8) 
type="TC2"
drawplot(idx,traits)
dev.off()

outblue_name <- paste0("Output/23HNXM3_TC2_datalist_bluelist.rds")

saveRDS(idx,file = outblue_name)

####HNXM3_P####
traits=c("YLD14",	"MST","PLOTWT","TWT","PHT","EHT")#

combined_outblue = list()

# 循环遍历不同的 type
for (i in 1) {
  # 根据type选择特定的数据子集
  if (i == 1) {
    
    type="Field5_P"
    subset_ph = ph[ph$Pass %in% 99:101 & ph$Range %in% 2:22, ]
  } 
  # 运行程序并生成 outblue
  outblue.tmp <-  programRun(dat = subset_ph, traits)
  
  # 将当前的 outblue 添加到 combined_outblue 列表中
  combined_outblue = c(combined_outblue, list(outblue.tmp))
  
}

idx=merge.fun(combined_outblue,traits)

count=analysis(idx,traits)
count$Type="P"
write.table(count,"Output/count.csv",sep=",",row.names=FALSE,col.names=FALSE,append=TRUE)

#####绘图

pdf("Output/HNXM3_plot_P.pdf")

par(mfrow=c(2,3),mar = c(5, 4, 4, 2) + 0.1,cex.main = 0.8) 
type="P"
drawplot(idx,traits)
dev.off()

outblue_name <- paste0("Output/23HNXM3_P_datalist_bluelist.rds")

saveRDS(idx,file = outblue_name)





# HNXM4_P -----------------------------------------------------------------
cau.data <- phenos %>% dplyr::filter(Location %in% "HNXM4")
cau.data$AOA="SCSU"
cau.data$Year=2023


target_cols <- c( "HarvNote")#"苗期综合评价",,"收获前综合评价"

# 搜索包含字符的行，并返回逻辑向量
has_char <- apply(cau.data[target_cols], 1, function(x) any(grepl("[[:alnum:][:punct:]]+", x)))


discard=cau.data[has_char,]
discard[,c("HarvNote")]
#discard[,c(27,28,87)]
ph_filtered=cau.data[!has_char,]
N_discard=nrow(ph_filtered[ph_filtered$PlotDiscarded %in% c("Yes"), ])

ph_filtered <- ph_filtered[!(ph_filtered$PlotDiscarded %in% c("Yes")), ]
names(ph_filtered)
ph=ph_filtered



traits=c("YLD14",	"MST","PLOTWT","TWT","PHT","EHT")#

combined_outblue = list()

# 循环遍历不同的 type
for (i in 1:2) {
  # 根据type选择特定的数据子集
  if (i == 1) {
    
    type="Field1_2_P1"
    subset_ph = ph[ph$Pass %in% 1:36 & ph$Range %in% 1:10, ]
  } else if(i==2){
    type="Field3_P1"
    subset_ph = ph[ph$Pass %in% 37:48 & ph$Range %in% 1:7, ]
  }
  # 运行程序并生成 outblue
  outblue.tmp <-  programRun(dat = subset_ph, traits)
  
  # 将当前的 outblue 添加到 combined_outblue 列表中
  combined_outblue = c(combined_outblue, list(outblue.tmp))
  
}

idx=merge.fun(combined_outblue,traits)

count=analysis(idx,traits)
count$Type="P"
write.table(count,"Output/count.csv",sep=",",row.names=FALSE,col.names=FALSE,append=TRUE)

#####绘图

pdf("Output/HNXM4_plot_P.pdf")

par(mfrow=c(2,3),mar = c(5, 4, 4, 2) + 0.1,cex.main = 0.8) 
type="P"
drawplot(idx,traits)
dev.off()

outblue_name <- paste0("Output/23HNXM4_P_datalist_bluelist.rds")

saveRDS(idx,file = outblue_name)






# HNXM5 -----------------------------------------------------------------
cau.data <- phenos %>% dplyr::filter(Location %in% "HNXM5")
cau.data$AOA="SCSU"
cau.data$Year=2023


target_cols <- c( "HarvNote")#"苗期综合评价",,"收获前综合评价"

# 搜索包含字符的行，并返回逻辑向量
has_char <- apply(cau.data[target_cols], 1, function(x) any(grepl("[[:alnum:][:punct:]]+", x)))


discard=cau.data[has_char,]
discard[,c("HarvNote")]
#discard[,c(27,28,87)]
ph_filtered=cau.data[!has_char,]
N_discard=nrow(ph_filtered[ph_filtered$PlotDiscarded %in% c("Yes"), ])

ph_filtered <- ph_filtered[!(ph_filtered$PlotDiscarded %in% c("Yes")), ]
names(ph_filtered)
ph=ph_filtered

####HNXM5_P####
traits=c("YLD14",	"MST","PLOTWT","TWT","PHT","EHT")#

combined_outblue = list()

# 循环遍历不同的 type
for (i in 1:3) {
  # 根据type选择特定的数据子集
  if (i == 1) {
    
    type="Field1_P1"
    subset_ph = ph[ph$Pass %in% 1:16 & ph$Range %in% 1:6, ]
    subset_ph = subset_ph[substr(subset_ph$Type,1,1)=="P",]
  } else if(i==2){
    type="Field2_P1"
    subset_ph = ph[ph$Pass %in% 17:20 & ph$Range %in% 1:4, ]
    subset_ph = subset_ph[substr(subset_ph$Type,1,1)=="P",]
  } else if(i==3){
    type="Field3_P1"
    subset_ph = ph[ph$Pass %in% 21:28 & ph$Range %in% 1:2, ]
    subset_ph = subset_ph[substr(subset_ph$Type,1,1)=="P",]
  }
  # 运行程序并生成 outblue
  outblue.tmp <-  programRun(dat = subset_ph, traits)
  
  # 将当前的 outblue 添加到 combined_outblue 列表中
  combined_outblue = c(combined_outblue, list(outblue.tmp))
  
}

idx=merge.fun(combined_outblue,traits)

count=analysis(idx,traits)
count$Type="P"
write.table(count,"Output/count.csv",sep=",",row.names=FALSE,col.names=FALSE,append=TRUE)

#####绘图

pdf("Output/HNXM5_plot_P.pdf")

par(mfrow=c(2,3),mar = c(5, 4, 4, 2) + 0.1,cex.main = 0.8) 
type="P"
drawplot(idx,traits)
dev.off()

outblue_name <- paste0("Output/23HNXM5_P_datalist_bluelist.rds")

saveRDS(idx,file = outblue_name)



####HNXM5_TC2####
traits=c("YLD14",	"MST","PLOTWT","TWT","PHT","EHT")#

combined_outblue = list()

# 循环遍历不同的 type
for (i in 1:3) {
  # 根据type选择特定的数据子集
  if (i == 1) {
    
    type="Field1_P1"
    subset_ph = ph[ph$Pass %in% 1:16 & ph$Range %in% 1:6, ]
    subset_ph = subset_ph[substr(subset_ph$Type,1,3)=="TC2",]
  } else if(i==2){
    type="Field2_P1"
    subset_ph = ph[ph$Pass %in% 17:20 & ph$Range %in% 1:4, ]
    subset_ph = subset_ph[substr(subset_ph$Type,1,3)=="TC2",]
  } else if(i==3){
    type="Field3_P1"
    subset_ph = ph[ph$Pass %in% 21:28 & ph$Range %in% 1:2, ]
    subset_ph = subset_ph[substr(subset_ph$Type,1,3)=="TC2",]
  }
  # 运行程序并生成 outblue
  outblue.tmp <-  programRun(dat = subset_ph, traits)
  
  # 将当前的 outblue 添加到 combined_outblue 列表中
  combined_outblue = c(combined_outblue, list(outblue.tmp))
  
}

idx=merge.fun(combined_outblue,traits)

count=analysis(idx,traits)
#count$Type="P"
write.table(count,"Output/count.csv",sep=",",row.names=FALSE,col.names=FALSE,append=TRUE)

#####绘图

pdf("Output/HNXM5_plot_TC2.pdf")

par(mfrow=c(2,3),mar = c(5, 4, 4, 2) + 0.1,cex.main = 0.8) 
type="TC2"
drawplot(idx,traits)
dev.off()

outblue_name <- paste0("Output/23HNXM5_TC2_datalist_bluelist.rds")

saveRDS(idx,file = outblue_name)








# SDCX --------------------------------------------------------------------

cau.data <- phenos %>% dplyr::filter(Location %in% "SDCX")
cau.data$AOA="MCSU"
cau.data$Year=2023


target_cols <- c( "HarvNote")#"苗期综合评价",,"收获前综合评价"

# 搜索包含字符的行，并返回逻辑向量
has_char <- apply(cau.data[target_cols], 1, function(x) any(grepl("[[:alnum:][:punct:]]+", x)))


discard=cau.data[has_char,]
discard[,c("HarvNote")]
#discard[,c(27,28,87)]
ph_filtered=cau.data[!has_char,]
nrow(ph_filtered[ph_filtered$PlotDiscarded %in% c("Yes"), ])

ph_filtered <- ph_filtered[!(ph_filtered$PlotDiscarded %in% c("Yes")), ]
names(ph_filtered)
ph=ph_filtered

####SDCX_P####
traits=c("YLD14",	"MST","PLOTWT","PHT","EHT")#,"TWT"

combined_outblue = list()

# 循环遍历不同的 type
for (i in 1:5) {
  # 根据type选择特定的数据子集
  if (i == 1) {
    
    type="Field1_P3_P4"
    subset_ph = ph[ph$Pass %in% 1:23 & ph$Range %in% 8:34, ]
  } else if (i==2){
    type="Field1_P1"
    subset_ph = ph[ph$Pass %in% 1:24 & ph$Range %in% 35:36, ]
  } else if (i==3){
    type="Field2_P1"
    subset_ph = ph[ph$Pass %in% 25:48 & ph$Range %in% 8:36, ]
  } else if (i==4){
    type="Field3_4_5_P1"
    subset_ph = ph[ph$Pass %in% 49:120 & ph$Range %in% 2:36, ]
  } else if (i==5){
    type="Field6_P1"
    subset_ph = ph[ph$Pass %in% 121:144 & ph$Range %in% 2:20, ]
  }
  # 运行程序并生成 outblue
  outblue.tmp <-  programRun(dat = subset_ph, traits)
  
  # 将当前的 outblue 添加到 combined_outblue 列表中
  combined_outblue = c(combined_outblue, list(outblue.tmp))
  
}

idx=merge.fun(combined_outblue,traits)

count=analysis(idx,traits)
count$Type="P"
write.table(count,"Output/count.csv",sep=",",row.names=FALSE,col.names=FALSE,append=TRUE)

#####绘图

pdf("Output/SDCX_plot_P.pdf")

par(mfrow=c(2,3),mar = c(5, 4, 4, 2) + 0.1,cex.main = 0.8) 
type="P"
drawplot(idx,traits)
dev.off()

outblue_name <- paste0("Output/23SDCX_P_datalist_bluelist.rds")

saveRDS(idx,file = outblue_name)



####SDCX_TC2####
traits=c("YLD14",	"MST","PLOTWT","PHT","EHT")#,"TWT"

combined_outblue = list()

# 循环遍历不同的 type
for (i in 1:3) {
  # 根据type选择特定的数据子集
  if (i == 1) {
    
    type="Field6_TC2"
    subset_ph = ph[ph$Pass %in% 121:144 & ph$Range %in% 21:33, ]
  } else if (i==2){
    type="Field7_8_TC2"
    subset_ph = ph[ph$Pass %in% 145:192 & ph$Range %in% 2:32, ]
  } else if (i==3){
    type="Field9_TC2"
    subset_ph = ph[ph$Pass %in% 193:216 & ph$Range %in% 2:15, ]
  } 
  # 运行程序并生成 outblue
  outblue.tmp <-  programRun(dat = subset_ph, traits)
  
  # 将当前的 outblue 添加到 combined_outblue 列表中
  combined_outblue = c(combined_outblue, list(outblue.tmp))
  
}

idx=merge.fun(combined_outblue,traits)

count=analysis(idx,traits)
#count$Type="P"
write.table(count,"Output/count.csv",sep=",",row.names=FALSE,col.names=FALSE,append=TRUE)

#####绘图

pdf("Output/SDCX_plot_TC2.pdf")

par(mfrow=c(2,3),mar = c(5, 4, 4, 2) + 0.1,cex.main = 0.8) 
type="TC2"
drawplot(idx,traits)
dev.off()

outblue_name <- paste0("Output/23SDCX_TC2_datalist_bluelist.rds")

saveRDS(idx,file = outblue_name)



# HNYJ --------------------------------------------------------------------
cau.data <- phenos %>% dplyr::filter(Location %in% "HNYJ")
cau.data$AOA="MCSU"
cau.data$Year=2023


target_cols <- c("HarvNote")#"收获前综合评价", ,"苗期综合评价"

# 搜索包含字符的行，并返回逻辑向量
has_char <- apply(cau.data[target_cols], 1, function(x) any(grepl("[[:alnum:][:punct:]]+", x)))


discard=cau.data[has_char,]
discard[,c("HarvNote")]

ph_filtered=cau.data[!has_char,]
N_discard=nrow(ph_filtered[ph_filtered$PlotDiscarded %in% c("Yes"), ])
N_discard
ph_filtered <- ph_filtered[!(ph_filtered$PlotDiscarded %in% c("Yes")), ]
names(ph_filtered)
ph=ph_filtered


####HNYJ_P####
traits=c("YLD14",	"MST","PLOTWT")#	,"TWT","PHT","EHT"

combined_outblue <- list()

# 循环遍历不同的 type
for (i in 1:4) {
  # 根据type选择特定的数据子集
  if (i == 1) {
    
    type="Field2_P1"
    subset_ph = ph[ph$Pass %in% 25:48 & ph$Range %in% 37:54, ]
  } else if (i == 2) {
    type="Field3_P1"
    subset_ph = ph[ph$Pass %in% 49:72 & ph$Range %in% 2:54, ]
  } else if (i == 3) {
    type="Field4_P1"
    subset_ph = ph[ph$Pass %in% 73:96 & ph$Range %in% 2:45, ]
  } else if (i == 4) {
    type="Field5_P1"
    subset_ph = ph[ph$Pass %in% 97:108 & ph$Range %in% 3:52, ]
  } 
  
  # 运行程序并生成 outblue
  outblue.tmp = programRun(dat = subset_ph, traits)
  
  # 将当前的 outblue 添加到 combined_outblue 列表中
  combined_outblue = append(combined_outblue, list(outblue.tmp))
  
}

idx=merge.fun(combined_outblue,traits)

count=analysis(idx,traits)
count$Type="P"
write.table(count,"Output/count.csv",sep=",",row.names=FALSE,col.names=FALSE,append=TRUE)
#####绘图

pdf("Output/HNYJ_plot_P.pdf")

par(mfrow=c(2,3),mar = c(5, 4, 4, 2) + 0.1,cex.main = 0.8) 
type="P"
drawplot(idx,traits)
dev.off()

outblue_name <- paste0("Output/23HNYJ_P_datalist_bluelist.rds")

saveRDS(idx,file = outblue_name )


####HNYJ_TC2####
traits=c("YLD14",	"MST","PLOTWT")#	,"TWT","PHT","EHT"

combined_outblue <- list()

# 循环遍历不同的 type
for (i in 1:2) {
  # 根据type选择特定的数据子集
  if (i == 1) {
    
    type="Field1_TC2"
    subset_ph = ph[ph$Pass %in% 1:24 & ph$Range %in% 2:54, ]
  } else if (i == 2) {
    type="Field2_TC2"
    subset_ph = ph[ph$Pass %in% 25:48 & ph$Range %in% 2:36, ]
  } 
  # 运行程序并生成 outblue
  outblue.tmp = programRun(dat = subset_ph, traits)
  
  # 将当前的 outblue 添加到 combined_outblue 列表中
  combined_outblue = append(combined_outblue, list(outblue.tmp))
  
}

idx=merge.fun(combined_outblue,traits)

count=analysis(idx,traits)
#count$Type="P"
write.table(count,"Output/count.csv",sep=",",row.names=FALSE,col.names=FALSE,append=TRUE)
#####绘图

pdf("Output/HNYJ_plot_TC2.pdf")

par(mfrow=c(2,3),mar = c(5, 4, 4, 2) + 0.1,cex.main = 0.8) 
type="TC2"
drawplot(idx,traits)
dev.off()

outblue_name <- paste0("Output/23HNYJ_TC2_datalist_bluelist.rds")

saveRDS(idx,file = outblue_name )



# HNXP --------------------------------------------------------------------
cau.data <- phenos %>% dplyr::filter(Location %in% "HNXP")
cau.data$AOA="SCSU"
cau.data$Year=2023


target_cols <- c("HarvNote")#"收获前综合评价", ,"苗期综合评价"

# 搜索包含字符的行，并返回逻辑向量
has_char <- apply(cau.data[target_cols], 1, function(x) any(grepl("[[:alnum:][:punct:]]+", x)))


discard=cau.data[has_char,]
discard[,c("HarvNote")]

ph_filtered=cau.data[!has_char,]
N_discard=nrow(ph_filtered[ph_filtered$PlotDiscarded %in% c("Yes"), ])
N_discard
ph_filtered <- ph_filtered[!(ph_filtered$PlotDiscarded %in% c("Yes")), ]
names(ph_filtered)
ph=ph_filtered

####HNXP_TC1####
traits=c("YLD14",	"MST","PLOTWT","TWT")#	,"PHT","EHT"

combined_outblue <- list()

# 循环遍历不同的 type
for (i in 1) {
  # 根据type选择特定的数据子集
  if (i == 1) {
    
    type="Field_TC1"
    subset_ph = ph[ph$Pass %in% 1:168 & ph$Range %in% 2:32, ]
  } 
  # 运行程序并生成 outblue
  outblue.tmp = programRun(dat = subset_ph, traits)
  
  # 将当前的 outblue 添加到 combined_outblue 列表中
  combined_outblue = append(combined_outblue, list(outblue.tmp))
  
}

idx=merge.fun(combined_outblue,traits)

count=analysis(idx,traits)
#count$Type="P"
write.table(count,"Output/count.csv",sep=",",row.names=FALSE,col.names=FALSE,append=TRUE)
#####绘图

pdf("Output/HNXP_plot_TC1.pdf")

par(mfrow=c(2,3),mar = c(5, 4, 4, 2) + 0.1,cex.main = 0.8) 
type="TC1"
drawplot(idx,traits)
dev.off()

outblue_name <- paste0("Output/23HNXP_TC1_datalist_bluelist.rds")

saveRDS(idx,file = outblue_name )




# HBCX_P--------------------------------------------------------------------
cau.data <- phenos %>% dplyr::filter(Location %in% "HBCX")
cau.data$AOA="MCSU"
cau.data$Year=2023


target_cols <- c("HarvNote")#"收获前综合评价", ,"苗期综合评价"

# 搜索包含字符的行，并返回逻辑向量
has_char <- apply(cau.data[target_cols], 1, function(x) any(grepl("[[:alnum:][:punct:]]+", x)))


discard=cau.data[has_char,]
discard[,c("HarvNote")]

ph_filtered=cau.data[!has_char,]
N_discard=nrow(ph_filtered[ph_filtered$PlotDiscarded %in% c("Yes"), ])
N_discard
ph_filtered <- ph_filtered[!(ph_filtered$PlotDiscarded %in% c("Yes")), ]
names(ph_filtered)
ph=ph_filtered


traits=c("YLD14",	"MST","PLOTWT","PHT","EHT")#	,"TWT"

combined_outblue <- list()

# 循环遍历不同的 type
for (i in 1:4) {
  # 根据type选择特定的数据子集
  if (i == 1) {
    
    type="Field1_P3_P4"
    subset_ph = ph[ph$Pass %in% 1:23 & ph$Range %in% 2:27, ]
  } else if(i==2){
    type="Field2_P3_P4"
    subset_ph = ph[ph$Pass %in% 25:47 & ph$Range %in% 2:6, ]
  } else if(i==3){
    type="Field2_P1"
    subset_ph = ph[ph$Pass %in% 25:48 & ph$Range %in% 7:27, ]
  } else if(i==4){
    type="Field3_4_5_6_7_P1"
    subset_ph = ph[ph$Pass %in% 49:168 & ph$Range %in% 2:27, ]
  }
  # 运行程序并生成 outblue
  outblue.tmp = programRun(dat = subset_ph, traits)
  
  # 将当前的 outblue 添加到 combined_outblue 列表中
  combined_outblue = append(combined_outblue, list(outblue.tmp))
  
}

idx=merge.fun(combined_outblue,traits)

count=analysis(idx,traits)
count$Type="P"
write.table(count,"Output/count.csv",sep=",",row.names=FALSE,col.names=FALSE,append=TRUE)
#####绘图

pdf("Output/HNCX_plot_P.pdf")

par(mfrow=c(2,3),mar = c(5, 4, 4, 2) + 0.1,cex.main = 0.8) 
type="P"
drawplot(idx,traits)
dev.off()

outblue_name <- paste0("Output/23HNCX_P_datalist_bluelist.rds")

saveRDS(idx,file = outblue_name )









# SDCQ1_P --------------------------------------------------------------------

cau.data <- phenos %>% dplyr::filter(Location %in% "SDCQ1")
cau.data$AOA="MCSU"
cau.data$Year=2023


target_cols <- c( "HarvNote")#"苗期综合评价",,"收获前综合评价"

# 搜索包含字符的行，并返回逻辑向量
has_char <- apply(cau.data[target_cols], 1, function(x) any(grepl("[[:alnum:][:punct:]]+", x)))


discard=cau.data[has_char,]
discard[,c("HarvNote")]
#discard[,c(27,28,87)]
ph_filtered=cau.data[!has_char,]
N_discard=nrow(ph_filtered[ph_filtered$PlotDiscarded %in% c("Yes"), ])

ph_filtered <- ph_filtered[!(ph_filtered$PlotDiscarded %in% c("Yes")), ]
names(ph_filtered)
ph=ph_filtered


traits=c("YLD14",	"MST","PLOTWT","PHT","EHT")#,"TWT"

combined_outblue = list()

# 循环遍历不同的 type
for (i in 1:3) {
  # 根据type选择特定的数据子集
  if (i == 1) {
    
    type="Field1_P3_P4"
    subset_ph = ph[ph$Pass %in% 1:23 & ph$Range %in% 2:32, ]
  } else if(i==2){
    type="Field1_P1"
    subset_ph = ph[ph$Pass %in% 1:24 & ph$Range %in% 33:34, ]
  } else if(i==3){
    type="Field2_3_4_P1"
    subset_ph = ph[ph$Pass %in% 25:80 & ph$Range %in% 2:34, ]
  } 
  # 运行程序并生成 outblue
  outblue.tmp <-  programRun(dat = subset_ph, traits)
  
  # 将当前的 outblue 添加到 combined_outblue 列表中
  combined_outblue = c(combined_outblue, list(outblue.tmp))
  
}

idx=merge.fun(combined_outblue,traits)

count=analysis(idx,traits)
count$Type="P"
write.table(count,"Output/count.csv",sep=",",row.names=FALSE,col.names=FALSE,append=TRUE)

#####绘图

pdf("Output/SDCQ1_plot_P.pdf")

par(mfrow=c(2,3),mar = c(5, 4, 4, 2) + 0.1,cex.main = 0.8) 
type="P"
drawplot(idx,traits)
dev.off()

outblue_name <- paste0("Output/23SDCQ1_P_datalist_bluelist.rds")

saveRDS(idx,file = outblue_name)


# SDCQ2_P --------------------------------------------------------------------

cau.data <- phenos %>% dplyr::filter(Location %in% "SDCQ2")
cau.data$AOA="MCSU"
cau.data$Year=2023


target_cols <- c( "HarvNote")#"苗期综合评价",,"收获前综合评价"

# 搜索包含字符的行，并返回逻辑向量
has_char <- apply(cau.data[target_cols], 1, function(x) any(grepl("[[:alnum:][:punct:]]+", x)))


discard=cau.data[has_char,]
discard[,c("HarvNote")]
#discard[,c(27,28,87)]
ph_filtered=cau.data[!has_char,]
N_discard=nrow(ph_filtered[ph_filtered$PlotDiscarded %in% c("Yes"), ])

ph_filtered <- ph_filtered[!(ph_filtered$PlotDiscarded %in% c("Yes")), ]
names(ph_filtered)
ph=ph_filtered


traits=c("YLD14",	"MST","PLOTWT","PHT","EHT")#,"TWT"

combined_outblue = list()

# 循环遍历不同的 type
for (i in 1) {
  # 根据type选择特定的数据子集
  if (i == 1) {
    
    type="Field1_2_3_P1"
    subset_ph = ph[ph$Pass %in% 1:50 & ph$Range %in% 2:36, ]
  } 
  # 运行程序并生成 outblue
  outblue.tmp <-  programRun(dat = subset_ph, traits)
  
  # 将当前的 outblue 添加到 combined_outblue 列表中
  combined_outblue = c(combined_outblue, list(outblue.tmp))
  
}

idx=merge.fun(combined_outblue,traits)

count=analysis(idx,traits)
count$Type="P"
write.table(count,"Output/count.csv",sep=",",row.names=FALSE,col.names=FALSE,append=TRUE)

#####绘图

pdf("Output/SDCQ2_plot_P.pdf")

par(mfrow=c(2,3),mar = c(5, 4, 4, 2) + 0.1,cex.main = 0.8) 
type="P"
drawplot(idx,traits)
dev.off()

outblue_name <- paste0("Output/23SDCQ2_P_datalist_bluelist.rds")

saveRDS(idx,file = outblue_name)








# HNSP --------------------------------------------------------------------

cau.data <- phenos %>% dplyr::filter(Location %in% "HNSP")
cau.data$AOA="SCSU"
cau.data$Year=2023


target_cols <- c( "HarvNote")#"苗期综合评价",,"收获前综合评价"

# 搜索包含字符的行，并返回逻辑向量
has_char <- apply(cau.data[target_cols], 1, function(x) any(grepl("[[:alnum:][:punct:]]+", x)))


discard=cau.data[has_char,]
discard[,c("HarvNote")]
#discard[,c(27,28,87)]
ph_filtered=cau.data[!has_char,]
N_discard=nrow(ph_filtered[ph_filtered$PlotDiscarded %in% c("Yes"), ])

ph_filtered <- ph_filtered[!(ph_filtered$PlotDiscarded %in% c("Yes")), ]
names(ph_filtered)
ph=ph_filtered

####HNSP_P####
traits=c("YLD14",	"MST","PLOTWT","TWT")#,"PHT","EHT"

combined_outblue = list()

# 循环遍历不同的 type
for (i in 1:3) {
  # 根据type选择特定的数据子集
  if (i == 1) {
    
    type="Field2_3_4_5_P1"
    subset_ph = ph[ph$Pass %in% 25:120 & ph$Range %in% 4:29, ]
  } else if(i==2) {
    type="Field6_P1"
    subset_ph = ph[ph$Pass %in% 121:144 & ph$Range %in% 3:28, ]
  } else if(i==3) {
    type="Field7_P1"
    subset_ph = ph[ph$Pass %in% 145:168 & ph$Range %in% 3:20, ]
  }
  # 运行程序并生成 outblue
  outblue.tmp <-  programRun(dat = subset_ph, traits)
  
  # 将当前的 outblue 添加到 combined_outblue 列表中
  combined_outblue = c(combined_outblue, list(outblue.tmp))
  
}

idx=merge.fun(combined_outblue,traits)

count=analysis(idx,traits)
count$Type="P"
write.table(count,"Output/count.csv",sep=",",row.names=FALSE,col.names=FALSE,append=TRUE)

#####绘图

pdf("Output/HNSP_plot_P.pdf")

par(mfrow=c(2,3),mar = c(5, 4, 4, 2) + 0.1,cex.main = 0.8) 
type="P"
drawplot(idx,traits)
dev.off()

outblue_name <- paste0("Output/23HNSP_P_datalist_bluelist.rds")

saveRDS(idx,file = outblue_name)









####HNSP_TC2####
traits=c("YLD14",	"MST","PLOTWT","TWT")#,"PHT","EHT"

combined_outblue = list()

# 循环遍历不同的 type
for (i in 1:4) {
  # 根据type选择特定的数据子集
  if (i == 1) {
    
    type="Field7_TC2"
    subset_ph = ph[ph$Pass %in% 145:168 & ph$Range %in% 21:28, ]
  } else if(i==2) {
    type="Field8_9_TC2"
    subset_ph = ph[ph$Pass %in% 169:216 & ph$Range %in% 3:28, ]
  } else if(i==3) {
    type="Field10_TC2"
    subset_ph = ph[ph$Pass %in% 217:240 & ph$Range %in% 2:27, ]
  } else if(i==4) {
    type="Field11_TC2"
    subset_ph = ph[ph$Pass %in% 241:246 & ph$Range %in% 2:21, ]
  }
  # 运行程序并生成 outblue
  outblue.tmp <-  programRun(dat = subset_ph, traits)
  
  # 将当前的 outblue 添加到 combined_outblue 列表中
  combined_outblue = c(combined_outblue, list(outblue.tmp))
  
}

idx=merge.fun(combined_outblue,traits)

count=analysis(idx,traits)
#count$Type="P"
write.table(count,"Output/count.csv",sep=",",row.names=FALSE,col.names=FALSE,append=TRUE)

#####绘图

pdf("Output/HNSP_plot_TC2.pdf")

par(mfrow=c(2,3),mar = c(5, 4, 4, 2) + 0.1,cex.main = 0.8) 
type="TC2"
drawplot(idx,traits)
dev.off()

outblue_name <- paste0("Output/23HNSP_TC2_datalist_bluelist.rds")

saveRDS(idx,file = outblue_name)











# SDST_P --------------------------------------------------------------------

cau.data <- phenos %>% dplyr::filter(Location %in% "SDST")
cau.data$AOA="MCSU"
cau.data$Year=2023


target_cols <- c( "HarvNote")#"苗期综合评价",,"收获前综合评价"

# 搜索包含字符的行，并返回逻辑向量
has_char <- apply(cau.data[target_cols], 1, function(x) any(grepl("[[:alnum:][:punct:]]+", x)))


discard=cau.data[has_char,]
discard[,c("HarvNote")]
#discard[,c(27,28,87)]
ph_filtered=cau.data[!has_char,]
N_discard=nrow(ph_filtered[ph_filtered$PlotDiscarded %in% c("Yes"), ])

ph_filtered <- ph_filtered[!(ph_filtered$PlotDiscarded %in% c("Yes")), ]
names(ph_filtered)
ph=ph_filtered

traits=c("YLD14",	"MST","PLOTWT")#,"PHT","EHT","TWT"

combined_outblue = list()

# 循环遍历不同的 type
for (i in 1:5) {
  # 根据type选择特定的数据子集
  if (i == 1) {
    
    type="Field1_P3_P4"
    subset_ph = ph[ph$Pass %in% 1:23 & ph$Range %in% 2:28, ]
  } else if(i==2) {
    type="Field1_P1"
    subset_ph = ph[ph$Pass %in% 1:24 & ph$Range %in% 29:44, ]
  } else if(i==3) {
    type="Field2_3_4_P1"
    subset_ph = ph[ph$Pass %in% 25:76 & ph$Range %in% 2:44, ]
  } else if(i==4) {
    type="Field5_P1"
    subset_ph = ph[ph$Pass %in% 83:96 & ph$Range %in% 2:44, ]
  } else if(i==5) {
    type="Field6_P1"
    subset_ph = ph[ph$Pass %in% 97:108 & ph$Range %in% 2:34, ]
  }
  # 运行程序并生成 outblue
  outblue.tmp <-  programRun(dat = subset_ph, traits)
  
  # 将当前的 outblue 添加到 combined_outblue 列表中
  combined_outblue = c(combined_outblue, list(outblue.tmp))
  
}

idx=merge.fun(combined_outblue,traits)

count=analysis(idx,traits)
count$Type="P"
write.table(count,"Output/count.csv",sep=",",row.names=FALSE,col.names=FALSE,append=TRUE)

#####绘图

pdf("Output/SDST_plot_P.pdf")

par(mfrow=c(2,3),mar = c(5, 4, 4, 2) + 0.1,cex.main = 0.8) 
type="P"
drawplot(idx,traits)
dev.off()

outblue_name <- paste0("Output/23SDST_P_datalist_bluelist.rds")

saveRDS(idx,file = outblue_name)












# HNXH_P --------------------------------------------------------------------

cau.data <- phenos %>% dplyr::filter(Location %in% "HNXH")
cau.data$AOA="SCSU"
cau.data$Year=2023


target_cols <- c( "HarvNote")#"苗期综合评价",,"收获前综合评价"

# 搜索包含字符的行，并返回逻辑向量
has_char <- apply(cau.data[target_cols], 1, function(x) any(grepl("[[:alnum:][:punct:]]+", x)))


discard=cau.data[has_char,]
discard[,c("HarvNote")]
#discard[,c(27,28,87)]
ph_filtered=cau.data[!has_char,]
N_discard=nrow(ph_filtered[ph_filtered$PlotDiscarded %in% c("Yes"), ])

ph_filtered <- ph_filtered[!(ph_filtered$PlotDiscarded %in% c("Yes")), ]
names(ph_filtered)
ph=ph_filtered

traits=c("YLD14",	"MST","PLOTWT")#,"PHT","EHT","TWT"

combined_outblue = list()

# 循环遍历不同的 type
for (i in 1:3) {
  # 根据type选择特定的数据子集
  if (i == 1) {
    
    type="Field1_P3_P4"
    subset_ph = ph[ph$Pass %in% 1:23 & ph$Range %in% 2:28, ]
  } else if(i==2) {
    type="Field1_P1"
    subset_ph = ph[ph$Pass %in% 1:24 & ph$Range %in% 29:67, ]
  } else if(i==3) {
    type="Field2_3_4_P1"
    subset_ph = ph[ph$Pass %in% 25:64 & ph$Range %in% 2:67, ]
  } 
  # 运行程序并生成 outblue
  outblue.tmp <-  programRun(dat = subset_ph, traits)
  
  # 将当前的 outblue 添加到 combined_outblue 列表中
  combined_outblue = c(combined_outblue, list(outblue.tmp))
  
}

idx=merge.fun(combined_outblue,traits)

count=analysis(idx,traits)
count$Type="P"
write.table(count,"Output/count.csv",sep=",",row.names=FALSE,col.names=FALSE,append=TRUE)

#####绘图

pdf("Output/HNXH_plot_P.pdf")

par(mfrow=c(2,3),mar = c(5, 4, 4, 2) + 0.1,cex.main = 0.8) 
type="P"
drawplot(idx,traits)
dev.off()

outblue_name <- paste0("Output/23HNXH_P_datalist_bluelist.rds")

saveRDS(idx,file = outblue_name)











# SDZQ1 --------------------------------------------------------------------

cau.data <- phenos %>% dplyr::filter(Location %in% "SDZQ1")
cau.data$AOA="MCSU"
cau.data$Year=2023


target_cols <- c( "HarvNote")#"苗期综合评价",,"收获前综合评价"

# 搜索包含字符的行，并返回逻辑向量
has_char <- apply(cau.data[target_cols], 1, function(x) any(grepl("[[:alnum:][:punct:]]+", x)))


discard=cau.data[has_char,]
discard[,c("HarvNote")]
#discard[,c(27,28,87)]
ph_filtered=cau.data[!has_char,]
N_discard=nrow(ph_filtered[ph_filtered$PlotDiscarded %in% c("Yes"), ])

ph_filtered <- ph_filtered[!(ph_filtered$PlotDiscarded %in% c("Yes")), ]
names(ph_filtered)
ph=ph_filtered


####SDZQ1_TC2####
traits=c("YLD14",	"MST","PLOTWT","PHT","EHT")#,"TWT"

combined_outblue = list()

# 循环遍历不同的 type
for (i in 1:3) {
  # 根据type选择特定的数据子集
  if (i == 1) {
    
    type="Field5_TC2"
    subset_ph = ph[ph$Pass %in% 97:120 & ph$Range %in% 33:51, ]
  } else if(i==2) {
    type="Field6_7_TC2"
    subset_ph = ph[ph$Pass %in% 121:168 & ph$Range %in% 2:51, ]
  } else if(i==3) {
    type="Field8_TC2"
    subset_ph = ph[ph$Pass %in% 169:192 & ph$Range %in% 2:25, ]
  } else if(i==4) {
    type="Field12_TC2"
    subset_ph = ph[ph$Pass %in% 249:253 & ph$Range %in% 2:51, ]
    subset_ph = subset_ph[substr(subset_ph$Name,1,7)=="23TC2XM",]
  } 
  
  
  # 运行程序并生成 outblue
  outblue.tmp <-  programRun(dat = subset_ph, traits)
  
  # 将当前的 outblue 添加到 combined_outblue 列表中
  combined_outblue = c(combined_outblue, list(outblue.tmp))
  
}

idx=merge.fun(combined_outblue,traits)

count=analysis(idx,traits)
#count$Type="P"
write.table(count,"Output/count.csv",sep=",",row.names=FALSE,col.names=FALSE,append=TRUE)

#####绘图

pdf("Output/SDZQ1_plot_TC2.pdf")

par(mfrow=c(2,3),mar = c(5, 4, 4, 2) + 0.1,cex.main = 0.8) 
type="TC2"
drawplot(idx,traits)
dev.off()

outblue_name <- paste0("Output/23SDZQ1_TC2_datalist_bluelist.rds")

saveRDS(idx,file = outblue_name)


####SDZQ1_P####
traits=c("YLD14",	"MST","PLOTWT","PHT","EHT")#,"TWT"

combined_outblue = list()

# 循环遍历不同的 type
for (i in 1:6) {
  # 根据type选择特定的数据子集
  if (i == 1) {
    
    type="Field1_P3_P4"
    subset_ph = ph[ph$Pass %in% 1:23 & ph$Range %in% 2:32, ]
  } else if(i==2) {
    type="Field1_P1"
    subset_ph = ph[ph$Pass %in% 1:24 & ph$Range %in% 33:47, ]
  } else if(i==3) {
    type="Field2_P1"
    subset_ph = ph[ph$Pass %in% 25:48 & ph$Range %in% 2:47, ]
  } else if(i==4) {
    type="Field3_4_P1"
    subset_ph = ph[ph$Pass %in% 49:96 & ph$Range %in% 2:51, ]
  } else if(i==5) {
    type="Field5_P1"
    subset_ph = ph[ph$Pass %in% 97:120 & ph$Range %in% 2:32, ]
  } else if(i==6) {
    type="Field12_P1"
    subset_ph = ph[ph$Pass %in% 249:253 & ph$Range %in% 2:51, ]
    subset_ph = subset_ph[substr(subset_ph$Type,1,1)=="P",]
  } 
  
  
  # 运行程序并生成 outblue
  outblue.tmp <-  programRun(dat = subset_ph, traits)
  
  # 将当前的 outblue 添加到 combined_outblue 列表中
  combined_outblue = c(combined_outblue, list(outblue.tmp))
  
}

idx=merge.fun(combined_outblue,traits)

count=analysis(idx,traits)
count$Type="P"
write.table(count,"Output/count.csv",sep=",",row.names=FALSE,col.names=FALSE,append=TRUE)

#####绘图

pdf("Output/SDZQ1_plot_P.pdf")

par(mfrow=c(2,3),mar = c(5, 4, 4, 2) + 0.1,cex.main = 0.8) 
type="P"
drawplot(idx,traits)
dev.off()

outblue_name <- paste0("Output/23SDZQ1_P_datalist_bluelist.rds")

saveRDS(idx,file = outblue_name)



####SDZQ1_TC1####
traits=c("YLD14",	"MST","PLOTWT","PHT","EHT")#,"TWT"

combined_outblue = list()

# 循环遍历不同的 type
for (i in 1:3) {
  # 根据type选择特定的数据子集
  if (i == 1) {
    
    type="Field8_TC1"
    subset_ph = ph[ph$Pass %in% 169:192 & ph$Range %in% 26:51, ]
  } else if(i==2) {
    type="Field9_10_11_TC1"
    subset_ph = ph[ph$Pass %in% 193:248 & ph$Range %in% 2:51, ]
  } else if(i==3) {
    type="Field12_TC1"
    subset_ph = ph[ph$Pass %in% 249:253 & ph$Range %in% 2:51, ]
    subset_ph = subset_ph[substr(subset_ph$Name,1,7)=="23TC1XM",]
  } 
  
  
  # 运行程序并生成 outblue
  outblue.tmp <-  programRun(dat = subset_ph, traits)
  
  # 将当前的 outblue 添加到 combined_outblue 列表中
  combined_outblue = c(combined_outblue, list(outblue.tmp))
  
}

idx=merge.fun(combined_outblue,traits)

count=analysis(idx,traits)
#count$Type="P"
write.table(count,"Output/count.csv",sep=",",row.names=FALSE,col.names=FALSE,append=TRUE)

#####绘图

pdf("Output/SDZQ1_plot_TC1.pdf")

par(mfrow=c(2,3),mar = c(5, 4, 4, 2) + 0.1,cex.main = 0.8) 
type="TC1"
drawplot(idx,traits)
dev.off()

outblue_name <- paste0("Output/23SDZQ1_TC1_datalist_bluelist.rds")

saveRDS(idx,file = outblue_name)

# SDZQ2_TC1 --------------------------------------------------------------------

cau.data <- phenos %>% dplyr::filter(Location %in% "SDZQ2")
cau.data$AOA="MCSU"
cau.data$Year=2023


target_cols <- c( "HarvNote")#"苗期综合评价",,"收获前综合评价"

# 搜索包含字符的行，并返回逻辑向量
has_char <- apply(cau.data[target_cols], 1, function(x) any(grepl("[[:alnum:][:punct:]]+", x)))


discard=cau.data[has_char,]
discard[,c("HarvNote")]
#discard[,c(27,28,87)]
ph_filtered=cau.data[!has_char,]
N_discard=nrow(ph_filtered[ph_filtered$PlotDiscarded %in% c("Yes"), ])

ph_filtered <- ph_filtered[!(ph_filtered$PlotDiscarded %in% c("Yes")), ]
names(ph_filtered)
ph=ph_filtered




traits=c("YLD14",	"MST","PLOTWT","PHT","EHT")#,"TWT"

combined_outblue = list()

# 循环遍历不同的 type
for (i in 1) {
  # 根据type选择特定的数据子集
  if (i == 1) {
    
    type="Field1_2_3_TC1"
    subset_ph = ph[ph$Pass %in% 1:64 & ph$Range %in% 2:50, ]
  } 
  
  
  # 运行程序并生成 outblue
  outblue.tmp <-  programRun(dat = subset_ph, traits)
  
  # 将当前的 outblue 添加到 combined_outblue 列表中
  combined_outblue = c(combined_outblue, list(outblue.tmp))
  
}

idx=merge.fun(combined_outblue,traits)

count=analysis(idx,traits)
#count$Type="P"
write.table(count,"Output/count.csv",sep=",",row.names=FALSE,col.names=FALSE,append=TRUE)

#####绘图

pdf("Output/SDZQ2_plot_TC1.pdf")

par(mfrow=c(2,3),mar = c(5, 4, 4, 2) + 0.1,cex.main = 0.8) 
type="TC1"
drawplot(idx,traits)
dev.off()

outblue_name <- paste0("Output/23SDZQ2_TC1_datalist_bluelist.rds")

saveRDS(idx,file = outblue_name)





# HBRZ_P --------------------------------------------------------------------

cau.data <- phenos %>% dplyr::filter(Location %in% "HBRZ")
cau.data$AOA="MCSU"
cau.data$Year=2023


target_cols <- c( "HarvNote")#"苗期综合评价",,"收获前综合评价"

# 搜索包含字符的行，并返回逻辑向量
has_char <- apply(cau.data[target_cols], 1, function(x) any(grepl("[[:alnum:][:punct:]]+", x)))


discard=cau.data[has_char,]
discard[,c("HarvNote")]
#discard[,c(27,28,87)]
ph_filtered=cau.data[!has_char,]
N_discard=nrow(ph_filtered[ph_filtered$PlotDiscarded %in% c("Yes"), ])

ph_filtered <- ph_filtered[!(ph_filtered$PlotDiscarded %in% c("Yes")), ]
names(ph_filtered)
ph=ph_filtered

traits=c("YLD14",	"MST","PLOTWT","PHT","EHT")#,"TWT"

combined_outblue = list()

# 循环遍历不同的 type
for (i in 1:4) {
  # 根据type选择特定的数据子集
  if (i == 1) {
    
    type="Field1_P3_P4"
    subset_ph = ph[ph$Pass %in% 1:23 & ph$Range %in% 2:32, ]
  } else if(i==2) {
    type="Field1_P1"
    subset_ph = ph[ph$Pass %in% 1:24 & ph$Range %in% 33:39, ]
  } else if(i==3) {
    type="Field2_3_4_5_P1"
    subset_ph = ph[ph$Pass %in% 25:108 & ph$Range %in% 2:39, ]
  } else if(i==4) {
    type="Field6_P1"
    subset_ph = ph[ph$Pass %in% 109:116 & ph$Range %in% 2:34, ]
  } 
  # 运行程序并生成 outblue
  outblue.tmp <-  programRun(dat = subset_ph, traits)
  
  # 将当前的 outblue 添加到 combined_outblue 列表中
  combined_outblue = c(combined_outblue, list(outblue.tmp))
  
}

idx=merge.fun(combined_outblue,traits)

count=analysis(idx,traits)
count$Type="P"
write.table(count,"Output/count.csv",sep=",",row.names=FALSE,col.names=FALSE,append=TRUE)

#####绘图

pdf("Output/HBRZ_plot_P.pdf")

par(mfrow=c(2,3),mar = c(5, 4, 4, 2) + 0.1,cex.main = 0.8) 
type="P"
drawplot(idx,traits)
dev.off()

outblue_name <- paste0("Output/23HBRZ_P_datalist_bluelist.rds")

saveRDS(idx,file = outblue_name)





# AHSX2 --------------------------------------------------------------------

cau.data <- phenos %>% dplyr::filter(Location %in% "AHSX2")
cau.data$AOA="SCSU"
cau.data$Year=2023


target_cols <- c( "HarvNote")#"苗期综合评价",,"收获前综合评价"

# 搜索包含字符的行，并返回逻辑向量
has_char <- apply(cau.data[target_cols], 1, function(x) any(grepl("[[:alnum:][:punct:]]+", x)))


discard=cau.data[has_char,]
discard[,c("HarvNote")]
#discard[,c(27,28,87)]
ph_filtered=cau.data[!has_char,]
N_discard=nrow(ph_filtered[ph_filtered$PlotDiscarded %in% c("Yes"), ])

ph_filtered <- ph_filtered[!(ph_filtered$PlotDiscarded %in% c("Yes")), ]
names(ph_filtered)
ph=ph_filtered

####AHSX2_TC2####
traits=c("YLD14",	"MST","PLOTWT")#,"TWT","PHT","EHT"

combined_outblue = list()

# 循环遍历不同的 type
for (i in 1:2) {
  # 根据type选择特定的数据子集
    if(i==1) {
    type="Field4_5_TC2"
    subset_ph = ph[ph$Pass %in% 49:84 & ph$Range %in% 4:40, ]
  } else if(i==2) {
    type="Field6_TC2"
    subset_ph = ph[ph$Pass %in% 85:88 & ph$Range %in% 4:24, ]
  } 
  # 运行程序并生成 outblue
  outblue.tmp <-  programRun(dat = subset_ph, traits)
  
  # 将当前的 outblue 添加到 combined_outblue 列表中
  combined_outblue = c(combined_outblue, list(outblue.tmp))
  
}

idx=merge.fun(combined_outblue,traits)

count=analysis(idx,traits)
count$Type="TC2"
write.table(count,"Output/count.csv",sep=",",row.names=FALSE,col.names=FALSE,append=TRUE)

#####绘图

pdf("Output/AHSX2_plot_TC2.pdf")

par(mfrow=c(2,3),mar = c(5, 4, 4, 2) + 0.1,cex.main = 0.8) 
type="TC2"
drawplot(idx,traits)
dev.off()

outblue_name <- paste0("Output/23AHSX2_TC2_datalist_bluelist.rds")

saveRDS(idx,file = outblue_name)




# SDGX1_P --------------------------------------------------------------------

cau.data <- phenos %>% dplyr::filter(Location %in% "SDGX1")
cau.data$AOA="MCSU"
cau.data$Year=2023


target_cols <- c( "HarvNote")#"苗期综合评价",,"收获前综合评价"

# 搜索包含字符的行，并返回逻辑向量
has_char <- apply(cau.data[target_cols], 1, function(x) any(grepl("[[:alnum:][:punct:]]+", x)))


discard=cau.data[has_char,]
discard[,c("HarvNote")]
#discard[,c(27,28,87)]
ph_filtered=cau.data[!has_char,]
N_discard=nrow(ph_filtered[ph_filtered$PlotDiscarded %in% c("Yes"), ])

ph_filtered <- ph_filtered[!(ph_filtered$PlotDiscarded %in% c("Yes")), ]
names(ph_filtered)
ph=ph_filtered


traits=c("YLD14",	"MST","PLOTWT")#,"TWT","PHT","EHT"

combined_outblue = list()

# 循环遍历不同的 type
for (i in 1:4) {
  # 根据type选择特定的数据子集
  if (i == 1) {
    
    type="Field1_P3_P4"
    subset_ph = ph[ph$Pass %in% 1:23 & ph$Range %in% 2:32, ]
  } else if(i==2) {
    type="Field1_P1"
    subset_ph = ph[ph$Pass %in% 1:24 & ph$Range %in% 33:35, ]
  } else if(i==3) {
    type="Field2_3_4_5_P1"
    subset_ph = ph[ph$Pass %in% 25:120 & ph$Range %in% 2:35, ]
  } else if(i==4) {
    type="Field6_P1"
    subset_ph = ph[ph$Pass %in% 121:128 & ph$Range %in% 1:36, ]
  } 
  # 运行程序并生成 outblue
  outblue.tmp <-  programRun(dat = subset_ph, traits)
  
  # 将当前的 outblue 添加到 combined_outblue 列表中
  combined_outblue = c(combined_outblue, list(outblue.tmp))
  
}

idx=merge.fun(combined_outblue,traits)

count=analysis(idx,traits)
count$Type="P"
write.table(count,"Output/count.csv",sep=",",row.names=FALSE,col.names=FALSE,append=TRUE)

#####绘图

pdf("Output/SDGX1_plot_P.pdf")

par(mfrow=c(2,3),mar = c(5, 4, 4, 2) + 0.1,cex.main = 0.8) 
type="P"
drawplot(idx,traits)
dev.off()

outblue_name <- paste0("Output/23SDGX1_P_datalist_bluelist.rds")

saveRDS(idx,file = outblue_name)



# SDGX2_TC2 --------------------------------------------------------------------

cau.data <- phenos %>% dplyr::filter(Location %in% "SDGX2")
cau.data$AOA="MCSU"
cau.data$Year=2023


target_cols <- c( "HarvNote")#"苗期综合评价",,"收获前综合评价"

# 搜索包含字符的行，并返回逻辑向量
has_char <- apply(cau.data[target_cols], 1, function(x) any(grepl("[[:alnum:][:punct:]]+", x)))


discard=cau.data[has_char,]
discard[,c("HarvNote")]
#discard[,c(27,28,87)]
ph_filtered=cau.data[!has_char,]
N_discard=nrow(ph_filtered[ph_filtered$PlotDiscarded %in% c("Yes"), ])

ph_filtered <- ph_filtered[!(ph_filtered$PlotDiscarded %in% c("Yes")), ]
names(ph_filtered)
ph=ph_filtered


traits=c("YLD14",	"MST","PLOTWT")#,"TWT","PHT","EHT"

combined_outblue = list()

# 循环遍历不同的 type
for (i in 1:2) {
  # 根据type选择特定的数据子集
  if (i == 1) {
    
    type="Field1_2_TC2"
    subset_ph = ph[ph$Pass %in% 1:48 & ph$Range %in% 2:38, ]
  } else if(i==2) {
    type="Field3_TC2"
    subset_ph = ph[ph$Pass %in% 49:60 & ph$Range %in% 2:35, ]
  } 
  # 运行程序并生成 outblue
  outblue.tmp <-  programRun(dat = subset_ph, traits)
  
  # 将当前的 outblue 添加到 combined_outblue 列表中
  combined_outblue = c(combined_outblue, list(outblue.tmp))
  
}

idx=merge.fun(combined_outblue,traits)

count=analysis(idx,traits)
#count$Type="P"
write.table(count,"Output/count.csv",sep=",",row.names=FALSE,col.names=FALSE,append=TRUE)

#####绘图

pdf("Output/SDGX2_plot_TC2.pdf")

par(mfrow=c(2,3),mar = c(5, 4, 4, 2) + 0.1,cex.main = 0.8) 
type="TC2"
drawplot(idx,traits)
dev.off()

outblue_name <- paste0("Output/23SDGX2_TC2_datalist_bluelist.rds")

saveRDS(idx,file = outblue_name)




# SDWS1_P --------------------------------------------------------------------

cau.data <- phenos %>% dplyr::filter(Location %in% "SDWS1")
cau.data$AOA="MCSU"
cau.data$Year=2023


target_cols <- c( "HarvNote")#"苗期综合评价",,"收获前综合评价"

# 搜索包含字符的行，并返回逻辑向量
has_char <- apply(cau.data[target_cols], 1, function(x) any(grepl("[[:alnum:][:punct:]]+", x)))


discard=cau.data[has_char,]
discard[,c("HarvNote")]
#discard[,c(27,28,87)]
ph_filtered=cau.data[!has_char,]
N_discard=nrow(ph_filtered[ph_filtered$PlotDiscarded %in% c("Yes"), ])

ph_filtered <- ph_filtered[!(ph_filtered$PlotDiscarded %in% c("Yes")), ]
names(ph_filtered)
ph=ph_filtered


traits=c("YLD14",	"MST","PLOTWT")#,"TWT","PHT","EHT"

combined_outblue = list()

# 循环遍历不同的 type
for (i in 1:4) {
  # 根据type选择特定的数据子集
  if (i == 1) {
    
    type="Field1_P3_P4"
    subset_ph = ph[ph$Pass %in% 1:23 & ph$Range %in% 2:28, ]
  } else if(i==2) {
    type="Field1_P1"
    subset_ph = ph[ph$Pass %in% 1:24 & ph$Range %in% 29:40, ]
  } else if(i==3) {
    type="Field2_3_P1"
    subset_ph = ph[ph$Pass %in% 25:72 & ph$Range %in% 3:41, ]
  } else if(i==4) {
    type="Field4_P1"
    subset_ph = ph[ph$Pass %in% 73:96 & ph$Range %in% 4:42, ]
  } 
  # 运行程序并生成 outblue
  outblue.tmp <-  programRun(dat = subset_ph, traits)
  
  # 将当前的 outblue 添加到 combined_outblue 列表中
  combined_outblue = c(combined_outblue, list(outblue.tmp))
  
}

idx=merge.fun(combined_outblue,traits)

count=analysis(idx,traits)
count$Type="P"
write.table(count,"Output/count.csv",sep=",",row.names=FALSE,col.names=FALSE,append=TRUE)

#####绘图

pdf("Output/SDWS1_plot_P.pdf")

par(mfrow=c(2,3),mar = c(5, 4, 4, 2) + 0.1,cex.main = 0.8) 
type="P"
drawplot(idx,traits)
dev.off()

outblue_name <- paste0("Output/23SDWS1_P_datalist_bluelist.rds")

saveRDS(idx,file = outblue_name)






# SDWS2 --------------------------------------------------------------------

cau.data <- phenos %>% dplyr::filter(Location %in% "SDWS2")
cau.data$AOA="MCSU"
cau.data$Year=2023


target_cols <- c( "HarvNote")#"苗期综合评价",,"收获前综合评价"

# 搜索包含字符的行，并返回逻辑向量
has_char <- apply(cau.data[target_cols], 1, function(x) any(grepl("[[:alnum:][:punct:]]+", x)))


discard=cau.data[has_char,]
discard[,c("HarvNote")]
#discard[,c(27,28,87)]
ph_filtered=cau.data[!has_char,]
N_discard=nrow(ph_filtered[ph_filtered$PlotDiscarded %in% c("Yes"), ])

ph_filtered <- ph_filtered[!(ph_filtered$PlotDiscarded %in% c("Yes")), ]
names(ph_filtered)
ph=ph_filtered

####SDWS2_P####
traits=c("YLD14",	"MST","PLOTWT")#,"TWT","PHT","EHT"

combined_outblue = list()

# 循环遍历不同的 type
for (i in 1) {
  # 根据type选择特定的数据子集
  if (i == 1) {
    
    type="Field1_P1"
    subset_ph = ph[ph$Pass %in% 1:24 & ph$Range %in% 2:34, ]
    subset_ph = subset_ph[substr(subset_ph$Type,1,1)=="P",]
  } 
  # 运行程序并生成 outblue
  outblue.tmp <-  programRun(dat = subset_ph, traits)
  
  # 将当前的 outblue 添加到 combined_outblue 列表中
  combined_outblue = c(combined_outblue, list(outblue.tmp))
  
}

idx=merge.fun(combined_outblue,traits)

count=analysis(idx,traits)
count$Type="P"
write.table(count,"Output/count.csv",sep=",",row.names=FALSE,col.names=FALSE,append=TRUE)

#####绘图

pdf("Output/SDWS2_plot_P.pdf")

par(mfrow=c(2,3),mar = c(5, 4, 4, 2) + 0.1,cex.main = 0.8) 
type="P"
drawplot(idx,traits)
dev.off()

outblue_name <- paste0("Output/23SDWS2_P_datalist_bluelist.rds")

saveRDS(idx,file = outblue_name)

####SDWS2_TC2####
traits=c("YLD14",	"MST","PLOTWT")#,"TWT","PHT","EHT"

combined_outblue = list()

# 循环遍历不同的 type
for (i in 1:2) {
  # 根据type选择特定的数据子集
  if (i == 1) {
    
    type="Field1_TC2"
    subset_ph = ph[ph$Pass %in% 1:24 & ph$Range %in% 2:34, ]
    subset_ph = subset_ph[substr(subset_ph$Type,1,3)=="TC2",]
  } else if(i == 2) {
    
    type="Field2_3_4_TC2"
    subset_ph = ph[ph$Pass %in% 25:84 & ph$Range %in% 2:34, ]
    
  } 
  # 运行程序并生成 outblue
  outblue.tmp <-  programRun(dat = subset_ph, traits)
  
  # 将当前的 outblue 添加到 combined_outblue 列表中
  combined_outblue = c(combined_outblue, list(outblue.tmp))
  
}

idx=merge.fun(combined_outblue,traits)

count=analysis(idx,traits)
#count$Type="P"
write.table(count,"Output/count.csv",sep=",",row.names=FALSE,col.names=FALSE,append=TRUE)

#####绘图

pdf("Output/SDWS2_plot_TC2.pdf")

par(mfrow=c(2,3),mar = c(5, 4, 4, 2) + 0.1,cex.main = 0.8) 
type="TC2"
drawplot(idx,traits)
dev.off()

outblue_name <- paste0("Output/23SDWS2_TC2_datalist_bluelist.rds")

saveRDS(idx,file = outblue_name)







# HBNJ2_TC1 --------------------------------------------------------------------

cau.data <- phenos %>% dplyr::filter(Location %in% "HBNJ2")
cau.data$AOA="MCSU"
cau.data$Year=2023


target_cols <- c( "HarvNote")#"苗期综合评价",,"收获前综合评价"

# 搜索包含字符的行，并返回逻辑向量
has_char <- apply(cau.data[target_cols], 1, function(x) any(grepl("[[:alnum:][:punct:]]+", x)))


discard=cau.data[has_char,]
discard[,c("HarvNote")]
#discard[,c(27,28,87)]
ph_filtered=cau.data[!has_char,]
N_discard=nrow(ph_filtered[ph_filtered$PlotDiscarded %in% c("Yes"), ])

ph_filtered <- ph_filtered[!(ph_filtered$PlotDiscarded %in% c("Yes")), ]
names(ph_filtered)
ph=ph_filtered




traits=c("YLD14",	"MST","PLOTWT","PHT","EHT")#,"TWT"

combined_outblue = list()

# 循环遍历不同的 type
for (i in 1) {
  # 根据type选择特定的数据子集
  if (i == 1) {
    
    type="Field1_2_3_4_5_TC1"
    subset_ph = ph[ph$Pass %in% 1:104 & ph$Range %in% 2:49, ]
  } 
  
  
  # 运行程序并生成 outblue
  outblue.tmp <-  programRun(dat = subset_ph, traits)
  
  # 将当前的 outblue 添加到 combined_outblue 列表中
  combined_outblue = c(combined_outblue, list(outblue.tmp))
  
}

idx=merge.fun(combined_outblue,traits)

count=analysis(idx,traits)
#count$Type="P"
write.table(count,"Output/count.csv",sep=",",row.names=FALSE,col.names=FALSE,append=TRUE)

#####绘图

pdf("Output/HBNJ2_plot_TC1.pdf")

par(mfrow=c(2,3),mar = c(5, 4, 4, 2) + 0.1,cex.main = 0.8) 
type="TC1"
drawplot(idx,traits)
dev.off()

outblue_name <- paste0("Output/23HBNJ2_TC1_datalist_bluelist.rds")

saveRDS(idx,file = outblue_name)




# HBWQ1 --------------------------------------------------------------------

cau.data <- phenos %>% dplyr::filter(Location %in% "HBWQ1")
cau.data$AOA="NCSU"
cau.data$Year=2023


target_cols <- c( "HarvNote")#"苗期综合评价",,"收获前综合评价"

# 搜索包含字符的行，并返回逻辑向量
has_char <- apply(cau.data[target_cols], 1, function(x) any(grepl("[[:alnum:][:punct:]]+", x)))


discard=cau.data[has_char,]
discard[,c("HarvNote")]
#discard[,c(27,28,87)]
ph_filtered=cau.data[!has_char,]
N_discard=nrow(ph_filtered[ph_filtered$PlotDiscarded %in% c("Yes"), ])

ph_filtered <- ph_filtered[!(ph_filtered$PlotDiscarded %in% c("Yes")), ]
names(ph_filtered)
ph=ph_filtered



####HBWQ1_P####
traits=c("YLD14",	"MST","PLOTWT","PHT","EHT","TWT")#

combined_outblue = list()

# 循环遍历不同的 type
for (i in 1:5) {
  # 根据type选择特定的数据子集
  if (i == 1) {
    
    type="Field1_P3_P4"
    subset_ph = ph[ph$Pass %in% 1:23 & ph$Range %in% 2:32, ]
  } else if (i==2) {
    type="Field1_P1"
    subset_ph = ph[ph$Pass %in% 1:24 & ph$Range %in% 33:41, ]
  } else if (i==3) {
    type="Field2_3_4_P1"
    subset_ph = ph[ph$Pass %in% 25:96 & ph$Range %in% 2:41, ]
  } else if (i==4) {
    type="Field5_P1"
    subset_ph = ph[ph$Pass %in% 97:120 & ph$Range %in% 2:7, ]
  } else if (i==5) {
    type="Field7_P1"
    subset_ph = ph[ph$Pass %in% 145:168 & ph$Range %in% 13:28, ]
  }
  
  
  # 运行程序并生成 outblue
  outblue.tmp <-  programRun(dat = subset_ph, traits)
  
  # 将当前的 outblue 添加到 combined_outblue 列表中
  combined_outblue = c(combined_outblue, list(outblue.tmp))
  
}

idx=merge.fun(combined_outblue,traits)

count=analysis(idx,traits)
count$Type="P"
write.table(count,"Output/count.csv",sep=",",row.names=FALSE,col.names=FALSE,append=TRUE)

#####绘图

pdf("Output/HBWQ1_plot_P.pdf")

par(mfrow=c(2,3),mar = c(5, 4, 4, 2) + 0.1,cex.main = 0.8) 
type="P"
drawplot(idx,traits)
dev.off()

outblue_name <- paste0("Output/23HBWQ1_P_datalist_bluelist.rds")

saveRDS(idx,file = outblue_name)



####HBWQ1_TC2####
traits=c("YLD14",	"MST","PLOTWT","PHT","EHT","TWT")#

combined_outblue = list()

# 循环遍历不同的 type
for (i in 1:4) {
  # 根据type选择特定的数据子集
  if (i == 1) {
    
    type="Field5_TC2"
    subset_ph = ph[ph$Pass %in% 97:120 & ph$Range %in% 8:41, ]
  } else if (i==2) {
    type="Field6_TC2"
    subset_ph = ph[ph$Pass %in% 121:144 & ph$Range %in% 2:41, ]
  } else if (i==3) {
    type="Field7_TC2"
    subset_ph = ph[ph$Pass %in% 145:168 & ph$Range %in% 2:12, ]
  } else if (i==4) {
    type="Field7_TC2"
    subset_ph = ph[ph$Pass %in% 145:168 & ph$Range %in% 29:36, ]
  }
  
  
  # 运行程序并生成 outblue
  outblue.tmp <-  programRun(dat = subset_ph, traits)
  
  # 将当前的 outblue 添加到 combined_outblue 列表中
  combined_outblue = c(combined_outblue, list(outblue.tmp))
  
}

idx=merge.fun(combined_outblue,traits)

count=analysis(idx,traits)
#count$Type="P"
write.table(count,"Output/count.csv",sep=",",row.names=FALSE,col.names=FALSE,append=TRUE)

#####绘图

pdf("Output/HBWQ1_plot_TC2.pdf")

par(mfrow=c(2,3),mar = c(5, 4, 4, 2) + 0.1,cex.main = 0.8) 
type="TC2"
drawplot(idx,traits)
dev.off()

outblue_name <- paste0("Output/23HBWQ1_TC2_datalist_bluelist.rds")

saveRDS(idx,file = outblue_name)





# HBZX_P--------------------------------------------------------------------

cau.data <- phenos %>% dplyr::filter(Location %in% "HBZX")
cau.data$AOA="NCSU"
cau.data$Year=2023


target_cols <- c( "HarvNote")#"苗期综合评价",,"收获前综合评价"

# 搜索包含字符的行，并返回逻辑向量
has_char <- apply(cau.data[target_cols], 1, function(x) any(grepl("[[:alnum:][:punct:]]+", x)))


discard=cau.data[has_char,]
discard[,c("HarvNote")]
#discard[,c(27,28,87)]
ph_filtered=cau.data[!has_char,]
N_discard=nrow(ph_filtered[ph_filtered$PlotDiscarded %in% c("Yes"), ])

ph_filtered <- ph_filtered[!(ph_filtered$PlotDiscarded %in% c("Yes")), ]
names(ph_filtered)
ph=ph_filtered



traits=c("YLD14",	"MST","PLOTWT")#,"PHT","EHT","TWT"

combined_outblue = list()

# 循环遍历不同的 type
for (i in 1:5) {
  # 根据type选择特定的数据子集
  if (i == 1) {
    
    type="Field1_P3_P4"
    subset_ph = ph[ph$Pass %in% 1:23 & ph$Range %in% 2:32, ]
  } else if (i==2) {
    type="Field1_P1"
    subset_ph = ph[ph$Pass %in% 1:24 & ph$Range %in% 33:45, ]
  } else if (i==3) {
    type="Field2_3_4_P1"
    subset_ph = ph[ph$Pass %in% 25:96 & ph$Range %in% 2:45, ]
  } else if (i==4) {
    type="Field5_P1"
    subset_ph = ph[ph$Pass %in% 97:100 & ph$Range %in% 2:37, ]
  } else if (i==5) {
    type="Field6_P1"
    subset_ph = ph[ph$Pass %in% 101:102 & ph$Range %in% 1:8, ]
  }
  
  
  # 运行程序并生成 outblue
  outblue.tmp <-  programRun(dat = subset_ph, traits)
  
  # 将当前的 outblue 添加到 combined_outblue 列表中
  combined_outblue = c(combined_outblue, list(outblue.tmp))
  
}

idx=merge.fun(combined_outblue,traits)

count=analysis(idx,traits)
count$Type="P"
write.table(count,"Output/count.csv",sep=",",row.names=FALSE,col.names=FALSE,append=TRUE)

#####绘图

pdf("Output/HBZX_plot_P.pdf")

par(mfrow=c(2,3),mar = c(5, 4, 4, 2) + 0.1,cex.main = 0.8) 
type="P"
drawplot(idx,traits)
dev.off()

outblue_name <- paste0("Output/23HBZX_P_datalist_bluelist.rds")

saveRDS(idx,file = outblue_name)









# SDSG--------------------------------------------------------------------

cau.data <- phenos %>% dplyr::filter(Location %in% "SDSG")
cau.data$AOA="MCSU"
cau.data$Year=2023


target_cols <- c( "HarvNote")#"苗期综合评价",,"收获前综合评价"

# 搜索包含字符的行，并返回逻辑向量
has_char <- apply(cau.data[target_cols], 1, function(x) any(grepl("[[:alnum:][:punct:]]+", x)))


discard=cau.data[has_char,]
discard[,c("HarvNote")]
#discard[,c(27,28,87)]
ph_filtered=cau.data[!has_char,]
N_discard=nrow(ph_filtered[ph_filtered$PlotDiscarded %in% c("Yes"), ])

ph_filtered <- ph_filtered[!(ph_filtered$PlotDiscarded %in% c("Yes")), ]
names(ph_filtered)
ph=ph_filtered


####SDSG_TC2_HNSP####
traits=c("YLD14",	"MST","PLOTWT")#,"PHT","EHT","TWT"

combined_outblue = list()

# 循环遍历不同的 type
for (i in 1:3) {
  # 根据type选择特定的数据子集
  if (i == 1) {
    
    type="Field5_TC2"
    subset_ph = ph[ph$Pass %in% 97:120 & ph$Range %in% 11:44, ]
  } else if (i==2) {
    type="Field6_TC2"
    subset_ph = ph[ph$Pass %in% 121:144 & ph$Range %in% 2:44, ]
  } else if (i==3) {
    type="Field7_TC2"
    subset_ph = ph[ph$Pass %in% 145:150 & ph$Range %in% 2:45, ]
  } 
  
  
  # 运行程序并生成 outblue
  outblue.tmp <-  programRun(dat = subset_ph, traits)
  
  # 将当前的 outblue 添加到 combined_outblue 列表中
  combined_outblue = c(combined_outblue, list(outblue.tmp))
  
}

idx=merge.fun(combined_outblue,traits)

count=analysis(idx,traits)
#count$Type="P"
write.table(count,"Output/count.csv",sep=",",row.names=FALSE,col.names=FALSE,append=TRUE)

#####绘图

pdf("Output/SDSG_plot_TC2.pdf")

par(mfrow=c(2,3),mar = c(5, 4, 4, 2) + 0.1,cex.main = 0.8) 
type="TC2"
drawplot(idx,traits)
dev.off()

outblue_name <- paste0("Output/23SDSG_TC2_datalist_bluelist.rds")

saveRDS(idx,file = outblue_name)



####SDSG_P####
traits=c("YLD14",	"MST","PLOTWT")#,"PHT","EHT","TWT"

combined_outblue = list()

# 循环遍历不同的 type
for (i in 1:4) {
  # 根据type选择特定的数据子集
  if (i == 1) {
    
    type="Field1_P3_P4"
    subset_ph = ph[ph$Pass %in% 1:23 & ph$Range %in% 2:32, ]
  } else if (i==2) {
    type="Field1_P1"
    subset_ph = ph[ph$Pass %in% 1:24 & ph$Range %in% 33:44, ]
  } else if (i==3) {
    type="Field2_3_4_P1"
    subset_ph = ph[ph$Pass %in% 25:96 & ph$Range %in% 2:44, ]
  } else if (i==4) {
    type="Field5_P1"
    subset_ph = ph[ph$Pass %in% 97:120 & ph$Range %in% 2:10, ]
  } 
  
  
  # 运行程序并生成 outblue
  outblue.tmp <-  programRun(dat = subset_ph, traits)
  
  # 将当前的 outblue 添加到 combined_outblue 列表中
  combined_outblue = c(combined_outblue, list(outblue.tmp))
  
}

idx=merge.fun(combined_outblue,traits)

count=analysis(idx,traits)
count$Type="P"
write.table(count,"Output/count.csv",sep=",",row.names=FALSE,col.names=FALSE,append=TRUE)

#####绘图

pdf("Output/SDSG_plot_P.pdf")

par(mfrow=c(2,3),mar = c(5, 4, 4, 2) + 0.1,cex.main = 0.8) 
type="P"
drawplot(idx,traits)
dev.off()

outblue_name <- paste0("Output/23SDSG_P_datalist_bluelist.rds")

saveRDS(idx,file = outblue_name)







# HBNJ1 --------------------------------------------------------------------

cau.data <- phenos %>% dplyr::filter(Location %in% "HBNJ1")
cau.data$AOA="MCSU"
cau.data$Year=2023


target_cols <- c( "HarvNote")#"苗期综合评价",,"收获前综合评价"

# 搜索包含字符的行，并返回逻辑向量
has_char <- apply(cau.data[target_cols], 1, function(x) any(grepl("[[:alnum:][:punct:]]+", x)))


discard=cau.data[has_char,]
discard[,c("HarvNote")]
#discard[,c(27,28,87)]
ph_filtered=cau.data[!has_char,]
N_discard=nrow(ph_filtered[ph_filtered$PlotDiscarded %in% c("Yes"), ])

ph_filtered <- ph_filtered[!(ph_filtered$PlotDiscarded %in% c("Yes")), ]
names(ph_filtered)
ph=ph_filtered



####HBNJ1_p####
traits=c("YLD14",	"MST","PLOTWT","PHT","EHT")#,"TWT"

combined_outblue = list()

# 循环遍历不同的 type
for (i in 1:4) {
  # 根据type选择特定的数据子集
  if (i == 1) {
    
    type="Field1_P3_P4"
    subset_ph = ph[ph$Pass %in% 1:23 & ph$Range %in% 2:32, ]
  } else if(i==2){
    type="Field1_P1"
    subset_ph = ph[ph$Pass %in% 1:24 & ph$Range %in% 33:41, ]
  } else if(i==3){
    type="Field2_3_4_P1"
    subset_ph = ph[ph$Pass %in% 25:96 & ph$Range %in% 2:41, ]
  } else if(i==4){
    type="Field5_P1"
    subset_ph = ph[ph$Pass %in% 97:120 & ph$Range %in% 2:23, ]
  }
  
  
  # 运行程序并生成 outblue
  outblue.tmp <-  programRun(dat = subset_ph, traits)
  
  # 将当前的 outblue 添加到 combined_outblue 列表中
  combined_outblue = c(combined_outblue, list(outblue.tmp))
  
}

idx=merge.fun(combined_outblue,traits)

count=analysis(idx,traits)
count$Type="P"
write.table(count,"Output/count.csv",sep=",",row.names=FALSE,col.names=FALSE,append=TRUE)

#####绘图

pdf("Output/HBNJ1_plot_P.pdf")

par(mfrow=c(2,3),mar = c(5, 4, 4, 2) + 0.1,cex.main = 0.8) 
type="P"
drawplot(idx,traits)
dev.off()

outblue_name <- paste0("Output/23HBNJ1_P_datalist_bluelist.rds")

saveRDS(idx,file = outblue_name)


####HBNJ1_TC2####
traits=c("YLD14",	"MST","PLOTWT","PHT","EHT")#,"TWT"

combined_outblue = list()

# 循环遍历不同的 type
for (i in 1:3) {
  # 根据type选择特定的数据子集
  if (i == 1) {
    
    type="Field5_TC2"
    subset_ph = ph[ph$Pass %in% 97:120 & ph$Range %in% 24:41, ]
  } else if(i==2){
    type="Field6_TC2"
    subset_ph = ph[ph$Pass %in% 121:144 & ph$Range %in% 2:41, ]
  } else if(i==3){
    type="Field7_TC2"
    subset_ph = ph[ph$Pass %in% 145:168 & ph$Range %in% 2:34, ]
  } 
  
  # 运行程序并生成 outblue
  outblue.tmp <-  programRun(dat = subset_ph, traits)
  
  # 将当前的 outblue 添加到 combined_outblue 列表中
  combined_outblue = c(combined_outblue, list(outblue.tmp))
  
}

idx=merge.fun(combined_outblue,traits)

count=analysis(idx,traits)
#count$Type="P"
write.table(count,"Output/count.csv",sep=",",row.names=FALSE,col.names=FALSE,append=TRUE)

#####绘图

pdf("Output/HBNJ1_plot_TC2.pdf")

par(mfrow=c(2,3),mar = c(5, 4, 4, 2) + 0.1,cex.main = 0.8) 
type="TC2"
drawplot(idx,traits)
dev.off()

outblue_name <- paste0("Output/23HBNJ1_TC2_datalist_bluelist.rds")

saveRDS(idx,file = outblue_name)



# HBWQ2 --------------------------------------------------------------------

cau.data <- phenos %>% dplyr::filter(Location %in% "HBWQ2")
cau.data$AOA="NCSU"
cau.data$Year=2023


target_cols <- c( "HarvNote")#"苗期综合评价",,"收获前综合评价"

# 搜索包含字符的行，并返回逻辑向量
has_char <- apply(cau.data[target_cols], 1, function(x) any(grepl("[[:alnum:][:punct:]]+", x)))


discard=cau.data[has_char,]
discard[,c("HarvNote")]
#discard[,c(27,28,87)]
ph_filtered=cau.data[!has_char,]
N_discard=nrow(ph_filtered[ph_filtered$PlotDiscarded %in% c("Yes"), ])

ph_filtered <- ph_filtered[!(ph_filtered$PlotDiscarded %in% c("Yes")), ]
names(ph_filtered)
ph=ph_filtered

####HBWQ2_TC1####
traits=c("YLD14",	"MST","PLOTWT","PHT","EHT","TWT")#

combined_outblue = list()

# 循环遍历不同的 type
for (i in 1) {
  # 根据type选择特定的数据子集
  if (i == 1) {
    
    type="Field1_2_3_4_5_6_TC1"
    subset_ph = ph[ph$Pass %in% 1:124 & ph$Range %in% 2:42, ]
  }
  
  # 运行程序并生成 outblue
  outblue.tmp <-  programRun(dat = subset_ph, traits)
  
  # 将当前的 outblue 添加到 combined_outblue 列表中
  combined_outblue = c(combined_outblue, list(outblue.tmp))
  
}

idx=merge.fun(combined_outblue,traits)

count=analysis(idx,traits)
#count$Type="P"
write.table(count,"Output/count.csv",sep=",",row.names=FALSE,col.names=FALSE,append=TRUE)

#####绘图

pdf("Output/HBWQ2_plot_TC1.pdf")

par(mfrow=c(2,3),mar = c(5, 4, 4, 2) + 0.1,cex.main = 0.8) 
type="TC1"
drawplot(idx,traits)
dev.off()

outblue_name <- paste0("Output/23HBWQ2_TC1_datalist_bluelist.rds")

saveRDS(idx,file = outblue_name)



# HBJZ --------------------------------------------------------------------

cau.data <- phenos %>% dplyr::filter(Location %in% "HBJZ")
cau.data$AOA="NCSU"
cau.data$Year=2023


target_cols <- c( "HarvNote")#"苗期综合评价",,"收获前综合评价"

# 搜索包含字符的行，并返回逻辑向量
has_char <- apply(cau.data[target_cols], 1, function(x) any(grepl("[[:alnum:][:punct:]]+", x)))


discard=cau.data[has_char,]
discard[,c("HarvNote")]
#discard[,c(27,28,87)]
ph_filtered=cau.data[!has_char,]
N_discard=nrow(ph_filtered[ph_filtered$PlotDiscarded %in% c("Yes"), ])

ph_filtered <- ph_filtered[!(ph_filtered$PlotDiscarded %in% c("Yes")), ]
names(ph_filtered)
ph=ph_filtered



####HBJZ_p####
traits=c("YLD14",	"MST","PLOTWT")#,"TWT","PHT","EHT"

combined_outblue = list()

# 循环遍历不同的 type
for (i in 1:3) {
  # 根据type选择特定的数据子集
  if (i == 1) {
    
    type="Field1_P3_P4"
    subset_ph = ph[ph$Pass %in% 1:23 & ph$Range %in% 2:32, ]
  } else if(i==2){
    type="Field2_3_4_5_P1"
    subset_ph = ph[ph$Pass %in% 25:120 & ph$Range %in% 2:32, ]
  } else if(i==3){
    type="Field6_P1"
    subset_ph = ph[ph$Pass %in% 121:144 & ph$Range %in% 2:28, ]
  } 
  
  
  # 运行程序并生成 outblue
  outblue.tmp <-  programRun(dat = subset_ph, traits)
  
  # 将当前的 outblue 添加到 combined_outblue 列表中
  combined_outblue = c(combined_outblue, list(outblue.tmp))
  
}

idx=merge.fun(combined_outblue,traits)

count=analysis(idx,traits)
count$Type="P"
write.table(count,"Output/count.csv",sep=",",row.names=FALSE,col.names=FALSE,append=TRUE)

#####绘图

pdf("Output/HBJZ_plot_P.pdf")

par(mfrow=c(2,3),mar = c(5, 4, 4, 2) + 0.1,cex.main = 0.8) 
type="P"
drawplot(idx,traits)
dev.off()

outblue_name <- paste0("Output/23HBJZ_P_datalist_bluelist.rds")

saveRDS(idx,file = outblue_name)


####HBJZ_TC2####
traits=c("YLD14",	"MST","PLOTWT")#,"TWT","PHT","EHT"

combined_outblue = list()

# 循环遍历不同的 type
for (i in 1:3) {
  # 根据type选择特定的数据子集
  if (i == 1) {
    
    type="Field6_TC2"
    subset_ph = ph[ph$Pass %in% 121:144 & ph$Range %in% 29:32, ]
  } else if(i==2){
    type="Field7_8_9_10_TC2"
    subset_ph = ph[ph$Pass %in% 145:212 & ph$Range %in% 2:32, ]
  } else if(i==3){
    type="Field11_TC2"
    subset_ph = ph[ph$Pass %in% 213:214 & ph$Range %in% 5:7, ]
  } 
  
  # 运行程序并生成 outblue
  outblue.tmp <-  programRun(dat = subset_ph, traits)
  
  # 将当前的 outblue 添加到 combined_outblue 列表中
  combined_outblue = c(combined_outblue, list(outblue.tmp))
  
}

idx=merge.fun(combined_outblue,traits)

count=analysis(idx,traits)
#count$Type="P"
write.table(count,"Output/count.csv",sep=",",row.names=FALSE,col.names=FALSE,append=TRUE)

#####绘图

pdf("Output/HBJZ_plot_TC2.pdf")

par(mfrow=c(2,3),mar = c(5, 4, 4, 2) + 0.1,cex.main = 0.8) 
type="TC2"
drawplot(idx,traits)
dev.off()

outblue_name <- paste0("Output/23HBJZ_TC2_datalist_bluelist.rds")

saveRDS(idx,file = outblue_name)




#########################SWCN#######################

# SCCD --------------------------------------------------------------------

cau.data <- phenos %>% dplyr::filter(Location %in% "SCCD")
cau.data$AOA="SWCN"
cau.data$Year=2023


target_cols <- c( "HarvNote")#"苗期综合评价",,"收获前综合评价"

# 搜索包含字符的行，并返回逻辑向量
has_char <- apply(cau.data[target_cols], 1, function(x) any(grepl("[[:alnum:][:punct:]]+", x)))


discard=cau.data[has_char,]
discard[,c("HarvNote")]
#discard[,c(27,28,87)]
ph_filtered=cau.data[!has_char,]
N_discard=nrow(ph_filtered[ph_filtered$PlotDiscarded %in% c("Yes"), ])

ph_filtered <- ph_filtered[!(ph_filtered$PlotDiscarded %in% c("Yes")), ]
names(ph_filtered)
ph=ph_filtered

####SCCD_P####
traits=c("YLD14",	"MST","PLOTWT","PHT","EHT")#"TWT",

combined_outblue = list()

# 循环遍历不同的 type
for (i in 1:1) {
  # 根据type选择特定的数据子集
  if (i == 1) {
    
    type="Field1_P1"
    subset_ph = ph[ph$Pass %in% 1:24 & ph$Range %in% 1:8, ]
  } 
  # 运行程序并生成 outblue
  outblue.tmp <-  programRun(dat = subset_ph, traits)
  
  # 将当前的 outblue 添加到 combined_outblue 列表中
  combined_outblue = c(combined_outblue, list(outblue.tmp))
  
}

idx=merge.fun(combined_outblue,traits)

count=analysis(idx,traits)
count$Type="P"
write.table(count,"Output/count.csv",sep=",",row.names=FALSE,col.names=FALSE,append=TRUE)

#####绘图

pdf("Output/SCCD_plot_P.pdf")

par(mfrow=c(2,3),mar = c(5, 4, 4, 2) + 0.1,cex.main = 0.8) 
type="P"
drawplot(idx,traits)
dev.off()

outblue_name <- paste0("Output/23SCCD_P_datalist_bluelist.rds")

saveRDS(idx,file = outblue_name)




####SCCD_TC2####
traits=c("YLD14",	"MST","PLOTWT","PHT","EHT")#"TWT",

combined_outblue = list()

# 循环遍历不同的 type
for (i in 1:1) {
  # 根据type选择特定的数据子集
  if (i == 1) {
    
    type="Field2_TC2"
    subset_ph = ph[ph$Pass %in% 1:24 & ph$Range %in% 9:17, ]
  } 
  # 运行程序并生成 outblue
  outblue.tmp <-  programRun(dat = subset_ph, traits)
  
  # 将当前的 outblue 添加到 combined_outblue 列表中
  combined_outblue = c(combined_outblue, list(outblue.tmp))
  
}

idx=merge.fun(combined_outblue,traits)

count=analysis(idx,traits)
#count$Type="P"
write.table(count,"Output/count.csv",sep=",",row.names=FALSE,col.names=FALSE,append=TRUE)

#####绘图

pdf("Output/SCCD_plot_TC2.pdf")

par(mfrow=c(2,3),mar = c(5, 4, 4, 2) + 0.1,cex.main = 0.8) 
type="TC2"
drawplot(idx,traits)
dev.off()

outblue_name <- paste0("Output/23SCCD_TC2_datalist_bluelist.rds")

saveRDS(idx,file = outblue_name)





# SCST_P1 --------------------------------------------------------------------

cau.data <- phenos %>% dplyr::filter(Location %in% "SCST")
cau.data$AOA="SWCN"
cau.data$Year=2023


target_cols <- c( "HarvNote")#"苗期综合评价",,"收获前综合评价"

# 搜索包含字符的行，并返回逻辑向量
has_char <- apply(cau.data[target_cols], 1, function(x) any(grepl("[[:alnum:][:punct:]]+", x)))


discard=cau.data[has_char,]
discard[,c("HarvNote")]
#discard[,c(27,28,87)]
ph_filtered=cau.data[!has_char,]
N_discard=nrow(ph_filtered[ph_filtered$PlotDiscarded %in% c("Yes"), ])

ph_filtered <- ph_filtered[!(ph_filtered$PlotDiscarded %in% c("Yes")), ]
names(ph_filtered)
ph=ph_filtered

# names(ph)[names(ph)=="YLD14"]="YLD14_0.6"
# names(ph)[names(ph)=="YLD14_BE"]="YLD14"
traits=c("YLD14",	"MST","PLOTWT","PHT","EHT")#"TWT",

combined_outblue = list()

# 循环遍历不同的 type
for (i in 1) {
  # 根据type选择特定的数据子集
  if (i == 1) {
    
    type="Field1_P1"
    subset_ph = ph[ph$Pass %in% 1:24 & ph$Range %in% 1:8, ]
  } 
  # 运行程序并生成 outblue
  outblue.tmp <-  programRun(dat = subset_ph, traits)
  
  # 将当前的 outblue 添加到 combined_outblue 列表中
  combined_outblue = c(combined_outblue, list(outblue.tmp))
  
}

idx=merge.fun(combined_outblue,traits)

count=analysis(idx,traits)
count$Type="P"
write.table(count,"Output/count.csv",sep=",",row.names=FALSE,col.names=FALSE,append=TRUE)

#####绘图

pdf("Output/SCST_plot_P.pdf")

par(mfrow=c(2,3),mar = c(5, 4, 4, 2) + 0.1,cex.main = 0.8) 
type="P"
drawplot(idx,traits)
dev.off()

outblue_name <- paste0("Output/23SCST_P_datalist_bluelist.rds")

saveRDS(idx,file = outblue_name)






# GZAS --------------------------------------------------------------------

cau.data <- phenos %>% dplyr::filter(Location %in% "GZAS")
cau.data$AOA="SWCN"
cau.data$Year=2023


target_cols <- c( "HarvNote")#"苗期综合评价",,"收获前综合评价"

# 搜索包含字符的行，并返回逻辑向量
has_char <- apply(cau.data[target_cols], 1, function(x) any(grepl("[[:alnum:][:punct:]]+", x)))


discard=cau.data[has_char,]
discard[,c("HarvNote")]
#discard[,c(27,28,87)]
ph_filtered=cau.data[!has_char,]
N_discard=nrow(ph_filtered[ph_filtered$PlotDiscarded %in% c("Yes"), ])

ph_filtered <- ph_filtered[!(ph_filtered$PlotDiscarded %in% c("Yes")), ]
names(ph_filtered)
ph=ph_filtered

####GZAS_P####
traits=c("YLD14",	"MST","PLOTWT","PHT","EHT")#"TWT",

combined_outblue = list()
par(mfrow=c(2,3),mar = c(5, 4, 4, 2) + 0.1,cex.main = 0.8)
# 循环遍历不同的 type
for (i in 1) {
  # 根据type选择特定的数据子集
  if (i == 1) {
    
    type="Field1_P1"
    subset_ph = ph[ph$Pass %in% 1:24 & ph$Range %in% 1:8, ]
  } 
  # 运行程序并生成 outblue
  outblue.tmp <-  programRun(dat = subset_ph, traits)
  
  # 将当前的 outblue 添加到 combined_outblue 列表中
  combined_outblue = c(combined_outblue, list(outblue.tmp))
  
}

idx=merge.fun(combined_outblue,traits)

count=analysis(idx,traits)
count$Type="P"
write.table(count,"Output/count.csv",sep=",",row.names=FALSE,col.names=FALSE,append=TRUE)

#####绘图

pdf("Output/GZAS_plot_P.pdf")

par(mfrow=c(2,3),mar = c(5, 4, 4, 2) + 0.1,cex.main = 0.8) 
type="P"
drawplot(idx,traits)
dev.off()

outblue_name <- paste0("Output/23GZAS_P_datalist_bluelist.rds")

saveRDS(idx,file = outblue_name)

####GZAS_TC2####
traits=c("YLD14",	"MST","PLOTWT","PHT","EHT")#"TWT",

combined_outblue = list()
par(mfrow=c(2,3),mar = c(5, 4, 4, 2) + 0.1,cex.main = 0.8)
# 循环遍历不同的 type
for (i in 1) {
  # 根据type选择特定的数据子集
  if (i == 1) {
    
    type="Field2_TC2"
    subset_ph = ph[ph$Pass %in% 1:24 & ph$Range %in% 9:17, ]
  } 
  # 运行程序并生成 outblue
  outblue.tmp <-  programRun(dat = subset_ph, traits)
  
  # 将当前的 outblue 添加到 combined_outblue 列表中
  combined_outblue = c(combined_outblue, list(outblue.tmp))
  
}

idx=merge.fun(combined_outblue,traits)

count=analysis(idx,traits)
#count$Type="P"
write.table(count,"Output/count.csv",sep=",",row.names=FALSE,col.names=FALSE,append=TRUE)

#####绘图

pdf("Output/GZAS_plot_TC2.pdf")

par(mfrow=c(2,3),mar = c(5, 4, 4, 2) + 0.1,cex.main = 0.8) 
type="TC2"
drawplot(idx,traits)
dev.off()

outblue_name <- paste0("Output/23GZAS_TC2_datalist_bluelist.rds")

saveRDS(idx,file = outblue_name)


####GZAS_TC1####
traits=c("YLD14",	"MST","PLOTWT","PHT","EHT")#"TWT",

combined_outblue = list()
par(mfrow=c(2,3),mar = c(5, 4, 4, 2) + 0.1,cex.main = 0.8)
# 循环遍历不同的 type
for (i in 1) {
  # 根据type选择特定的数据子集
  if (i == 1) {
    
    type="Field3_TC1"
    subset_ph = ph[ph$Pass %in% 1:24 & ph$Range %in% 18:165, ]
  } 
  # 运行程序并生成 outblue
  outblue.tmp <-  programRun(dat = subset_ph, traits)
  
  # 将当前的 outblue 添加到 combined_outblue 列表中
  combined_outblue = c(combined_outblue, list(outblue.tmp))
  
}

idx=merge.fun(combined_outblue,traits)

count=analysis(idx,traits)
count$Type="TC1"
write.table(count,"Output/count.csv",sep=",",row.names=FALSE,col.names=FALSE,append=TRUE)

#####绘图

pdf("Output/GZAS_plot_TC1.pdf")

par(mfrow=c(2,3),mar = c(5, 4, 4, 2) + 0.1,cex.main = 0.8) 
type="TC1"
drawplot(idx,traits)
dev.off()

outblue_name <- paste0("Output/23GZAS_TC1_datalist_bluelist.rds")

saveRDS(idx,file = outblue_name)





# 

# YNFM_P1 --------------------------------------------------------------------

cau.data <- phenos %>% dplyr::filter(Location %in% "YNFM")
cau.data$AOA="SWCN"
cau.data$Year=2023


target_cols <- c( "HarvNote")#"苗期综合评价",,"收获前综合评价"

# 搜索包含字符的行，并返回逻辑向量
has_char <- apply(cau.data[target_cols], 1, function(x) any(grepl("[[:alnum:][:punct:]]+", x)))


discard=cau.data[has_char,]
discard[,c("HarvNote")]
#discard[,c(27,28,87)]
ph_filtered=cau.data[!has_char,]
N_discard=nrow(ph_filtered[ph_filtered$PlotDiscarded %in% c("Yes"), ])

ph_filtered <- ph_filtered[!(ph_filtered$PlotDiscarded %in% c("Yes")), ]
names(ph_filtered)
ph=ph_filtered

traits=c("YLD14",	"MST","PLOTWT","PHT","EHT")#"TWT",

combined_outblue = list()

# 循环遍历不同的 type
for (i in 1:1) {
  # 根据type选择特定的数据子集
  if (i == 1) {
    
    type="Field1_P1"
    subset_ph = ph
  } 
  # 运行程序并生成 outblue
  outblue.tmp <-  programRun(dat = subset_ph, traits)
  
  # 将当前的 outblue 添加到 combined_outblue 列表中
  combined_outblue = c(combined_outblue, list(outblue.tmp))
  
}

idx=merge.fun(combined_outblue,traits)

count=analysis(idx,traits)
count$Type="P"
write.table(count,"Output/count.csv",sep=",",row.names=FALSE,col.names=FALSE,append=TRUE)

#####绘图

pdf("Output/YNFM_plot_P.pdf")

par(mfrow=c(2,3),mar = c(5, 4, 4, 2) + 0.1,cex.main = 0.8) 
type="P"
drawplot(idx,traits)
dev.off()

outblue_name <- paste0("Output/23YNFM_P_datalist_bluelist.rds")

saveRDS(idx,file = outblue_name)









# YNZY --------------------------------------------------------------------

cau.data <- phenos %>% dplyr::filter(Location %in% "YNZY")
cau.data$AOA="SWCN"
cau.data$Year=2023


target_cols <- c( "HarvNote")#"苗期综合评价",,"收获前综合评价"

# 搜索包含字符的行，并返回逻辑向量
has_char <- apply(cau.data[target_cols], 1, function(x) any(grepl("[[:alnum:][:punct:]]+", x)))


discard=cau.data[has_char,]
discard[,c("HarvNote")]
#discard[,c(27,28,87)]
ph_filtered=cau.data[!has_char,]
N_discard=nrow(ph_filtered[ph_filtered$PlotDiscarded %in% c("Yes"), ])

ph_filtered <- ph_filtered[!(ph_filtered$PlotDiscarded %in% c("Yes")), ]
names(ph_filtered)
ph=ph_filtered

###YNZY_P###
traits=c("YLD14",	"MST","PLOTWT","PHT","EHT")#"TWT",

combined_outblue = list()

# 循环遍历不同的 type
for (i in 1:1) {
  # 根据type选择特定的数据子集
  if (i == 1) {
    
    type="Field1_P1"
    subset_ph = ph[ph$Pass %in% 1:24 & ph$Range %in% 1:8, ]
  } 
  # 运行程序并生成 outblue
  outblue.tmp <-  programRun(dat = subset_ph, traits)
  
  # 将当前的 outblue 添加到 combined_outblue 列表中
  combined_outblue = c(combined_outblue, list(outblue.tmp))
  
}

idx=merge.fun(combined_outblue,traits)

count=analysis(idx,traits)
count$Type="P"
write.table(count,"Output/count.csv",sep=",",row.names=FALSE,col.names=FALSE,append=TRUE)

#####绘图

pdf("Output/YNZY_plot_P.pdf")

par(mfrow=c(2,3),mar = c(5, 4, 4, 2) + 0.1,cex.main = 0.8) 
type="P"
drawplot(idx,traits)
dev.off()

outblue_name <- paste0("Output/23YNZY_P_datalist_bluelist.rds")

saveRDS(idx,file = outblue_name)

####YNZY_TC2
traits=c("YLD14",	"MST","PLOTWT","PHT","EHT")#"TWT",

combined_outblue = list()

# 循环遍历不同的 type
for (i in 1:1) {
  # 根据type选择特定的数据子集
  if (i == 1) {
    
    type="Field1_TC2"
    subset_ph = ph[ph$Pass %in% 1:24 & ph$Range %in% 9:16, ]
  } 
  # 运行程序并生成 outblue
  outblue.tmp <-  programRun(dat = subset_ph, traits)
  
  # 将当前的 outblue 添加到 combined_outblue 列表中
  combined_outblue = c(combined_outblue, list(outblue.tmp))
  
}

idx=merge.fun(combined_outblue,traits)

count=analysis(idx,traits)
count$Type="TC2"
write.table(count,"Output/count.csv",sep=",",row.names=FALSE,col.names=FALSE,append=TRUE)

#####绘图

pdf("Output/YNZY_plot_TC2.pdf")

par(mfrow=c(2,3),mar = c(5, 4, 4, 2) + 0.1,cex.main = 0.8) 
type="TC2"
drawplot(idx,traits)
dev.off()

outblue_name <- paste0("Output/23YNZY_TC2_datalist_bluelist.rds")

saveRDS(idx,file = outblue_name)

####YNZY_TC1
traits=c("YLD14",	"MST","PLOTWT","PHT","EHT")#"TWT",

combined_outblue = list()

# 循环遍历不同的 type
for (i in 1:1) {
  # 根据type选择特定的数据子集
  if (i == 1) {
    
    type="Field1_TC1"
    subset_ph = ph[ph$Pass %in% 1:24 & ph$Range %in% 17:113, ]
  } 
  # 运行程序并生成 outblue
  outblue.tmp <-  programRun(dat = subset_ph, traits)
  
  # 将当前的 outblue 添加到 combined_outblue 列表中
  combined_outblue = c(combined_outblue, list(outblue.tmp))
  
}

idx=merge.fun(combined_outblue,traits)

count=analysis(idx,traits)
count$Type="TC1"
write.table(count,"Output/count.csv",sep=",",row.names=FALSE,col.names=FALSE,append=TRUE)

#####绘图

pdf("Output/YNZY_plot_TC1.pdf")

par(mfrow=c(2,3),mar = c(5, 4, 4, 2) + 0.1,cex.main = 0.8) 
type="TC1"
drawplot(idx,traits)
dev.off()

outblue_name <- paste0("Output/23YNZY_TC1_datalist_bluelist.rds")

saveRDS(idx,file = outblue_name)


# ########DISCARD#####
# 
# # SXXD --------------------------------------------------------------------
# #山西小店数据情况不好，地势不均匀
# cau.data <- phenos %>% dplyr::filter(Location %in% "SXXD")
# cau.data$AOA="LMSP"
# cau.data$Year=2023
# 
# 
# target_cols <- c( "HarvNote")#"苗期综合评价",,"收获前综合评价"
# 
# # 搜索包含字符的行，并返回逻辑向量
# has_char <- apply(cau.data[target_cols], 1, function(x) any(grepl("[[:alnum:][:punct:]]+", x)))
# 
# 
# discard=cau.data[has_char,]
# discard[,c("HarvNote")]
# #discard[,c(27,28,87)]
# ph_filtered=cau.data[!has_char,]
# N_discard=nrow(ph_filtered[ph_filtered$PlotDiscarded %in% c("Yes"), ])
# 
# ph_filtered <- ph_filtered[!(ph_filtered$PlotDiscarded %in% c("Yes")), ]
# names(ph_filtered)
# ph=ph_filtered
# 
# 
# traits=c("YLD14",	"MST","PLOTWT","TWT","PHT","EHT")#
# 
# combined_outblue = list()
# par(mfrow=c(2,3),mar = c(5, 4, 4, 2) + 0.1,cex.main = 0.8)
# # 循环遍历不同的 type
# for (i in 1:3) {
#   # 根据type选择特定的数据子集
#   if (i == 1) {
#     
#     type="Field1_P3_P4"
#     subset_ph = ph[ph$Pass %in% 1:23 & ph$Range %in% 2:15, ]
#   } else if(i==2){
#     type="Field1_P1"
#     subset_ph = ph[ph$Pass %in% 1:24 & ph$Range %in% 16:63, ]
#   } else if(i==3){
#     type="Field2_3_P1"
#     subset_ph = ph[ph$Pass %in% 25:64 & ph$Range %in% 2:63, ]
#   }
#   # 运行程序并生成 outblue
#   outblue.tmp <-  programRun(dat = subset_ph, traits)
#   
#   # 将当前的 outblue 添加到 combined_outblue 列表中
#   combined_outblue = c(combined_outblue, list(outblue.tmp))
#   
# }
# 
# idx=merge.fun(combined_outblue,traits)
# 
# count=analysis(idx,traits)
# count$Type="P"
# write.table(count,"Output/count.csv",sep=",",row.names=FALSE,col.names=FALSE,append=TRUE)
# 
# #####绘图
# 
# pdf("Output/SXXD_plot_P.pdf")
# 
# par(mfrow=c(2,3),mar = c(5, 4, 4, 2) + 0.1,cex.main = 0.8) 
# type="P"
# drawplot(idx,traits)
# dev.off()
# 
# outblue_name <- paste0("Output/23SXXD_P_datalist_bluelist.rds")
# 
# saveRDS(idx,file = outblue_name)
# 
# ####SCCD_TC1####
# traits=c("YLD14",	"MST","PLOTWT","PHT","EHT")#"TWT",
# 
# combined_outblue = list()
# 
# # 循环遍历不同的 type
# for (i in 1:1) {
#   # 根据type选择特定的数据子集
#   if (i == 1) {
#     
#     type="Field3_TC1"
#     subset_ph = ph[ph$Pass %in% 1:24 & ph$Range %in% 18:104, ]
#   } 
#   # 运行程序并生成 outblue
#   outblue.tmp <-  programRun(dat = subset_ph, traits)
#   
#   # 将当前的 outblue 添加到 combined_outblue 列表中
#   combined_outblue = c(combined_outblue, list(outblue.tmp))
#   
# }
# 
# idx=merge.fun(combined_outblue,traits)
# 
# count=analysis(idx,traits)
# #count$Type="P"
# write.table(count,"Output/count.csv",sep=",",row.names=FALSE,col.names=FALSE,append=TRUE)
# 
# #####绘图
# 
# pdf("Output/SCCD_plot_TC1.pdf")
# 
# par(mfrow=c(2,3),mar = c(5, 4, 4, 2) + 0.1,cex.main = 0.8) 
# type="TC1"
# drawplot(idx,traits)
# dev.off()
# 
# outblue_name <- paste0("Output/23SCCD_TC1_datalist_bluelist.rds")
# 
# saveRDS(idx,file = outblue_name)
# 
# 
# # AHSX1_P --------------------------------------------------------------------
# 
# cau.data <- phenos %>% dplyr::filter(Location %in% "AHSX1")
# cau.data$AOA="SCSU"
# cau.data$Year=2023
# 
# 
# target_cols <- c( "HarvNote")#"苗期综合评价",,"收获前综合评价"
# 
# # 搜索包含字符的行，并返回逻辑向量
# has_char <- apply(cau.data[target_cols], 1, function(x) any(grepl("[[:alnum:][:punct:]]+", x)))
# 
# 
# discard=cau.data[has_char,]
# discard[,c("HarvNote")]
# #discard[,c(27,28,87)]
# ph_filtered=cau.data[!has_char,]
# N_discard=nrow(ph_filtered[ph_filtered$PlotDiscarded %in% c("Yes"), ])
# 
# ph_filtered <- ph_filtered[!(ph_filtered$PlotDiscarded %in% c("Yes")), ]
# names(ph_filtered)
# ph=ph_filtered
# 
# traits=c("YLD14",	"MST","PLOTWT","PHT","EHT")#,"TWT"
# 
# combined_outblue = list()
# 
# # 循环遍历不同的 type
# for (i in 1:3) {
#   # 根据type选择特定的数据子集
#   if (i == 1) {
#     
#     type="Field1_P3_P4"
#     subset_ph = ph[ph$Pass %in% 1:23 & ph$Range %in% 2:28, ]
#   } else if(i==2) {
#     type="Field1_P1"
#     subset_ph = ph[ph$Pass %in% 1:24 & ph$Range %in% 29:36, ]
#   } else if(i==3) {
#     type="Field2_3_4_5_P1"
#     subset_ph = ph[ph$Pass %in% 25:88 & ph$Range %in% 3:37, ]
#   } 
#   # 运行程序并生成 outblue
#   outblue.tmp <-  programRun(dat = subset_ph, traits)
#   
#   # 将当前的 outblue 添加到 combined_outblue 列表中
#   combined_outblue = c(combined_outblue, list(outblue.tmp))
#   
# }
# 
# idx=merge.fun(combined_outblue,traits)
# 
# count=analysis(idx,traits)
# count$Type="P"
# write.table(count,"Output/count.csv",sep=",",row.names=FALSE,col.names=FALSE,append=TRUE)
# 
# #####绘图
# 
# pdf("Output/AHSX1_plot_P.pdf")
# 
# par(mfrow=c(2,3),mar = c(5, 4, 4, 2) + 0.1,cex.main = 0.8) 
# type="P"
# drawplot(idx,traits)
# dev.off()
# 
# outblue_name <- paste0("Output/23AHSX1_P_datalist_bluelist.rds")
# 
# saveRDS(idx,file = outblue_name)
# 
# 
# 
# 
# 
# 
# 
# ####AHSX2_P####
# traits=c("YLD14",	"MST","PLOTWT","PHT","EHT")#,"TWT"
# 
# combined_outblue = list()
# 
# # 循环遍历不同的 type
# for (i in 1:2) {
#   # 根据type选择特定的数据子集
#   if (i == 1) {
#     
#     type="Field1_P1"
#     subset_ph = ph[ph$Pass %in% 1:24 & ph$Range %in% 2:38, ]
#   } else if(i==2) {
#     type="Field2_P1"
#     subset_ph = ph[ph$Pass %in% 25:48 & ph$Range %in% 3:13, ]
#   } 
#   # 运行程序并生成 outblue
#   outblue.tmp <-  programRun(dat = subset_ph, traits)
#   
#   # 将当前的 outblue 添加到 combined_outblue 列表中
#   combined_outblue = c(combined_outblue, list(outblue.tmp))
#   
# }
# 
# idx=merge.fun(combined_outblue,traits)
# 
# count=analysis(idx,traits)
# count$Type="P"
# write.table(count,"Output/count.csv",sep=",",row.names=FALSE,col.names=FALSE,append=TRUE)
# 
# #####绘图
# 
# pdf("Output/AHSX2_plot_P.pdf")
# 
# par(mfrow=c(2,3),mar = c(5, 4, 4, 2) + 0.1,cex.main = 0.8) 
# type="P"
# drawplot(idx,traits)
# dev.off()
# 
# outblue_name <- paste0("Output/23AHSX2_P_datalist_bluelist.rds")
# 
# saveRDS(idx,file = outblue_name)

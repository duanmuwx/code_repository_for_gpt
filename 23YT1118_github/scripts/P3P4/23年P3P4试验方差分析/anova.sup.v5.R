library(readxl)
library(writexl)
library(biosim)
library(tidyverse)
library(reshape2)
library(stringi)
library(asreml)
library(ggplot2)
library(ggpubr)
library(ggpmisc)


P4Extraction <- function(blue,p4,aoa=NULL){
    if (aoa=='SP'){
        aoa <- c('MMSP','LMSP')
    }else if(aoa=='SU'){
        aoa <- c('NCSU','MCSU','SCSU')
    }
    # blue <- merge(blue,ped,by='Name',all.x=T)
    table(blue$Type)
    blue$Name <- as.character(blue$Name)
    # Extracting P4s
    blue <- blue[blue$AOA %in% aoa,]
    dim(blue)
    idlist <- unique(blue$Name[blue$Name %in% p4$Name])
    #idlist2 <- unique(blue$Name[blue$Name %in% p4$Hybrid])

    idlist <- idlist[order(idlist)]
    length(idlist) #21
    head(idlist)
    blue4 <- blue[blue$Name %in% idlist,]
    dim(blue4) #1160
    blue4cc <- blue4[blue4$AOA %in% aoa,]
    dim(blue4cc) #603
    table(blue4cc$Type,blue4cc$Location)

    table(blue4cc$Type,blue4cc$Year)

    blue4cc <- blue4cc[blue4cc$Year == 2023, ]
    # blue4cc$Type[blue4cc$Year==2023] <- 'P4'
    # blue4cc$Type[blue4cc$Year==2022] <- 'P3'
    # blue4cc$Type[blue4cc$Year<2022] <- 'P1'

    table(blue4cc$Type,blue4cc$Year)

    #        2020 2021 2022 2023
    # P1       15   88    0    0
    # P3        0    0  169    0
    # P4        0    0    0  193

    return (blue4cc)
}

P3Extraction <- function(bluelist,p3,trait,aoa=NULL){
    if (aoa=='SP'){
        aoa <- c('MMSP','LMSP')
    }else if(aoa=='SU'){
        aoa <- c('NCSU','MCSU','SCSU')
    }
    blue=bluelist[[trait]]
    blue <- blue[blue$Name %in% p3 & blue$AOA %in% aoa,]
    blue <- blue[blue$Year == 2023,]
    # blue <- blue[blue$Year>2020,]
    # blue$Type[blue$Year<2023] <- 'P1'
    # blue$Type[blue$Year==2023] <- 'P3'

    dim(blue)
    table(blue$Type,blue$Year)
    #    2021 2022 2023
    # P1   25  337    0
    # P3    0    0  417
    unilength(blue$Name) #54

    locations <- unique(blue$Location)
    locations <- locations[order(locations)]
    locations

    return(list(blue=blue,idlist=unique(blue$Name),locations=locations))
}

hybridEvaluation <- function(blue,phen,trait,control,aoa=NULL,type=NULL){
    dim(blue)
    if (aoa=='EMSP' & type=='P3') {
        blue <- phen
        names(blue)[names(blue)==trait] <- 'predicted.value'
    }
    if (aoa=='SWCN' & type=='P3') {
        blue <- phen
        names(blue)[names(blue)==trait] <- 'predicted.value'
    }
    if (aoa=='NCSU' & type=='P3') {
        blue <- phen
        names(blue)[names(blue)==trait] <- 'predicted.value'
    }

    if (aoa=='EMSP' & type=='P4') {
        blue <- phen
        names(blue)[names(blue)==trait] <- 'predicted.value'
    }
    if (aoa=='SWCN' & type=='P4') {
        blue <- phen
        names(blue)[names(blue)==trait] <- 'predicted.value'
    }
    if (aoa=='NCSU' & type=='P4') {
        blue <- phen
        names(blue)[names(blue)==trait] <- 'predicted.value'
        # blue$Year <- "2023"
    }

    blue$Name <- as.factor(blue$Name)
    blue$Location <- as.factor(blue$Location)
    blue$Year <- as.factor(blue$Year)
    blue$Type <- as.factor(blue$Type)

    nyears    <- unilength(blue$Year)
    nlocation <- unilength(blue$Location)

    if (nyears==1){
        if(nlocation == 1){
            asr <- asreml(predicted.value~1,random = ~Name,
                          residual = ~idv(units),
                          data = blue, pworkspace="5Gb")
        }else{
            asr <- asreml(predicted.value~1,random = ~Name*Location,
                          residual = ~idv(units),
                          data = blue, pworkspace="5Gb")
        }

        # nyears <- 1
    } else {
        asr <- asreml(predicted.value~1,random = ~Name*Location*Year,
                      residual = ~idv(units),
                      data = blue, pworkspace="5Gb")
    }
    summary(asr)$varcomp


    cnt=data.frame(table(phen$Name,phen$Location,phen$Year))
    cnt$id <- paste0(cnt$Var1,'_',cnt$Var2,'_',cnt$Var3)
    cnt <- cnt[,c('id','Freq')]
    cnt=cnt[cnt$Freq!=0,]
    if (nyears==1){
        if(nlocation == 1){
            bvp      <- predict(asr, classify = "Name")$pvals
            bvp$Year <- unique(phen$Year)
            bvp$Location <- unique(phen$Location)
        }else{
            bvp      <- predict(asr, classify = "Name*Location")$pvals
            bvp$Year <- unique(phen$Year)
        }
    } else {
        bvp <- predict(asr, classify = "Name*Location*Year")$pvals
    }
    bvp$id <- paste0(bvp$Name,'_',bvp$Location,'_',bvp$Year)
    bvp <- merge(bvp,cnt,by='id',all.x = T)
    bvp$Freq <- ifelse(is.na(bvp$Freq),0,bvp$Freq)
    bvp <- bvp[,c(-1)]
    #plot(std.error~Freq,bvp)
    bvp$CKPCT <- NA

    years <- unique(bvp$Year); years <- years[order(years)]
    locations <- unique(bvp$Location); locations <- locations[order(locations)]
    for (i in 1:unilength(bvp$Year)){
        for (j in 1:unilength(bvp$Location)){
            w <- bvp[bvp$Year==years[i] & bvp$Location==locations[j],]
            if(length(w$predicted.value[w$Name==control]) != 0){
                bvp$CKPCT[bvp$Year==years[i] & bvp$Location==locations[j]] <-
                    round(w$predicted.value/w$predicted.value[w$Name==control]*100,0)
            }else{
                bvp$CKPCT[bvp$Year==years[i] & bvp$Location==locations[j]] <- NA
            }

        }
    }


    bvpred0 <- predict(asr, classify = "Name")$pvals

    cnt=data.frame(table(phen$Name))
    bvpred0 <- merge(bvpred0,cnt,by.x='Name',by.y='Var1',all.x = T)
    bvpred0$Location <- ''
    bvpred0$Year <- ''
    if(length(w$predicted.value[w$Name==control]) != 0){
        bvpred0$CKPCT <- round(bvpred0$predicted.value/bvpred0$predicted.value[bvpred0$Name==control]*100,2)
    }else{
        bvpred0$CKPCT <- NA
    }
    bvpred0 <- bvpred0[,names(bvp)]
    bvpred0$Location <- 'All'
    bvpred0$Year <- 'All'
    brpred0 <- bvpred0[order(-bvpred0$predicted.value),]
    bvpred <- rbind(bvpred0,bvp)
    names(bvpred)[names(bvpred)=='Freq'] <- 'nobs'
    bvpred <- bvpred[,!names(bvpred) %in% c('status')]
    bvpred <- bvpred[,c(1,3,2,4:7)]
    # bvpred <- bvpred[order(bvpred$Name,bvpred$Year,bvpred$Location),]
    bvpred$phen <- NA
    head(bvpred)
    dim(bvpred)
    #=======================================================
    # round.tmp <- round(bvpred$predicted.value,2)
    # bvpred[,2] <- round.tmp
    # bvpred$std.error <- round(bvpred$std.error,0)
    bvpred$Name <- as.character(bvpred$Name)
    blue$Name <- as.character(blue$Name)

    for (i in 1:nrow(bvpred)){
        if (bvpred$Location[i]!='All'){
            w <- phen[phen$Name==bvpred$Name[i] &
                          phen$Year==as.character(bvpred$Year[i]) &
                          as.character(phen$Location)==as.character(bvpred$Location[i]),]
        } else {
            w <- phen[as.character(phen$Name)==as.character(bvpred$Name[i]),]
            #print(w)
        }
        if(nrow(w)>0)bvpred$phen[i] <- round(mean(w[[trait]],na.rm=T),2)
    }
    head(bvpred,10)
    bvpred$id <- paste0(bvpred$Name,'_',bvpred$Year,'_',bvpred$Location)
    blue$Name <- as.factor(blue$Name)

    varcomp <- summary(asr)$varcomp
    varcomp$item <- row.names(varcomp)
    varcomp$trait <- trait
    varcomp <- varcomp[,c(6,7,1:3)]
    varcomp$component <- round(varcomp$component,2)
    varcomp$std.error <- round(varcomp$std.error,2)
    varcomp$z.ratio <- round(varcomp$z.ratio,2)
    varcomp <- varcomp[varcomp$item!='units!R',]

    return(list(ebv=bvpred,varcomp=varcomp))

}

P4Estimation <- function(p4,bluelist,phenlist,traits,aoa,blupEBV,control){
    outdir <- paste0("./output/Hybrids_Anova_Res/", aoa)
    if(!dir.exists(outdir)) dir.create(outdir)

    varcomp <- c()
    for (i in 1:length(traits)){
        trait <- traits[i]
        cat('P4Estimation:  ',trait,'\n')
        p4blue <- P4Extraction(blue=bluelist[[trait]],p4,aoa=aoa)
        p4phen <- P4Extraction(blue=phenlist[[trait]],p4,aoa=aoa)

        # if(aoa == "NCSU"){
        #     p4blue$Year <- "2023C"
        #     p4phen$Year <- "2023C"
        # }

        bvalue <- hybridEvaluation(blue=p4blue,
                                   phen=p4phen,
                                   trait=trait,
                                   control = control,
                                   aoa = aoa,
                                   type = "P4")
        bvalue$ebv <- bvalue$ebv[,c('Name','Year','Location','phen','predicted.value','CKPCT','nobs','id')]
        names(bvalue$ebv)[4:7] <- paste0(c('Phen_','PredVal_','CKpct_','Obs_'),trait)

        if (i==1){
            EBV <- bvalue$ebv
        } else {
            EBV <- merge(EBV,bvalue$ebv[,4:8],by='id',all=T)
        }
        varcomp <- rbind(varcomp,bvalue$varcomp)
    }
    varcomp
    EBV <- EBV[,c(-1)]
    head(EBV)

    outlist <- list(EstimatedValue=EBV[EBV$Obs_YLD14>0,],Varcomp=varcomp)

    writexl::write_xlsx(outlist,path=paste0(outdir, '/P4_Hybrids_Evaluation_',aoa,'.xlsx'))

    ebvall <- EBV[EBV$Year=='All',]
    ebvall <- EBV[EBV$Year=='All',c('Name',paste0('CKpct_',traits))]
    ebvall$color <- ifelse(ebvall$Name %in% c('DK159','XY335','XY1483','ZD958'),3,2)
    ebvall$color <- factor(ebvall$color)
    #===============================================================
    outfile  = paste0(outdir, '/P4_YLD-MST_comparison_',aoa,'.png')
    # png(file = outfile,
        # width = 30, height = 20, units = 'cm',res = 600, bg = 'transparent')
    mstp <- ggplot(ebvall,aes(x=CKpct_MST,y=CKpct_YLD14,group=color))+
        geom_point(aes(color=color,shape=color), size = 2) +
        scale_colour_hue(name="",labels=levels(factor(ebvall$color)),l=40)  +
        geom_text(label=ebvall$Name, nudge_x=0, nudge_y=0.2,check_overlap=T,size=3)+
        ggtitle(paste0('P4:Relationship between YLD and MST at ',aoa)) +
        theme_bw(base_size = 8) +
        theme(legend.key.width = unit(0.5, "cm"))+
        theme(legend.key.height = unit(6, "mm")) +
        theme(panel.background = element_rect(fill = "white"),
              plot.margin = margin(0.4, 0.4, 0.4, 0.4, "cm"),
              axis.text.x = element_text(angle = 0, hjust = 0.9,size=7),
              plot.background = element_rect(fill = "white",colour = "black",size = 0.2),
              legend.text = element_text(size = 8),
              legend.position = 'none',
              plot.title=element_text(size=20,  #大小
                                      hjust=0.5, #位置
                                      vjust=0.5))
    print(mstp)
    # dev.off()
    #===============================================================
    #===============================================================
    for (i in 1:length(traits)){
        trait <- traits[i]
        cat('Plotting for ',trait,'\n')
        bvall <- EBV[EBV$Year=='All',c('Name','Year','Location',paste0('CKpct_',trait),
                                       paste0('Obs_',trait))]
        bvyl <- EBV[,c('Name','Year','Location',paste0('CKpct_',trait),paste0('Obs_',trait))]
        names(bvall)[4] <- 'CKpct'
        names(bvyl)[4] <- 'CKpct'
        names(bvall)[5] <- 'Obs'
        names(bvyl)[5] <- 'Obs'

        bvall <- bvall[order(-bvall$CKpct),]
        bvyl$Name <- factor(bvyl$Name,levels = bvall$Name)
        bvyl <- bvyl[bvyl$Obs!=0,]

        outfile = paste0(outdir, '/P4_LocYear_Variation_',trait,'_',aoa,'.png')

        png(file = outfile,
            width = 40, height = 20, units = 'cm',res = 600, bg = 'transparent')

        p <- ggplot(bvyl, aes(x=Name,y=CKpct,group=interaction(Name,Year))) +
            geom_point(aes(x=Name,y=CKpct,shape=Year,color=Year),size=1, fill="white") +
            xlab("") +
            facet_wrap(.~Location,scales = 'free_y') +
            scale_colour_hue(name="Year",labels=levels(factor(bvyl$Year)),l=40) +
            scale_shape_manual(name='Year',values=1:unilength(bvyl$Year)) +
            ggtitle(paste0('P4:Hybrid evaluation for ',trait,' at ',aoa)) +
            theme_bw(base_size = 8) +
            theme(legend.key.width = unit(0.5, "cm"))+
            theme(legend.key.height = unit(6, "mm")) +
            theme(panel.background = element_rect(fill = "white"),
                  plot.margin = margin(0.4, 0.4, 0.4, 0.4, "cm"),
                  axis.text.x = element_text(angle = 90, hjust = 0.9,size=7),
                  plot.background = element_rect(fill = "white",colour = "black",size = 0.2),
                  legend.text = element_text(size = 8),
                  legend.position = 'right',
                  plot.title=element_text(size=20,  #大小
                                          hjust=0.5, #位置
                                          vjust=0.5))
        print(p)
        dev.off()

        # dev.print(tiff, file = paste0('images/hybrid_perform_comparison_P4_',trait,'_',aoa,'.tiff'),
        #           res = 600, units = 'cm',width = 30, height = 20,pointsize = 5, compression="lzw")
    }
    #===============================================================
    # Comparision with BLUP EBV

    # perf <- blupEBV#data.frame(read_excel('Output/P_allTraits1026_SP_diease_YLD(1).xlsx',sheet = 3))
    # cols <- c('Name','AOA',paste0('Perf_',traits))
    # perf <- perf[perf$AOA=='MLSP',cols]
    # perf <- perf[perf$Name %in% c(p4$Name,p4$Hybrid),]
    # perf <- perf[!is.na(perf[,3]),]
    # perf$Name[perf$Name=='D1980839'] <- 'QTY812'
    # head(perf)
    # bv <- EBV[EBV$Year=='All',names(EBV) %in% c('Name',paste0("PredVal_",traits))]
    # names(bv)[2:5] <- traits
    # names(perf)[3:6] <- traits
    # perfl <- melt(perf,id.vars = 'Name',measure.vars = 3:6,variable.name = 'Trait',value.name = 'Perf')
    # bvl <- melt(bv,id.vars = 'Name',measure.vars = 2:5,variable.name = 'Trait',value.name = 'Perf')
    # bvp <- merge(bvl,perfl,by=c('Name','Trait'))
    # names(bvp)[3:4] <- c('MLM','AOV')
    # head(bvp)

    # png(file = paste0(outdir, '/P4_Model_comparison_',aoa,'.png'),
    #     width = 30, height = 20, units = 'cm',res = 600, bg = 'transparent')
    #
    # my.formula <- y ~ x
    # library(ggpmisc)
    # gg <- ggplot(bvp, aes(x=MLM, y=AOV, color=Trait)) +
    #     geom_point() +
    #     geom_smooth(formula = y ~ x, method='lm',se=FALSE) + #layer 2
    #     stat_poly_eq(formula = my.formula,
    #                  aes(label = paste( ..eq.label..,..rr.label.., sep = "~~~")),
    #                  parse = TRUE) +
    #     scale_colour_hue(name="Trait",labels=levels(factor(bvp$Trait)),l=40) +
    #     labs(title="P4:Comparison between BLUP EBVs and Hybrid Evaluation",
    #          x="BLUP EBV", y="ANOVA") + # add axis lables and plot title.
    #     theme(plot.title=element_text(size=15, face="bold"),
    #           axis.text.x=element_text(size=8),
    #           axis.text.y=element_text(size=8),
    #           axis.title.x=element_text(size=10),
    #           axis.title.y=element_text(size=10)) +
    #     # add title and axis text, change legend title.
    #     #scale_color_discrete(name="Trait") +
    #     facet_wrap( ~ Trait, ncol=2,scales = 'free')
    # print(gg)
    #
    # dev.off()

    return(outlist$EstimatedValue)
    # dev.print(tiff, file = paste0('images/hybrid_perform_BLUP_comparison_P4_',aoa,'.tiff'),
    #           res = 600, units = 'cm',width = 30, height = 20,pointsize = 5, compression="lzw")

}

P3Estimation <- function(p3,bluelist,phenlist,traits,aoa,blupEBV,control){
    outdir <- paste0("./output/Hybrids_Anova_Res/", aoa)
    if(!dir.exists(outdir)) dir.create(outdir)

    varcomp <- c()
    for (i in 1:length(traits)){
        trait <- traits[i]
        cat('P3Estimation:  ',trait,'\n')
        p3blue <- P3Extraction(bluelist,p3,trait,aoa=aoa)$blue
        p3phen <- P3Extraction(phenlist,p3,trait,aoa=aoa)$blue
        #bvalue <- hybridEvaluation(blue=p3$blue,trait=trait)
        bvalue <- hybridEvaluation(blue=p3blue,
                                   phen=p3phen,
                                   trait=trait,
                                   control=control,
                                   aoa = aoa,
                                   type = 'P3')
        bvalue$ebv <- bvalue$ebv[,c('Name','Year','Location','phen','predicted.value','CKPCT','nobs','id')]
        names(bvalue$ebv)[4:7] <- paste0(c('Phen_','PredVal_','CKpct_','Obs_'),trait)
        # names(bvalue$ebv)[4:8] <- paste0(c('PredValue_','se_','Obs_','CKpct_','phen_'),trait)
        # bvalue$ebv <- bvalue$ebv[,c(1:4,6:9)]
        if (i==1){
            EBV <- bvalue$ebv
        } else {
            EBV <- merge(EBV,bvalue$ebv[,4:8],by='id',all=T)
        }
        varcomp <- rbind(varcomp,bvalue$varcomp)
    }
    varcomp
    EBV <- EBV[,c(-1)]
    head(EBV)

    outlist <- list(EstimatedValue=EBV[EBV$Obs_YLD14>0,],Varcomp=varcomp)

    writexl::write_xlsx(outlist,path=paste0(outdir, '/P3_Hybrids_Evaluation_',aoa,'.xlsx'))


    ebvall <- EBV[EBV$Year=='All',]
    ebvall <- EBV[EBV$Year=='All',c('Name',paste0('CKpct_',traits))]
    ebvall$color <- ifelse(ebvall$Name %in% c('DK159','XY335','XY1483','ZD958'
                                              ,'CD99','ZY335','DH605','XY1446','DMY1','DMY3','C1563'),3,2)
    ebvall$color <- factor(ebvall$color)

    outfile = paste0(outdir, '/P3_YLD-MST_comparison_',aoa,'.png')
    # png(file = outfile,
        # width = 30, height = 20, units = 'cm',res = 600, bg = 'transparent')

    mstp <- ggplot(ebvall,aes(x=CKpct_MST,y=CKpct_YLD14,group=color))+
        geom_point(aes(color=color,shape=color), size = 2) +
        scale_colour_hue(name="",labels=levels(factor(ebvall$color)),l=40)  +
        geom_text(label=ebvall$Name, nudge_x=0, nudge_y=0.1,check_overlap=T,size=2)+
        ggtitle(paste0('P3:Relationship between YLD and MST at ',aoa)) +
        theme_bw(base_size = 8) +
        theme(legend.key.width = unit(0.5, "cm"))+
        theme(legend.key.height = unit(6, "mm")) +
        theme(panel.background = element_rect(fill = "white"),
              plot.margin = margin(0.4, 0.4, 0.4, 0.4, "cm"),
              axis.text.x = element_text(angle = 0, hjust = 0.9,size=7),
              plot.background = element_rect(fill = "white",colour = "black",size = 0.2),
              legend.text = element_text(size = 8),
              legend.position = 'none',
              plot.title=element_text(size=20,  #大小
                                      hjust=0.5, #位置
                                      vjust=0.5))
    print(mstp)
    # dev.off()

    # dev.print(tiff, file = paste0('images/hybrid_perform_comparison_P3_',aoa,'.tiff'), res = 600, units = 'cm',
    #           width = 30, height = 20,pointsize = 5, compression="lzw")
    #

    for (i in 1:length(traits)){
        trait <- traits[i]
        cat('Plotting for ',trait,'\n')

        bvall <- EBV[EBV$Year=='All',c('Name','Year','Location',paste0('CKpct_',trait),
                                       paste0('Obs_',trait))]
        bvyl <- EBV[,c('Name','Year','Location',paste0('CKpct_',trait),paste0('Obs_',trait))]
        names(bvall)[4] <- 'CKpct'
        names(bvyl)[4] <- 'CKpct'
        names(bvall)[5] <- 'Obs'
        names(bvyl)[5] <- 'Obs'

        bvall <- bvall[order(-bvall$CKpct),]
        bvyl$Name <- factor(bvyl$Name,levels = bvall$Name)
        bvyl <- bvyl[bvyl$Obs!=0,]

        outfile = paste0(outdir, '/P3_LocYear_Variation_',trait,'_',aoa,'.png')

        png(file = outfile,
            width = 40, height = 20, units = 'cm',res = 600, bg = 'transparent')

        p <- ggplot(bvyl, aes(x=Name,y=CKpct,group=interaction(Name,Year))) +
            geom_point(aes(x=Name,y=CKpct,shape=Year,color=Year),size=2, fill="white") +
            xlab("") +
            facet_wrap(.~Location,scales = 'free_y') +
            scale_colour_hue(name="Year",labels=levels(factor(bvyl$Year)),l=40) +
            scale_shape_manual(name='Year',values=1:unilength(bvyl$Year)) +
            ggtitle(paste0('P3:Hybrid evaluation for ',trait,' at ',aoa)) +
            theme_bw(base_size = 8) +
            theme(legend.key.width = unit(0.5, "cm"))+
            theme(legend.key.height = unit(6, "mm")) +
            theme(panel.background = element_rect(fill = "white"),
                  plot.margin = margin(0.4, 0.4, 0.4, 0.4, "cm"),
                  axis.text.x = element_text(angle = 90, hjust = 0.9,size=7),
                  plot.background = element_rect(fill = "white",colour = "black",size = 0.2),
                  legend.text = element_text(size = 8),
                  legend.position = 'right',
                  plot.title=element_text(size=20,  #大小
                                          hjust=0.5, #位置
                                          vjust=0.5))
        print(p)
        dev.off()
        # dev.print(tiff, file = paste0('images/hybrid_perform_comparison_P3_',trait,'_',aoa,'.tiff'),
        #           res = 600, units = 'cm',width = 30, height = 20,pointsize = 5, compression="lzw")
    }



    # perf <- blupEBV
    # cols <- c('Name','AOA',paste0('Perf_',traits))
    # perf <- perf[perf$AOA=='MLSP',cols]
    # perf <- perf[perf$Name %in% p3,]
    # perf <- perf[!is.na(perf[,3]),]
    #
    # head(perf)
    # bv <- EBV[EBV$Year=='All',names(EBV) %in% c('Name',paste0("PredVal_",traits))]
    # names(bv)[2:5] <- traits
    # names(perf)[3:6] <- traits
    # perfl <- melt(perf,id.vars = 'Name',measure.vars = 3:6,variable.name = 'Trait',value.name = 'Perf')
    # bvl <- melt(bv,id.vars = 'Name',measure.vars = 2:5,variable.name = 'Trait',value.name = 'Perf')
    # bvp <- merge(bvl,perfl,by=c('Name','Trait'))
    # names(bvp)[3:4] <- c('MLM','AOV')
    # head(bvp)
    #
    # png(file = paste0(outdir, '/P3_Model_comparison_',aoa,'.png'),
    #     width = 30, height = 20, units = 'cm',res = 600, bg = 'transparent')
    #
    # my.formula <- y ~ x
    # gg <- ggplot(bvp, aes(x=MLM, y=AOV, color=Trait)) +
    #     geom_point() +
    #     geom_smooth(formula = y ~ x, method='lm',se=FALSE) + #layer 2
    #     stat_poly_eq(formula = my.formula,
    #                  aes(label = paste( ..eq.label..,..rr.label.., sep = "~~~")),
    #                  parse = TRUE) +
    #     scale_colour_hue(name="Trait",labels=levels(factor(bvp$Trait)),l=40) +
    #     labs(title="P3:Comparison between BLUP EBVs and Hybrid Evaluation",
    #          x="BLUP EBV", y="ANOVA") + # add axis lables and plot title.
    #     theme(plot.title=element_text(size=15, face="bold"),
    #           axis.text.x=element_text(size=8),
    #           axis.text.y=element_text(size=8),
    #           axis.title.x=element_text(size=10),
    #           axis.title.y=element_text(size=10)) +
    #     # add title and axis text, change legend title.
    #     #scale_color_discrete(name="Trait") +
    #     facet_wrap( ~ Trait, ncol=2,scales = 'free')
    # print(gg)
    #
    # dev.off()

    return(outlist$EstimatedValue)
    # dev.print(tiff, file = paste0('images/hybrid_perform_BLUP_comparison_P3_',aoa,'.tiff'),
    #           res = 600, units = 'cm',width = 30, height = 20,pointsize = 5, compression="lzw")

}

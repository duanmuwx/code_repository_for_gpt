testfun <- function(){
    #=================================================================================
    #测试所有已经导出的函数文档说明
    ?dataOutliers_Identifying
    ?dataOutliers.for

    ?missingAdding
    ?replicationCheck
    ?blueEstimation

    ?extracting
    ?GCA_SCA_Estimation
    ?programRun

    ?mgxeEstimation
    ?mulTraitEstimation
    #=================================================================================
    #测试数据过滤
    #Test data outlier
    rm(list = ls())
    data("phenotypes_P_alltrails")
    w <- phdemo$SpML

    outliers.Res <- dataOutliers_Identifying(data=w,trait ="MST",model='3SD') #|v| Test example

    outliers.Res <- dataOutliers_Identifying(data=su,trait ="MST",model='3SD')  #|v| Test su data
    res <- dataOutliers.for(data=su,traits =c('MST',"YLD","PHT","EHT"),model='3SD')

    outliers.Res <- dataOutliers_Identifying(data=su,trait ="MST",model='3SD') #|v| Test ml data
    res <- dataOutliers.for(data=su,traits =c('MST',"YLD","PHT","EHT"),model='3SD')
    #=================================================================================
    #测试空间校正（BLUE分析）
    rm(list = ls())
    data("phenotypes_P_alltrails")
    trait  <- "MST"
    outliers.Res <- dataOutliers_Identifying(data = phdemo$SpML,
                                             trait = trait,
                                             model = '3SD')
    dat    <- outliers.Res$datclean
    rm.idx <- unique(unlist(apply(dat, 2, FUN = function(x){return(which(is.na(x)))})))
    dat    <- dat[-rm.idx, ]
    #空间分析的基本单位是小区，需要各选择1种类型的Location，Field，Type
    dat    <- dat[dat$Location==unique(dat$Location)[1],]
    dat    <- dat[dat$Field==unique(dat$Field)[1],]
    datdf  <- dat[dat$Type==unique(dat$Type)[1],]
    print(table(datdf$Range,datdf$Pass))
    names(datdf)[names(datdf)==trait]<-"trait"
    repct<-replicationCheck(df=datdf)
    repct
    datdf <- missingAdding(w=datdf)
    datdf <- datdf[order(datdf$Range,datdf$Pass),]
    print(table(datdf$Range,datdf$Pass))
    blueout<-blueEstimation(datdf,repct,trait)
    #=================================================================================
    #测试BLUP分析
    rm(list = ls())
    data("TestAinv")
    data("TestPed")
    data("TestBlueList")
    names(ped)[names(ped)=='FGenoID'] <- 'DAM'
    names(ped)[names(ped)=='MGenoID'] <- 'SIRE'
    head(ped)
    bluedat <- lapply(bluelist, function(x) x$blue)
    regions <- names(bluedat)
    regions
    #[1] "EM" "ML" "SU" "SW"
    ebv <- list()
    for (i in 1:4){
        reg <- regions[i]
        if (reg=='EM') traits <- c('MST', 'YLD', 'PHT', 'EHT')#"LDG","NCLB"
        if (reg=='SpML') traits <- c('MST', 'TWT', 'YLD', 'PHT', 'EHT')#"GLS","PMD"
        if (reg=='SU') traits <- c('MST', 'TWT', 'YLD', 'PHT', 'EHT')#"PMD"
        if (reg=='SW') traits <- c('MST', 'YLD', 'PHT', 'EHT')#"PMD","ER"
        print(traits)
        ebv[[reg]] <- programRun(bluelist=bluedat[[reg]], ainv, ped, traits)
        print(names(ebv[[reg]]))
        #writexl::write_xlsx(ebv[[reg]], path=paste0('Output/blup_EBV_', reg, '.xlsx'))
    }
    #saveRDS(ebv, 'Output/blup_breedingValues_allRegions.rds')
    #=================================================================================
    #测试多性状联合分析,表型遗传相关
    rm(list = ls())
    data("TestAinv")
    data("TestPed")
    data("TestBlueList")
    filter.res <- dataOutliers.for(data   = phdemo$SpML,
                            traits = c('MST',"YLD","PHT","EHT"),
                            model  = '3SD')
    testCor <- mulTraitEstimation(datdf = filter.res$result_df,
                                  trait.vec = c("YLD", "MST.x", "MST.y", "PHT", "EHT"),
                                  ainv = ainv, ped = ped)
    testCor
    #=================================================================================
    #多点分析（GxE），需要将不同的Location进行分组
    rm(list = ls())
    data("TestAinv")
    data("TestPed")
    data("TestBlueList")
    filter.res <- dataOutliers.for(data   = phdemo$SpML,
                                   traits = c('MST',"YLD","PHT","EHT"),
                                   model  = '3SD')
    grp.by  <- list(set1 = c("6HLJSC", "6JLDF", "6JLDH"),
                    set2 = c("6JLGZL", "7LNTL", "7NMTL"))
    for (tt in c("YLD","PHT", "EHT")) {
        message(paste0("==================== ", tt, " ===================="))
        testGxE <- mgxeEstimation(datdf = filter.res$result_df,
                                  trait = tt,
                                  ped   = ped,
                                  ainv  = ainv,
                                  loc.grp.by = grp.by, verbose = T)
    }
    plot(testGxE$mgxe$F.EBV.set1$EBV, testGxE$mgxe$F.EBV.set2$EBV)
    plot(testGxE$mgxe$M.EBV.set1$EBV, testGxE$mgxe$M.EBV.set2$EBV)
    return(NULL)
}

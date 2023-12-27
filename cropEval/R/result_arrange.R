#' Title
#'
#' @param parent a data.frame of parents
#'
#' @return result a data.frame of number of hybrid observations
#' @export
#'
#'
#'
Hybrid_Observations <- function(parent){

  #pheno <- readRDS("datafiles/pheno.rds")  ###
  result <- list()

  for (i in 1:length(pheno)) {
    sub_list <- pheno[[i]]
    region_name = names(pheno)[i]

    pheno1 <- sub_list[!(sub_list$IsCK == "CK" & !is.na(sub_list$IsCK)), ]
    #---pheno2<- pheno1[!is.na(pheno1$FGenoID),]
    pheno2<- pheno1[!is.na(pheno1[[parent]]),]
    Geno <- unique(pheno2[[parent]])
    Geno_count <- data.frame(GenoID = character(),
                             HyridCount = numeric())

    for (j in seq_along(Geno)) {
      name=Geno[j]
      ZMN<-pheno2[pheno2[[parent]]==name,]
      ZMN<-ZMN[!is.na(ZMN$Name),]

      hyrid<-length(unique(ZMN$Name))

      temp <- data.frame(GenoID = name, HyridCount = hyrid)

      Geno_count <- rbind(Geno_count, temp)
      #print(FGeno_count)
      #cat("\n")
    }

    result[[region_name]] <- Geno_count
  }

  #writexl::write_xlsx(result, paste0("Output/", parent, "_HyridCount.xlsx"))
  return(result)


}


#' Collation of genetic evaluation results
#'
#' @param popmean Population average
#' @param region  Maize varieties
#' @param traits  phenotypic traits
#' @param ctrlnsample  Reference male parent
#' @param ctrlssample Reference female parent
#' @param adv advanced individuals in 2022
#' @param hybobs No. of hybrid observations
#' @param blupebv breeding value
#' @param Aov Variance test table
#'
#' @return outebv A list of data frames
#'
#' @examples
#' \donttest{
#' rm(list = ls())
#' library(readxl)
#' library(writexl)
#' library(biosim)
#' library(tidyverse)
#' library(reshape2)
#' library(stringi)
#' library(ggplot2)
#' library(ggpmisc)
#' region = 'SU'
#' traits <- c('YLD','MST','TWT','PHT','EHT','PMD','LDG')
#' ctrlnsample = 'Chang72'  # Reference male parent
#' ctrlssample = 'Zheng58'  # Reference female parent
#' del_result_list <- del_result(popmean, region, traits, ctrlnsample, ctrlssample, adv, hybobs, blupebv, Aov)
#' print(names(del_result_list))
#' }
#'
#' @export
#'
del_result <- function(popmean, region, traits, ctrlnsample, ctrlssample, adv, hybobs, blupebv, Aov) {

  traits_1 = unlist(lapply(traits, function(x){c(x, str_c(x, 'acc'), str_c(x, 'CKpct'))}))

  popmean <- data.frame(popmean)

  #popmean <- data.frame(read.table(textConnection("Trait  ML   SU
    #MST  29.31  23.76
    #TWT  73.12  71.94
    #YLD  774.76 606.97
    #PHT  337.58 288.06
    #EHT  131.50  110.34
    #PMD  6.628126  7.62
    #LDG  8.239709  7.56
    #NCLB 3.276354  NA
    #GLS  5.742301  NA"),
    #header = T, sep = ""))

  if('ML' %in% names(popmean)){
    popmean <- popmean %>% rename(SpML = ML)
  }

  ssebv <- Hybrid_Observations(parent = 'FGenoID')[[region]]  # 母本
  nsebv  <- Hybrid_Observations(parent = 'MGenoID')[[region]] # 父本

  # breeds
  # ---------------------------------------------------------------
  #--- traits <- c("YLD", "MST", "TWT", "PHT", "EHT", "PMD", "LDG", "NCLB", "GLS")

  #--names(blupebv)[1] <- "ML"

  # Extract GCA EBVs for NS and SS inbred lines for ML breeds
  #----  ns <- blupebv$ML$NS_GCA
  #----- ss <- blupebv$ML$SS_GCA
  ns <- blupebv[[region]][['NS_GCA']]
  ss <- blupebv[[region]][['SS_GCA']]

  names(ss)[1] <- "GenoID"
  names(ns)[1] <- "GenoID"

  #mlaov <- aov$SpML
  #suaov <- aov$SU
  mlaov <- Aov[[region]]

  nsd <- mlaov$Mebv
  ssd <- mlaov$Febv
  names(ssd)[1] <- "GenoID"
  names(nsd)[1] <- "GenoID"
  #---- head(nsd)

  nsd <- nsd[, !names(nsd) %in% paste0(traits, "se")]
  ssd <- ssd[, !names(ssd) %in% paste0(traits, "se")]

  #for (i in 2:9) {
  #nsd[, i] <- round(nsd[, i], 3)
  #ssd[, i] <- round(ssd[, i], 3)
  #}


  ssd <- ssd %>% mutate(across(where(is.numeric), ~round(.x, 3)))
  nsd <- nsd %>% mutate(across(where(is.numeric), ~round(.x, 3)))

  traitsdisease <- traits[traits %in% c("PMD", "LDG", "NCLB", "GLS")]
  #---ctrln <- nsd[nsd$GenoID == "ZMN00545", c("GenoID", traits[6:9])]
  #---ctrls <- ssd[ssd$GenoID == "ZMN00080", c("GenoID", traits[6:9])]
  ctrln <- nsd[nsd$GenoID == ctrlnsample, c("GenoID", traitsdisease)]
  ctrls <- ssd[ssd$GenoID == ctrlssample, c("GenoID", traitsdisease)]


  #---for (i in 6:length(traits)) {
  for (i in 1:length(traitsdisease)) {
    #---mu <- popmean$ML[popmean$Trait == traits[i]]
    mu <- popmean[[region]][popmean$Trait == traitsdisease[i]]
    mu <- as.numeric(mu)
    nsd[[paste0(traitsdisease[i], "CKpct")]] <- round((nsd[[traitsdisease[i]]] + mu)/(ctrln[[traitsdisease[i]]] + mu) * 100, 0)
    ssd[[paste0(traitsdisease[i], "CKpct")]] <- round((ssd[[traitsdisease[i]]] + mu)/(ctrls[[traitsdisease[i]]] + mu) * 100, 0)
  }


  scad <- mlaov$SCAebv
  scad <- scad[, !names(scad) %in% paste0(traits, "se")]
  scad <- scad[, !names(scad) %in% paste0(traits, "acc")]

  #---for (i in 5:8) {
  #---scad[, i] <- round(scad[, i], 3)
  #--}

  for (i in traitsdisease) {
    scad[, i] <- round(scad[, i], 3)
  }


  names(scad)[names(scad) == "SIRE"] <- "MGenoID" # 父本
  names(scad)[names(scad) == "DAM"] <- "FGenoID"  # 母本
  names(scad)[names(scad) == "Fam"] <- "GenoID"
  scad$GenoID <- paste0(scad$FGenoID, "/", scad$MGenoID)
  scad <- scad[scad$GenoID != "0/0", ]


  #---ctrl <- scad[scad$GenoID == "ZMN00080/ZMN00545", c(traits[6:9])]
  ctrl <- scad[scad$GenoID == paste0(ctrlssample, '/', ctrlnsample), traitsdisease]
  ctrl

  #---for (i in 6:length(traits)) {
  for (i in 1:length(traitsdisease)) {
    #---mu <- popmean$ML[popmean$Trait == traits[i]]
    mu <- popmean[[region]][popmean$Trait == traitsdisease[i]]
    mu <- as.numeric(mu)
    scad[[paste0("CKpct_", traitsdisease[i])]] <- round((scad[[traitsdisease[i]]] +
                                                           mu)/(ctrl[[traitsdisease[i]]] + mu) * 100, 0)
  }

  #---for (i in 6:length(traits)) {
  for (i in 1:length(traitsdisease)) {
    names(scad)[names(scad) == traitsdisease[i]] <- paste0(traitsdisease[i], "_SCA")
  }

  #---nsd_gca <- nsd[, c("GenoID", traits[6:9])]
  nsd_gca <- nsd[, c("GenoID", traitsdisease)]

  #---names(nsd_gca)[2:5] <- paste0(names(nsd_gca)[2:5], "_GCA_MPar")
  names(nsd_gca)[2:(length(traitsdisease) + 1)] <- paste0(names(nsd_gca)[2: (length(traitsdisease) + 1)], "_GCA_MPar")

  #--ssd_gca <- ssd[, c("GenoID", traits[6:9])]
  ssd_gca <- ssd[, c("GenoID", traitsdisease)]

  #---names(ssd_gca)[2:5] <- paste0(names(ssd_gca)[2:5], "_GCA_FPar")
  names(ssd_gca)[2: (length(traitsdisease) + 1)] <- paste0(names(ssd_gca)[2: (length(traitsdisease) + 1)], "_GCA_FPar")

  scad <- merge(scad, ssd_gca, by.x = "FGenoID", by.y = "GenoID", all.x = T)
  scad <- merge(scad, nsd_gca, by.x = "MGenoID", by.y = "GenoID", all.x = T)
  scad <- scad[, !names(scad) %in% c("FamNum", "FGenoID", "MGenoID")]

  #for (i in 6:9) {
  for(i in 1:length(traitsdisease)){
    names(scad)[names(scad) == paste0(traitsdisease[i], "_SCA")] <- paste0("SCA_",
                                                                           traitsdisease[i])
    names(scad)[names(scad) == paste0(traitsdisease[i], "_GCA_FPar")] <- paste0("FGenoID_GCA_",
                                                                                traitsdisease[i])
    names(scad)[names(scad) == paste0(traitsdisease[i], "_GCA_MPar")] <- paste0("MGenoID_GCA_",
                                                                                traitsdisease[i])
    names(scad)[names(scad) == paste0(traitsdisease[i], "_perf")] <- paste0("Perf_",
                                                                            traitsdisease[i])
    names(scad)[names(scad) == paste0(traitsdisease[i], "_CKpct")] <- paste0("CKpct_",
                                                                             traitsdisease[i])
  }


  #----for (i in 6:length(traits)) {
  for(i in 1:length(traitsdisease)){
    scad[paste0("Perf_", traitsdisease[i])] <- scad[paste0("MGenoID_GCA_",
                                                           traitsdisease[i])] + scad[paste0("FGenoID_GCA_", traitsdisease[i])] + scad[paste0("SCA_",
                                                                                                                                             traitsdisease[i])]
  }


  # Extract GCA EBVs for NS and SS inbred lines for ML breeds

  #--ns <- merge(ns, nsebv[, c(1, 5)], by = "GenoID", all.x = T)
  ns <- merge(ns, nsebv[, c(1, 2)], by = "GenoID", all.x = T)

  names(ns)[names(ns)=='HyridCount'] <- "nHybrids"

  #--names(ns)[17] <- c("nHybrids")

  #--ss <- merge(ss, ssebv[, c(1, 5)], by = "GenoID", all.x = T)
  ss <- merge(ss, ssebv[, c(1, 2)], by = "GenoID", all.x = T)

  #--names(ss)[17] <- c("nHybrids")
  names(ss)[names(ss)=='HyridCount'] <- "nHybrids"

  names(ns)[names(ns) == "YLD_acc"] <- "YLDacc"
  names(ss)[names(ss) == "YLD_acc"] <- "YLDacc"

  #---ss <- ss[, c(1, 17, 2:7, 14:16, 8:13)]
  ss <- ss[ ,c('GenoID' ,'nHybrids' ,traits_1[traits_1 %in% names(ss)[2:(length(names(ss))-1)]])]
  #---ns <- ns[, c(1, 17, 2:7, 14:16, 8:13)]
  ns <- ns[ ,c('GenoID' ,'nHybrids' ,traits_1[traits_1 %in% names(ns)[2:(length(names(ns))-1)]])]

  ns <- merge(ns, nsd, by = "GenoID")
  ss <- merge(ss, ssd, by = "GenoID")

  for (i in 1:length(traits)) {
    #--mu <- popmean$ML[popmean$Trait == traits[i]]
    mu <- popmean[[region]][popmean$Trait == traits[i]]
    mu <- as.numeric(mu)
    ns[[traits[i]]] <- ns[[traits[i]]] + mu
    ss[[traits[i]]] <- ss[[traits[i]]] + mu
  }

  ns <- ns[order(-ns$YLDCKpct), ]
  ss <- ss[order(-ss$YLDCKpct), ]

  names(ns)[names(ns) == "YLD_acc"] <- "YLDacc"
  names(ss)[names(ss) == "YLD_acc"] <- "YLDacc"

  colname <- c("GenoID", "nHybrids")
  cn <- c("", "acc", "CKpct")

  #---for (i in 1:9) {
  for(i in 1:length(traits)){
    colname <- c(colname, paste0(traits[i], cn))
  }
  colname

  ns <- ns[, colname]
  ss <- ss[, colname]


  # Extract SCA breeding values for ML breeds
  #---sca <- blupebv$ML$Cross_SCA  #ebv[,c('Fam','SIRE','DAM',traits,paste0(traits,'acc'))]
  sca <- blupebv[[region]][['Cross_SCA']]
  names(sca)[names(sca) == "Name"] <- c("Hybrid")
  names(sca)[names(sca) == "Fam"] <- c("GenoID")
  sca$GenoID <- paste0(sca$FGenoID, "/", sca$MGenoID)
  dim(sca)  #1867
  sca <- unique(sca)
  dim(sca)

  sca <- sca[order(-sca$YLD_CKpct), ]

  # Merge SCA Bvs and No of hybrid observations
  #---hobs <- hybobs$SpML
  hobs <- hybobs[[region]]
  hobs <- hobs[, c("Hybrid", traits)]
  names(hobs)
  #---names(hobs)[2:10] <- paste0("nHybObs_", traits)
  names(hobs)[2:(length(traits) + 1)] <- paste0("nHybObs_", traits)
  dim(hobs)  #1859
  sca <- merge(sca, hobs, by = "Hybrid", all.x = T, all.y = T)

  for (i in 1:length(traits)) {
    names(sca)[names(sca) == paste0(traits[i], "_SCA")] <- paste0("SCA_",
                                                                  traits[i])
    names(sca)[names(sca) == paste0(traits[i], "_GCA_FPar")] <- paste0("FGenoID_GCA_",
                                                                       traits[i])
    names(sca)[names(sca) == paste0(traits[i], "_GCA_MPar")] <- paste0("MGenoID_GCA_",
                                                                       traits[i])
    names(sca)[names(sca) == paste0(traits[i], "_perf")] <- paste0("Perf_",
                                                                   traits[i])
    names(sca)[names(sca) == paste0(traits[i], "_CKpct")] <- paste0("CKpct_",
                                                                    traits[i])
  }

  sca <- merge(sca, scad, by = "GenoID")
  # Re-arrange columns
  colname <- c()
  cn <- c("FGenoID_GCA_", "MGenoID_GCA_", "SCA_", "Perf_", "CKpct_",
          "nHybObs_")
  for (i in 1:length(traits)) {
    colname <- c(colname, paste0(cn, traits[i]))
  }

  sca <- sca[, c("Hybrid", "GenoID", "FGenoID", "MGenoID", colname)]

  sca <- merge(sca, adv[, 1:2], by = "Hybrid", all.x = T)   ###
  sca$Adv22 <- ifelse(!is.na(sca$ZMN_ID), "YES", "")
  sca <- sca[order(-sca$CKpct_YLD), ]

  #---h2 <- rbind(blupebv$ML$h2, mlaov$h2)
  #---d2 <- rbind(blupebv$ML$d2, mlaov$d2)
  #---varcomp <- rbind(blupebv$ML$varcomp, mlaov$varcomp)
  h2 <- rbind(blupebv[[region]]$h2, mlaov$h2)
  d2 <- rbind(blupebv[[region]]$d2, mlaov$d2)
  varcomp <- rbind(blupebv[[region]]$varcomp, mlaov$varcomp)
  h2$model <- c(rep("BLUP", length(traits) - length(traitsdisease)), rep("ANOVA", length(traitsdisease)))
  d2$model <- c(rep("BLUP", length(traits) - length(traitsdisease)), rep("ANOVA", length(traitsdisease)))

  outebv <- list(SS_GCA = ss, NS_GCA = ns, Cross_SCA = sca, popMean = popmean,
                 h2 = h2, d2 = d2, varcomp = varcomp, Adv22 = adv)

  return (outebv)
  #writexl::write_xlsx(outebv, path = "Output/BLUPbreedingValues_GCA-SCA_SpML_allTraits_20230801.xlsx")

}



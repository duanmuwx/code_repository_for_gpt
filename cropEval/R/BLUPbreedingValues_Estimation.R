#' Extracting some important statistics
#'
#' The function extracts some important statistics from the model,
#' such as variance component, standard error, etc
#' and returns a data frame containing some important statistics.
#'
#' @param asr  Established model.
#' @param var  Variance component corresponding to factor.
#' @param individual  Samples corresponding to factors.
#' @param term  The name of each factor displayed in the model results.
#' @param gxe   Identifier.
#'
#' @return A data frame with important statistics.
#'
#' @export extracting

extracting <- function(asr, var, individual, term, gxe=FALSE) {
  # Extract some important statistics from the model,
  # such as variance component, standard error, etc.
  ebv <- summary(asr, coef=T)$coef.random
  ebv <- as.data.frame(ebv)
  bvalue <- as.data.frame(ebv[grepl(term,row.names(ebv)),])

  head(bvalue)
  dim(bvalue)
  bvalue$id <- 0
  for (i in 1:dim(bvalue)[1]){
    s <- strsplit(row.names(bvalue)[i],"[_]")[[1]]
    if (gxe){
      bvalue$id[i] <- s[3]
    } else {
      bvalue$id[i] <- s[2]
    }
  }
  v <- (1-bvalue$'std.error'^2/var)
  bvalue$rel <- ifelse(v<=0,0,v)
  bvalue$acc <- sqrt(bvalue$rel)

  bvalue <- bvalue[bvalue$id %in% individual,]
  bvalue <- bvalue[,c(4,1,2,5,6)]
  row.names(bvalue) <- seq(1:nrow(bvalue))
  names(bvalue) <- c("ID","EBV","se","rel","acc")
  print(head(bvalue))

  return(bvalue)

}


#' Establishing multiple mixed linear models
#'
#' This function generates additive effect, non-additive effect,
#' generalized heritability and so on according to the blue value,
#' pedigree file and inverse matrix of genetic relationship.
#'
#' @param bluehb Data frame with blue value.
#' @param ainv Inverse matrix of kinship.
#' @param ped Genealogical file of the sample.
#'
#' @return A list with model results.
#'
#' @export GCA_SCA_Estimation

GCA_SCA_Estimation <- function(bluehb, ainv, ped){

    bluehb <- merge(bluehb, ped, by='Name',all.x=T) # Merge Blue files and genealogy files by using left join

    bluehb$SIRE <- as.character(bluehb$SIRE)
    bluehb$DAM <- as.character(bluehb$DAM)
    bluehb$SIRE[is.na(bluehb$SIRE)] <- '0' # Missing values are represented by 0
    bluehb$DAM[is.na(bluehb$DAM)] <- '0'   # Missing values are represented by 0
    bluehb$SIRE <- as.factor(bluehb$SIRE)
    bluehb$DAM <- as.factor(bluehb$DAM)


    bluehb$Fam <- paste0(bluehb$SIRE, '_', bluehb$DAM) # Connecting parent sample and parent sample.
    bluehb$Fam[bluehb$Fam=='NA_NA'] <- '0_0'           # Missing values are represented by 0_0.

    # Group according to paternal group and maternal group, and mark with numbers.
    fam <- unique(bluehb[,c('Fam','SIRE','DAM')])
    fam$FamNum <- 1:nrow(fam)
    bluehb <- merge(bluehb, fam[,c('Fam','FamNum')], by='Fam')
    head(bluehb)
    #writecsv(bluehb,'temp.csv')

    bluehb$FamNum <- as.factor(bluehb$FamNum)
    bluehb$Location  <- as.factor(bluehb$Location )
    bluehb$Field  <- as.factor(bluehb$Field)
    bluehb$LF <- paste0(bluehb$Location,'_',bluehb$Field) # Combine Location and Field as one factor to get LF.
    bluehb$LF  <- as.factor(bluehb$LF)

    # LF as a fixed factor, maternal and paternal and maternal_paternal as a random factor,
    # establish a linear model
    asr2 <- asreml(predicted.value ~ LF,
                   random = ~ vm(SIRE,ainv) + vm(DAM,ainv) + idv(FamNum),
                   residual = ~idv(units),
                   data = bluehb)

    # Extracting variance components from the model.
    fixed<-asr2$coefficients$fixed
    popmean <- fixed[row.names(fixed)=='(Intercept)',]

    # The generalized heritability(H2) is calculated according to the variance component of each factor.
    summary(asr2)$varcomp
    vs <- summary(asr2)$varcomp$component[2]
    vd <- summary(asr2)$varcomp$component[3]
    vf <- summary(asr2)$varcomp$component[1]
    h2 <- vpredict(asr2,h2~(V2+V3)/(V1+V2+V3+V4))
    print(h2)
    d2 <- vpredict(asr2,h2~V1/(V1+V2+V3+V4))
    print(d2)

    # The maximum likelihood value is obtained by establishing the model,
    # and whether the factor(FamNum) is significant is tested.
    loglik2 <- asr2$loglik
    asr3 <- asreml(predicted.value ~ LF,
                   random = ~ vm(SIRE,ainv) + vm(DAM,ainv),# + idv(FamNum),
                   residual = ~idv(units),
                   data = bluehb)
    loglik3 <- asr3$loglik
    d2$LRTvalue <- (loglik2-loglik3)*2


    bvalue <- summary(asr2,coef=T)$coef.random
    head(bvalue)
    #write.csv(bvalue,'temp.csv')
    ind <- unique(c(bluehb$SIRE,bluehb$DAM))
    ebv <- extracting(asr2,var=vs,individual = ind,term='vm(Name, ainv)*')
    head(ebv)

    # Extracting breeding value from the model, that is, non-additive effect
    febv <- extracting(asr2,var=vf,individual = fam$FamNum,term='FamNum*')
    febv <- merge(fam,febv,by.x='FamNum',by.y='ID')
    febv <- febv[order(-febv$EBV),]
    head(febv)
    #write.csv(febv,'SCA.csv')

    # The breeding value is extracted from the model, that is, the additive effect of the maternal line.
    ind <- unique(c(bluehb$SIRE))
    sebv <- extracting(asr2,var=vs,individual = ind,term='vm(SIRE, ainv)*')
    sebv <- sebv[sebv$EBV!=0,]
    sebv$par='SIRE'
    #sebv <- sebv[1:315,]

    # The breeding value is extracted from the model, that is, the additive effect of the paternal line.
    ind <- unique(c(bluehb$DAM))
    debv <- extracting(asr=asr2,var=vd,individual = ind,term='vm(DAM, ainv)*')
    debv <- debv[debv$EBV!=0,]
    debv$par='DAM'

    sebv <- sebv[!substr(row.names(sebv),1,6)=='vm(DAM',]
    debv <- debv[!substr(row.names(debv),1,7)=='vm(SIRE',]

    outlist <- list(h2=h2,d2=d2,varcomp=summary(asr2)$varcomp,Mebv=sebv,Febv=debv,
                    SCAebv=febv,popmean=popmean)
    return(outlist)
}

#' Calculating the statistics of all characters to be processed
#'
#' This function processes each trait in turn, returns their respective
#' statistical values and summarizes them in a list.
#'
#' @param bluelist List with blue value.
#' @param ainv Inverse matrix of kinship.
#' @param ped Genealogical file of the sample.
#' @param traits All characters to be treated.
#'
#' @return A list containing statistics of all traits.
#'
#' @examples
#' \donttest{
#' library(biosim)
#' library(matrixcalc)
#' library(asreml)
#' library(MCMCglmm)
#' names(ped)[names(ped)=='FGenoID'] <- 'DAM'
#' names(ped)[names(ped)=='MGenoID'] <- 'SIRE'
#' print(head(ped))
#' bluedat <- lapply(bluedat, function(x) x$blue)
#' traits <- c('MST', 'TWT', 'YLD', 'PHT', 'EHT')
#' ebv <- programRun(bluedat[['SpML']], ainv, ped, traits)
#' print(names(ebv))
#' }
#'
#' @export programRun

programRun <- function(bluelist, ainv, ped, traits){
    h2 <- c()
    d2 <- c()
    varc <- c()
    SCAebv <- c()
    popMean <- data.frame(trait=traits,popmean=0)

    Mebv <- c()
    Febv <- c

    # Processing each character in turn and calculating the breeding value
    # Modify the name, merge the data of each character and output it to a list.
    for (i in 1:length(traits)){
        trait <- traits[i]
        blue <- bluelist[[traits[i]]]
        asrout <- GCA_SCA_Estimation(bluehb=blue,ainv=ainv,ped=ped)

        popMean$popmean[i] <- asrout$popmean

        asrout$h2
        asrout$d2
        asrout$varcomp
        if (i==1){
            SCAebv <- asrout$SCAebv[,c(1:6,8)]
            Mebv <- asrout$Mebv[,c(1:3,5)]
            Febv <- asrout$Febv[,c(1:3,5)]
            names(SCAebv)[names(SCAebv)=='EBV'] <- trait
            names(Mebv)[names(Mebv)=='EBV'] <- trait
            names(Febv)[names(Febv)=='EBV'] <- trait
            names(SCAebv)[names(SCAebv)=='se'] <- paste0(trait,'se')
            names(Mebv)[names(Mebv)=='se'] <- paste0(trait,'se')
            names(Febv)[names(Febv)=='se'] <- paste0(trait,'se')

            names(SCAebv)[names(SCAebv)=='acc'] <- paste0(trait,'acc')
            names(Mebv)[names(Mebv)=='acc'] <- paste0(trait,'acc')
            names(Febv)[names(Febv)=='acc'] <- paste0(trait,'acc')

        } else {
            names(asrout$SCAebv)[names(asrout$SCAebv)=='EBV'] <- trait
            names(asrout$Mebv)[names(asrout$Mebv)=='EBV'] <- trait
            names(asrout$Febv)[names(asrout$Febv)=='EBV'] <- trait

            names(asrout$SCAebv)[names(asrout$SCAebv)=='se'] <- paste0(trait,'se')
            names(asrout$Mebv)[names(asrout$Mebv)=='se'] <- paste0(trait,'se')
            names(asrout$Febv)[names(asrout$Febv)=='se'] <- paste0(trait,'se')

            names(asrout$SCAebv)[names(asrout$SCAebv)=='acc'] <- paste0(trait,'acc')
            names(asrout$Mebv)[names(asrout$Mebv)=='acc'] <- paste0(trait,'acc')
            names(asrout$Febv)[names(asrout$Febv)=='acc'] <- paste0(trait,'acc')

            SCAebv <- merge(SCAebv,asrout$SCAebv[,c(1,5,6,8)],by='FamNum')

            Mebv<-merge(Mebv,asrout$Mebv[,c(1:3,5)],by='ID')
            Febv<-merge(Febv,asrout$Febv[,c(1:3,5)],by='ID')
        }

        asrout$h2$trait <- trait
        asrout$d2$trait <- trait
        asrout$varcomp$trait <- trait
        #filename=paste0("Febv_",trait,".csv")
        #write.csv(asrout$Febv,filename)
        #filename=paste0("Mebv_",trait,".csv")
        #write.csv(asrout$Mebv,filename)

        h2 <- rbind(h2,asrout$h2)
        d2 <- rbind(d2,asrout$d2)
        varc <- rbind(varc,asrout$varcomp)
    }

    #dim(ebv)
    h2
    d2
    #head(febv)
    varc$item <- row.names(varc)
    varc <- varc[,c(6,7,1:5)]

    varc
    exout <- list(h2=h2,d2=d2,varcomp=varc,Mebv=Mebv,
                  Febv=Febv,SCAebv=SCAebv,popMean=popMean)

    return(exout)
}

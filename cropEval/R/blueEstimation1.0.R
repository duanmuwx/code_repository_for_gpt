
#' Adding Missing positions for spatial analysis
#'
#' @param w A data frame
#'
#' @return A data frame with missing value added
#' @export
#'
#' @examples
#' \donttest{
#'  df <- data.frame(read.table(textConnection(
#'    '     Range Pass Set Rep     Name  YLD
#'   1    1  24   2  ID1   901
#'   1    3  24   2  ID2   852
#'   1    4  24   2  ID3   753
#'   2    1  23   2  ID4   765
#'   2    2  23   2  ID5   964
#'   2    3  23   2  ID6   963
#'   2    4  23   2  ID7   737
#'   3    1  22   2  ID8   803
#'   3    2  22   2  ID9   963
#'   3    3  22   2 ID10   769
#'   3    4  22   2 ID11   852
#'   4    1  21   2 ID12   778
#'   4    2  21   2 ID13   1021
#'   4    3  21   2 ID14   667
#'   4    4  21   2 ID15   852'
#'   ),header=T,sep=''))
#'
#' print(table(df$Range,df$Pass))
#' dat <- missingAdding(w=df)
#' dat <- dat[order(dat$Range,dat$Pass),]
#' print(table(dat$Range,dat$Pass))
#' dat
#' }


missingAdding <- function(w){
  m <- table(w$Range,w$Pass) #create a contingency table "m" using the table() function, based on the variables "Range" and "Pass" from the data frame "w"
  m <- as.matrix(m) #The contingency table "m" is then converted into a matrix using as.matrix()
  colnames(m) #retrieve the column names of the matrix "m"
  row.names(m) #retrieve the row names
  nr=nrow(m)
  nc=ncol(m) #The number of rows and columns in the matrix "m" are stored in the variables "nr" and "nc", respectively
  ndiff <- nc*nr-sum(m) #calculate the difference between the product of "nc" and "nr" and the sum of all elements in the matrix "m"

  if (ndiff>0){ #If "ndiff" is greater than zero (which means there are missing values), the code proceeds with the following actions
    missing <- c() #create an empty vector
    for (i in 1:nr){
      for (j in 1:nc){ #A nested loop iterates through each element in the matrix "m"
        if (m[i,j]==0){
          rn <- rownames(m)[i]
          cn <- colnames(m)[j] #If an element is equal to zero, the corresponding row and column names are stored in the variables "rn" and "cn", respectively

          #The row and column names are then added to the "missing" vector using rbind()
          missing <- rbind(missing, c(as.numeric(as.character(rn)),

                                      as.numeric(as.character(cn))))
        }
      }
    }
    #A new data frame named "zeros" is created with NA values, and its column names are set to match the column names of the data frame "w"
    zeros <- data.frame(matrix(NA,nrow(missing),ncol(w)))
    names(zeros) <- names(w)
    for (i in 1:nrow(missing)){
      #The "Pass" and "Range" columns of the "zeros" data frame are filled with the values from the "missing" vector
      zeros$Pass[i] <- missing[i,2]
      zeros$Range[i] <- missing[i,1]
      #Additional columns ("Name", "Field", and "Location") are also populated with specific values from the data frame "w"
      zeros$Name[i] <- paste0('DUMMY',i)
      zeros$Field[i] <- w$Field[1]
      zeros$Location[i] <- w$Location[1]
    }
    #The "zeros" data frame is appended to the original data frame "w" using rbind(), and the resulting data frame is stored in the variable "ww"
    ww <- rbind(w,zeros)
    table(ww$Range,ww$Pass)
  } else {
    #If there are no missing values (i.e., "ndiff" is not greater than zero), the code assigns the original data frame "w" to the variable "ww"
    ww <- w
  }


  return (ww)

}

#' the ratio of the number of unique values occurring more than once to the total number of unique values in the "Name" column
#'
#' @param df a data.frame
#' @return a number
#' @examples
#' \donttest{
#'  df <- data.frame(read.table(textConnection(
#'    '   Name Range Pass Set Rep YLD
#'   1    1  24   2  ID1   901
#'   2    3  24   2  ID2   852
#'   2    4  24   2  ID3   753
#'   2    1  23   2  ID4   765
#'   2    2  23   2  ID5   964
#'   2    3  23   2  ID6   963
#'   2    4  23   2  ID7   737
#'   3    1  22   2  ID8   803
#'   3    2  22   2  ID9   963
#'   3    3  22   2 ID10   769
#'   3    4  22   2 ID11   852
#'   4    1  21   2 ID12   778
#'   4    2  21   2 ID13   1021
#'   4    3  21   2 ID14   667
#'   4    4  21   2 ID15   852'
#'   ),header=T,sep=''))
#'
#' print(table(df$Name))
#' dat <- replicationCheck(df)
#' dat
#' }
#' @export replicationCheck

replicationCheck <- function(df){
  #create a data frame "q" using the table() function. It counts the frequency of occurrences of each unique value in the "Name" column of the data frame "df"
  q <- data.frame(table(df$Name))
  #"q$Freq > 1" filters the rows in the data frame "q" where the frequency is greater than 1.
  #nrow(q[q$Freq > 1, ]) calculates the number of rows in the filtered data frame, representing the number of unique values in the "Name" column that occur more than once
  #nrow(q) calculates the total number of unique values in the "Name" column.
  #returns the ratio of the number of unique values occurring more than once to the total number of unique values in the "Name" column
  return (nrow(q[q$Freq>1,])/nrow(q))
}

#' BLUE value estimation
#'
#' @param datdf a data.frame
#' @param repct a number
#' @param trait trait name
#' @return A list
#' @examples
#' \donttest{
#' names(phdemo)
#' dat<-phdemo$SU
#' dat <- dat[!is.na(dat$Range) & !is.na(dat$Pass),]
#' dat <- dat[!is.na(dat$Name),]
#' dat <- dat[!is.na(dat$Location),]
#' dat <- dat[!is.na(dat$Field),]
#' dat <- dat[!is.na(dat$Type),]
#' dat<-dat[dat$Location==unique(dat$Location)[1],]
#' dat<-dat[dat$Field==unique(dat$Field)[1],]
#' datdf<-dat[dat$Type==unique(dat$Type)[1],]
#'
#' print(table(datdf$Range,datdf$Pass))
#'
#' trait<-"YLD"
#' names(datdf)[names(datdf)=="YLD"]<-"trait"
#'
#' repct<-replicationCheck(df=datdf)
#' repct
#' datdf<-missingAdding(w=datdf)
#' datdf <- datdf[order(datdf$Range,datdf$Pass),]
#' print(table(datdf$Range,datdf$Pass))
#' library(asreml)
#' blueout<-blueEstimation(datdf,repct,trait)
#' }
#' @export blueEstimation

blueEstimation <- function(datdf,repct,trait){
#checks if the column names of the data frame datdf contain 'Name', 'Range', 'Pass', 'Location', 'Field', 'Set', 'Type', 'trait'.
#If these column names are present, it converts the Range, Pass, Name, Set, and Location columns of datdf into factor variables.
      if (all(c('Name','Range','Pass','Location','Field','Set',
              'Type','trait') %in% names(datdf))){
        datdf$Range <- as.factor(datdf$Range)
        datdf$Pass <- as.factor(datdf$Pass)
        datdf$Name <- as.factor(datdf$Name)
        datdf$Set <- as.factor(datdf$Set)
        datdf$Location <- as.factor(datdf$Location)

          #if the number of rows in datdf that have non-missing values in the trait column is less than 130
          #create a linear mixed-effects model using the asreml() function. This model assumes that the trait variable is determined by the Range and Pass.
          #The random effects are specified by the idv(Name) term, which suggests that the variable Name represents individual-level random effects.
          #The residual variance is modeled using the idv(units) term. The data for model fitting is provided through the data = datdf argument.
          if (nrow(datdf[!is.na(datdf$trait),])<130){
          cat('asr <- asreml(trait ~ Range + Pass,random = ~idv(Name),
                              residual = ~idv(units), data = datdf)\n')
          asr <- asreml(trait ~ Range + Pass,random = ~idv(Name),
                        residual = ~idv(units), data = datdf)

          #if the number of rows in datdf that have non-missing values in the trait column is less than 150
          #create a linear mixed-effects model using the asreml() function. This model assumes that the trait variable is determined only by the intercept (constant term).
          #The random effects are specified by the idv(Name) term, which suggests that the variable Name represents individual-level random effects.
          #The residual variance is modeled using the autoregressive structure ar1(Range):ar1(Pass). The data for model fitting is provided through the data = datdf argument.
        } else if (nrow(datdf[!is.na(datdf$trait),])<150){
          cat('asr <- asreml(trait ~ 1,random = ~idv(Name),
                              residual = ~ar1(Range):ar1(Pass), data = datdf)\n')
          asr <- asreml(trait ~ 1,random = ~idv(Name),
                        residual = ~ar1(Range):ar1(Pass), data = datdf)

          #if the variable trait is equal to "TWT"
          #create a linear mixed-effects model using the asreml() function. This model assumes that the trait variable depends on the Pass variable.
          #The random effects are specified by the term idv(Name), indicating that the variable Name represents individual-level random effects.
          #The residual variance is modeled using the autoregressive structure ar1(Range):ar1(Pass). The data for model fitting is provided through the data = datdf argument.
        } else if(trait=='TWT'){
            cat('asr <- asreml(trait ~ Pass,random = ~idv(Name),
                          residual = ~ar1(Range):ar1(Pass), data = datdf)\n')
            asr <- asreml(trait ~ Pass,random = ~idv(Name),
                          residual = ~ar1(Range):ar1(Pass), data = datdf)

          #When all the above conditions are not met, execute the following model.
        } else {
            cat('asr <- asreml(trait ~ Range + Pass,random = ~idv(Name),
                          residual = ~ar1(Range):ar1(Pass), data = datdf)\n')
            asr <- asreml(trait ~ Range + Pass,random = ~idv(Name),
                          residual = ~ar1(Range):ar1(Pass), data = datdf)

        }

        #iteratively update the model until convergence is achieved
        asr <- update(asr)
        if(!asr$converge)asr <- update(asr)
        if(!asr$converge)asr <- update(asr)
        if(!asr$converge)asr <- update(asr)

        #Viewing the significance of fixed effects
        print(wald(asr))

        #Viewing variance components
        summary(asr)$varcomp

        #Calculating heritability and standard errors,H2 is a data.frame
        H2 <- vpredict(asr,H2~V1/(V1+V2))

        #assign the unique values from specific columns of the datdf dataset to properties (Location, Field, Type, and N) of the variable H2.
        H2$Location <- unique(datdf$Location)[1]
        H2$Field <- unique(datdf$Field)[1]
        H2$Type <- unique(datdf$Type)[1]
        H2$N <- nrow(datdf)
        H2


        #Predict means for individuals
        blue <- predict(asr, classify = "Name",pworkspace="800mb",
                        data=datdf)$pvals

          #Check if the number of rows in blue (excluding rows where predicted.value is not NA) is not equal to zero
          if(nrow(blue[!is.na(blue$predicted.value),]!=0)){
          #Removes rows from blue where the Name column is equal to '0
          blue <- blue[blue$Name!='0',]
          #Sets the Location, Field, Type, and N properties of blue using the unique values from the corresponding columns in the datdf dataset
          blue$Location <- unique(datdf$Location)[1]
          blue$Field <- unique(datdf$Field)[1]
          blue$Type <- unique(datdf$Type)[1]
          blue$N <- nrow(datdf)
          blue$converge <- asr$converge
        } else {
          #if the all predicted.value is NA
          #extract asr$coefficients$random as blue
          blue <- data.frame(asr$coefficients$random)
          #Modify the content in the blue data.frame
          names(blue)[1] <- 'predicted.value'
          newcol <- stringsplit(row.names(blue),'_')
          blue$Name <- newcol$X2
          blue <- blue[blue$Name!='0',]
          blue <- blue[,c(2,1)]
          blue$std.error <- 'NA'
          blue$status <- 'NA'
          blue$Location <- unique(datdf$Location)[1]
          blue$Field <- unique(datdf$Field)[1]
          blue$N <- nrow(datdf)
          blue$converge <- asr$converge

        }
        head(blue)
        #Draw a QQ chart for the predicted.value in blue to check whether it conforms to the Normal distribution
        title <- paste(trait,'BLUE:  Location =',
                       unique(datdf$Location)[1],
                       ' Field =',unique(datdf$Field)[1],
                       ' H2 =',round(H2[1,1],2),'N =',nrow(datdf),
                       '\n')
        qqnorm(blue$predicted.value,main=title)
        qqline(blue$predicted.value)
        hist(blue$predicted.value,main = title)
    } else {
        cat('Not all necessary columns are included in data frame\n')
        head(datdf)
    }
   #create a list named outlist and assigns the variables H2, blue, and asr to its corresponding elements.
    outlist <- list(H2=H2,blue=blue,asr=asr)

    return (outlist)

}

extracting <- function(asr,var,individual,term) {
  ebv <- summary(asr,coef=T)$coef.random
  ebv <- as.data.frame(ebv)
  bvalue <- as.data.frame(ebv[grepl(term,row.names(ebv)),])

  head(bvalue)
  dim(bvalue)
  #print(dim(bvalue))
  bvalue$id <- 0
  bvalue$id <- 0
  for (i in 1:dim(bvalue)[1]){
    #print(i)
    s <- strsplit(row.names(bvalue)[i],"[_]")[[1]]
    bvalue$id[i] <- s[2]
  }
  v <- (1-bvalue$'std.error'^2/var)
  bvalue$rel <- ifelse(v<=0,0,v)
  bvalue$acc <- sqrt(bvalue$rel)


  bvalue <- bvalue[bvalue$id %in% individual,]

  bvalue <- bvalue[,c(4,1,2,5,6)]
  names(bvalue) <- c("ID","EBV","se","rel","acc")
  #print(head(bvalue))

  print(dim(bvalue))
  return(bvalue)

}

GCA_SCA_Estimation <- function(bluehb,ainv,ped){

  bluehb <- merge(bluehb,ped,by='Name',all.x=T)

  bluehb$SIRE <- as.character(bluehb$SIRE)
  bluehb$DAM <- as.character(bluehb$DAM)
  bluehb$SIRE[is.na(bluehb$SIRE)] <- '0'
  bluehb$DAM[is.na(bluehb$DAM)] <- '0'
  bluehb$SIRE <- as.factor(bluehb$SIRE)
  bluehb$DAM <- as.factor(bluehb$DAM)


  bluehb$Fam <- paste0(bluehb$SIRE,'_',bluehb$DAM)

  bluehb$Fam[bluehb$Fam=='NA_NA'] <- '0_0'

  fam <- unique(bluehb[,c('Fam','SIRE','DAM')])
  fam$FamNum <- 1:nrow(fam)
  bluehb <- merge(bluehb,fam[,c('Fam','FamNum')],by='Fam')
  head(bluehb)
  #writecsv(bluehb,'temp.csv')
  bluehb$FamNum <- as.factor(bluehb$FamNum)
  bluehb$Location  <- as.factor(bluehb$Location )
  bluehb$Field  <- as.factor(bluehb$Field)
  bluehb$LF <- paste0(bluehb$Location,'_',bluehb$Field)
  bluehb$LF  <- as.factor(bluehb$LF)

  asr2 <- asreml(predicted.value ~ LF,
                 random = ~ vm(SIRE,ainv) + vm(DAM,ainv) + idv(FamNum),
                 residual = ~idv(units),
                 data = bluehb)
  fixed<-asr2$coefficients$fixed
  popmean <- fixed[row.names(fixed)=='(Intercept)',]

  summary(asr2)$varcomp
  vs <- summary(asr2)$varcomp$component[2]
  vd <- summary(asr2)$varcomp$component[3]
  vf <- summary(asr2)$varcomp$component[1]

  h2 <- vpredict(asr2,h2~(V2+V3)/(V1+V2+V3+V4))
  print(h2)

  d2 <- vpredict(asr2,h2~V1/(V1+V2+V3+V4))
  print(d2)
  loglik2 <- asr2$loglik
  asr3 <- asreml(predicted.value ~ LF,
                 random = ~ vm(SIRE,ainv) + vm(DAM,ainv),# + idv(FamNum),
                 residual = ~idv(units),
                 data = bluehb)
  loglik3 <- asr3$loglik

  d2$LRTvalue <- abs(loglik2-loglik3)*2

  bvalue <- summary(asr2,coef=T)$coef.random
  head(bvalue)
  #write.csv(bvalue,'temp.csv')
  ind <- unique(c(bluehb$SIRE,bluehb$DAM))
  ebv <- extracting(asr2,var=vs,individual = ind,term='vm(Name, ainv)*')
  head(ebv)
  febv <- extracting(asr2,var=vf,individual = fam$FamNum,term='FamNum*')
  febv <- merge(fam,febv,by.x='FamNum',by.y='ID')
  febv <- febv[order(-febv$EBV),]
  head(febv)
  #write.csv(febv,'SCA.csv')

  ind <- unique(c(bluehb$SIRE))
  sebv <- extracting(asr2,var=vs,individual = ind,term='vm(SIRE, ainv)*')
  sebv <- sebv[sebv$EBV!=0,]
  sebv$par='SIRE'
  #sebv <- sebv[1:315,]
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


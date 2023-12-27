library(ggplot2)

#' Identify and process outliers in the data for a single trait
#'
#' This function identifies and processes outliers in the data for a single trait.
#'
#' @param data The input data as a data.frame.
#' @param trait The name of the trait to be processed.
#' @param model The method to be used for identifying outliers, either '3SD' or 'boxplot'.
#'
#' @return A list containing processed data for the given trait.
#'
#'
#' @export
#'
#' @examples
#' w <- phdemo$SpML
#' mst.Outliers<-dataOutliers_Identifying(data=w,trait ="MST",model='3SD')
dataOutliers_Identifying <- function(data, trait, model) {
  #处理缺失值与排序
  ww <- missingAdding(data)
  ww <- ww[order(ww$Range, ww$Pass),]

  #提取需要处理的数据
  if (trait == 'YLD') {
    datw <- ww[, c('Name', 'Range', 'Pass', 'Location',
                   'Field', 'Type', 'Set', trait, 'MST')]
  } else {
    datw <- ww[, c('Name', 'Range', 'Pass', 'Location',
                   'Field', 'Type', 'Set', trait)]
  }
  names(datw)[8] <- 'trait'

  # 3SD方法处理异常值
  if (model == '3SD') {
    std <- sd(datw$trait, na.rm = TRUE)
    mn <- mean(datw$trait, na.rm = TRUE)
    datw$abnormal <- ifelse(datw$trait < mn - std * 3 | datw$trait > mn + std * 3, 2, 1)
  }

  #boxplot方法处理异常值
  if (model == 'boxplot') {
    bx <- boxplot(datw$YLD)
    min <- bx$stats[1, 1]
    max <- bx$stats[5, 1]
    c(min, max)
    datw$abnormal <- ifelse(datw$trait < min | datw$trait > max, 2, 1)
  }

  #其他常识性设定范围处理异常值
  if (trait == 'MST') datw$abnormal[datw$trait >= 43] <- 2
  if (trait == 'YLD') datw$abnormal[datw$trait < 200] <- 2
  if (trait == 'YLD') datw$abnormal[datw$MST >= 43] <- 2
  if (trait == 'LDG') datw$abnormal <- 1
  if (trait == 'GLS') datw$abnormal <- 1
  if (trait == 'NCLB') datw$abnormal <- 1
  if (trait == 'ER') datw$abnormal <- 1
  if (trait == 'PMD') datw$abnormal <- 1

  #将标记为异常值的部分处理为NA
  outliers <- datw[datw$abnormal == 2,]
  raw <- datw
  datw$trait[datw$abnormal == 2] <- NA

  #将处理后数据的被处理列名称改回原名称
  names(datw)[8] <- trait

  # 将结果存储到列表中
  outlist <- list(datclean = datw, outliers = outliers, rawdata = raw)
  return(outlist)
}


#' Identify and process outliers in the data for multiple traits
#'
#' This function identifies and processes outliers in the data for multiple traits.
#'
#' @param data The input data as a data.frame.
#' @param traits A character vector containing the names of the traits to be processed.
#' @param model The method to be used for identifying outliers, either '3SD' or 'boxplot'.
#'
#' @return A list containing processed data for each trait separately.
#'
#' @export
#'
#' @examples
#' w <- phdemo$SpML
#' yld <- dataOutliers.for(data=w,traits =c('MST',"YLD","PHT","EHT"),model='3SD')
#'
dataOutliers.for <- function(data, traits, model) {
  # 创建一个空列表用于存储结果
  result_list <- list()
  # 获取原始数据的ID
  ID <- paste0(data$Name, "_", data$Range, "_", data$Pass, "_", data$Location,
               "_", data$Field, "_", data$Type, "_", data$Set, "_")
  data$ID <- ID
  result_df <- cbind(data[, 1:16], ID = data$ID)

  for (trait in traits) {
    # 调用 dataOutliers_Identifying 函数
    result <- dataOutliers_Identifying(data, trait, model)

    # 获取MST的已处理数据的ID
    ID2 <- paste0(result$datclean$Name, "_", result$datclean$Range,
                  "_", result$datclean$Pass, "_", result$datclean$Location,
                  "_", result$datclean$Field, "_", result$datclean$Type,
                  "_", result$datclean$Set, "_")
    result$datclean$ID <- ID2

    # 提取 trait 列，并且设置其正确的名称
    new_trait_name <- names(result$datclean)[8]
    mst <- data.frame(result$datclean[, 8], ID = result$datclean$ID)
    names(mst)[1] <- new_trait_name

    # 合并到 result_df 中
    result_df <- merge(result_df, mst, by = "ID", all = TRUE)

    # 绘制处理前的QQ图
    par(mfrow = c(1, 2))
    qqnorm(result$rawdata[, 8], main = paste("Rawdata QQ -", trait), col = result$rawdata$abnormal)
    qqline(result$rawdata[, 8])

    # 绘制处理后的QQ图
    qqnorm(result$datclean[, 8], main = paste("Dataclean QQ -", trait), col = result$datclean$abnormal)
    qqline(result$datclean[, 8])

    # 将处理后的 result_list 保存到列表中
    result_list[[trait]] <- result
  }

  # 从result_df中删除ID列
  result_df <- result_df[, -grep("^ID$", names(result_df))]

  # 将所有结果进行返回，包括result_list
  return(list(result_list = result_list, result_df = result_df))
}

#用包中的SpML数据进行试运行试运行
# w <- phdemo$SpML
# yld <- dataOutliers.for(data=w,traits =c('MST',"YLD","PHT","EHT"),model='3SD')





















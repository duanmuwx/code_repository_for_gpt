
library(openxlsx)
library(readxl)
library(writexl)
library(biosim)
library(tidyverse)
library(reshape2)
library(stringi) 
library(asreml)

options("openxlsx.borderColour" = "black")
options("openxlsx.borderStyle" = "thin")
options("openxlsx.dateFormat" = "mm/dd/yyyy")
options("openxlsx.datetimeFormat" = "yyyy-mm-dd hh:mm:ss")
options("openxlsx.numFmt" = NULL)

## Change the default border colour to #4F81BD
options("openxlsx.borderColour" = "#4F81BD")

d <- data.frame(read_excel('Output/2023年-P1试验-杂交种产量及病害分析-20231103.xlsx',
                           sheet = 'EMSP'))

wb <- createWorkbook()
addWorksheet(wb, "EMSP")

hs1 <- createStyle(
    fontColour = "black", fgFill = "azure",
    halign = "center", valign = "center", textDecoration = "bold",
    border = "TopBottomLeftRight"
)

writeData(wb, "EMSP", d[,1:5],
          rowNames = FALSE, 
          startCol = 1, startRow = 1,
          borders = "surrounding", borderColour = "black",
          headerStyle = hs1, borderStyle = "double"
) ## black border

setColWidths(wb, "EMSP", cols = 1:5, widths = 5)
for (i in 1:18){
    istart <- (i-1)*4+6
    iend <- (i-1)*4+9
    writeData(wb, "EMSP", d[,istart:iend],
              rowNames = FALSE, 
              startCol = istart, startRow = 1,
              borders = "surrounding", borderColour = "black",
              headerStyle = hs1, borderStyle = "double"
    ) ## black border
    
    setColWidths(wb, "EMSP", cols = istart:iend, widths = 5)
    
}

writeDataTable(wb, 2, df, startRow = 23, startCol = 2, tableStyle = "TableStyleMedium21")
setColWidths(wb, "EMSP", cols = 1:5, widths = 5)
# for (i in 1:18){
#     istart <- (i-1)*4+6
#     iend <- (i-1)*4+9
#     writeDataTable(wb, "EMSP", d[,istart:iend],
#               rowNames = FALSE, 
#               startCol = istart, startRow = 1,
#               borders = "surrounding", borderColour = "black",
#               headerStyle = hs1, borderStyle = "double",
#               tableStyle = "TableStyleMedium21"
#     ) ## black border
#     
#     setColWidths(wb, "EMSP", cols = istart:iend, widths = 5)
#     
# }



saveWorkbook(wb, "temp2.xlsx", overwrite = TRUE)

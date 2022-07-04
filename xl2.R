#!/usr/bin/env Rscript
# from
# http://www.sthda.com/english/wiki/r-xlsx-package-a-quick-start-guide-to-manipulate-excel-files-in-r
library(ggplot2)
library(Cairo)
library(xlsx)

# create a new workbook for outputs
# possible values for type are : "xls" and "xlsx"
wb <- createWorkbook(type="xlsx")
sheet <- createSheet(wb, sheetName = "US State Facts")

# Styles for the data table row/column names
TABLE_ROWNAMES_STYLE <- CellStyle(wb) + Font(wb, isBold=TRUE)
TABLE_COLNAMES_STYLE <- CellStyle(wb) + Font(wb, isBold=TRUE) +
  Alignment(wrapText=TRUE, horizontal="ALIGN_CENTER") +
  Border(color="black", position="BOTTOM", pen="BORDER_THICK") 

# Create a new sheet in the workbook

# Add a table into a worksheet
#++++++++++++++++++++++++++++++++++++
addDataFrame(state.x77, sheet, startRow=1, startColumn=1, 
             colnamesStyle = TABLE_COLNAMES_STYLE,
             rownamesStyle = TABLE_ROWNAMES_STYLE)

# Change column width
setColumnWidth(sheet, colIndex=c(1:ncol(state.x77)), colWidth=11)
cb <- CellBlock(sheet, 10, 1, 1, 5, create=F) # the create=F is vital if you don't want value overwriting.
# cf <- Font(wb, color="#FF0000", heightInPoints=12)
cf <- Font(wb, color="#FF0000", name="Liberation Sans", heightInPoints=10) # these are my system defaults, but they had to be first detected and then set here.
CB.setFont(cb, cf, 1, 1:5)
cf <- Fill(backgroundColor="#AA5522")
CB.setFill(cb, cf, 1, 2:5)

# Save it all.
saveWorkbook(wb, "r-xlsx-report-example.xlsx")

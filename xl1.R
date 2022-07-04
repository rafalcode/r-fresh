#!/usr/bin/env Rscript
# from
# http://www.sthda.com/english/wiki/r-xlsx-package-a-quick-start-guide-to-manipulate-excel-files-in-r
library(ggplot2)
library(Cairo)
library(xlsx)

# create a new workbook for outputs
# possible values for type are : "xls" and "xlsx"
wb <- createWorkbook(type="xlsx")

# Define some cell styles
# Title and sub title styles
TITLE_STYLE <- CellStyle(wb)+ Font(wb,  heightInPoints=16, color="blue", isBold=TRUE, underline=1)
SUB_TITLE_STYLE <- CellStyle(wb) + Font(wb,  heightInPoints=14, isItalic=TRUE, isBold=FALSE)

# Styles for the data table row/column names
# COLROW <- CellStyle(wb) + Fill(foregroundColor="blue")
TABLE_ROWNAMES_STYLE <- CellStyle(wb) + Font(wb, isBold=TRUE)
TABLE_COLNAMES_STYLE <- CellStyle(wb) + Font(wb, isBold=TRUE) +
  Alignment(wrapText=TRUE, horizontal="ALIGN_CENTER") +
  Border(color="black", position=c("TOP", "BOTTOM"), pen=c("BORDER_THIN", "BORDER_THICK")) 

# Create a new sheet in the workbook
sheet <- createSheet(wb, sheetName = "US State Facts")

#++++++++++++++++++++++++
# Helper function to add titles
#++++++++++++++++++++++++
# - sheet : sheet object to contain the title
# - rowIndex : numeric value indicating the row to 
  #contain the title
# - title : the text to use as title
# - titleStyle : style object to use for title
xlsx.addTitle<-function(sheet, rowIndex, title, titleStyle){
  rows <-createRow(sheet,rowIndex=rowIndex)
  sheetTitle <-createCell(rows, colIndex=1)
  setCellValue(sheetTitle[[1,1]], title)
  setCellStyle(sheetTitle[[1,1]], titleStyle)
}

# Add title and sub title into a worksheet
#++++++++++++++++++++++++++++++++++++
# Add title
xlsx.addTitle(sheet, rowIndex=1, title="US State Facts", titleStyle = TITLE_STYLE)
# Add sub title
xlsx.addTitle(sheet, rowIndex=2, title="Data sets related to the 50 states of USA.", titleStyle = SUB_TITLE_STYLE)

# Add a table into a worksheet
#++++++++++++++++++++++++++++++++++++
addDataFrame(state.x77, sheet, startRow=3, startColumn=1, 
             colnamesStyle = TABLE_COLNAMES_STYLE,
             rownamesStyle = TABLE_ROWNAMES_STYLE)

# Change column width
setColumnWidth(sheet, colIndex=c(1:ncol(state.x77)), colWidth=11)
cb <- CellBlock(sheet, 10, 1, 1, 8, create=F) # the create=F is vital if you don't want value overwriting.
# cf <- Font(wb, color="#FF0000", heightInPoints=12)
cf <- Font(wb, color="#FF0000", name="Liberation Sans", heightInPoints=10) # these are my system defaults, but they had to be first detected and then set here.
CB.setFont(cb, cf, 1, 1:8)
cf <- Fill(backgroundColor="#AA5522")
CB.setFill(cb, cf, 1, 2:8)

# Save it all.
saveWorkbook(wb, "r-xlsx-report-example.xlsx")

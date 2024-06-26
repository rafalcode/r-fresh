#!/usr/bin/env Rscript
# https://rpubs.com/kylewbrown/r-colors
# actually this proved difficult to work out
# it's old type R graphics.
library(Cairo)

# 1.Define R Color Data ----
# RGB codes
color.rgb <- t(col2rgb(colors()))
# Hexadecimal codes
color.hex <- rgb(color.rgb[,1], color.rgb[,2], color.rgb[,3], maxColorValue = 255)
# Text highlighting
color.text <- ifelse(apply(color.rgb, 1, mean) > 127, "black", "white")
# Consolidate
color.df <- data.frame(name = colors(),
                       red = color.rgb[, "red"],
                       green = color.rgb[, "green"],
                       blue = color.rgb[, "blue"],
                       hex = color.hex,
                       text = color.text)

# 2.Plot R Colors By Name ----
# configure graphical device
CairoPNG("colr0.png", 1000, 1000)
n.col <- 11
n.row <- 60
par(pin = c(11.692, 6.267), mai=c(0, 0, 0, 0))
# create plot
plot(c(0, n.col), c(0, n.row), 
     type = "n", 
     bty = "n", 
     ylab = "", 
     xlab = "", 
     axes = FALSE)

jj <- 0
for(i in 1:n.col){
  color.count <- (i-1) * n.row
  color.mod <- length(colors()) - color.count
  y.val <- ifelse(color.mod < n.row, n.row - color.mod + 1, 1)
  color.names <- as(color.df[color.count + 1:n.row, "name"], "character")
  lab.color.names <- paste0(color.names, "(", jj,")")
  rect(i - 1, y.val - 0.5, i, n.row:y.val + 0.5, border = "black", col = color.names)
  text.color <- as(color.df[color.count + 1:n.row, "text"], "character")
  text(i-0.5, n.row:y.val, labels = lab.color.names, cex = 0.75, col = text.color)
}

# Cairo image template
# put plot command here
dev.off()

CairoPNG("colr01.png", 1000, 1000)
par(pin = c(11.692, 6.267), mai=c(0, 0, 0, 0))

# create plot
plot(c(0, n.col), c(0, n.row), 
     type = "n", 
     bty = "n", 
     ylab = "", 
     xlab = "", 
     axes = FALSE)

for(i in 1:n.col){
  color.count <- (i-1) * n.row
  color.mod <- length(colors()) - color.count
  y.val <- ifelse(color.mod < n.row, n.row-color.mod + 1, 1)
  color.names <- as(color.df[color.count + 1:n.row, "hex"], "character")
  rect(i - 1, y.val - 0.5, i, n.row:y.val + 0.5, border = "black", col = color.names)
  text.color <- as(color.df[color.count + 1:n.row, "text"], "character")
  text(i-0.5, n.row:y.val, labels = color.names, cex = 0.75, col = text.color)
}
dev.off()

# from rbloggers
# ref. https://www.r-bloggers.com/2024/02/demystifying-the-melt-function-in-r/
# protoype of the metl() function
# as you know!) this conevrt wide form (i.e. usualy table form
# to long form which usually 3 but can be two, columns variable (the colname repeated for nrows)
# and the value column
library(reshape)

scores <- data.frame(
  student = c("Alice", "Bob", "Charlie"),
  math = c(90, 80, 85),
  english = c(85, 90, 80))

lfscores <- melt(scores, id.vars = "student", measure.vars = c("math", "english"))

# You can source this script and compare scores and lfscroes.

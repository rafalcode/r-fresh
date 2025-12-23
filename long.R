# the why of long formats
# all because of Hadley, it's his concept fo tidy data which (i you read most of it) is touched upon
# in:
# https://stackoverflow.com/questions/64390747/ggplot-why-do-i-have-to-transform-the-data-into-the-long-format
library(tidyverse)

# Data in wide format
  df_wide <- data.frame(
   Horizons = seq(1,10,1),
   Country1 = c(2.5, 2.3, 2.2, 2.2, 2.1, 2.0, 1.7, 1.8, 1.7, 1.6),
   Country2 = c(3.5, 3.3, 3.2, 3.2, 3.1, 3.0, 3.7, 3.8, 3.7, 3.6),
   Country3 = c(1.5, 1.3, 1.2, 1.2, 1.1, 1.0, 0.7, 0.8, 0.7, 0.6)
   )

# Convert to long format
  df_long <- df_wide %>%
   gather(key = "variable", value = "value", -Horizons)

# Plot the lines
#  plotstov <- ggplot(df_long, aes(x = Horizons, y = value)) +
#   geom_line(aes(colour = variable, group = variable))+
#   theme_bw()


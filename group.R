#!/usr/bin/env Rscript
# this script does what?  library(dplyr)

# You have a df a column representing groups
# it could also easily be repeats of the principal label.
# i.e. genes that are diff methyl ... ther will be several CpGs corresponding to one gene
# tso tht gene will be repeated.

# We can use "aggegrate()" from the stats package (quite native of course)
# but as usaly dplyr has a method too.

# this example from https://statisticsglobe.com/mean-by-group-in-r

#  the famous iris
# has no obsevration names, but you can use the indexed rownames for that. 150 obs.
# now, the df gets a "Speices" which groups each foot eh 150.
# so it's ideal to try group on.
# in these two cases we concentrate on 


ag <- aggregate(x = iris$Sepal.Length,                # Specify data column
          by = list(iris$Species),              # Specify group indicator
          FUN = mean)                           # Specify function (i.e. mean)
 
#    Group.1     x
#     setosa 5.006
# versicolor 5.936
#  virginica 6.588



ag2 <- iris %>%                                        # Specify data frame
  group_by(Species) %>%                         # Specify group indicator
  summarise_at(vars(Sepal.Length),              # Specify column
               list(name = mean))               # Specify function
 
# A tibble: 3 x 2
# Species    Sepal.Length
#              
# setosa             5.01
# versicolor         5.94
# virginica          6.59o

# NOTE: this is a sort of old invocation of dplyr

# HW says
# ‘vars()’ is superseded because it is only needed for the scoped verbs (i.e. ‘mutate_at()’, ‘summarise_at()’, and friends), whic have been been superseded in favour of ‘across()’. See                                                             │ 34   group_by(Species) %>%                         # Specify group indicator ‘vignette("colwise")’ for details.


# dplyr manual


tb <- tibble(g = c(1, 1, 2, 3), v1 = 10:13, v2 = 20:23) # soo like a dataframe
# gdf <- tibble(g = c(1, 1, 2, 3), v1 = 10:13, v2 = 20:23) %>% group_by(g)
gdf <- tb %>% group_by(g)
# note how group_by does not actaully transofrm the table, but the tittle is now aware that one of the columns is a grouper.


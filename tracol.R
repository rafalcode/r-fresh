# grDevices::color()
# will give you the 600 or so colors tht have names in R
# many of them are similar to the html colors but there's quite a few exceptions too.

# in any case you won't get opacity with the names though
# you have to convert via col2rgb()

# here is a ref and function
# ref. https://www.dataanalytics.org.uk/make-transparent-colors-in-r/
## Transparent colors
## Mark Gardener 2015
## www.dataanalytics.org.uk

# my version
trcol <- function(color, alpha = .5)
{
    rgb.val <- col2rgb(color)

    ## Make new color using input color as base and alpha set by transparency
    tracol <- rgb(rgb.val[1], rgb.val[2], rgb.val[3], max = 255, alpha = 255*alpha)

    ## Save the color
    return(tracol)
}

t_col <- function(color, percent = 50, name = NULL)
{
      #      color = color name
      #    percent = % transparency
      #       name = an optional name for the color

    ## Get RGB values for named color
    rgb.val <- col2rgb(color)

    ## Make new color using input color as base and alpha set by transparency
    t.col <- rgb(rgb.val[1], rgb.val[2], rgb.val[3],
                              max = 255,
                                           alpha = (100 - percent) * 255 / 100,
                                           names = name)

    ## Save the color
    invisible(t.col)
}

tc <- trcol("seagreen", 0.8)

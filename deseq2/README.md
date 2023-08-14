
# DESeq styles

We can look at the simulation function to gain insights into DESeq and its terminology.

  makeExampleDESeqDataSet
                             Make a simulated DESeqDataSet
Description
    Constructs a simulated dataset of Negative Binomial data from two conditions. By default, there
    are no fold changes between the two conditions, but this can be adjusted with the betaSD argument.
Usage
    makeExampleDESeqDataSet(
       n = 1000,
       m = 12,
       betaSD = 0,
       interceptMean = 4,
     interceptSD = 2,
     dispMeanRel = function(x) 4/x + 0.1,
     sizeFactors = rep(1, m)
   )
Arguments
   n              number of rows
   m              number of columns
   betaSD         the standard deviation for non-intercept betas, i.e. beta ~ N(0,betaSD)
   interceptMean  the mean of the intercept betas (log2 scale)
   interceptSD    the standard deviation of the intercept betas (log2 scale)
   dispMeanRel    a function specifying the relationship of the dispersions on 2^trueIntercept
   sizeFactors    multiplicative factors for each sample

What can we tell? 
Well, sizefactors are a per sample thing, so probably the full meaning is library sizes or the number of reads a sample wa sable to provoke.

# xampds.R
means Example Dataset, here I daable with the log2FC value.

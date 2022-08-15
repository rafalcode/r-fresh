# Rf's version of Sebastian Sauer's p-value simulation
# note it's not an artcile Rf likes very much but I'm curious as to ho he's doing it.
library(ggplot2)
library(dplyr)
library(Cairo)

# Original R code by Norman Markgraf and Sebastian Sauer

stipros <- function(n_stipros = 1000, mean = 100, sd = 15, n = 30, distribution = "normal")
{
    # returns histogram of simulated samples
    # arguments:
    # n_stipros: number of samples to be simulated
    # mean: mean of distribution
    # sd: sd of distribtuion
    # distribution: type of distribution. Either "normal" (default) or "uniform"

    stopifnot(distribution %in% c("normal", "uniform"))

    if (distribution == "normal") {
        result <- replicate(n_stipros, mean(rnorm(n = n, mean = mean, sd = sd)))
        hist(result, main = "Histogramm zu Stichproben \naus einer Normalverteilung")
    }
    if (distribution == "uniform") {
        result <- replicate(n_stipros, mean(runif(n = n, min=0, max=1)))
        hist(result, main = "Histogramm zu Stichproben \naus einer Gleichverteilung")
    }
}  # end stipros

simu_p <- function(n_samples = 1000, mean = 100, sd = 15, n = 30, distribution = "normal", sample_mean = 107)
{
    # returns: histogram of simulated samples including shaded area under the curve according to p-value as per simulation
    # arguments:
    # n_stipros: number of samples to be simulated
    # mean: mean of distribution
    # sd: sd of distribution
    # n: sample size
    # distribution: type of distribution. Either "normal" (default) or "uniform"
    # sample_mean: mean of sample

    stopifnot(distribution %in% c("normal", "uniform"))

    if (distribution == "normal") {
        # Watch this, the means is calcuated for a small (30) sample from the normal,
        # and it's done n_samples amount of times.
        samples <- replicate(n_samples, mean(rnorm(n = n, mean = mean, sd = sd)))
        # so really it should be called a set of samplemeans.
        # note as we're given a sample_mean, it's possible to perform a t-test
        #although he doesn't do that here.
        df <- data.frame(samples = samples)

        df %>% mutate(perc_num = percent_rank(samples),
                      max_5perc = ifelse(perc_num > (trunc(.95*n)/n), 1, 0),
                      greater_than_sample = ifelse(samples > sample_mean, 1, 0)) -> df2
        # so it's recorded wich are over a certain mean
        p_value <- round(mean(df2$greater_than_sample), 3)
        browser()

        CairoPNG("clt0.png", 800, 800)
        ggplot(df2) + aes(x = samples) + geom_histogram() +
            labs(title = paste("Histogram of ", n_samples, " samples", "\n from a normal distribution", sep = ""),
                 caption = paste("mean-pop=", mean, ", sd-pop=",sd, sep = "", ", mean in sample=", sample_mean), x = "sample means") +
        geom_histogram(data = filter(df2, perc_num > .95), fill = "red") +
        theme(plot.title = element_text(hjust = 0.5)) +
        geom_vline(xintercept = sample_mean, linetype = "dashed", color = "grey40") +
        annotate("text", x = Inf, y = Inf, label = paste("p =", p_value), hjust = 1, vjust = 1)
        dev.off()
    }

    if (distribution == "uniform") {
        samples <- replicate(n_samples, mean(runif(n = n, min=0, max=1)))
        df <- data.frame(samples = samples)

        if(sample_mean > 1)
            sample_mean <- .99

        df %>% mutate(perc_num = percent_rank(samples),
                      max_5perc = ifelse(perc_num > (trunc(.95*n)/n), 1, 0),
                      greater_than_sample = ifelse(samples > sample_mean, 1, 0)) -> df2

        p_value <- round(mean(df2$greater_than_sample), 3)

        CairoPNG("clt0.png", 800, 800)
        ggplot(df2) + aes(x = samples) + geom_histogram() +
            labs(title = paste("Histogram of ", n_samples, " samples", "\n from a uniform distribution", sep = ""),
                 caption = paste("sample mean =", sample_mean, ", min-pop = 0, max-pop = 1"), x = "sample means") +
        geom_histogram(data = filter(df2, perc_num > .95), fill = "red") +
        theme(plot.title = element_text(hjust = 0.5)) +
        geom_vline(xintercept = sample_mean, linetype = "dashed", color = "grey40") +
        annotate("text", x = Inf, y = Inf, label = paste("p =",p_value), hjust = 1, vjust = 1)
        dev.off()
    }
}

simu_p()

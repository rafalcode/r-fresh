library(ggplot2)
library(Cairo)

# df <- iris[1:4]
df <- iris[c(10:12, 60:62, 110:112),]
df$idx <- rownames(df)
rownames(df) <- 1:nrow(df) # restart row indexing.

pca_res <- prcomp(df[1:4], scale. = TRUE)
percentVar <- round(100*pca_res$sdev^2 / sum( pca_res$sdev^2 ), 1)
df_proj <- as.data.frame(pca_res$x)[1:2]
colnames(df_proj) <- paste0(colnames(df_proj), " (", percentVar[1:2], "%)")
df_proj$species <- df$Species # this can be done, because thr prcomp retains the df2 order - not apparent but I checked and yes
rownames(df_proj) <- 1:nrow(df_proj) # not really necessary, but it reminds us how brushPoints sees the reactive df.
df_proj$idx <- df$idx

CairoPNG("ggpca.png", 800, 800)
# ggplot(df_proj, aes(x = .data[[1]], y = .data[[2]], label=idx)) +
# ggplot(df_proj, aes(x = df_proj[[1]], y = df_proj[[2]], label=idx)) +
# ggplot(df_proj, aes_string(x = names(df_proj)[1], y = names(df_proj)[2], label=names(df_proj)[4])) +
ggplot(df_proj, aes_string(x = get(names(df_proj)[1]), y = get(names(df_proj)[2]))) +
        geom_point(aes(col=df_proj$species), size=3.5) +
        xlim(-3, 3) +
        ylim(-2.5, 2.5)
dev.off()

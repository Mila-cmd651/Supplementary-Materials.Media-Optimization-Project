setwd("~/ETHZ/Autumn 2025/Research project/Heatmap")
heatmap <- read.table("heatmap.txt", sep="\t", dec=",", header=TRUE, na.strings = c(""))
heatmap$initial.pH <- NULL
heatmap$Supplement.mix <- NULL
heatmap$Donor <- NULL

# https://www.youtube.com/watch?v=Z0yOHBK4-2o

library(ggplot2)
library(reshape2)
library(dplyr)

cor_matrix2 <- data.frame(cor(heatmap[,c(6,7,8,9)], method="kendall"))
cor_matrix2

## Find the p-values for each correlation

p_matrix2 <- matrix(NA,
                   nrow = nrow(cor_matrix2),
                   ncol = ncol(cor_matrix2),
                   dimnames = dimnames(cor_matrix2))

for (i in rownames(cor_matrix2)) {
  for (j in colnames(cor_matrix2)) {
    test <- cor.test(heatmap[[i]], heatmap[[j]], method="kendall")
    p_matrix2[i, j] <- test$p.value
  }
}

p_to_stars <- function(p) {
  if (p <= 0.0001) {
    return("****")
  } else if (p <= 0.001) {
    return("***")
  } else if (p <= 0.01) {
    return("**")
  } else if (p <= 0.05) {
    return("*")
  } else if (p > 0.05) {
    return("")
  }
}

# Convert our correlation matrix into a long format that can be used by ggplot2

cor_long2 <- melt(as.matrix(cor_matrix2))
p_long2 <- melt(p_matrix2, varnames = c("Var1", "Var2"), value.name = "pvalue")
p_long2$Significance <- sapply(p_long2$pvalue, p_to_stars)

cor_final2 <- merge(cor_long2, p_long2, by = c("Var1", "Var2"))
cor_final2


# Plotting the heatmap
ggplot(data = cor_final2, aes(x = Var1, y = Var2, fill = value)) +
  geom_tile() +
  ylab("") +
  xlab("") +
  geom_text(aes(label = paste0(round(value, 2), Significance)), color = "black", size = 3) +
  scale_fill_gradient2(low = "darkred", high = "navy", mid = "white", midpoint = 0,
                       limit = c(-1,1), space = "Lab", name = "Correlation") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))


## Do the exact same thing but only with 5.8
heatmap <- read.table("heatmap.txt", sep="\t", dec=",", header=TRUE, na.strings = c(""))
heatmap_filt2 <- heatmap[heatmap$initial.pH == "5.8", ]
heatmap_filt2$initial.pH <- NULL
heatmap_filt2$Supplement.mix <- NULL
heatmap_filt2$Donor <- NULL

cor_matrix_filt2 <- data.frame(cor(heatmap_filt2[,c(6,7,8,9)], method="kendall"))


## Find the p-values for each correlation

p_matrix_filt2 <- matrix(NA,
                        nrow = nrow(cor_matrix_filt2),
                        ncol = ncol(cor_matrix_filt2),
                        dimnames = dimnames(cor_matrix_filt2))

for (i in rownames(cor_matrix_filt2)) {
  for (j in colnames(cor_matrix_filt2)) {
    test <- cor.test(heatmap_filt2[[i]], heatmap_filt2[[j]], method="kendall")
    p_matrix_filt2[i, j] <- test$p.value
  }
}


# Convert our correlation matrix into a long format that can be used by ggplot2

cor_long_filt2 <- melt(as.matrix(cor_matrix_filt2))
p_long_filt2 <- melt(p_matrix_filt2, varnames = c("Var1", "Var2"), value.name = "pvalue")
p_long_filt2$Significance <- sapply(p_long_filt2$pvalue, p_to_stars)

cor_final_filt2 <- merge(cor_long_filt2, p_long_filt2, by = c("Var1", "Var2"))
cor_final_filt2

cor_triangle <- cor_final_filt2 %>%
  dplyr::filter(as.numeric(Var1) <= as.numeric(Var2))


# Plotting the heatmap
ggplot(data = cor_triangle, aes(x = Var1, y = Var2, fill = value)) +
  geom_tile() +
  ylab("") +
  xlab("") +
  geom_text(aes(label = paste0(round(value, 2), Significance)), color = "black", size = 3.5) +
  scale_fill_gradient2(low = "darkred", high = "navy", mid = "white", midpoint = 0,
                       limit = c(-1,1), space = "Lab", name = "Correlation") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))

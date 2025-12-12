setwd("~/ETHZ/Autumn 2025/Research project/Heatmap")
heatmap <- read.table("heatmap.txt", sep="\t", dec=",", header=TRUE, na.strings = c(""))
heatmap$initial.pH <- NULL
heatmap$Supplement.mix <- NULL
heatmap$Donor <- NULL

# https://www.youtube.com/watch?v=Z0yOHBK4-2o

library(ggplot2)
library(reshape2)
sessionInfo()

cor_matrix <- data.frame(cor(heatmap, method="kendall"))

## Convert correlation matrix to a have the information in the shape I want
cor_matrix$Mucin <- NULL
cor_matrix$GOS <- NULL
cor_matrix$X6c <- NULL
cor_matrix$HMOs <- NULL
cor_matrix$Fuc.Lac <- NULL
cor_matrix$Mean.Corrected.OD <- NULL
cor_matrix$DNA.conc. <- NULL
cor_matrix <- cor_matrix[-c(6, 7, 8, 9), ]

## Find the p-values for each correlation

p_matrix <- matrix(NA,
                   nrow = nrow(cor_matrix),
                   ncol = ncol(cor_matrix),
                   dimnames = dimnames(cor_matrix))

for (i in rownames(cor_matrix)) {
  for (j in colnames(cor_matrix)) {
    test <- cor.test(heatmap[[i]], heatmap[[j]], method="kendall")
    p_matrix[i, j] <- test$p.value
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

cor_long <- melt(as.matrix(cor_matrix))
p_long <- melt(p_matrix, varnames = c("Var1", "Var2"), value.name = "pvalue")
p_long$Significance <- sapply(p_long$pvalue, p_to_stars)

cor_final <- merge(cor_long, p_long, by = c("Var1", "Var2"))
cor_final


# Plotting the heatmap
ggplot(data = cor_final, aes(x = Var1, y = Var2, fill = value)) +
  geom_tile() +
  ylab("") +
  xlab("") +
  geom_text(aes(label = paste0(round(value, 2), Significance)), color = "black", size = 3) +
  scale_fill_gradient2(low = "darkred", high = "navy", mid = "white", midpoint = 0,
                      limit = c(-1,1), space = "Lab", name = "Correlation") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))


## Do the exact same thing but only with 5.8, only targets + OD and DNA concentration
heatmap <- read.table("heatmap.txt", sep="\t", dec=",", header=TRUE, na.strings = c(""))
heatmap_filt <- heatmap[heatmap$initial.pH == "5.8", ]
heatmap_filt$initial.pH <- NULL
heatmap_filt$Supplement.mix <- NULL
heatmap_filt$Donor <- NULL

cor_matrix_filt <- data.frame(cor(heatmap_filt, method="kendall"))

## Convert correlation matrix to a have the information in the shape I want
cor_matrix_filt$Mucin <- NULL
cor_matrix_filt$GOS <- NULL
cor_matrix_filt$X6c <- NULL
cor_matrix_filt$HMOs <- NULL
cor_matrix_filt$Fuc.Lac <- NULL
cor_matrix_filt$Mean.Corrected.OD <- NULL
cor_matrix_filt$DNA.conc. <- NULL
cor_matrix_filt <- cor_matrix_filt[-c(1, 2, 3, 4, 5, 6, 7, 8, 9), ]

## Find the p-values for each correlation

p_matrix_filt <- matrix(NA,
                   nrow = nrow(cor_matrix_filt),
                   ncol = ncol(cor_matrix_filt),
                   dimnames = dimnames(cor_matrix_filt))

for (i in rownames(cor_matrix_filt)) {
  for (j in colnames(cor_matrix_filt)) {
    test <- cor.test(heatmap_filt[[i]], heatmap_filt[[j]], method="kendall")
    p_matrix_filt[i, j] <- test$p.value
  }
}


# Convert our correlation matrix into a long format that can be used by ggplot2

cor_long_filt <- melt(as.matrix(cor_matrix_filt))
p_long_filt <- melt(p_matrix_filt, varnames = c("Var1", "Var2"), value.name = "pvalue")
p_long_filt$Significance <- sapply(p_long_filt$pvalue, p_to_stars)

cor_final_filt <- merge(cor_long_filt, p_long_filt, by = c("Var1", "Var2"))
cor_final_filt


# Plotting the heatmap
ggplot(data = cor_final_filt, aes(x = Var1, y = Var2, fill = value)) +
  geom_tile() +
  ylab("") +
  xlab("") +
  geom_text(aes(label = paste0(round(value, 2), Significance)), color = "black", size = 4) +
  scale_fill_gradient2(low = "darkred", high = "navy", mid = "white", midpoint = 0,
                       limit = c(-1,1), space = "Lab", name = "Correlation") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))


## Do the exact same thing but only with 5.8, only targets + supplements
heatmap <- read.table("heatmap.txt", sep="\t", dec=",", header=TRUE, na.strings = c(""))
heatmap_filt <- heatmap[heatmap$initial.pH == "5.8", ]
heatmap_filt$initial.pH <- NULL
heatmap_filt$Supplement.mix <- NULL
heatmap_filt$Donor <- NULL

cor_matrix_filt2 <- data.frame(cor(heatmap_filt, method="kendall"))

## Convert correlation matrix to a have the information in the shape I want
cor_matrix_filt2$Mucin <- NULL
cor_matrix_filt2$GOS <- NULL
cor_matrix_filt2$X6c <- NULL
cor_matrix_filt2$HMOs <- NULL
cor_matrix_filt2$Fuc.Lac <- NULL
cor_matrix_filt2$Mean.Corrected.OD <- NULL
cor_matrix_filt2$DNA.conc. <- NULL
cor_matrix_filt2 <- cor_matrix_filt2[-c(6, 7, 8, 9, 10, 11), ]

## Find the p-values for each correlation

p_matrix_filt2 <- matrix(NA,
                        nrow = nrow(cor_matrix_filt2),
                        ncol = ncol(cor_matrix_filt2),
                        dimnames = dimnames(cor_matrix_filt2))

for (i in rownames(cor_matrix_filt2)) {
  for (j in colnames(cor_matrix_filt2)) {
    test <- cor.test(heatmap_filt[[i]], heatmap_filt[[j]], method="kendall")
    p_matrix_filt2[i, j] <- test$p.value
  }
}


# Convert our correlation matrix into a long format that can be used by ggplot2

cor_long_filt2 <- melt(as.matrix(cor_matrix_filt2))
p_long_filt2 <- melt(p_matrix_filt2, varnames = c("Var1", "Var2"), value.name = "pvalue")
p_long_filt2$Significance <- sapply(p_long_filt2$pvalue, p_to_stars)

cor_final_filt2 <- merge(cor_long_filt2, p_long_filt2, by = c("Var1", "Var2"))
cor_final_filt2


# Plotting the heatmap
ggplot(data = cor_final_filt2, aes(x = Var1, y = Var2, fill = value)) +
  geom_tile() +
  ylab("") +
  xlab("") +
  geom_text(aes(label = paste0(round(value, 2), Significance)), color = "black", size = 4) +
  scale_fill_gradient2(low = "darkred", high = "navy", mid = "white", midpoint = 0,
                       limit = c(-1,1), space = "Lab", name = "Correlation") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))



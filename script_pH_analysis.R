setwd("~/ETHZ/Autumn 2025/Research project/OD_pH_measurements")
df <- read.table("OD_pH.txt", sep="\t", dec=",", header=TRUE, na.strings = c(""))
df$initial.pH <- factor(df$initial.pH)
df$Supplement.mix <- factor(df$Supplement.mix, levels = c("H2O", "HMOs", "GOS", "HMOs+GOS", "HMOs+GOS+muc", "HMOs+fuc+lac+muc", "HMOs+fuc+lac+muc+GOS", "HMOs+fuc+lac+muc+GOS+6c", "muc+6c", "HMOs+muc+6c"))
df$Mucin <- factor(df$Mucin)
df$GOS <- factor(df$GOS)
df$X6c <- factor(df$X6c)
df$HMOs <- factor(df$HMOs)
df$Fuc.Lac <- factor(df$Fuc.Lac)

library(ggplot2)
library(ggpubr)
sessionInfo()

## Removing the one where the pH meter broke

which(is.na(df$X48h))
df2 <- df[-168, ]
which(df$pH.diff > 1)

## Global difference of pH change between pH

ggplot(df2, aes(x = initial.pH, y = pH.diff, fill = initial.pH)) +
  geom_boxplot(position = position_dodge(width = 0.8), alpha = 0.8) +
  theme_minimal() +
  ylab("pH difference") +
  xlab("") +
  labs(fill = "Initial pH") +
  scale_fill_manual(values = c("4.8" = "tomato2", "5.8" = "brown4")) +
  stat_compare_means(
    aes(group = 1),
    comparisons = list(c("4.8", "5.8")),
    label = "p.signif",
    method = "wilcox.test",
    paired = FALSE,
    hide.ns = TRUE)

#In order to get the stats of this comparison:
wilcox.test(df2$pH.diff[df2$initial.pH=="4.8"], df2$pH.diff[df2$initial.pH=="5.8"])

#To put in my report:
median(df2$pH.diff[df2$initial.pH=="4.8"])
wilcox.test(df2$pH.diff[df2$initial.pH=="4.8"], mu=0)

median(df2$pH.diff[df2$initial.pH=="5.8"])
wilcox.test(df2$pH.diff[df2$initial.pH=="5.8"])

med_58_pH <- median(df2$pH.diff[df2$initial.pH == "5.8"])
med_48_pH <- median(df2$pH.diff[df2$initial.pH == "4.8"])
fold_change_pH <- med_58_pH / med_48_pH
fold_change_pH

#To see if the final pH was different according to the inital pH
wilcox.test(df2$X48h.pH[df2$initial.pH=="4.8"], df2$X48h.pH[df2$initial.pH=="5.8"])
median(df2$X48h.pH[df2$initial.pH=="4.8"])
median(df2$X48h.pH[df2$initial.pH=="5.8"])

## pH change according to Supplement mix

# No statistical tests were done, because difficult with only n=3. So only visuals here

p1 <- ggplot(df2, aes(x = Supplement.mix, y = pH.diff, fill = Supplement.mix)) +
  geom_hline(yintercept = 0, color = "darkred", linewidth = 0.8, alpha = 0.8) +
  geom_boxplot(position = position_dodge(width = 0.8), alpha = 0.8) +
  theme_minimal() + 
  ylab("pH difference") +
  xlab("") +
  labs(fill = "Supplement mix") +
  facet_wrap("initial.pH") +
  theme(panel.border = element_rect(color = "grey70", fill = NA, linewidth = 0.4), strip.text = element_text(size = 12)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1), plot.title = element_text(hjust = 0.5)) +
  theme(plot.margin = unit(c(0, 0.5, 0, 0), "inches")) +
  scale_fill_manual(values = c("H2O" = "firebrick3", 
                               "HMOs" = "sienna3", 
                               "GOS" = "gold3",
                               "HMOs+GOS" = "olivedrab4",
                               "HMOs+GOS+muc" = "aquamarine4",
                               "HMOs+fuc+lac+muc" = "deepskyblue3",
                               "HMOs+fuc+lac+muc+GOS" = "royalblue4",
                               "HMOs+fuc+lac+muc+GOS+6c" = "mediumpurple4",
                               "muc+6c" = "hotpink3",
                               "HMOs+muc+6c" = "palevioletred4"))


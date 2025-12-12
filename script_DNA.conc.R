setwd("~/ETHZ/Autumn 2025/Research project/OD_DNA.conc")
DNA.conc <- read.table("OD_DNA.conc.txt", sep="\t", dec=",", header=TRUE, na.strings = c(""))
DNA.conc$Supplement.mix <- factor(DNA.conc$Supplement.mix, levels = c("H2O", "HMOs", "GOS", "HMOs+GOS", "HMOs+GOS+muc", "HMOs+fuc+lac+muc", "HMOs+fuc+lac+muc+GOS", "HMOs+fuc+lac+muc+GOS+6c", "muc+6c", "HMOs+muc+6c"))
DNA.conc$initial.pH <- factor(DNA.conc$initial.pH)

library(ggplot2)
library(ggpubr)

## Global difference of corrected OD between pH

ggplot(DNA.conc, aes(x = initial.pH, y = DNA.conc., fill = initial.pH)) +
  geom_boxplot(position = position_dodge(width = 0.8), alpha = 0.8) +
  theme_minimal() +
  ylab("DNA concentrationi [ng/uL]") +
  xlab("") +
  labs(fill = "Initial pH") +
  scale_fill_manual(values = c("4.8" = "tomato2", "5.8" = "brown4")) +
  stat_compare_means(
    aes(group = 1),
    comparisons = list(c("4.8", "5.8")),
    label = "p.signif",
    method = "t.test",
    paired = FALSE,
    hide.ns = TRUE)

t.test(DNA.conc$DNA.conc.[DNA.conc$initial.pH=="4.8"], DNA.conc$DNA.conc.[DNA.conc$initial.pH=="5.8"], var.equal = FALSE)

mean(DNA.conc$DNA.conc.[DNA.conc$initial.pH=="4.8"])
sd(DNA.conc$DNA.conc.[DNA.conc$initial.pH=="4.8"])
wilcox.test(DNA.conc$DNA.conc.[DNA.conc$initial.pH=="4.8"], mu=0)

mean(DNA.conc$DNA.conc.[DNA.conc$initial.pH=="5.8"])
sd(DNA.conc$DNA.conc.[DNA.conc$initial.pH=="5.8"])
wilcox.test(DNA.conc$DNA.conc.[DNA.conc$initial.pH=="5.8"])

mean_58_DNA.conc <- mean(DNA.conc$DNA.conc.[DNA.conc$initial.pH == "5.8"])
mean_48_DNA.conc <- mean(DNA.conc$DNA.conc.[DNA.conc$initial.pH == "4.8"])
fold_change_DNA.conc <- mean_58_DNA.conc / mean_48_DNA.conc
fold_change_DNA.conc

## OD change according to Supplement mix

# No statistical tests here because not enough observations per group (n=3)

ggplot(DNA.conc, aes(x = Supplement.mix, y = DNA.conc., fill = Supplement.mix)) +
  geom_boxplot(position = position_dodge(width = 0.8), alpha = 0.8) +
  theme_minimal() + 
  ylab("DNA concentrationi [ng/uL]") +
  xlab("") +
  labs(fill = "Supplement mix") +
  facet_wrap("initial.pH") +
  theme(panel.border = element_rect(color = "grey70", fill = NA, linewidth = 0.4), strip.text = element_text(size = 12)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1), plot.title = element_text(hjust = 0.5)) +
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


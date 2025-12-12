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

## Global difference of corrected OD between pH

ggplot(df, aes(x = initial.pH, y = Corrected.OD, fill = initial.pH)) +
  geom_boxplot(position = position_dodge(width = 0.8), alpha = 0.8) +
  theme_minimal() +
  ylab("Corrected OD value") +
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
#The data seems normally distributed, but variance is not the same. So Welch

#In order to get the stats of this comparison:
wilcox.test(df2$Corrected.OD[df2$initial.pH=="4.8"], df2$Corrected.OD[df2$initial.pH=="5.8"])

#To put in my report:
median(df$Corrected.OD[df$initial.pH=="4.8"])
wilcox.test(df$Corrected.OD[df$initial.pH=="4.8"], mu=0)

median(df$Corrected.OD[df$initial.pH=="5.8"])
wilcox.test(df$Corrected.OD[df$initial.pH=="5.8"])

med_58_OD <- median(df$Corrected.OD[df$initial.pH == "5.8"])
med_48_OD <- median(df$Corrected.OD[df$initial.pH == "4.8"])
fold_change_OD <- med_58_OD / med_48_OD
fold_change_OD

# Max OD 

## OD change according to Supplement mix

# No statistical tests here because not enough observations per group (n=3)

p2 <- ggplot(df, aes(x = Supplement.mix, y = Corrected.OD, fill = Supplement.mix)) +
  geom_hline(yintercept = 0, color = "darkred", linewidth = 0.8, alpha = 0.8) +
  geom_boxplot(position = position_dodge(width = 0.8), alpha = 0.8) +
  theme_minimal() + 
  ylab("Corrected OD value") +
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


# To arrange the plots of OD and pH change according to supplement media in one plot (note you have to run the script_pH_analysis.R before to have p1):
ggarrange(p1, p2, 
          labels = c("A", "B"), 
          common.legend = TRUE, 
          legend = "right",
          ncol = 2, nrow = 1)




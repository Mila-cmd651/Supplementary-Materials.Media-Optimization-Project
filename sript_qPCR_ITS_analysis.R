setwd("~/ETHZ/Autumn 2025/Research project/qPCR")
ITS <- read.table("qPCR_All.txt", sep="\t", dec=",", header=TRUE, na.strings = c(""))
ITS$Supplement.mix <- factor(ITS$Supplement.mix, levels = c("H2O", "HMOs", "GOS", "HMOs+GOS", "HMOs+GOS+muc", "HMOs+fuc+lac+muc", "HMOs+fuc+lac+muc+GOS", "HMOs+fuc+lac+muc+GOS+6c", "muc+6c", "HMOs+muc+6c", "Raw feces"))
ITS$pH <- factor(ITS$pH)
ITS$Mucin <- factor(ITS$Mucin)
ITS$GOS <- factor(ITS$GOS)
ITS$X6c <- factor(ITS$X6c)
ITS$HMOs <- factor(ITS$HMOs)
ITS$Fuc.Lac <- factor(ITS$Fuc.Lac)

library(ggplot2)

ITS_filtered <- ITS[ITS$pH == "5,8" | ITS$pH == "Raw feces", ]

#Condition, separated by pH

ggplot(ITS, aes(x = Supplement.mix, y = Copies.ITS.mL, fill=pH)) +
  geom_boxplot(position = position_dodge(width = 0.8), alpha = 0.8) +
  theme_minimal() +  
  coord_cartesian(ylim = c(-1e7, 2.5e7)) +
  ylab("ITS copies per mL") +
  xlab("") +
  labs(fill = "pH") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1), plot.title = element_text(hjust = 0.5)) +
  scale_fill_manual(values = c("4,8" = "tomato2", "5,8" = "brown4", "Raw feces" = "ivory2"))

#Condition, only 5.8
ITS_cond <- setdiff(unique(ITS_filtered$Supplement.mix), "Raw feces")
ITS_p_values = c()
ITS_sig = c()

for (n in ITS_cond) {
  x <-ITS_filtered$Copies.ITS.mL[ITS_filtered$Supplement.mix == n]
  y <- ITS_filtered$Copies.ITS.mL[ITS_filtered$Supplement.mix == "Raw feces"]
  wilc <- wilcox.test(x, y, paired=FALSE)
  p_val <- wilc$p.value
  ITS_p_values = c(ITS_p_values, p_val)
  if (p_val > 0.05) {
    ITS_sig <- c(ITS_sig, " ")
  } else if (p_val <= 0.05 & p_val > 0.01) {
    ITS_sig <- c(ITS_sig, "*")
  } else if (p_val <= 0.01 & p_val > 0.001) {
    ITS_sig <- c(ITS_sig, "**")
  } else if (p_val <= 0.001 & p_val > 0.0001) {
    ITS_sig <- c(ITS_sig, "***")
  } else if (p_val <= 0.0001) {
    ITS_sig <- c(ITS_sig, "****")
  }
}

ITS_stats <- data.frame("Supplement.mix" = ITS_cond, "Significance" = ITS_sig, "P values" = ITS_p_values)

ggplot(ITS_filtered, aes(x = Supplement.mix, y = Copies.ITS.mL, fill=Supplement.mix)) +
  geom_boxplot(position = position_dodge(width = 0.8), alpha = 0.8) +
  geom_jitter(aes(shape = Donor),
              width = 0.275, height = 0, size = 1.5, color = "black") +
  scale_shape_manual(values = c("D6" = 16, "D7" = 17, "D8" = 15)) +
  theme_minimal() +  
  scale_y_continuous(breaks = seq(0, 2.3e7, by = 5e6)) +
  coord_cartesian(ylim = c(0, 2.3e7)) +
  ylab("Copies per mL") +
  xlab("") +
  labs(fill = "Supplement mix") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1), legend.position = "right", legend.justification = c(0.5, 1)) +
  scale_fill_manual(values = c("H2O" = "firebrick3", 
                               "HMOs" = "sienna3", 
                               "GOS" = "gold3",
                               "HMOs+GOS" = "olivedrab4",
                               "HMOs+GOS+muc" = "aquamarine4",
                               "HMOs+fuc+lac+muc" = "deepskyblue3",
                               "HMOs+fuc+lac+muc+GOS" = "royalblue4",
                               "HMOs+fuc+lac+muc+GOS+6c" = "mediumpurple4",
                               "muc+6c" = "hotpink3",
                               "HMOs+muc+6c" = "palevioletred4",
                               "Raw feces" = "ivory2")) +
  geom_text(data = ITS_stats, aes(x = Supplement.mix, y = 2.3e7, label = ITS_sig))


ITS_med_copies <- c()
ITS_name <- c()
for (n in unique(ITS_filtered$Supplement.mix)) {
  median <- formatC(median(ITS_filtered$Copies.ITS.mL[ITS_filtered$Supplement.mix==n]), format = "e", digit = 2)
  ITS_med_copies <- c(ITS_med_copies, median)
  ITS_name <- c(ITS_name, n)
}

ITS.amounts <- data.frame("Supplement mix" = ITS_name, "Median copies" = ITS_med_copies)

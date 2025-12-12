setwd("~/ETHZ/Autumn 2025/Research project/qPCR")
Bact <- read.table("qPCR_All.txt", sep="\t", dec=",", header=TRUE, na.strings = c(""))
Bact$Supplement.mix <- factor(Bact$Supplement.mix, levels = c("H2O", "HMOs", "GOS", "HMOs+GOS", "HMOs+GOS+muc", "HMOs+fuc+lac+muc", "HMOs+fuc+lac+muc+GOS", "HMOs+fuc+lac+muc+GOS+6c", "muc+6c", "HMOs+muc+6c", "Raw feces"))
Bact$pH <- factor(Bact$pH)
Bact$Mucin <- factor(Bact$Mucin)
Bact$GOS <- factor(Bact$GOS)
Bact$X6c <- factor(Bact$X6c)
Bact$HMOs <- factor(Bact$HMOs)
Bact$Fuc.Lac <- factor(Bact$Fuc.Lac)

library(ggplot2)

Bact_filtered <- Bact[ Bact$pH == "5,8" | Bact$pH == "Raw feces", ]

#Condition, separated by pH

ggplot(Bact, aes(x = Supplement.mix, y = Copies.16S.mL, fill=pH)) +
  geom_boxplot(position = position_dodge(width = 0.8), alpha = 0.7) +
  theme_minimal() +
  scale_y_continuous(limit = c(-1e11, 1.6e12), breaks = seq(0, 1.6e12, by = 2e11)) +
  ylab("16S copies per mL") +
  xlab("") +
  labs(fill = "Initial pH") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1), plot.title = element_text(hjust = 0.5)) +
  scale_fill_manual(values = c("4,8" = "tomato2", "5,8" = "brown4", "Raw feces" = "ivory2"))


#Condition, only 5.8
unique(Bact_filtered$Supplement.mix)
Bact_cond <- setdiff(unique(Bact_filtered$Supplement.mix), "Raw feces")
Bact_cond
Bact_p_values = c()
Bact_sig = c()


for (n in Bact_cond) {
  x <-Bact_filtered$Copies.16S.mL[Bact_filtered$Supplement.mix == n]
  y <- Bact_filtered$Copies.16S.mL[Bact_filtered$Supplement.mix == "Raw feces"]
  wilc <- wilcox.test(x, y, paired=FALSE)
  p_val <- wilc$p.value
  Bact_p_values = c(Bact_p_values, p_val)
  if (p_val > 0.05) {
    Bact_sig <- c(Bact_sig, " ")
  } else if (p_val <= 0.05 & p_val > 0.01) {
    Bact_sig <- c(Bact_sig, "*")
  } else if (p_val <= 0.01 & p_val > 0.001) {
    Bact_sig <- c(Bact_sig, "**")
  } else if (p_val <= 0.001 & p_val > 0.0001) {
    Bact_sig <- c(Bact_sig, "***")
  } else if (p_val <= 0.0001) {
    Bact_sig <- c(Bact_sig, "****")
  }
}

Bact_stats <- data.frame("Supplement.mix" = Bact_cond, "Significance" = Bact_sig, "P values" = Bact_p_values)

ggplot(Bact_filtered, aes(x = Supplement.mix, y = Copies.16S.mL, fill=Supplement.mix)) +
  geom_boxplot(position = position_dodge(width = 0.8), alpha = 0.8) +
  geom_jitter(aes(shape = Donor),
              width = 0.275, height = 0, size = 1.5, color = "black") +
  scale_shape_manual(values = c("D6" = 16, "D7" = 17, "D8" = 15)) +
  theme_minimal() +  
  scale_y_continuous(breaks = seq(0, 1.6e12, by = 2e11)) +
  coord_cartesian(ylim = c(0, 1.7e12)) +
  ylab("16S copies per mL") +
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
  geom_text(data = Bact_stats, aes(x = Supplement.mix, y = 1.7e12, label = Bact_sig))

Bact_med_copies <- c()
Bact_name <- c()
for (n in unique(Bact_filtered$Supplement.mix)) {
  median <- formatC(median(Bact_filtered$Copies.16S.mL[Bact_filtered$Supplement.mix==n]), format = "e", digits = 2)
  Bact_med_copies <- c(Bact_med_copies, median)
  Bact_name <- c(Bact_name, n)
}

Bact.amounts <- data.frame("Supplement mix" = Bact_name, "Median copies" = Bact_med_copies)



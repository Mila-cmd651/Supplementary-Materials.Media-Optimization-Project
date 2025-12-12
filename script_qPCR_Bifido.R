setwd("~/ETHZ/Autumn 2025/Research project/qPCR")
Bifido <- read.table("qPCR_All.txt", sep="\t", dec=",", header=TRUE, na.strings = c(""))
Bifido$Supplement.mix <- factor(Bifido$Supplement.mix, levels = c("H2O", "HMOs", "GOS", "HMOs+GOS", "HMOs+GOS+muc", "HMOs+fuc+lac+muc", "HMOs+fuc+lac+muc+GOS", "HMOs+fuc+lac+muc+GOS+6c", "muc+6c", "HMOs+muc+6c", "Raw feces"))
Bifido$pH <- factor(Bifido$pH)
Bifido$Mucin <- factor(Bifido$Mucin)
Bifido$GOS <- factor(Bifido$GOS)
Bifido$X6c <- factor(Bifido$X6c)
Bifido$HMOs <- factor(Bifido$HMOs)
Bifido$Fuc.Lac <- factor(Bifido$Fuc.Lac)

library(ggplot2)

Bifido_filtered <- Bifido[Bifido$pH == "5,8" | Bifido$pH == "Raw feces", ]

#Condition, separated by pH

ggplot(Bifido, aes(x = Supplement.mix, y = Bifido.mL, fill=pH)) +
  geom_boxplot(position = position_dodge(width = 0.8), alpha = 0.8) +
  theme_minimal() +  
  scale_y_continuous(breaks = seq(0, 2.25e9, by = 5e8)) +
  coord_cartesian(ylim = c(0, 2.3e9)) +
  ylab("Bacteria number per mL") +
  xlab("") +
  labs(fill = "pH") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1), plot.title = element_text(hjust = 0.5)) +
  scale_fill_manual(values = c("4,8" = "tomato2", "5,8" = "brown4", "Raw feces" = "ivory2"))

#Condition, only 5.8
Bifido_cond <- setdiff(unique(Bifido_filtered$Supplement.mix), "Raw feces")
Bifido_p_values = c()
Bifido_sig = c()

for (n in Bifido_cond) {
  x <-Bifido_filtered$Bifido.mL[Bifido_filtered$Supplement.mix == n]
  y <- Bifido_filtered$Bifido.mL[Bifido_filtered$Supplement.mix == "Raw feces"]
  wilc <- wilcox.test(x, y, paired=FALSE)
  p_val <- wilc$p.value
  Bifido_p_values = c(Bifido_p_values, p_val)
  if (p_val > 0.05) {
    Bifido_sig <- c(Bifido_sig, " ")
  } else if (p_val <= 0.05 & p_val > 0.01) {
    Bifido_sig <- c(Bifido_sig, "*")
  } else if (p_val <= 0.01 & p_val > 0.001) {
    Bifido_sig <- c(Bifido_sig, "**")
  } else if (p_val <= 0.001 & p_val > 0.0001) {
    Bifido_sig <- c(Bifido_sig, "***")
  } else if (p_val <= 0.0001) {
    Bifido_sig <- c(Bifido_sig, "****")
  }
}

Bifido_stats <- data.frame("Supplement.mix" = Bifido_cond, "Significance" = Bifido_sig, "P values" = Bifido_p_values)

ggplot(Bifido_filtered, aes(x = Supplement.mix, y = Bifido.mL, fill=Supplement.mix)) +
  geom_boxplot(position = position_dodge(width = 0.8), alpha = 0.8) +
  geom_jitter(aes(shape = Donor),
              width = 0.275, height = 0, size = 1.5, color = "black") +
  scale_shape_manual(values = c("D6" = 16, "D7" = 17, "D8" = 15)) +
  theme_minimal() +  
  scale_y_continuous(breaks = seq(0, 2.25e9, by = 5e8)) +
  coord_cartesian(ylim = c(0, 2.3e9)) +
  ylab("Bacteria number per mL") +
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
  geom_text(data = Bifido_stats, aes(x = Supplement.mix, y = 2.32e9, label = Bifido_sig))

#Creating a table with the mean and sd for each mix

Bifido_med_copies <- c()
Bifido_name <- c()
for (n in unique(Bifido_filtered$Supplement.mix)) {
  median <- formatC(median(Bifido_filtered$Bifido.mL[Bifido_filtered$Supplement.mix==n]), format = "e", digit = 2)
  Bifido_med_copies <- c(Bifido_med_copies, median)
  Bifido_name <- c(Bifido_name, n)
}

Bifido.amounts <- data.frame("Supplement mix" = Bifido_name, "Median copies" = Bifido_med_copies)
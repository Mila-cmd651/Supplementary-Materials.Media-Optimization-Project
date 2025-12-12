setwd("~/ETHZ/Autumn 2025/Research project/qPCR")
Entero <- read.table("qPCR_All.txt", sep="\t", dec=",", header=TRUE, na.strings = c(""))
Entero$Supplement.mix <- factor(Entero$Supplement.mix, levels = c("H2O", "HMOs", "GOS", "HMOs+GOS", "HMOs+GOS+muc", "HMOs+fuc+lac+muc", "HMOs+fuc+lac+muc+GOS", "HMOs+fuc+lac+muc+GOS+6c", "muc+6c", "HMOs+muc+6c", "Raw feces"))
Entero$pH <- factor(Entero$pH)
Entero$Mucin <- factor(Entero$Mucin)
Entero$GOS <- factor(Entero$GOS)
Entero$X6c <- factor(Entero$X6c)
Entero$HMOs <- factor(Entero$HMOs)
Entero$Fuc.Lac <- factor(Entero$Fuc.Lac)

library(ggplot2)

Entero_filtered <- Entero[Entero$pH == "5,8" | Entero$pH == "Raw feces", ]

#Condition, separated by pH

ggplot(Entero, aes(x = Supplement.mix, y = Entero.mL, fill=pH)) +
  geom_boxplot(position = position_dodge(width = 0.8), alpha = 0.8) +
  theme_minimal() +
  scale_y_continuous(breaks = seq(0, 7.5e9, by = 1e9)) +
  coord_cartesian(ylim = c(0, 7.55e9)) +
  ylab("Bacteria number per mL") +
  xlab("") +
  labs(fill = "pH") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1), plot.title = element_text(hjust = 0.5)) +
  scale_fill_manual(values = c("4,8" = "tomato2", "5,8" = "brown4", "Raw feces" = "ivory2"))

#Condition, only 5.8
Entero_cond <- setdiff(unique(Entero_filtered$Supplement.mix), "Raw feces")
Entero_p_values = c()
Entero_sig = c()

for (n in Entero_cond) {
  x <-Entero_filtered$Entero.mL[Entero_filtered$Supplement.mix == n]
  y <- Entero_filtered$Entero.mL[Entero_filtered$Supplement.mix == "Raw feces"]
  wilc <- wilcox.test(x, y, paired=FALSE)
  p_val <- wilc$p.value
  Entero_p_values = c(Entero_p_values, p_val)
  if (p_val > 0.05) {
    Entero_sig <- c(Entero_sig, " ")
  } else if (p_val <= 0.05 & p_val > 0.01) {
    Entero_sig <- c(Entero_sig, "*")
  } else if (p_val <= 0.01 & p_val > 0.001) {
    Entero_sig <- c(Entero_sig, "**")
  } else if (p_val <= 0.001 & p_val > 0.0001) {
    Entero_sig <- c(Entero_sig, "***")
  } else if (p_val <= 0.0001) {
    Entero_sig <- c(Entero_sig, "****")
  }
}

Entero_stats <- data.frame("Supplement.mix" = Entero_cond, "Significance" = Entero_sig, "P values" = Entero_p_values)

ggplot(Entero_filtered, aes(x = Supplement.mix, y = Entero.mL, fill=Supplement.mix)) +
  geom_boxplot(position = position_dodge(width = 0.8), alpha = 0.8) +
  geom_jitter(aes(shape = Donor),
              width = 0.275, height = 0, size = 1.5, color = "black") +
  scale_shape_manual(values = c("D6" = 16, "D7" = 17, "D8" = 15)) +
  theme_minimal() +  
  scale_y_continuous(breaks = seq(0, 7.5e9, by = 1e9)) +
  coord_cartesian(ylim = c(0, 7.55e9)) +
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
  geom_text(data = Entero_stats, aes(x = Supplement.mix, y = 7.57e9, label = Entero_sig))

#Creating a table with the mean and sd for each mix

Entero_med_copies <- c()
Entero_name <- c()
for (n in unique(Entero_filtered$Supplement.mix)) {
  median <- formatC(median(Entero_filtered$Entero.mL[Entero_filtered$Supplement.mix==n]), format = "e", digit = 2)
  Entero_med_copies <- c(Entero_med_copies, median)
  Entero_name <- c(Entero_name, n)
}

Entero.amounts <- data.frame("Supplement mix" = Entero_name, "Median copies" = Entero_med_copies)

setwd("~/ETHZ/Autumn 2025/Research project/OD_DNA.conc")
OD_DNA <- read.table("OD_DNA.conc.txt", sep="\t", dec=",", header=TRUE, na.strings = c(""))
OD_DNA$Supplement.mix <- factor(OD_DNA$Supplement.mix, levels = c("H2O", "HMOs", "GOS", "HMOs+GOS", "HMOs+GOS+muc", "HMOs+fuc+lac+muc", "HMOs+fuc+lac+muc+GOS", "HMOs+fuc+lac+muc+GOS+6c", "muc+6c", "HMOs+muc+6c"))
OD_DNA$initial.pH <- factor(OD_DNA$initial.pH)

library(ggplot2)

#Testing normality of the data

qqnorm(OD_DNA$Mean.Corrected.OD)
qqline(OD_DNA$Mean.Corrected.OD)

qqnorm(OD_DNA$DNA.conc.)
qqline(OD_DNA$DNA.conc.)

plot(OD_DNA$Mean.Corrected.OD, OD_DNA$DNA.conc.)

# --> Not normal, but no ex-aequos so let's take spearman

cor.test(OD_DNA$Mean.Corrected.OD, OD_DNA$DNA.conc., method = "spearman")

ggplot(OD_DNA, aes(x = DNA.conc., y = Mean.Corrected.OD, color = Supplement.mix, shape=initial.pH)) +
  geom_point(size = 2) +
  theme_minimal() +
  xlab("DNA concentration (ng/uL)") +
  ylab("Average corrected OD value") +
  labs(color = "Supplement mix", shape = "Initial pH") +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_color_manual(values = c("H2O" = "firebrick3", 
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
  scale_shape_manual(values = c("4.8" = 16, "5.8" = 15))


# Let's do the same but only with the 5.8

OD_DNA_filtered <- OD_DNA[OD_DNA$initial.pH == 5.8, ]
length(OD_DNA_filtered$Mean.Corrected.OD)
sum(OD_DNA_filtered$Mean.Corrected.OD > 0.5)
sum(OD_DNA_filtered$Mean.Corrected.OD > 1)

cor.test(OD_DNA_filtered$Mean.Corrected.OD, OD_DNA_filtered$DNA.conc., method = "spearman")

ggplot(OD_DNA_filtered, aes(x = DNA.conc., y = Mean.Corrected.OD, color = Supplement.mix)) +
  geom_point(size = 2, shape = 16) +
  theme_minimal() +
  xlab("DNA concentration (ng/uL)") +
  ylab("Corrected OD value") +
  labs(color = "Supplement mix") +
  scale_color_manual(values = c("H2O" = "firebrick3", 
                                "HMOs" = "sienna3", 
                                "GOS" = "gold3",
                                "HMOs+GOS" = "olivedrab4",
                                "HMOs+GOS+muc" = "aquamarine4",
                                "HMOs+fuc+lac+muc" = "deepskyblue3",
                                "HMOs+fuc+lac+muc+GOS" = "royalblue4",
                                "HMOs+fuc+lac+muc+GOS+6c" = "mediumpurple4",
                                "muc+6c" = "hotpink3",
                                "HMOs+muc+6c" = "palevioletred4",
                                "Raw feces" = "ivory2"))


# If we try removing the two outliers that have a DNA.conc. > 100
which(OD_DNA_filtered$DNA.conc. > 100)
OD_DNA_test <- OD_DNA_filtered[-c(29,30), ]

cor.test(OD_DNA_test$Mean.Corrected.OD, OD_DNA_test$DNA.conc., method = "spearman")
# --> removing the outliers doesn't help

# Let's do the same but only with the 4.8

OD_DNA_filtered2 <- OD_DNA[OD_DNA$initial.pH == 4.8, ]
max(OD_DNA_filtered2$Mean.Corrected.OD)

cor.test(OD_DNA_filtered2$Mean.Corrected.OD, OD_DNA_filtered2$DNA.conc., method = "spearman")

ggplot(OD_DNA_filtered2, aes(x = DNA.conc., y = Mean.Corrected.OD, color = Supplement.mix)) +
  geom_point(size = 2) +
  theme_minimal() +
  xlab("DNA concentration (ng/uL)") +
  ylab("Corrected OD value") +
  labs(color = "Supplement mix") +
  scale_color_manual(values = c("H2O" = "firebrick3", 
                                "HMOs" = "sienna3", 
                                "GOS" = "gold3",
                                "HMOs+GOS" = "olivedrab4",
                                "HMOs+GOS+muc" = "aquamarine4",
                                "HMOs+fuc+lac+muc" = "deepskyblue3",
                                "HMOs+fuc+lac+muc+GOS" = "royalblue4",
                                "HMOs+fuc+lac+muc+GOS+6c" = "mediumpurple4",
                                "muc+6c" = "hotpink3",
                                "HMOs+muc+6c" = "palevioletred4",
                                "Raw feces" = "ivory2"))
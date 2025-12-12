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

df_filt_58 <- df[df$initial.pH=="5.8", ]
which(is.na(df_filt_58$X48h))
df_filt_582 <- df_filt_58[-78, ]

# Data normality

qqnorm(df_filt_582$pH.diff)
qqline(df_filt_582$pH.diff)

qqnorm(df_filt_582$Corrected.OD)
qqline(df_filt_582$Corrected.OD)

plot(df_filt_582$Corrected.OD, df_filt_582$pH.diff)

# --> not normal, Spearman

cor.test(df_filt_582$pH.diff, df_filt_582$Corrected.OD, method="spearman")

# --> got the warning that we have overlaps, so use Kendall

cor.test(df_filt_582$pH.diff, df_filt_582$Corrected.OD, method="kendall")


ggplot(df_filt_582, aes(x = Corrected.OD, y = pH.diff, color = Supplement.mix)) +
  geom_point(size = 2) +
  theme_minimal() +
  xlab("Corrected OD value") +
  ylab("pH difference") +
  labs(color = "Supplement mix", shape = "Initial pH") +
  scale_color_manual(values = c("H2O" = "firebrick3", 
                                "HMOs" = "sienna3", 
                                "GOS" = "gold3",
                                "HMOs+GOS" = "olivedrab4",
                                "HMOs+GOS+muc" = "aquamarine4",
                                "HMOs+fuc+lac+muc" = "deepskyblue3",
                                "HMOs+fuc+lac+muc+GOS" = "royalblue4",
                                "HMOs+fuc+lac+muc+GOS+6c" = "mediumpurple4",
                                "muc+6c" = "hotpink3",
                                "HMOs+muc+6c" = "palevioletred4"))

#### Let's do the same but with all values

which(is.na(df$X48h))
df_nafilt <- df[-168, ]

# Data normality

qqnorm(df_nafilt$pH.diff)
qqline(df_nafilt$pH.diff)

qqnorm(df_nafilt$Corrected.OD)
qqline(df_nafilt$Corrected.OD)

plot(df_nafilt$Corrected.OD, df_nafilt$pH.diff)

# Still not normal, let's use Kendall

cor.test(df_nafilt$pH.diff, df_nafilt$Corrected.OD, method="kendall")

ggplot(df_nafilt, aes(x = Corrected.OD, y = pH.diff, color = Supplement.mix, shape = initial.pH)) +
  geom_point(size = 2) +
  theme_minimal() +
  xlab("Corrected OD value") +
  ylab("pH difference") +
  labs(color = "Supplement mix", shape = "Initial pH") +
  scale_color_manual(values = c("H2O" = "firebrick3", 
                                "HMOs" = "sienna3", 
                                "GOS" = "gold3",
                                "HMOs+GOS" = "olivedrab4",
                                "HMOs+GOS+muc" = "aquamarine4",
                                "HMOs+fuc+lac+muc" = "deepskyblue3",
                                "HMOs+fuc+lac+muc+GOS" = "royalblue4",
                                "HMOs+fuc+lac+muc+GOS+6c" = "mediumpurple4",
                                "muc+6c" = "hotpink3",
                                "HMOs+muc+6c" = "palevioletred4")) +
  scale_shape_manual(values = c("4.8" = 16, "5.8" = 15))


ggplot(df_nafilt, aes(x = Corrected.OD, y = pH.diff)) +
  geom_point(aes(color = Supplement.mix, shape = initial.pH), size = 2) +
  stat_ellipse(data = df_nafilt,
               aes(x = Corrected.OD, y = pH.diff, color = initial.pH, group = initial.pH),
               type = "euclid",
               linetype = 2,
               inherit.aes = FALSE) +
  theme_minimal() +
  xlab("Corrected OD value") +
  ylab("pH difference") +
  labs(color = "Supplement mix and initial pH", shape = "Initial pH") +
  scale_color_manual(
    values = c(
      # Colors for points
      "H2O" = "firebrick3", 
      "HMOs" = "sienna3", 
      "GOS" = "gold3",
      "HMOs+GOS" = "olivedrab4",
      "HMOs+GOS+muc" = "aquamarine4",
      "HMOs+fuc+lac+muc" = "deepskyblue3",
      "HMOs+fuc+lac+muc+GOS" = "royalblue4",
      "HMOs+fuc+lac+muc+GOS+6c" = "mediumpurple4",
      "muc+6c" = "hotpink3",
      "HMOs+muc+6c" = "palevioletred4",
      # Colors for ellipses
      "4.8" = "tomato2",
      "5.8" = "brown4"
    )
  ) +
  scale_shape_manual(values = c("4.8" = 16, "5.8" = 15))

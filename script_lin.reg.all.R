setwd("~/ETHZ/Autumn 2025/Research project/qPCR")
all <- read.table("qPCR_All.txt", sep="\t", dec=",", header=TRUE, na.strings = c(""))
all$Supplement.mix <- factor(all$Supplement.mix, levels = c("H2O", "HMOs", "GOS", "HMOs+GOS", "HMOs+GOS+muc", "HMOs+fuc+lac+muc", "HMOs+fuc+lac+muc+GOS", "HMOs+fuc+lac+muc+GOS+6c", "muc+6c", "HMOs+muc+6c", "Raw feces"))
all$pH <- factor(all$pH)
all$Mucin <- factor(all$Mucin)
all$GOS <- factor(all$GOS)
all$X6c <- factor(all$X6c)
all$HMOs <- factor(all$HMOs)
all$Fuc.Lac <- factor(all$Fuc.Lac)

library(ggplot2)
sessionInfo()

all_filtered <- all[ all$pH == "5,8", ]


## Do a regression for ITS, with supplements and donor as explication variables

# Testing normality

m.ITS.All <- lm(Copies.ITS.mL ~ Mucin + GOS + X6c + HMOs + Fuc.Lac + Donor, data = all_filtered)
m.ITS.Donors <- lm(Copies.ITS.mL ~ Donor, data = all_filtered)
m.ITS.Supp <- lm(Copies.ITS.mL ~ Mucin + GOS + X6c + HMOs + Fuc.Lac, data = all_filtered)
par(mfrow=c(2,2))
plot(m.ITS.All)
plot(m.ITS.Donors)
plot(m.ITS.Supp)
# --> yes it's alright

summary(m.ITS.All)
summary(m.ITS.Donors)
summary(m.ITS.Supp)

# Comparing the effect on the impact of Donors and Supplements
anova(m.ITS.Donors, m.ITS.All)
anova(m.ITS.Supp, m.ITS.All)


## Same but for 16S
m.Bact.All <- lm(Copies.16S.mL ~ Mucin + GOS + X6c + HMOs + Fuc.Lac + Donor, data = all_filtered)
m.Bact.Donors <- lm(Copies.16S.mL ~ Donor, data = all_filtered)
m.Bact.Supp <- lm(Copies.16S.mL ~ Mucin + GOS + X6c + HMOs + Fuc.Lac, data = all_filtered)
par(mfrow=c(2,2))
plot(m.Bact.All)
plot(m.Bact.Donors)
plot(m.Bact.Supp)
# --> yes it's alright

summary(m.Bact.All)
summary(m.Bact.Donors)
summary(m.Bact.Supp)

# Comparing the effect on the impact of Donors and Supplements
anova(m.Bact.Donors, m.Bact.All)
anova(m.Bact.Supp, m.Bact.All)


## Same but for Bifido
m.Bifido.All <- lm(Bifido.mL ~ Mucin + GOS + X6c + HMOs + Fuc.Lac + Donor, data = all_filtered)
m.Bifido.Donors <- lm(Bifido.mL ~ Donor, data = all_filtered)
m.Bifido.Supp <- lm(Bifido.mL ~ Mucin + GOS + X6c + HMOs + Fuc.Lac, data = all_filtered)
par(mfrow=c(2,2))
plot(m.Bifido.All)
plot(m.Bifido.Donors)
plot(m.Bifido.Supp)
# --> yes it's alright

summary(m.Bifido.All)
summary(m.Bifido.Donors)
summary(m.Bifido.Supp)

# Comparing the effect on the impact of Donors and Supplements
anova(m.Bifido.Donors, m.Bifido.All)
anova(m.Bifido.Supp, m.Bifido.All)


## Same but for Entero
m.Entero.All <- lm(Entero.mL ~ Mucin + GOS + X6c + HMOs + Fuc.Lac + Donor, data = all_filtered)
m.Entero.Donors <- lm(Entero.mL ~ Donor, data = all_filtered)
m.Entero.Supp <- lm(Entero.mL ~ Mucin + GOS + X6c + HMOs + Fuc.Lac, data = all_filtered)
par(mfrow=c(2,2))
plot(m.Entero.All)
plot(m.Entero.Donors)
plot(m.Entero.Supp)
# --> yes it's alright

summary(m.Entero.All)
summary(m.Entero.Donors)
summary(m.Entero.Supp)

# Comparing the effect on the impact of Donors and Supplements
anova(m.Entero.Donors, m.Entero.All)
anova(m.Entero.Supp, m.Entero.All)

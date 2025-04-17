#################################
## Partie 1 modèles statistiques ##
#################################

library(phyloseq)
library(lme4)
library(lmerTest)
library(nlme)
library(dplyr)
library(lattice)

## Recover the phylum rank from our phyloseq object ## 
physeq_phylum <- tax_glom(physeq_normalize, taxrank = "Phylum")
physeq_family <- tax_glom(physeq_normalize, taxrank = "Family")
physeq_genus <- tax_glom(physeq_normalize, taxrank = "Genus")

# Melt phyloseq object
phylum <- psmelt(physeq_phylum)
family <- psmelt(physeq_family)
genus <- psmelt(physeq_genus)

# View the updated dataset
head(phylum)
head(family)
head(genus)

# Ensure variables are factors
phylum$Condition <- factor(phylum$Condition, levels = c("Untreated", "Treated"))
phylum$Calves <- as_factor(phylum$Calves)

family$Condition <- factor(family$Condition, levels = c("Untreated", "Treated"))
family$Calves <- as_factor(family$Calves)

genus$Condition <- factor(genus$Condition, levels = c("Untreated", "Treated"))
genus$Calves <- as_factor(genus$Calves)

## Time_continuous : time as a continuous variable and not a factor for further analysis ## 
# Create a mapping for Time levels
time_mapping <- c("T0" = 0, "T1" = 6, "T2" = 30, "T3" = 60)
# Map Time factor to continuous values
phylum$Time_continuous <- time_mapping[as.character(phylum$Time)]
phylum <-phylum[phylum$Time_continuous !=60,]

family$Time_continuous <- time_mapping[as.character(family$Time)]
family <-family[family$Time_continuous !=60,]

genus$Time_continuous <- time_mapping[as.character(genus$Time)]
genus <-genus[genus$Time_continuous !=60,]

# Split data by Phylum
phylum_list <- split(phylum, phylum$Phylum)

names(phylum_list)

# Split data by Family
family_list <- split(family, family$Family)
names(family_list)

# Split data by Genus
genus_list <- split(genus, genus$Genus)
names(genus_list)



## Choose "Bacillota", "Bacteroidota" "Actinomycetota" "Pseudomonadota" for the Phylum
# Access a specific phylum by name
bacillota_phylum <- phylum_list[["Bacillota"]]
bacteroidota_phylum <- phylum_list[["Bacteroidota"]]
actinomycetota_phylum <- phylum_list[["Actinomycetota"]]
pseudomonadota_phylum <- phylum_list[["Pseudomonadota"]]

## BACILLOTA ## 
# Create bar plot
ggplot(bacillota_phylum, aes(x = Time_continuous, y = Abundance)) +
  geom_line(aes(group = Calves), size = 1, color = "blue") +  # Line per calf
  geom_point(size = 2) +  # Points for observations
  facet_grid(Condition ~ Calves, scales = "free_y") +  # Separate rows by Condition and columns by Calf
  labs(
    x = "Time",
    y = "Abundance",
    title = "Abundance of Bacillota Over Time (Treated vs Untreated)"
  ) +
  theme_bw() +
  theme(
    strip.background = element_blank(),
    strip.text.x = element_text(size = 10),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

## Modèle 
## PHYLUM 
#BACILLOTA
model_bacillota <- lme(Abundance ~ Time_continuous*Condition+ I(Time_continuous^2)*Condition, random = ~1+Time_continuous|Calves, weights = varIdent(~1|Condition), data = bacillota_phylum,control=list(msMaxIter=2000,opt="optim",msVerbose=TRUE))
summary(model_bacillota)
plot(model_bacillota, resid(.)~Time_continuous)
qq <- qqnorm(resid(model_bacillota), main = "QQ Plot with Sample Labels ModeL", xlab = "Theoretical Quantiles", ylab = "Sample Quantiles")
qqline(resid(model_bacillota), col = "red")  # Add reference line

# Plot des effets fixes
pdata <- expand.grid(Time_continuous=c(0,6,30), Condition=c('Untreated', 'Treated'))
pdata$Abundance <- predict(model_bacillota, pdata, level=0)
plot(pdata$Time_continuous, pdata$Abundance, type='n', xlab='Time', ylab=' Relative Abundance')
points(pdata$Time_continuous[1:3], pdata$Abundance[1:3], type='b', pch=19, lwd=2)
points(pdata$Time_continuous[4:6], pdata$Abundance[4:6], type='b', pch=22, lwd=2,col="hotpink2") ## For Treated Calves 

##BACTEROIDOTA
model_bacteroidota <- lme(Abundance ~ Time_continuous*Condition+ I(Time_continuous^2)*Condition, random = ~1+Time_continuous|Calves, weights = varIdent(~1|Condition), data = bacteroidota_phylum,control=list(msMaxIter=2000,opt="optim",msVerbose=TRUE))
summary(model_bacteroidota)
plot(model_bacteroidota, resid(.)~Time_continuous)
qq <- qqnorm(resid(model_bacteroidota), main = "QQ Plot with Sample Labels Model", xlab = "Theoretical Quantiles", ylab = "Sample Quantiles")
qqline(resid(model_bacteroidota), col = "red")  # Add reference line

# Plot des effets fixes
pdata <- expand.grid(Time_continuous=c(0,6,30), Condition=c('Untreated', 'Treated'))
pdata$Abundance <- predict(model_bacteroidota, pdata, level=0)
plot(pdata$Time_continuous, pdata$Abundance, type='n', xlab='Time', ylab=' Relative Abundance')
points(pdata$Time_continuous[1:3], pdata$Abundance[1:3], type='b', pch=19, lwd=2)
points(pdata$Time_continuous[4:6], pdata$Abundance[4:6], type='b', pch=22, lwd=2,col="hotpink2") ## For Treated Calves 

##ACTINOMYCETOTA
model_actinomycetota <- lme(Abundance ~ Time_continuous*Condition+ I(Time_continuous^2)*Condition, random = ~1+Time_continuous|Calves, weights = varIdent(~1|Condition), data = actinomycetota_phylum,control=list(msMaxIter=2000,opt="optim",msVerbose=TRUE))
summary(model_actinomycetota)
plot(model_actinomycetota, resid(.)~Time_continuous)
qq <- qqnorm(resid(model_actinomycetota), main = "QQ Plot with Sample Labels Model", xlab = "Theoretical Quantiles", ylab = "Sample Quantiles")
qqline(resid(model_actinomycetota), col = "red")  # Add reference line

# Plot des effets fixes
pdata <- expand.grid(Time_continuous=c(0,6,30), Condition=c('Untreated', 'Treated'))
pdata$Abundance <- predict(model_actinomycetota, pdata, level=0)
plot(pdata$Time_continuous, pdata$Abundance, type='n', xlab='Time', ylab=' Relative Abundance')
points(pdata$Time_continuous[1:3], pdata$Abundance[1:3], type='b', pch=19, lwd=2)
points(pdata$Time_continuous[4:6], pdata$Abundance[4:6], type='b', pch=22, lwd=2,col="hotpink2") ## For Treated Calves 

##PSEUDOMONADOTA
model_pseudomonadota <- lme(Abundance ~ Time_continuous*Condition+ I(Time_continuous^2)*Condition, random = ~1+Time_continuous|Calves, weights = varIdent(~1|Condition), data = pseudomonadota_phylum,control=list(msMaxIter=2000,opt="optim",msVerbose=TRUE))
summary(model_pseudomonadota)
plot(model_pseudomonadota, resid(.)~Time_continuous)
qq <- qqnorm(resid(model_pseudomonadota), main = "QQ Plot with Sample Labels Model", xlab = "Theoretical Quantiles", ylab = "Sample Quantiles")
qqline(resid(model_pseudomonadota), col = "red")  # Add reference line

# Plot des effets fixes
pdata <- expand.grid(Time_continuous=c(0,6,30), Condition=c('Untreated', 'Treated'))
pdata$Abundance <- predict(model_pseudomonadota, pdata, level=0)
plot(pdata$Time_continuous, pdata$Abundance, type='n', xlab='Time', ylab=' Relative Abundance')
points(pdata$Time_continuous[1:3], pdata$Abundance[1:3], type='b', pch=19, lwd=2)
points(pdata$Time_continuous[4:6], pdata$Abundance[4:6], type='b', pch=22, lwd=2,col="hotpink2") ## For Treated Calves 

## FAMILY
## Choose 
# Access a specific family by name
acidaminococcaceae_family <- family_list[["Acidaminococcaceae"]]
anaerotignaceae_family <- family_list[["Anaerotignaceae"]]
bacteroidaceae_family <- family_list[["Bacteroidaceae"]]
clostridiaceae_family <- family_list[["Clostridiaceae"]]
enterobacteriaceae_family <- family_list[["Enterobacteriaceae"]]
enterococcaceae_family <- family_list[["Enterococcaceae"]]
eubacteriales_incertae_sedis_family <- family_list[["Eubacteriales_Incertae_sedis"]]
lachnospiraceae_family <- family_list[["Lachnospiraceae"]]
lactobacillaceae_family <- family_list[["Lactobacillaceae"]]
oscillospiraceae_family <- family_list[["Oscillospiraceae"]]
pasteurellaceae_family <- family_list[["Pasteurellaceae"]]
peptostreptococcaceae_family <- family_list[["Peptostreptococcaceae"]]
prevotellaceae_family <- family_list[["Prevotellaceae"]]
selenomonadaceae_family <- family_list[["Selenomonadaceae"]]
streptococcaceae_family <- family_list[["Streptococcaceae"]]


#Acidaminococcaceae
model_acidaminococcaceae <- lme(Abundance ~ Time_continuous*Condition+ I(Time_continuous^2)*Condition, random = ~1+Time_continuous|Calves, weights = varIdent(~1|Condition), data = acidaminococcaceae_family,control=list(msMaxIter=2000,opt="optim",msVerbose=TRUE))
summary(model_acidaminococcaceae)
plot(model_acidaminococcaceae, resid(.)~Time_continuous)
qq <- qqnorm(resid(model_acidaminococcaceae), main = "QQ Plot with Sample Labels Model", xlab = "Theoretical Quantiles", ylab = "Sample Quantiles")
qqline(resid(model_acidaminococcaceae), col = "red")  # Add reference line

# Plot des effets fixes
pdata <- expand.grid(Time_continuous=c(0,6,30), Condition=c('Untreated', 'Treated'))
pdata$Abundance <- predict(model_acidaminococcaceae, pdata, level=0)
plot(pdata$Time_continuous, pdata$Abundance, type='n', xlab='Time', ylab=' Relative Abundance')
points(pdata$Time_continuous[1:3], pdata$Abundance[1:3], type='b', pch=19, lwd=2)
points(pdata$Time_continuous[4:6], pdata$Abundance[4:6], type='b', pch=22, lwd=2,col="hotpink2") ## For Treated Calves 

#Anaerotignaceae
model_anaerotignaceae <- lme(Abundance ~ Time_continuous*Condition+ I(Time_continuous^2)*Condition, random = ~1+Time_continuous|Calves, weights = varIdent(~1|Condition), data = anaerotignaceae_family,control=list(msMaxIter=2000,opt="optim",msVerbose=TRUE))
summary(model_anaerotignaceae)
plot(model_anaerotignaceae, resid(.)~Time_continuous)
qq <- qqnorm(resid(model_anaerotignaceae), main = "QQ Plot with Sample Labels Model", xlab = "Theoretical Quantiles", ylab = "Sample Quantiles")
qqline(resid(model_anaerotignaceae), col = "red")  # Add reference line

# Plot des effets fixes
pdata <- expand.grid(Time_continuous=c(0,6,30), Condition=c('Untreated', 'Treated'))
pdata$Abundance <- predict(model_anaerotignaceae, pdata, level=0)
plot(pdata$Time_continuous, pdata$Abundance, type='n', xlab='Time', ylab=' Relative Abundance')
points(pdata$Time_continuous[1:3], pdata$Abundance[1:3], type='b', pch=19, lwd=2)
points(pdata$Time_continuous[4:6], pdata$Abundance[4:6], type='b', pch=22, lwd=2,col="hotpink2") ## For Treated Calves 

#Bacteroidaceae
model_bacteroidaceae <- lme(Abundance ~ Time_continuous*Condition+ I(Time_continuous^2)*Condition, random = ~1+Time_continuous|Calves, weights = varIdent(~1|Condition), data = bacteroidaceae_family,control=list(msMaxIter=2000,opt="optim",msVerbose=TRUE))
summary(model_bacteroidaceae)
plot(model_bacteroidaceae, resid(.)~Time_continuous)
qq <- qqnorm(resid(model_bacteroidaceae), main = "QQ Plot with Sample Labels Model", xlab = "Theoretical Quantiles", ylab = "Sample Quantiles")
qqline(resid(model_bacteroidaceae), col = "red")  # Add reference line

# Plot des effets fixes
pdata <- expand.grid(Time_continuous=c(0,6,30), Condition=c('Untreated', 'Treated'))
pdata$Abundance <- predict(model_bacteroidaceae, pdata, level=0)
plot(pdata$Time_continuous, pdata$Abundance, type='n', xlab='Time', ylab=' Relative Abundance')
points(pdata$Time_continuous[1:3], pdata$Abundance[1:3], type='b', pch=19, lwd=2)
points(pdata$Time_continuous[4:6], pdata$Abundance[4:6], type='b', pch=22, lwd=2,col="hotpink2") ## For Treated Calves 

#Clostridiaceae
model_clostridiaceae <- lme(Abundance ~ Time_continuous*Condition+ I(Time_continuous^2)*Condition, random = ~1+Time_continuous|Calves, weights = varIdent(~1|Condition), data = clostridiaceae_family,control=list(msMaxIter=2000,opt="optim",msVerbose=TRUE))
summary(model_clostridiaceae)
plot(model_clostridiaceae, resid(.)~Time_continuous)
qq <- qqnorm(resid(model_clostridiaceae), main = "QQ Plot with Sample Labels Model", xlab = "Theoretical Quantiles", ylab = "Sample Quantiles")
qqline(resid(model_clostridiaceae), col = "red")  # Add reference line

# Plot des effets fixes
pdata <- expand.grid(Time_continuous=c(0,6,30), Condition=c('Untreated', 'Treated'))
pdata$Abundance <- predict(model_clostridiaceae, pdata, level=0)
plot(pdata$Time_continuous, pdata$Abundance, type='n', xlab='Time', ylab=' Relative Abundance')
points(pdata$Time_continuous[1:3], pdata$Abundance[1:3], type='b', pch=19, lwd=2)
points(pdata$Time_continuous[4:6], pdata$Abundance[4:6], type='b', pch=22, lwd=2,col="hotpink2") ## For Treated Calves 

#Enterobacteriaceae
model_enterobacteriaceae <- lme(Abundance ~ Time_continuous*Condition+ I(Time_continuous^2)*Condition, random = ~1+Time_continuous|Calves, weights = varIdent(~1|Condition), data = enterobacteriaceae_family,control=list(msMaxIter=2000,opt="optim",msVerbose=TRUE))
summary(model_enterobacteriaceae)
plot(model_enterobacteriaceae, resid(.)~Time_continuous)
qq <- qqnorm(resid(model_enterobacteriaceae), main = "QQ Plot with Sample Labels Model", xlab = "Theoretical Quantiles", ylab = "Sample Quantiles")
qqline(resid(model_enterobacteriaceae), col = "red")  # Add reference line

# Plot des effets fixes
pdata <- expand.grid(Time_continuous=c(0,6,30), Condition=c('Untreated', 'Treated'))
pdata$Abundance <- predict(model_enterobacteriaceae, pdata, level=0)
plot(pdata$Time_continuous, pdata$Abundance, type='n', xlab='Time', ylab=' Relative Abundance')
points(pdata$Time_continuous[1:3], pdata$Abundance[1:3], type='b', pch=19, lwd=2)
points(pdata$Time_continuous[4:6], pdata$Abundance[4:6], type='b', pch=22, lwd=2,col="hotpink2") ## For Treated Calves 

#Enterococcaceae
model_enterococcaceae <- lme(Abundance ~ Time_continuous*Condition+ I(Time_continuous^2)*Condition, random = ~1+Time_continuous|Calves, weights = varIdent(~1|Condition), data = enterococcaceae_family,control=list(msMaxIter=2000,opt="optim",msVerbose=TRUE))
summary(model_enterococcaceae)
plot(model_enterococcaceae, resid(.)~Time_continuous)
qq <- qqnorm(resid(model_enterococcaceae), main = "QQ Plot with Sample Labels Model", xlab = "Theoretical Quantiles", ylab = "Sample Quantiles")
qqline(resid(model_enterococcaceae), col = "red")  # Add reference line

# Plot des effets fixes
pdata <- expand.grid(Time_continuous=c(0,6,30), Condition=c('Untreated', 'Treated'))
pdata$Abundance <- predict(model_enterococcaceae, pdata, level=0)
plot(pdata$Time_continuous, pdata$Abundance, type='n', xlab='Time', ylab=' Relative Abundance')
points(pdata$Time_continuous[1:3], pdata$Abundance[1:3], type='b', pch=19, lwd=2)
points(pdata$Time_continuous[4:6], pdata$Abundance[4:6], type='b', pch=22, lwd=2,col="hotpink2") ## For Treated Calves 

#Eubacteriales_incertae_sedis
model_eubacteriales_incertae_sedis <- lme(Abundance ~ Time_continuous*Condition+ I(Time_continuous^2)*Condition, random = ~1+Time_continuous|Calves, weights = varIdent(~1|Condition), data = eubacteriales_incertae_sedis_family,control=list(msMaxIter=2000,opt="optim",msVerbose=TRUE))
summary(model_eubacteriales_incertae_sedis)
plot(model_eubacteriales_incertae_sedis, resid(.)~Time_continuous)
qq <- qqnorm(resid(model_eubacteriales_incertae_sedis), main = "QQ Plot with Sample Labels Model", xlab = "Theoretical Quantiles", ylab = "Sample Quantiles")
qqline(resid(model_eubacteriales_incertae_sedis), col = "red")  # Add reference line

# Plot des effets fixes
pdata <- expand.grid(Time_continuous=c(0,6,30), Condition=c('Untreated', 'Treated'))
pdata$Abundance <- predict(model_eubacteriales_incertae_sedis, pdata, level=0)
plot(pdata$Time_continuous, pdata$Abundance, type='n', xlab='Time', ylab=' Relative Abundance')
points(pdata$Time_continuous[1:3], pdata$Abundance[1:3], type='b', pch=19, lwd=2)
points(pdata$Time_continuous[4:6], pdata$Abundance[4:6], type='b', pch=22, lwd=2,col="hotpink2") ## For Treated Calves 

#Lachnospiraceae
model_lachnospiraceae <- lme(Abundance ~ Time_continuous*Condition+ I(Time_continuous^2)*Condition, random = ~1+Time_continuous|Calves, weights = varIdent(~1|Condition), data = lachnospiraceae_family,control=list(msMaxIter=2000,opt="optim",msVerbose=TRUE))
summary(model_lachnospiraceae)
plot(model_lachnospiraceae, resid(.)~Time_continuous)
qq <- qqnorm(resid(model_lachnospiraceae), main = "QQ Plot with Sample Labels Model", xlab = "Theoretical Quantiles", ylab = "Sample Quantiles")
qqline(resid(model_lachnospiraceae), col = "red")  # Add reference line

# Plot des effets fixes
pdata <- expand.grid(Time_continuous=c(0,6,30), Condition=c('Untreated', 'Treated'))
pdata$Abundance <- predict(model_lachnospiraceae, pdata, level=0)
plot(pdata$Time_continuous, pdata$Abundance, type='n', xlab='Time', ylab=' Relative Abundance')
points(pdata$Time_continuous[1:3], pdata$Abundance[1:3], type='b', pch=19, lwd=2)
points(pdata$Time_continuous[4:6], pdata$Abundance[4:6], type='b', pch=22, lwd=2,col="hotpink2") ## For Treated Calves 


#Lactobacillaceae
model_lactobacillaceae <- lme(Abundance ~ Time_continuous*Condition+ I(Time_continuous^2)*Condition, random = ~1+Time_continuous|Calves, weights = varIdent(~1|Condition), data = lactobacillaceae_family,control=list(msMaxIter=2000,opt="optim",msVerbose=TRUE))
summary(model_lactobacillaceae)
plot(model_lactobacillaceae, resid(.)~Time_continuous)
qq <- qqnorm(resid(model_lactobacillaceae), main = "QQ Plot with Sample Labels Model", xlab = "Theoretical Quantiles", ylab = "Sample Quantiles")
qqline(resid(model_lactobacillaceae), col = "red")  # Add reference line

# Plot des effets fixes
pdata <- expand.grid(Time_continuous=c(0,6,30), Condition=c('Untreated', 'Treated'))
pdata$Abundance <- predict(model_lactobacillaceae, pdata, level=0)
plot(pdata$Time_continuous, pdata$Abundance, type='n', xlab='Time', ylab=' Relative Abundance')
points(pdata$Time_continuous[1:3], pdata$Abundance[1:3], type='b', pch=19, lwd=2)
points(pdata$Time_continuous[4:6], pdata$Abundance[4:6], type='b', pch=22, lwd=2,col="hotpink2") ## For Treated Calves 

#Oscillospiraceae
model_oscillospiraceae <- lme(Abundance ~ Time_continuous*Condition+ I(Time_continuous^2)*Condition, random = ~1+Time_continuous|Calves, weights = varIdent(~1|Condition), data = oscillospiraceae_family,control=list(msMaxIter=2000,opt="optim",msVerbose=TRUE))
summary(model_oscillospiraceae)
plot(model_oscillospiraceae, resid(.)~Time_continuous)
qq <- qqnorm(resid(model_oscillospiraceae), main = "QQ Plot with Sample Labels Model", xlab = "Theoretical Quantiles", ylab = "Sample Quantiles")
qqline(resid(model_oscillospiraceae), col = "red")  # Add reference line

# Plot des effets fixes
pdata <- expand.grid(Time_continuous=c(0,6,30), Condition=c('Untreated', 'Treated'))
pdata$Abundance <- predict(model_oscillospiraceae, pdata, level=0)
plot(pdata$Time_continuous, pdata$Abundance, type='n', xlab='Time', ylab=' Relative Abundance')
points(pdata$Time_continuous[1:3], pdata$Abundance[1:3], type='b', pch=19, lwd=2)
points(pdata$Time_continuous[4:6], pdata$Abundance[4:6], type='b', pch=22, lwd=2,col="hotpink2") ## For Treated Calves 

#Pasteurellaceae
model_pasteurellaceae <- lme(Abundance ~ Time_continuous*Condition+ I(Time_continuous^2)*Condition, random = ~1+Time_continuous|Calves, weights = varIdent(~1|Condition), data = pasteurellaceae_family,control=list(msMaxIter=2000,opt="optim",msVerbose=TRUE))
summary(model_pasteurellaceae)
plot(model_pasteurellaceae, resid(.)~Time_continuous)
qq <- qqnorm(resid(model_pasteurellaceae), main = "QQ Plot with Sample Labels Model", xlab = "Theoretical Quantiles", ylab = "Sample Quantiles")
qqline(resid(model_pasteurellaceae), col = "red")  # Add reference line

# Plot des effets fixes
pdata <- expand.grid(Time_continuous=c(0,6,30), Condition=c('Untreated', 'Treated'))
pdata$Abundance <- predict(model_pasteurellaceae, pdata, level=0)
plot(pdata$Time_continuous, pdata$Abundance, type='n', xlab='Time', ylab=' Relative Abundance')
points(pdata$Time_continuous[1:3], pdata$Abundance[1:3], type='b', pch=19, lwd=2)
points(pdata$Time_continuous[4:6], pdata$Abundance[4:6], type='b', pch=22, lwd=2,col="hotpink2") ## For Treated Calves 

#Peptostreptococcaceae
model_peptostreptococcaceae <- lme(Abundance ~ Time_continuous*Condition+ I(Time_continuous^2)*Condition, random = ~1+Time_continuous|Calves, weights = varIdent(~1|Condition), data = peptostreptococcaceae_family,control=list(msMaxIter=2000,opt="optim",msVerbose=TRUE))
summary(model_peptostreptococcaceae)
plot(model_peptostreptococcaceae, resid(.)~Time_continuous)
qq <- qqnorm(resid(model_peptostreptococcaceae), main = "QQ Plot with Sample Labels Model", xlab = "Theoretical Quantiles", ylab = "Sample Quantiles")
qqline(resid(model_peptostreptococcaceae), col = "red")  # Add reference line

# Plot des effets fixes
pdata <- expand.grid(Time_continuous=c(0,6,30), Condition=c('Untreated', 'Treated'))
pdata$Abundance <- predict(model_peptostreptococcaceae, pdata, level=0)
plot(pdata$Time_continuous, pdata$Abundance, type='n', xlab='Time', ylab=' Relative Abundance')
points(pdata$Time_continuous[1:3], pdata$Abundance[1:3], type='b', pch=19, lwd=2)
points(pdata$Time_continuous[4:6], pdata$Abundance[4:6], type='b', pch=22, lwd=2,col="hotpink2") ## For Treated Calves 

#Prevotellaceae
model_prevotellaceae <- lme(Abundance ~ Time_continuous*Condition+ I(Time_continuous^2)*Condition, random = ~1+Time_continuous|Calves, weights = varIdent(~1|Condition), data = prevotellaceae_family,control=list(msMaxIter=2000,opt="optim",msVerbose=TRUE))
summary(model_prevotellaceae)
plot(model_prevotellaceae, resid(.)~Time_continuous)
qq <- qqnorm(resid(model_prevotellaceae), main = "QQ Plot with Sample Labels Model", xlab = "Theoretical Quantiles", ylab = "Sample Quantiles")
qqline(resid(model_prevotellaceae), col = "red")  # Add reference line

# Plot des effets fixes
pdata <- expand.grid(Time_continuous=c(0,6,30), Condition=c('Untreated', 'Treated'))
pdata$Abundance <- predict(model_prevotellaceae, pdata, level=0)
plot(pdata$Time_continuous, pdata$Abundance, type='n', xlab='Time', ylab=' Relative Abundance')
points(pdata$Time_continuous[1:3], pdata$Abundance[1:3], type='b', pch=19, lwd=2)
points(pdata$Time_continuous[4:6], pdata$Abundance[4:6], type='b', pch=22, lwd=2,col="hotpink2") ## For Treated Calves 

#Selenomonadaceae
model_selenomonadaceae <- lme(Abundance ~ Time_continuous*Condition+ I(Time_continuous^2)*Condition, random = ~1+Time_continuous|Calves, weights = varIdent(~1|Condition), data = selenomonadaceae_family,control=list(msMaxIter=2000,opt="optim",msVerbose=TRUE))
summary(model_selenomonadaceae)
plot(model_selenomonadaceae, resid(.)~Time_continuous)
qq <- qqnorm(resid(model_selenomonadaceae), main = "QQ Plot with Sample Labels Model", xlab = "Theoretical Quantiles", ylab = "Sample Quantiles")
qqline(resid(model_selenomonadaceae), col = "red")  # Add reference line

# Plot des effets fixes
pdata <- expand.grid(Time_continuous=c(0,6,30), Condition=c('Untreated', 'Treated'))
pdata$Abundance <- predict(model_selenomonadaceae, pdata, level=0)
plot(pdata$Time_continuous, pdata$Abundance, type='n', xlab='Time', ylab=' Relative Abundance')
points(pdata$Time_continuous[1:3], pdata$Abundance[1:3], type='b', pch=19, lwd=2)
points(pdata$Time_continuous[4:6], pdata$Abundance[4:6], type='b', pch=22, lwd=2,col="hotpink2") ## For Treated Calves 

#Streptococcaceae
model_streptococcaceae <- lme(Abundance ~ Time_continuous*Condition+ I(Time_continuous^2)*Condition, random = ~1+Time_continuous|Calves, weights = varIdent(~1|Condition), data = streptococcaceae_family,control=list(msMaxIter=2000,opt="optim",msVerbose=TRUE))
summary(model_streptococcaceae)
plot(model_streptococcaceae, resid(.)~Time_continuous)
qq <- qqnorm(resid(model_streptococcaceae), main = "QQ Plot with Sample Labels Model", xlab = "Theoretical Quantiles", ylab = "Sample Quantiles")
qqline(resid(model_streptococcaceae), col = "red")  # Add reference line

# Plot des effets fixes
pdata <- expand.grid(Time_continuous=c(0,6,30), Condition=c('Untreated', 'Treated'))
pdata$Abundance <- predict(model_streptococcaceae, pdata, level=0)
plot(pdata$Time_continuous, pdata$Abundance, type='n', xlab='Time', ylab=' Relative Abundance')
points(pdata$Time_continuous[1:3], pdata$Abundance[1:3], type='b', pch=19, lwd=2)
points(pdata$Time_continuous[4:6], pdata$Abundance[4:6], type='b', pch=22, lwd=2,col="hotpink2") ## For Treated Calves 

## GENUS
## Choose 
# Access a specific family by name
anaerostipes_genus <- genus_list[["Anaerostipes"]]
blautia_genus <- genus_list[["Blautia"]]
clostridium_genus <- genus_list[["Clostridium"]]
coprococcus_genus <- genus_list[["Coprococcus"]]
enterococcus_genus <- genus_list[["Enterococcus"]]
escherichia_genus <- genus_list[["Escherichia"]]
faecalibacterium_genus <- genus_list[["Faecalibacterium"]]
lachnoclostridium_genus <- genus_list[["Lachnoclostridium"]]
lactobacillus_genus <- genus_list[["Lactobacillus"]]
ligilactobacillus_genus <- genus_list[["Ligilactobacillus"]]
mediterraneibacter_genus <- genus_list[["Mediterraneibacter"]]
megamonas_genus <- genus_list[["Megamonas"]]
romboutsia_genus <- genus_list[["Romboutsia"]]
ruminococcus_genus <- genus_list[["Ruminococcus"]]
streptococcus_genus <- genus_list[["Streptococcus"]]

#Anaerostipes
model_anaerostipes <- lme(Abundance ~ Time_continuous*Condition+ I(Time_continuous^2)*Condition, random = ~1+Time_continuous|Calves, weights = varIdent(~1|Condition), data = anaerostipes_genus,control=list(msMaxIter=2000,opt="optim",msVerbose=TRUE))
summary(model_anaerostipes)
plot(model_anaerostipes, resid(.)~Time_continuous)
qq <- qqnorm(resid(model_anaerostipes), main = "QQ Plot with Sample Labels Model", xlab = "Theoretical Quantiles", ylab = "Sample Quantiles")
qqline(resid(model_anaerostipes), col = "red")  # Add reference line

# Plot des effets fixes
pdata <- expand.grid(Time_continuous=c(0,6,30), Condition=c('Untreated', 'Treated'))
pdata$Abundance <- predict(model_anaerostipes, pdata, level=0)
plot(pdata$Time_continuous, pdata$Abundance, type='n', xlab='Time', ylab=' Relative Abundance')
points(pdata$Time_continuous[1:3], pdata$Abundance[1:3], type='b', pch=19, lwd=2)
points(pdata$Time_continuous[4:6], pdata$Abundance[4:6], type='b', pch=22, lwd=2,col="hotpink2") ## For Treated Calves 

#Blautia
model_blautia <- lme(Abundance ~ Time_continuous*Condition+ I(Time_continuous^2)*Condition, random = ~1+Time_continuous|Calves, weights = varIdent(~1|Condition), data = blautia_genus,control=list(msMaxIter=4000,opt="optim",msVerbose=TRUE))
summary(model_blautia)
plot(model_blautia, resid(.)~Time_continuous)
qq <- qqnorm(resid(model_blautia), main = "QQ Plot with Sample Labels Model", xlab = "Theoretical Quantiles", ylab = "Sample Quantiles")
qqline(resid(model_blautia), col = "red")  # Add reference line

# Plot des effets fixes
pdata <- expand.grid(Time_continuous=c(0,6,30), Condition=c('Untreated', 'Treated'))
pdata$Abundance <- predict(model_blautia, pdata, level=0)
plot(pdata$Time_continuous, pdata$Abundance, type='n', xlab='Time', ylab=' Relative Abundance')
points(pdata$Time_continuous[1:3], pdata$Abundance[1:3], type='b', pch=19, lwd=2)
points(pdata$Time_continuous[4:6], pdata$Abundance[4:6], type='b', pch=22, lwd=2,col="hotpink2") ## For Treated Calves 

#Clostridium
model_clostridium <- lme(Abundance ~ Time_continuous*Condition+ I(Time_continuous^2)*Condition, random = ~1+Time_continuous|Calves, weights = varIdent(~1|Condition), data = clostridium_genus,control=list(msMaxIter=4000,opt="optim",msVerbose=TRUE))
summary(model_clostridium)
plot(model_clostridium, resid(.)~Time_continuous)
qq <- qqnorm(resid(model_clostridium), main = "QQ Plot with Sample Labels Model", xlab = "Theoretical Quantiles", ylab = "Sample Quantiles")
qqline(resid(model_clostridium), col = "red")  # Add reference line

# Plot des effets fixes
pdata <- expand.grid(Time_continuous=c(0,6,30), Condition=c('Untreated', 'Treated'))
pdata$Abundance <- predict(model_clostridium, pdata, level=0)
plot(pdata$Time_continuous, pdata$Abundance, type='n', xlab='Time', ylab=' Relative Abundance')
points(pdata$Time_continuous[1:3], pdata$Abundance[1:3], type='b', pch=19, lwd=2)
points(pdata$Time_continuous[4:6], pdata$Abundance[4:6], type='b', pch=22, lwd=2,col="hotpink2") ## For Treated Calves 

#Coprococcus
model_coprococcus<- lme(Abundance ~ Time_continuous*Condition+ I(Time_continuous^2)*Condition, random = ~1+Time_continuous|Calves, weights = varIdent(~1|Condition), data = coprococcus_genus,control=list(msMaxIter=4000,opt="optim",msVerbose=TRUE))
summary(model_coprococcus)
plot(model_coprococcus, resid(.)~Time_continuous)
qq <- qqnorm(resid(model_coprococcus), main = "QQ Plot with Sample Labels Model", xlab = "Theoretical Quantiles", ylab = "Sample Quantiles")
qqline(resid(model_coprococcus), col = "red")  # Add reference line

# Plot des effets fixes
pdata <- expand.grid(Time_continuous=c(0,6,30), Condition=c('Untreated', 'Treated'))
pdata$Abundance <- predict(model_coprococcus, pdata, level=0)
plot(pdata$Time_continuous, pdata$Abundance, type='n', xlab='Time', ylab=' Relative Abundance')
points(pdata$Time_continuous[1:3], pdata$Abundance[1:3], type='b', pch=19, lwd=2)
points(pdata$Time_continuous[4:6], pdata$Abundance[4:6], type='b', pch=22, lwd=2,col="hotpink2") ## For Treated Calves 

#Enterococcus
model_enterococcus <- lme(Abundance ~ Time_continuous*Condition+ I(Time_continuous^2)*Condition, random = ~1+Time_continuous|Calves, weights = varIdent(~1|Condition), data = enterococcus_genus,control=list(msMaxIter=4000,opt="optim",msVerbose=TRUE))
summary(model_enterococcus)
plot(model_enterococcus, resid(.)~Time_continuous)
qq <- qqnorm(resid(model_enterococcus), main = "QQ Plot with Sample Labels Model", xlab = "Theoretical Quantiles", ylab = "Sample Quantiles")
qqline(resid(model_enterococcus), col = "red")  # Add reference line

# Plot des effets fixes
pdata <- expand.grid(Time_continuous=c(0,6,30), Condition=c('Untreated', 'Treated'))
pdata$Abundance <- predict(model_enterococcus, pdata, level=0)
plot(pdata$Time_continuous, pdata$Abundance, type='n', xlab='Time', ylab=' Relative Abundance')
points(pdata$Time_continuous[1:3], pdata$Abundance[1:3], type='b', pch=19, lwd=2)
points(pdata$Time_continuous[4:6], pdata$Abundance[4:6], type='b', pch=22, lwd=2,col="hotpink2") ## For Treated Calves 

#Escherichia
model_escherichia <- lme(Abundance ~ Time_continuous*Condition+ I(Time_continuous^2)*Condition, random = ~1+Time_continuous|Calves, weights = varIdent(~1|Condition), data = escherichia_genus,control=list(msMaxIter=4000,opt="optim",msVerbose=TRUE))
summary(model_escherichia)
plot(model_escherichia, resid(.)~Time_continuous)
qq <- qqnorm(resid(model_escherichia), main = "QQ Plot with Sample Labels Model", xlab = "Theoretical Quantiles", ylab = "Sample Quantiles")
qqline(resid(model_escherichia), col = "red")  # Add reference line

# Plot des effets fixes
pdata <- expand.grid(Time_continuous=c(0,6,30), Condition=c('Untreated', 'Treated'))
pdata$Abundance <- predict(model_escherichia, pdata, level=0)
plot(pdata$Time_continuous, pdata$Abundance, type='n', xlab='Time', ylab=' Relative Abundance')
points(pdata$Time_continuous[1:3], pdata$Abundance[1:3], type='b', pch=19, lwd=2)
points(pdata$Time_continuous[4:6], pdata$Abundance[4:6], type='b', pch=22, lwd=2,col="hotpink2") ## For Treated Calves 

#Faecalibacterium
model_faecalibacterium <- lme(Abundance ~ Time_continuous*Condition+ I(Time_continuous^2)*Condition, random = ~1+Time_continuous|Calves, weights = varIdent(~1|Condition), data = faecalibacterium_genus,control=list(msMaxIter=4000,opt="optim",msVerbose=TRUE))
summary(model_faecalibacterium)
plot(model_faecalibacterium, resid(.)~Time_continuous)
qq <- qqnorm(resid(model_faecalibacterium), main = "QQ Plot with Sample Labels Model", xlab = "Theoretical Quantiles", ylab = "Sample Quantiles")
qqline(resid(model_faecalibacterium), col = "red")  # Add reference line

# Plot des effets fixes
pdata <- expand.grid(Time_continuous=c(0,6,30), Condition=c('Untreated', 'Treated'))
pdata$Abundance <- predict(model_faecalibacterium, pdata, level=0)
plot(pdata$Time_continuous, pdata$Abundance, type='n', xlab='Time', ylab=' Relative Abundance')
points(pdata$Time_continuous[1:3], pdata$Abundance[1:3], type='b', pch=19, lwd=2)
points(pdata$Time_continuous[4:6], pdata$Abundance[4:6], type='b', pch=22, lwd=2,col="hotpink2") ## For Treated Calves 

#Lachnoclostridium
model_lachnoclostridium <- lme(Abundance ~ Time_continuous*Condition+ I(Time_continuous^2)*Condition, random = ~1+Time_continuous|Calves, weights = varIdent(~1|Condition), data = lachnoclostridium_genus,control=list(msMaxIter=4000,opt="optim",msVerbose=TRUE))
summary(model_lachnoclostridium)
plot(model_lachnoclostridium, resid(.)~Time_continuous)
qq <- qqnorm(resid(model_lachnoclostridium), main = "QQ Plot with Sample Labels Model", xlab = "Theoretical Quantiles", ylab = "Sample Quantiles")
qqline(resid(model_lachnoclostridium), col = "red")  # Add reference line

# Plot des effets fixes
pdata <- expand.grid(Time_continuous=c(0,6,30), Condition=c('Untreated', 'Treated'))
pdata$Abundance <- predict(model_lachnoclostridium, pdata, level=0)
plot(pdata$Time_continuous, pdata$Abundance, type='n', xlab='Time', ylab=' Relative Abundance')
points(pdata$Time_continuous[1:3], pdata$Abundance[1:3], type='b', pch=19, lwd=2)
points(pdata$Time_continuous[4:6], pdata$Abundance[4:6], type='b', pch=22, lwd=2,col="hotpink2") ## For Treated Calves 


#Lactobacillus
model_lactobacillus <- lme(Abundance ~ Time_continuous*Condition+ I(Time_continuous^2)*Condition, random = ~1+Time_continuous|Calves, weights = varIdent(~1|Condition), data = lactobacillus_genus,control=list(msMaxIter=4000,opt="optim",msVerbose=TRUE))
summary(model_lactobacillus)
plot(model_lactobacillus, resid(.)~Time_continuous)
qq <- qqnorm(resid(model_lactobacillus), main = "QQ Plot with Sample Labels Model", xlab = "Theoretical Quantiles", ylab = "Sample Quantiles")
qqline(resid(model_lactobacillus), col = "red")  # Add reference line

# Plot des effets fixes
pdata <- expand.grid(Time_continuous=c(0,6,30), Condition=c('Untreated', 'Treated'))
pdata$Abundance <- predict(model_lactobacillus, pdata, level=0)
plot(pdata$Time_continuous, pdata$Abundance, type='n', xlab='Time', ylab=' Relative Abundance')
points(pdata$Time_continuous[1:3], pdata$Abundance[1:3], type='b', pch=19, lwd=2)
points(pdata$Time_continuous[4:6], pdata$Abundance[4:6], type='b', pch=22, lwd=2,col="hotpink2") ## For Treated Calves 

#Ligilactobacillus
model_ligilactobacillus <- lme(Abundance ~ Time_continuous*Condition+ I(Time_continuous^2)*Condition, random = ~1+Time_continuous|Calves, weights = varIdent(~1|Condition), data = ligilactobacillus_genus,control=list(msMaxIter=4000,opt="optim",msVerbose=TRUE))
summary(model_ligilactobacillus)
plot(model_ligilactobacillus, resid(.)~Time_continuous)
qq <- qqnorm(resid(model_ligilactobacillus), main = "QQ Plot with Sample Labels Model", xlab = "Theoretical Quantiles", ylab = "Sample Quantiles")
qqline(resid(model_ligilactobacillus), col = "red")  # Add reference line

# Plot des effets fixes
pdata <- expand.grid(Time_continuous=c(0,6,30), Condition=c('Untreated', 'Treated'))
pdata$Abundance <- predict(model_ligilactobacillus, pdata, level=0)
plot(pdata$Time_continuous, pdata$Abundance, type='n', xlab='Time', ylab=' Relative Abundance')
points(pdata$Time_continuous[1:3], pdata$Abundance[1:3], type='b', pch=19, lwd=2)
points(pdata$Time_continuous[4:6], pdata$Abundance[4:6], type='b', pch=22, lwd=2,col="hotpink2") ## For Treated Calves 


#Mediterraneibacter
model_mediterraneibacter <- lme(Abundance ~ Time_continuous*Condition+ I(Time_continuous^2)*Condition, random = ~1+Time_continuous|Calves, weights = varIdent(~1|Condition), data = mediterraneibacter_genus,control=list(msMaxIter=4000,opt="optim",msVerbose=TRUE))
summary(model_mediterraneibacter)
plot(model_mediterraneibacter, resid(.)~Time_continuous)
qq <- qqnorm(resid(model_mediterraneibacter), main = "QQ Plot with Sample Labels Model", xlab = "Theoretical Quantiles", ylab = "Sample Quantiles")
qqline(resid(model_mediterraneibacter), col = "red")  # Add reference line

# Plot des effets fixes
pdata <- expand.grid(Time_continuous=c(0,6,30), Condition=c('Untreated', 'Treated'))
pdata$Abundance <- predict(model_mediterraneibacter, pdata, level=0)
plot(pdata$Time_continuous, pdata$Abundance, type='n', xlab='Time', ylab=' Relative Abundance')
points(pdata$Time_continuous[1:3], pdata$Abundance[1:3], type='b', pch=19, lwd=2)
points(pdata$Time_continuous[4:6], pdata$Abundance[4:6], type='b', pch=22, lwd=2,col="hotpink2") ## For Treated Calves 

#Megamonas
model_megamonas <- lme(Abundance ~ Time_continuous*Condition+ I(Time_continuous^2)*Condition, random = ~1+Time_continuous|Calves, weights = varIdent(~1|Condition), data = megamonas_genus,control=list(msMaxIter=4000,opt="optim",msVerbose=TRUE))
summary(model_megamonas)
plot(model_megamonas, resid(.)~Time_continuous)
qq <- qqnorm(resid(model_megamonas), main = "QQ Plot with Sample Labels Model", xlab = "Theoretical Quantiles", ylab = "Sample Quantiles")
qqline(resid(model_megamonas), col = "red")  # Add reference line

# Plot des effets fixes
pdata <- expand.grid(Time_continuous=c(0,6,30), Condition=c('Untreated', 'Treated'))
pdata$Abundance <- predict(model_megamonas, pdata, level=0)
plot(pdata$Time_continuous, pdata$Abundance, type='n', xlab='Time', ylab=' Relative Abundance')
points(pdata$Time_continuous[1:3], pdata$Abundance[1:3], type='b', pch=19, lwd=2)
points(pdata$Time_continuous[4:6], pdata$Abundance[4:6], type='b', pch=22, lwd=2,col="hotpink2") ## For Treated Calves 

#Romboutsia
model_romboutsia <- lme(Abundance ~ Time_continuous*Condition+ I(Time_continuous^2)*Condition, random = ~1+Time_continuous|Calves, weights = varIdent(~1|Condition), data = romboutsia_genus,control=list(msMaxIter=4000,opt="optim",msVerbose=TRUE))
summary(model_romboutsia)
plot(model_romboutsia, resid(.)~Time_continuous)
qq <- qqnorm(resid(model_romboutsia), main = "QQ Plot with Sample Labels Model", xlab = "Theoretical Quantiles", ylab = "Sample Quantiles")
qqline(resid(model_romboutsia), col = "red")  # Add reference line

# Plot des effets fixes
pdata <- expand.grid(Time_continuous=c(0,6,30), Condition=c('Untreated', 'Treated'))
pdata$Abundance <- predict(model_romboutsia, pdata, level=0)
plot(pdata$Time_continuous, pdata$Abundance, type='n', xlab='Time', ylab=' Relative Abundance')
points(pdata$Time_continuous[1:3], pdata$Abundance[1:3], type='b', pch=19, lwd=2)
points(pdata$Time_continuous[4:6], pdata$Abundance[4:6], type='b', pch=22, lwd=2,col="hotpink2") ## For Treated Calves 

#Ruminococcus
model_ruminococcus <- lme(Abundance ~ Time_continuous*Condition+ I(Time_continuous^2)*Condition, random = ~1+Time_continuous|Calves, weights = varIdent(~1|Condition), data = ruminococcus_genus,control=list(msMaxIter=4000,opt="optim",msVerbose=TRUE))
summary(model_ruminococcus)
plot(model_ruminococcus, resid(.)~Time_continuous)
qq <- qqnorm(resid(model_ruminococcus), main = "QQ Plot with Sample Labels Model", xlab = "Theoretical Quantiles", ylab = "Sample Quantiles")
qqline(resid(model_ruminococcus), col = "red")  # Add reference line

# Plot des effets fixes
pdata <- expand.grid(Time_continuous=c(0,6,30), Condition=c('Untreated', 'Treated'))
pdata$Abundance <- predict(model_ruminococcus, pdata, level=0)
plot(pdata$Time_continuous, pdata$Abundance, type='n', xlab='Time', ylab=' Relative Abundance')
points(pdata$Time_continuous[1:3], pdata$Abundance[1:3], type='b', pch=19, lwd=2)
points(pdata$Time_continuous[4:6], pdata$Abundance[4:6], type='b', pch=22, lwd=2,col="hotpink2") ## For Treated Calves 

#Streptococcus
model_streptococcus <- lme(Abundance ~ Time_continuous*Condition+ I(Time_continuous^2)*Condition, random = ~1+Time_continuous|Calves, weights = varIdent(~1|Condition), data = streptococcus_genus,control=list(msMaxIter=4000,opt="optim",msVerbose=TRUE))
summary(model_streptococcus)
plot(model_streptococcus, resid(.)~Time_continuous)
qq <- qqnorm(resid(model_streptococcus), main = "QQ Plot with Sample Labels Model", xlab = "Theoretical Quantiles", ylab = "Sample Quantiles")
qqline(resid(model_streptococcus), col = "red")  # Add reference line

# Plot des effets fixes
pdata <- expand.grid(Time_continuous=c(0,6,30), Condition=c('Untreated', 'Treated'))
pdata$Abundance <- predict(model_streptococcus, pdata, level=0)
plot(pdata$Time_continuous, pdata$Abundance, type='n', xlab='Time', ylab=' Relative Abundance')
points(pdata$Time_continuous[1:3], pdata$Abundance[1:3], type='b', pch=19, lwd=2)
points(pdata$Time_continuous[4:6], pdata$Abundance[4:6], type='b', pch=22, lwd=2,col="hotpink2") ## For Treated Calves 







#####################
## Modèle linéaire ##
#####################



library(phyloseq)
library(lme4)
library(lmerTest)
library(nlme)

install.packages("Matrix", type = "binary")
install.packages("lme4", type = "binary")

physeq_phylum <- tax_glom(physeq_rarefy, taxrank = "Phylum")

# Melt phyloseq object
phylum <- psmelt(physeq_phylum)

# Ensure variables are factors
phylum$Condition <- factor(phylum$Condition, levels = c("Untreated", "Treated"))
phylum$Time <- factor(phylum$Time, levels = c("T0", "T1", "T2", "T3"))
phylum$Calves <- as_factor(phylum$Calves)

# Split data by Phylum
phylum_list <- split(phylum, phylum$Phylum)

names(phylum_list)

## Choose "Bacillota", "Bacteroidota" "Actinomycetota" "Pseudomonadota"
# Access a specific phylum by name
bacillota_phylum <- phylum_list[["Bacillota"]]
bacteroidota_phylum <- phylum_list[["Bacteroidota"]]
actinomycetota_phylum <- phylum_list[["Actinomycetota"]]
pseudomonadota_phylum <- phylum_list[["Pseudomonadota"]]

# Fit the model
#model_bacillota <- lmer(Abundance ~ Condition:Time + (1|Calves), data = bacillota_phylum)
model_2 <- lme(fixed = Abundance ~ Condition*Time, random = ~1|Calves, data = bacillota_phylum)

# View the summary
summary(model_bacillot)
summary(model_2)

# Residual diagnostics
plot(model_2)  # Residuals vs fitted values
qqnorm(resid(model_2))  # Q-Q plot for residuals
qqline(resid(model_2))
anova(model_2)

library(emmeans)
posthoc <- emmeans(model_2, ~ Condition * Time)
pairs(posthoc)  # Pairwise comparisons

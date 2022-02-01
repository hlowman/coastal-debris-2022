# Coastal Fund Project Linear Mixed Effects Models
# January 7, 2022
# Heili Lowman

# The following script will create a series of linear mixed effect models with the newly
# curated dataset combining both lignin and pyC measures.

#### Setup and Load-in ####

# Load packages
library(tidyverse)
library(here)
library(lubridate)
library(calecopal)
library(patchwork)
library(naniar)
library(nlme)
library(multcomp)
library(emmeans)
library(GGally)
library(groupdata2)
library(MuMIn)

# Load data
dat <- read_csv("data_raw/lignin_pyc_results.csv")

#### Data Prep ####

# Tidy data
dat_ed <- dat %>%
  mutate(Date = mdy(Collection_Date)) %>% # properly format date
  mutate(Environment_f = factor(case_when(Location %in% c("Montecito Debris",
                                                          "Goleta Beach",
                                                          "Goleta Slough",
                                                          "Debris Deposition Site",
                                                          "Clay Debris") ~ "Beach",
                                          Location %in% c("Goleta Bay West",
                                                          "Goleta Bay East") ~ "Ocean"))) # create environment column

# Create beach dataset for use in the models below

# List to combine MD and DD sampling sites into a single name
list1 <- c("MD", "DD")

# Models for terrestrial sampling will be able to examine the effect of time
dat_beach <- dat_ed %>%
  filter(Environment_f == "Beach") %>%
  filter(Location != "Clay Debris") %>% # Removing the clay debris layer, because it was only sampled once
  mutate(Location_f = factor(case_when(Location_id %in% list1 ~ "DD",
                                       TRUE ~ Location_id))) %>% # level locations
  mutate(Date_f = factor(case_when(Date == "2018-02-02" ~ "A",
                            Date == "2018-02-23" ~ "B",
                            Date == "2018-04-24" ~ "C"))) %>% # add new column for date as factor
  mutate(Rep_f = factor(Replicate_id)) # creating column of replicate as a factor

# Create the marine dataset for use in the models below
# Models for marine sampling will NOT be able to examine the effect of time
dat_marine <- dat_ed %>%
  filter(Environment_f == "Ocean")

#### LMEM #1: Beach Lambda ####

# model structure:
# Lambda ~ Date + 1 | Site / Replicate

# Upon reflection, and running the below model including "Site" as a fixed effect,
# I've removed it and used it only as a random effect, since the Goleta Slough had
# disappearing Lambda values towards the end, which yielded p = NaN in the final results
# for "Site".

# Examine data
plot(Lambda ~ Date_f, data = dat_beach) # yeah, there's definitely a difference here
plot(Lambda ~ Location_f, data = dat_beach) # same
plot(Lambda ~ Rep_f, data = dat_beach) # more overlap
hist(dat_beach$Lambda)

dat_beach <- dat_beach %>%
  mutate(logLambda = log10(Lambda))

hist(dat_beach$logLambda) # much better

##STEP 1: Create a linear regression and check residuals.

# Create initial linear model with only fixed effects.
a1 <- lm(logLambda ~ Date_f, data = dat_beach) 

# Assigns standardized residuals to ra1.
ra1 <- data.frame(rstandard(a1))

# Because there's missing data, here's a workaround to compare residuals to data.
ra1 <- ra1 %>%
  rownames_to_column("record")
ra1_ed <- ra1 %>%
  mutate(rn = as.numeric(record))
db1 <- dat_beach %>%
  mutate(RECORD = seq(1,18)) %>%
  left_join(ra1_ed, by = c("RECORD" = "rn"))

# Plot said residuals.
ggplot(data = db1, aes(x = Date_f, y = rstandard.a1.)) +
  geom_point()

ggplot(data = db1, aes(x = Location_f, y = rstandard.a1.)) +
  geom_point() # doesn't appear to be a significant trend in either

##STEP 2: Fit the lm() with GLS and compare to lme().

# We cannot proceed with missing data in the dataset
db_Lambda_rmna <- dat_beach %>%
  drop_na(logLambda)

a2 <- gls(logLambda ~ Date_f, data = db_Lambda_rmna) # Effectively a linear regression
a3 <- lme(logLambda ~ Date_f, random =~1 | Location_f/Rep_f, data = db_Lambda_rmna) # Creates the first LMEM with random term.
anova(a2, a3) # Compares the two models. a2 preferred with AIC value of 28.51, but I am going to keep
# the random effect term in due to repeated sampling.

##STEP 3: Decide on a variance structure (aka random terms).

plot(a3, col=1) # Check the residuals before jumping right to applying a variance transformation.

qqnorm(a3) # This actually looks pretty good considering how few samples there are - so I'm going to skip the added variance structure. Keep in mind, log-transforming from the start is more powerful than adding in a variance component on the back end.

##STEP 4: Fit the lme().

# Using a3 <- lme(logLambda ~ Date_f, random =~1 | Location_f/Rep_f, data = db_Lambda_rmna)

##STEP 5: Compare the lm() and lme().

# See Step 2.

##STEP 6: Everything ok? Check residuals.

# See Step 3.

##STEP 7/8: Step-wise Optimal Fixed Structure

# See Step 4.

##STEP 9: Refit with REML

afinal <- lme(logLambda ~ Date_f, 
              random =~1 | Location_f/Rep_f,
              method = "REML", 
              data = db_Lambda_rmna)

# Output of the model.
# Note, summary() function looks at contrasts between singular effects.
summary(afinal)

# Checking residuals.
plot(afinal, col=1) # No pattern.
qqnorm(afinal) # Looks pretty good.

# Final results.
# anova() function looks at contrasts across all effects.
anova(afinal)

##STEP 10: What does this mean in WORDS?

# I applied a linear mixed effect modelling approach because the data were collected at each site (3) on multiple occasions (3). My model suggests there is a significant effect of date on log(Lambda) values of beach sediment samples. Random intercepts by site and replicate core were added.

# Equation: log(Lambda) = -0.32 - 0.42[Feb23] - 1.12[Apr24] + random
# Fixed Effect (Date): F(2, 5) = 6.90, p = 0.04

# Post-hoc:
aHSD <- glht(afinal, linfct=mcp(Date_f="Tukey")) # Run a Tukey's post hoc analysis on Date factor.
summary(aHSD) # February 2 & April 24 significantly different from one another (p < 0.001).

#### LMEM #2: Beach PyC ####

# model structure:
# Lambda ~ Date + 1 | Site / Replicate

# Similar to model above, using Site as a fixed effect yielded NaNs due to low sample size/power, so I've removed it
# from the model below, and focused instead on changes through time.

# Examine data
plot(pyC ~ Date_f, data = dat_beach) # yeah, there's definitely a difference here
plot(pyC ~ Location_f, data = dat_beach) # look similar to Date
plot(pyC ~ Rep_f, data = dat_beach) # serious overlap
hist(dat_beach$pyC) # high in the initial debris, so log-transforming also

dat_beach <- dat_beach %>%
  mutate(logpyC = log10(pyC))

hist(dat_beach$logpyC) # somewhat better

##STEP 1: Create a linear regression and check residuals.

# Create initial linear model with only fixed effects.
b1 <- lm(logpyC ~ Date_f, data = dat_beach) 

# Assigns standardized residuals to ra1.
rb1 <- data.frame(rstandard(b1))

db2 <- dat_beach %>%
  mutate(rb1 = rb1$rstandard.b1.)

# Plot said residuals.
ggplot(data = db2, aes(x = Date_f, y = rb1)) +
  geom_point()

##STEP 2: Fit the lm() with GLS and compare to lme().

# No missing pyC values, so proceeding as is.
b2 <- gls(logpyC ~ Date_f, data = dat_beach) # Effectively a linear regression
b3 <- lme(logpyC ~ Date_f, random =~1 | Location_f/Rep_f, data = dat_beach) # Creates the first LMEM with random term.
anova(b2, b3) # Compares the two models. b2 preferred with AIC value of 27.95, but I am going to keep
# the random effect term in due to repeated sampling.

##STEP 3: Decide on a variance structure (aka random terms).

plot(b3, col=1) # Check the residuals before jumping right to applying a variance transformation.

qqnorm(b3) # These aren't looking great, so I'm going to add a variance structure by Location,
# based on plots above.

b4 <- lme(logpyC ~ Date_f, 
          random =~1 | Location_f/Rep_f, 
          data = dat_beach,
          weights = varIdent(form =~1 | Location_f))

anova(b3, b4) # Compares the two models. b4 preferred with AIC value of 9.79.

plot(b4, col = 1)
qqnorm(b4) # Both plots look WAY better.

##STEP 4: Fit the lme().

# Using b4 <- lme(logpyC ~ Date_f + Location_f, random =~1 | Location_f/Rep_f, data = dat_beach, weights = varIdent(form =~1 | Location_f))

##STEP 5: Compare the lm() and lme().

anova(b2, b4) # Compares the two models. b4 preferred with AIC value of 9.79.

##STEP 6: Everything ok? Check residuals.

# See Step 3.

##STEP 7/8: Step-wise Optimal Fixed Structure

# See Step 4.

# STEP 9: Refit with REML

bfinal <- lme(logpyC ~ Date_f, 
              random =~1 | Location_f/Rep_f,
              method = "REML", 
              data = dat_beach,
              weights = varIdent(form =~1 | Location_f))

# Output of the model.
# Note, summary() function looks at contrasts between singular effects.
summary(bfinal)

# Checking residuals.
plot(bfinal, col=1) 
qqnorm(bfinal) # Looks same as above.

# Final results.
# anova() function looks at contrasts across all effects.
anova(bfinal)

##STEP 10: What does this mean in WORDS?

# My model suggests there is no significant effect of date on log(pyC) values of beach sediment samples. Random intercepts by site and replicate core were added as well as a variance term by sampling site.

# Equation: log(pyC) = -0.42 + 0.03[Feb23] + 0.04[Apr24] + random + var(Site)
# Fixed Effect (Date): F(2, 10) = 0.54, p = 0.60

#### LMs: Marine Analytes ####

# Running multiple linear regressions since no repeated sampling took place for the following analytes in marine sediment:
# Lambda
# PyC
# S/V
# C/V
# P/V+S
# 3,5 Bd/V

# comparisons by water depth - 5, 10, 20m; site - GOLBW, GOLBE; core section - 0-10, 10-20

## Lambda
hist(dat_marine$Lambda)

dat_marine <- dat_marine %>%
  mutate(logLambda = log10(Lambda)) %>% # log-transform Lambda values, as we did above
  mutate(Location_f = factor(Location_id)) %>% # make site a factor
  mutate(Water_Depth_f = factor(Water_Depth)) %>% # make depth a factor for tukey's post hoc use
  mutate(Core_Section_f = factor(Core_Section)) # make core section a factor as well

hist(dat_marine$logLambda) # better

lm1 <- lm(logLambda ~ Location_f + Water_Depth_f + Core_Section_f, data = dat_marine)

summary(lm1)
# Water_Depth_f20: p = 0.004

anova(lm1)
# Water_Depth_f: p = 0.01

# tukey's post-hoc
# null hypothesis for each of the pair-wise comparisons is still no difference in means
posthoc1 <- glht(lm1, linfct = mcp(Water_Depth_f = 'Tukey'))
summary(posthoc1)
# 20 - 5: p = 0.01

## PyC
hist(dat_marine$pyC)

dat_marine <- dat_marine %>%
  mutate(logpyC = log10(pyC)) # log-transform pyC values, as we did above

hist(dat_marine$logpyC) # better

lm2 <- lm(logpyC ~ Location_f + Water_Depth_f + Core_Section_f, data = dat_marine)

summary(lm2)
# Location_fGOBW: p = 0.0004
# Water_Depth_f20: p = 0.002

anova(lm2)
# Location_f: p = 0.003
# Water_Depth_f: p = 0.005

# tukey's post-hoc
# null hypothesis for each of the pair-wise comparisons is still no difference in means
posthoc2 <- glht(lm2, linfct = mcp(Location_f = 'Tukey'))
summary(posthoc2)
# GOBE - GOBW: p = 0.0004

posthoc3 <- glht(lm2, linfct = mcp(Water_Depth_f = 'Tukey'))
summary(posthoc3)
# 20 - 5: p = 0.006

## S/V

hist(dat_marine$SV)

lm3 <- lm(SV ~ Location_f + Water_Depth_f + Core_Section_f, data = dat_marine)

summary(lm3)

anova(lm3)
# None are significant

## C/V

hist(dat_marine$CV)

dat_marine <- dat_marine %>%
  mutate(logCV = log10(CV)) # log-transform CV values

hist(dat_marine$logCV) # yes, better

lm4 <- lm(logCV ~ Location_f + Water_Depth_f + Core_Section_f, data = dat_marine)

summary(lm4)

anova(lm4)
# None are significant

## P/V+S

hist(dat_marine$PVS)

dat_marine <- dat_marine %>%
  mutate(logPVS = log10(PVS)) # log-transform PVS values

hist(dat_marine$logPVS) # eh, i'm going to stick with the raw values

lm5 <- lm(logPVS ~ Location_f + Water_Depth_f + Core_Section_f, data = dat_marine)

summary(lm5)

anova(lm5)
# None are significant

## 3,5 Bd/V

hist(dat_marine$BdV)

dat_marine <- dat_marine %>%
  mutate(logBdV = log10(BdV)) # log-transform BdV values

hist(dat_marine$logBdV) # also better

lm6 <- lm(logBdV ~ Location_f + Water_Depth_f + Core_Section_f, data = dat_marine)

summary(lm6)
# Water_Depth_f10: p = 0.02
# Water_Depth_f20: p = 0. 002

anova(lm6)
# Water_Depth_f: p = 0.005

# tukey's post-hoc
# null hypothesis for each of the pair-wise comparisons is still no difference in means
posthoc4 <- glht(lm6, linfct = mcp(Water_Depth_f = 'Tukey'))
summary(posthoc4)
# 10 - 5: p = 0.049
# 20 - 5: p = 0.004

# End of script.
  
# Coastal Fund Project Linear Mixed Effects Models
# January 7, 2022
# Heili Lowman

# The following script will create a series of linear mixed effect modesl with the newly
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

# Similar to model above, using Site as a fixed effect yielded NaNs, so I've removed it
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

#### ANOVAs: Marine Analytes ####

# Running ANOVAs since no repeated sampling took place for the following analytes in marine sediment:
# Lambda
# PyC
# S/V
# C/V
# P/V+S
# 3,5 Bd/V

# anova by water depth - 5, 10, 20m

## Lambda
hist(dat_marine$Lambda)

dat_marine <- dat_marine %>%
  mutate(logLambda = log10(Lambda)) %>% # log-transform Lambda values, as we did above
  mutate(Water_Depth_f = factor(Water_Depth)) # make depth a factor for tukey's post hoc use

hist(dat_marine$logLambda) # better

# bartlett's test
# testing null hypothesis that variances across all groups are equal
vartest1 <- bartlett.test(logLambda ~ Water_Depth_f, data = dat_marine) # p = 0.31, so assumption met

aov1 <- aov(logLambda ~ Water_Depth_f, data = dat_marine)
summary(aov1) # p = 0.02

# tukey's post-hoc
# null hypothesis for each of the pair-wise comparisons is still no difference in means
posthoc1 <- TukeyHSD(aov1)
posthoc1 # 20-5 p = 0.02

## PyC
hist(dat_marine$pyC)

dat_marine <- dat_marine %>%
  mutate(logpyC = log10(pyC)) # log-transform pyC values, as we did above

hist(dat_marine$logpyC) # better

# bartlett's test
# testing null hypothesis that variances across all groups are equal
vartest2 <- bartlett.test(logpyC ~ Water_Depth_f, data = dat_marine) # p = 0.85, so assumption met

aov2 <- aov(logpyC ~ Water_Depth_f, data = dat_marine)
summary(aov2) # p = 0.2

## S/V

hist(dat_marine$SV)

# bartlett's test
# testing null hypothesis that variances across all groups are equal
vartest3 <- bartlett.test(SV ~ Water_Depth_f, data = dat_marine) # p = 0.42, so assumption met

aov3 <- aov(SV ~ Water_Depth_f, data = dat_marine)
summary(aov3) # p = 0.17

## C/V

hist(dat_marine$CV)

dat_marine <- dat_marine %>%
  mutate(logCV = log10(CV)) # log-transform CV values

hist(dat_marine$logCV) # yes, better

# bartlett's test
# testing null hypothesis that variances across all groups are equal
vartest4 <- bartlett.test(logCV ~ Water_Depth_f, data = dat_marine) # p = 0.05, so assumption met

aov4 <- aov(logCV ~ Water_Depth_f, data = dat_marine)
summary(aov4) # p = 0.81

## P/V+S

hist(dat_marine$PVS)

dat_marine <- dat_marine %>%
  mutate(logPVS = log10(PVS)) # log-transform PVS values

hist(dat_marine$logPVS) # eh, i'm going to stick with the raw values

# bartlett's test
# testing null hypothesis that variances across all groups are equal
vartest5 <- bartlett.test(PVS ~ Water_Depth_f, data = dat_marine) # p = 0.52, so assumption met

aov5 <- aov(PVS ~ Water_Depth_f, data = dat_marine)
summary(aov5) # p = 0.26

## 3,5 Bd/V

hist(dat_marine$BdV)

dat_marine <- dat_marine %>%
  mutate(logBdV = log10(BdV)) # log-transform BdV values

hist(dat_marine$logBdV) # also better

# bartlett's test
# testing null hypothesis that variances across all groups are equal
vartest6 <- bartlett.test(logBdV ~ Water_Depth_f, data = dat_marine) # p = 0.24, so assumption met

aov6 <- aov(logBdV ~ Water_Depth_f, data = dat_marine)
summary(aov6) # p = 0.006

# tukey's post-hoc
# null hypothesis for each of the pair-wise comparisons is still no difference in means
posthoc6 <- TukeyHSD(aov6)
posthoc6 # 20-5 p = 0.004

#### T-tests: Marine Analytes ####

# Running t-tests since no repeated sampling took place for the following analytes in marine sediment:
# Lambda
# PyC
# S/V
# C/V
# P/V+S
# 3,5 Bd/V

# t-test by site - GOLBW, GOLBE
# t-test by core section - 0-10, 10-20

dat_marine <- dat_marine %>%
  mutate(Core_Section_f = factor(Core_Section))

## Lambda

hist(dat_marine$logLambda)
plot(dat_marine$Location_f, dat_marine$logLambda)
plot(dat_marine$Core_Section_f, dat_marine$logLambda)

vt1 <- dat_marine %>%
  dplyr::select(UF_Sample_id, Location_f, logLambda) %>%
  pivot_wider(names_from = Location_f, values_from = logLambda)

# testing null hypothesis that variances across all groups are equal
var1 <- var.test(vt1$GOBW, vt1$GOBE) # p = 0.04, so null hypothesis disproven

# need to use Mann-Whitney U test instead, since assumptions were not met
wt1 <- wilcox.test(logLambda ~ Location_f, data = dat_marine)
wt1 # p = 0.94

vt2 <- dat_marine %>%
  dplyr::select(UF_Sample_id, Core_Section_f, logLambda) %>%
  pivot_wider(names_from = Core_Section_f, values_from = logLambda)

# testing null hypothesis that variances across all groups are equal
var2 <- var.test(vt2$`0-10`, vt2$`10-20`) # p = 0.69, so assumption met

tt2 <- t.test(logLambda ~ Core_Section, data = dat_marine)
tt2 # p = 0.75

## PyC

hist(dat_marine$logpyC)
plot(dat_marine$Location_f, dat_marine$logpyC)
plot(dat_marine$Core_Section_f, dat_marine$logpyC)

vt3 <- dat_marine %>%
  dplyr::select(UF_Sample_id, Location_f, logpyC) %>%
  pivot_wider(names_from = Location_f, values_from = logpyC)

# testing null hypothesis that variances across all groups are equal
var3 <- var.test(vt3$GOBW, vt3$GOBE) # p = 0.23, so assumption met

tt3 <- t.test(logpyC ~ Location_f, data = dat_marine)
tt3 # p = 0.009

vt4 <- dat_marine %>%
  dplyr::select(UF_Sample_id, Core_Section_f, logpyC) %>%
  pivot_wider(names_from = Core_Section_f, values_from = logpyC)

# testing null hypothesis that variances across all groups are equal
var4 <- var.test(vt4$`0-10`, vt4$`10-20`) # p = 0.29, so assumption met

tt4 <- t.test(logpyC ~ Core_Section, data = dat_marine)
tt4 # p = 0.33

## S/V

hist(dat_marine$SV)
plot(dat_marine$Location_f, dat_marine$SV)
plot(dat_marine$Core_Section_f, dat_marine$SV)

vt5 <- dat_marine %>%
  dplyr::select(UF_Sample_id, Location_f, SV) %>%
  pivot_wider(names_from = Location_f, values_from = SV)

# testing null hypothesis that variances across all groups are equal
var5 <- var.test(vt5$GOBW, vt5$GOBE) # p = 0.85, so assumption met

tt5 <- t.test(SV ~ Location_f, data = dat_marine)
tt5 # p = 0.92

vt6 <- dat_marine %>%
  dplyr::select(UF_Sample_id, Core_Section_f, SV) %>%
  pivot_wider(names_from = Core_Section_f, values_from = SV)

# testing null hypothesis that variances across all groups are equal
var6 <- var.test(vt6$`0-10`, vt6$`10-20`) # p = 0.77, so assumption met

tt6 <- t.test(SV ~ Core_Section, data = dat_marine)
tt6 # p = 0.95

## C/V

hist(dat_marine$logCV)
plot(dat_marine$Location_f, dat_marine$logCV)
plot(dat_marine$Core_Section_f, dat_marine$logCV)

vt7 <- dat_marine %>%
  dplyr::select(UF_Sample_id, Location_f, logCV) %>%
  pivot_wider(names_from = Location_f, values_from = logCV)

# testing null hypothesis that variances across all groups are equal
var7 <- var.test(vt7$GOBW, vt7$GOBE) # p = 0.14, so assumption met

tt7 <- t.test(logCV ~ Location_f, data = dat_marine)
tt7 # p = 0.36

vt8 <- dat_marine %>%
  dplyr::select(UF_Sample_id, Core_Section_f, logCV) %>%
  pivot_wider(names_from = Core_Section_f, values_from = logCV)

# testing null hypothesis that variances across all groups are equal
var8 <- var.test(vt8$`0-10`, vt8$`10-20`) # p = 0.97, so assumption met

tt8 <- t.test(logCV ~ Core_Section, data = dat_marine)
tt8 # p = 0.23

## P/V+S

hist(dat_marine$logPVS)
plot(dat_marine$Location_f, dat_marine$logPVS)
plot(dat_marine$Core_Section_f, dat_marine$logPVS)

vt9 <- dat_marine %>%
  dplyr::select(UF_Sample_id, Location_f, logPVS) %>%
  pivot_wider(names_from = Location_f, values_from = logPVS)

# testing null hypothesis that variances across all groups are equal
var9 <- var.test(vt9$GOBW, vt9$GOBE) # p = 0.14, so assumption met

tt9 <- t.test(logPVS ~ Location_f, data = dat_marine)
tt9 # p = 0.40

vt10 <- dat_marine %>%
  dplyr::select(UF_Sample_id, Core_Section_f, logPVS) %>%
  pivot_wider(names_from = Core_Section_f, values_from = logPVS)

# testing null hypothesis that variances across all groups are equal
var10 <- var.test(vt10$`0-10`, vt10$`10-20`) # p = 0.33, so assumption met

tt10 <- t.test(logPVS ~ Core_Section, data = dat_marine)
tt10 # p = 0.79

## 3,5Bd/V

hist(dat_marine$logBdV)
plot(dat_marine$Location_f, dat_marine$logBdV)
plot(dat_marine$Core_Section_f, dat_marine$logBdV)

vt11 <- dat_marine %>%
  dplyr::select(UF_Sample_id, Location_f, logBdV) %>%
  pivot_wider(names_from = Location_f, values_from = logBdV)

# testing null hypothesis that variances across all groups are equal
var11 <- var.test(vt11$GOBW, vt11$GOBE) # p = 0.10, so assumption met

tt11 <- t.test(logBdV ~ Location_f, data = dat_marine)
tt11 # p = 0.82

vt12 <- dat_marine %>%
  dplyr::select(UF_Sample_id, Core_Section_f, logBdV) %>%
  pivot_wider(names_from = Core_Section_f, values_from = logBdV)

# testing null hypothesis that variances across all groups are equal
var12 <- var.test(vt12$`0-10`, vt12$`10-20`) # p = 0.91, so assumption met

tt12 <- t.test(logBdV ~ Core_Section, data = dat_marine)
tt12 # p = 0.73

# End of script.
  
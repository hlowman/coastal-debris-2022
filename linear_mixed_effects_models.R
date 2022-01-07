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
# Lambda ~ Date + Location + 1 | Site / Replicate

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
a1 <- lm(logLambda ~ Date_f + Location_f, data = dat_beach) 

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

# After all that work above, it seems we cannot proceed with missing data in the dataset so...
db_Lambda_rmna <- dat_beach %>%
  drop_na(logLambda)

a2 <- gls(logLambda ~ Date_f + Location_f, data = db_Lambda_rmna) # Effectively a linear regression
a3 <- lme(logLambda ~ Date_f + Location_f, random =~1 | Location_f/Rep_f, data = db_Lambda_rmna) # Creates the first LMEM with random term.
anova(a2, a3) # Compares the two models. a2 preferred with AIC value of 26.06, but I am going to keep
# the random effect term in due to repeated sampling.

##STEP 3: Decide on a variance structure (aka random terms).

plot(a3, col=1) # Check the residuals before jumping right to applying a variance transformation.

qqnorm(a3) # This actually looks pretty good considering how few samples there are - so I'm going to skip the added variance structure. Keep in mind, log-transforming from the start is more powerful than adding in a variance component on the back end.

##STEP 4: Fit the lme().

# Using a3 <- lme(logLambda ~ Date_f + Location_f, random =~1 | Location_f/Rep_f, data = db_Lambda_rmna)

##STEP 5: Compare the lm() and lme().

# See Step 2.

##STEP 6: Everything ok? Check residuals.

# See Step 3.

##STEP 7/8: Step-wise Optimal Fixed Structure

a3_ml <- lme(logLambda ~ Date_f + Location_f, 
             random =~1 | Location_f/Rep_f,
             method = "ML", 
             data = db_Lambda_rmna) # Switch over to ML for fixed component editing portion.

a4 <- lme(logLambda ~ Location_f, 
          random =~1 | Location_f/Rep_f,
          method = "ML", 
          data = db_Lambda_rmna) # Remove "Date_f" as a fixed factor.

a5 <- lme(logLambda ~ Date_f, 
          random =~1 | Location_f/Rep_f,
          method = "ML", 
          data = db_Lambda_rmna) # Remove "Location_f" as a fixed factor.

anova(a3_ml, a4, a5) # Compare the three models. a3_ml preferred with AIC value of 23.18.

# STEP 9: Refit with REML

afinal <- lme(logLambda ~ Date_f + Location_f, 
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

#### LMEM #2: Beach PyC ####

# model structure:
# Lambda ~ Date + Site + 1 | Site / Replicate

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
a1 <- lm(logLambda ~ Date_f + Location_f, data = dat_beach) 

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

# After all that work above, it seems we cannot proceed with missing data in the dataset so...
db_Lambda_rmna <- dat_beach %>%
  drop_na(logLambda)

a2 <- gls(logLambda ~ Date_f + Location_f, data = db_Lambda_rmna) # Effectively a linear regression
a3 <- lme(logLambda ~ Date_f + Location_f, random =~1 | Location_f/Rep_f, data = db_Lambda_rmna) # Creates the first LMEM with random term.
anova(a2, a3) # Compares the two models. a2 preferred with AIC value of 26.06, but I am going to keep
# the random effect term in due to repeated sampling.

##STEP 3: Decide on a variance structure (aka random terms).

plot(a3, col=1) # Check the residuals before jumping right to applying a variance transformation.

qqnorm(a3) # This actually looks pretty good considering how few samples there are - so I'm going to skip the added variance structure. Keep in mind, log-transforming from the start is more powerful than adding in a variance component on the back end.

##STEP 4: Fit the lme().

# Using a3 <- lme(logLambda ~ Date_f + Location_f, random =~1 | Location_f/Rep_f, data = db_Lambda_rmna)

##STEP 5: Compare the lm() and lme().

# See Step 2.

##STEP 6: Everything ok? Check residuals.

# See Step 3.

##STEP 7/8: Step-wise Optimal Fixed Structure

a3_ml <- lme(logLambda ~ Date_f + Location_f, 
             random =~1 | Location_f/Rep_f,
             method = "ML", 
             data = db_Lambda_rmna) # Switch over to ML for fixed component editing portion.

a4 <- lme(logLambda ~ Location_f, 
          random =~1 | Location_f/Rep_f,
          method = "ML", 
          data = db_Lambda_rmna) # Remove "Date_f" as a fixed factor.

a5 <- lme(logLambda ~ Date_f, 
          random =~1 | Location_f/Rep_f,
          method = "ML", 
          data = db_Lambda_rmna) # Remove "Location_f" as a fixed factor.

anova(a3_ml, a4, a5) # Compare the three models. a3_ml preferred with AIC value of 23.18.

# STEP 9: Refit with REML

afinal <- lme(logLambda ~ Date_f + Location_f, 
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
  
# End of script.
  
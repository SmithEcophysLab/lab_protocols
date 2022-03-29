###############################################################################
## Load libraries
###############################################################################
library(dplyr)
library(tidyr)

###############################################################################
## Read in .csv with all plate reader data for ammonium assays
###############################################################################
df <- read.csv("NH4_NO3_data_example.csv")

###############################################################################
## Remove any incorrect or unreasonable blank values
###############################################################################
#df$abs.540[2] <- NA

###############################################################################
## Calculate mean blank values for each plate, then subtract mean blank values
## from abs.650 values to determine blank-corrected absorbance
###############################################################################
df.no3 <- df %>%
  filter(analysis.type == "no3") %>%
  select(plate, id, type, stand.conc, abs.540) %>%
  group_by(plate) %>%
  filter(id == "blank") %>%
  summarize(mean.blank = mean(abs.540, na.rm = TRUE)) %>%
  full_join(df) %>%
  filter(analysis.type == "no3" & id != "blank") %>%
  select(plate, id, type, stand.conc, abs.540, mean.blank) %>%
  mutate(abs.540.blank = as.numeric(abs.540) - as.numeric(mean.blank))
df.no3

###############################################################################
## Linear regression to calculate R^2 and linear model coefficients. R^2 should
## be greater than 0.98. Can remove points to improve model fit
###############################################################################
# Regression for plate 1
summary(lm(abs.540.blank ~ as.numeric(stand.conc),
           data = subset(df.no3, plate == 1)))

# Intercept and slope of regression equation for plate 1
coef(lm(abs.540.blank ~ as.numeric(stand.conc),
        data = subset(df.no3, plate == 1)))

## Regression for plate 2
summary(lm(abs.540.blank ~ as.numeric(stand.conc),
           data = subset(df.no3, plate == 2)))

# Intercept and slope of regression equation for plate 2
coef(lm(abs.540.blank ~ as.numeric(stand.conc),
        data = subset(df.no3, plate == 2)))

###############################################################################
## Add slope and intercept column to df.nh4. Slope and intercept will be used
## to calculate [NH4-N] for positive controls, standard curve, and unknowns
###############################################################################
df.no3$slope <- ifelse(df.no3$plate == 1,
                       coef(lm(abs.540.blank ~ as.numeric(stand.conc),
                               data = subset(df.no3, plate == 1)))[2],
                       coef(lm(abs.540.blank ~ as.numeric(stand.conc),
                               data = subset(df.no3, plate == 2)))[2])
df.no3$intercept <- ifelse(df.no3$plate == 1,
                           coef(lm(abs.540.blank ~ as.numeric(stand.conc),
                                   data = subset(df.no3,
                                                 plate == 1)))[1],
                           coef(lm(abs.540.blank ~ as.numeric(stand.conc),
                                   data = subset(df.no3,
                                                 plate == 2)))[1])

###############################################################################
## Use slope and intercept from regressions to calculate NH4 concentrations,
## remove any rows with NA
###############################################################################
no3.bag <- df.no3 %>%
  mutate(no3.conc = (abs.540.blank - intercept) / slope,
         no3.conc = ifelse(no3.conc < 0.4, NA, no3.conc)) %>%
  filter(complete.cases(no3.conc))
no3.bag

###############################################################################
## Calculate coefficient of variation for positive controls. Acceptable
## controls are less than 10%
###############################################################################
cv <- no3.bag %>%
  filter(type == "POS") %>%
  group_by(plate, stand.conc) %>%
  summarize(pc.sd = sd(no3.conc, na.rm = TRUE),
            pc.sd = ifelse(is.na(pc.sd) == TRUE, 0, pc.sd),
            pc.mean = mean(no3.conc, na.rm = TRUE),
            cv = pc.sd / pc.mean * 100)
cv
###############################################################################
## If coefficient of variation is less than 10%, calculate [NH4-N] mean of 
## duplicate samples
###############################################################################
no3.mean <- no3.bag %>%
  filter(type == "SPL") %>%
  group_by(id) %>%
  summarize(no3.conc = mean(no3.conc, na.rm = TRUE),
            no3.conc = ifelse(no3.conc < 0.4, NA, no3.conc)) %>%
  data.frame()
no3.mean

###############################################################################
## Write .csv with id and [NH4-N]. Be sure to save this file to maintain
## a reproducible workflow
###############################################################################
#write.csv(nh4.mean, "[insert path here]")
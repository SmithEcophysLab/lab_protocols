###############################################################################
## Load libraries
###############################################################################
library(dplyr)
library(tidyr)

###############################################################################
## Read in .csv with all plate reader data for ammonium assays
###############################################################################
df <- read.csv("NH4_NO3_data_example.csv")
head(df)

###############################################################################
## Remove any incorrect or unreliable blank values
###############################################################################
df$abs.650[2] <- NA

###############################################################################
## Calculate mean blank values for each plate, then subtract mean blank values
## from abs.650 values to determine blank-corrected absorbance
###############################################################################
df.nh4 <- df %>%
  filter(analysis.type == "nh4") %>%
  select(plate, id, type, stand.conc, abs.650) %>%
  group_by(plate) %>%
  filter(id == "blank") %>%
  summarize(mean.blank = mean(abs.650, na.rm = TRUE)) %>%
  full_join(df) %>%
  filter(analysis.type == "nh4" & id != "blank") %>%
  select(plate, id, type, stand.conc, abs.650, mean.blank) %>%
  mutate(abs.650.blank = as.numeric(abs.650) - as.numeric(mean.blank))
df.nh4

###############################################################################
## Linear regression to calculate R^2 and linear model coefficients. R^2 should
## be greater than 0.98. Can remove points to improve model fit
###############################################################################
# Regression for plate 1
summary(lm(abs.650.blank ~ as.numeric(stand.conc),
           data = subset(df.nh4, plate == 1 & type == "STD")))

# Intercept and slope of regression equation for plate 1
coef(lm(abs.650.blank ~ as.numeric(stand.conc),
        data = subset(df.nh4, plate == 1 & type == "STD")))

# Visualize standard curve for plate 1
plot(x = subset(df.nh4, plate == 1 & type == "STD")$stand.conc,
     y = subset(df.nh4, plate == 1 & type == "STD")$abs.650.blank)
abline(a = coef(lm(abs.650.blank ~ as.numeric(stand.conc),
                   data = subset(df.nh4, 
                                 plate == 1 & 
                                   type == "STD")))[[1]], 
       b = coef(lm(abs.650.blank ~ as.numeric(stand.conc),
                   data = subset(df.nh4, 
                                 plate == 1 & 
                                   type == "STD")))[[2]])

## Regression for plate 2
summary(lm(abs.650.blank ~ as.numeric(stand.conc),
           data = subset(df.nh4, plate == 2 & type == "STD")))

# Intercept and slope of regression equation for plate 2
coef(lm(abs.650.blank ~ as.numeric(stand.conc),
        data = subset(df.nh4, plate == 2 & type == "STD")))

# Visualize standard curve for plate 1
plot(x = subset(df.nh4, plate == 2 & type == "STD")$stand.conc,
     y = subset(df.nh4, plate == 2 & type == "STD")$abs.650.blank)
abline(a = coef(lm(abs.650.blank ~ as.numeric(stand.conc),
                   data = subset(df.nh4, 
                                 plate == 2 & 
                                   type == "STD")))[[1]], 
       b = coef(lm(abs.650.blank ~ as.numeric(stand.conc),
                   data = subset(df.nh4, 
                                 plate == 2 & 
                                   type == "STD")))[[2]])

###############################################################################
## Add slope and intercept column to df.nh4. Slope and intercept will be used
## to calculate [NH4-N] for positive controls, standard curve, and unknowns
###############################################################################
df.nh4$slope <- ifelse(df.nh4$plate == 1,
                       coef(lm(abs.650.blank ~ as.numeric(stand.conc),
                               data = subset(df.nh4, 
                                             plate == 1 & type == "STD")))[2],
                       coef(lm(abs.650.blank ~ as.numeric(stand.conc),
                               data = subset(df.nh4, plate == 2)))[2])
df.nh4$int <- ifelse(df.nh4$plate == 1,
                     coef(lm(abs.650.blank ~ as.numeric(stand.conc),
                             data = subset(df.nh4,
                                           plate == 1  & type == "STD")))[1],
                     coef(lm(abs.650.blank ~ as.numeric(stand.conc),
                             data = subset(df.nh4,
                                           plate == 2  & type == "STD")))[1])

head(df.nh4)
###############################################################################
## Use slope and intercept from regressions to calculate NH4 concentrations,
## remove any rows with NA
###############################################################################
nh4.bag <- df.nh4 %>%
  mutate(nh4.conc = (abs.650.blank - int) / slope,
         nh4.conc = ifelse(nh4.conc < 0.4, NA, nh4.conc)) %>%
  filter(complete.cases(nh4.conc))
nh4.bag

###############################################################################
## Calculate coefficient of variation for positive controls. Acceptable
## controls are less than 10% (cv is on percent basis)
###############################################################################
cv <- nh4.bag %>%
  filter(type == "POS") %>%
  group_by(plate, stand.conc) %>%
  summarize(pc.sd = sd(nh4.conc, na.rm = TRUE),
            pc.sd = ifelse(is.na(pc.sd) == TRUE, 0, pc.sd),
            pc.mean = mean(nh4.conc, na.rm = TRUE),
            cv = pc.sd / pc.mean * 100)
cv

###############################################################################
## If coefficient of variation is less than 10%, calculate [NH4-N] mean of 
## duplicate samples
###############################################################################
nh4.mean <- nh4.bag %>%
  filter(type == "SPL") %>%
  group_by(id) %>%
  summarize(nh4.conc = mean(nh4.conc, na.rm = TRUE),
            nh4.conc = ifelse(nh4.conc < 0.4, NA, nh4.conc)) %>%
  data.frame()
nh4.mean

###############################################################################
## Write .csv with id and [NH4-N]. Be sure to save this file to maintain
## a reproducible workflow
###############################################################################
## write.csv(nh4.mean, "./NH4_analysis_data.csv")




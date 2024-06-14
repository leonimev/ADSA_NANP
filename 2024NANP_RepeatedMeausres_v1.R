################### Mixed Model Repeated Measures Analysis
### Created by: Leoni F. Martins | The Pennsylvania State University
### June 13, 2024

################# Warning !!!
### All equations should be carefully evaluated before code is used for analyses.

rm(list=ls()) # to clean environment tab
graphics.off() # to delete open graphics

getwd() # to get working directory
setwd("~/Desktop/Leoni/Lectures/NANP/2024") # to set wd

# Install necessary packages, if not installed yet
#install.packages("lme4")
#install.packages("nlme") # package for repeated measures
#install.packages("emmeans") # For least-squares means and contrasts
#install.packages("readxl")
#install.packages("dplyr")
#install.packages("lubridate") # For preparing date and time data
#install.packages("ggplot2")
#install.packages("forcats")

# Load necessary libraries
library(lme4)
library(nlme) # package for repeated measures
library(emmeans) # For least-squares means and contrasts
library(readxl)
library(dplyr)
library(lubridate) # For preparing date and time data
library(ggplot2)
library(forcats)


# Read the data
data <- read_excel("2024NANP_EntericGasData.xlsx", sheet = "Study 1")

ls(data) # for column names

sapply(data, class) # check class of all variables

num = c("CO2", "CH4", "H2")
data[num] <- lapply(data[num], as.numeric) # for numeric

fact = c("Cow", "Period", "Parity", "Time", "Square", "Treatment")
data[fact] <- lapply(data[fact], as.factor) # for factor

data$Time <- fct_reorder(data$Time,as.integer(data$Time)) # for integer and reorder

summary(data)

sum(is.na(data$CH4)) # Check for missing values in the "CH4" variable

data1 <- na.omit(data) # Remove missing data

# Process the date and time data (this code split columns and round time base don sampling timepoints)

#data1$`Start Time` <- ymd_hms(data1$`Start Time`)
#str(data1)

# Extract components using dplyr and lubridate
#data1 <- data1 %>%
 # mutate(
  #  Year = year(`Start Time`),
   # Month = month(`Start Time`),
    # Day = day(`Start Time`),
    # TimeR = format(`Start Time`, format = "%H:%M:%S")
  #)

# Display the modified data frame
#print(data1)

#sampling_timepoints <- c(1, 3, 5, 8, 11, 13, 15, 17, 19, 21, 23) # Definition of sampling timepoints

#find_closest_timepoint <- function(datetime, timepoints) {
 # closest_timepoint <- timepoints[which.min(abs(datetime - timepoints))]
  #return(closest_timepoint)
#}

#data$RoundedTime <- sapply(data$`Start Time`, function(x) {
 # find_closest_timepoint(hour(x), sampling_timepoints)
#})



# Fit the mixed model using lmer:
# can provide Satterthwaite or Kenward-Roger approximations
# more suitable for small sample sizes or complex designs
# assumes independent and identically distributed residuals

#help(lmer) # learn more about lmer

lmod_mix <-lmer(CH4 ~ Treatment*Time + Period + (1|Cow:Square), data1)
anova(lmod_mix, ddf = "Satterthwaite") # Type 3 ANOVA with Satterthwaite's method # other option "Kenward-Roger"
joint_tests(lmod_mix) # Type 3 ANOVA Wald tests (not adjusted on df)

if(requireNamespace("pbkrtest", quietly = TRUE))
anova(lmod_mix, type=2, ddf="Kenward-Roger") # Type 2 with Kenward-Roger method # you can also use Satterthwaite's method

summary(lmod_mix)
#print(lmod_mix, correlation=TRUE)
#vcov(lmod_mix)

#ranova(lmod_mix)

#(vc.lm <-VarCorr(lmod_mix))
#print(vc.lm,comp=c("Variance"))

lsmeans_Trt= emmeans(lmod_mix,"Treatment") # overall marginal means for Treatment
summary(lsmeans_Trt)  

pairs(lsmeans_Trt, adj = "none") #pairwise comparisons between marginal means
multcomp::cld(lsmeans_Trt, by = NULL, Letters = "abcdefg", alpha = .05, adj = "none")

# Treatment x time interaction
TrtTime = emmeans(lmod_mix,"Treatment",by="Time") # slice by time
summary(TrtTime)
pairs(TrtTime,adjust="none")
plot(TrtTime,by = "Time")

TrtTimedf <- data.frame(TrtTime)

custom_colors <- c("Con" = "navy", "SF1" = "burlywood", "SF2" = "azure4") # define treatment colors

ggplot(data = TrtTimedf, aes(x = Time, y = emmean, fill = Treatment)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL),
                position = position_dodge(width = 0.9), width = 0.25) +
  labs(x = "Time", y = "Estimated Mean") +
  scale_fill_manual(values = custom_colors) +  # Adjust colors as needed
  theme_minimal()

emmip(lmod_mix,Treatment~Time,main="interaction means plot",ylab="Mean score")
multcomp::cld(TrtTime, by = "Time", Letters = "abcdefg", alpha = .05, adjust="none")


# Fit the mixed model using lme:
# provides functionality for modeling correlated errors and heteroscedasticity

# Compound symetry

fit.cs <- lme(CH4 ~ Treatment * Time + Period, random = ~ 1|Cow, data = data1,
               corr = corCompSymm(, form= Time ~ 1 |Cow/Square), method="REML")

summary(fit.cs)
anova(lmod_mix, ddf = "Satterthwaite") 
joint_tests(fit.cs)

# Ar(1)

fit.ar1 <- lme(CH4 ~ Treatment * Time + Period, random = ~ 1|Cow, data = data1,
               corr = corAR1(, form= Time ~ 1 |Cow/Square), method="REML")

summary(fit.ar1)
joint_tests(fit.ar1) #Type 3 ANOVA

# Compare models based on BIC
BIC_compare<-BIC(fit.cs, fit.ar1)
BIC_compare<-anova(fit.cs, fit.ar1)
BIC_compare = BIC_compare[order(BIC_compare$BIC),]

BIC_compare

# Compound symetry

fit.cs <- lme(CH4 ~ Treatment * Time + Period, random = ~ 1|Cow, data = data1,
              corr = corCompSymm(, form=  ~ 1 |Cow/Square), method="REML")

summary(fit.cs)
anova(lmod_mix, ddf = "Satterthwaite") 
joint_tests(fit.cs) #Type 3 ANOVA

lsmeans_Trt= emmeans(fit.cs,"Treatment") # overall marginal means for Treatment
summary(lsmeans_Trt)  

pairs(lsmeans_Trt, adj = "none") #pairwise comparisons between marginal means (default adjustment is tukey)
multcomp::cld(lsmeans_Trt, by = NULL, Letters = "abcdefg", alpha = .05, adj = "none")

TrtTime = emmeans(fit.cs,"Treatment",by="Time")
summary(TrtTime)
pairs(TrtTime,adjust="none")
plot(TrtTime,by = "Time")

# Plot TrtTime by Time

TrtTimedf <- data.frame(TrtTime)

custom_colors <- c("Con" = "navy", "SF1" = "burlywood", "SF2" = "azure4") # define treatment colors

ggplot(data = TrtTimedf, aes(x = Time, y = emmean, fill = Treatment)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL),
                position = position_dodge(width = 0.9), width = 0.25) +
  labs(x = "Time", y = "Estimated Mean") +
  scale_fill_manual(values = custom_colors) +  # Adjust colors as needed
  theme_minimal()

emmip(fit.cs,Treatment~Time,main="interaction means plot",ylab="Mean score")
multcomp::cld(TrtTime, by = "Time", Letters = "abcdefg", alpha = .05, adjust="none")

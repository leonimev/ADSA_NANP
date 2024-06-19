################### Circadian rhythm analysis
### Created by: Leoni F. Martins | The Pennsylvania State University
### June 13, 2024


################# Warning !!!
### All equations should be carefully evaluated before code is used for analyses. 
### Codes were created based on Bourdon et al. (1995), Seltman (1997), Niu et al. (2014)


rm(list=ls()) # to clean environment tab
graphics.off() # to delete open graphics

# Install necessary packages, if not installed yet
#install.packages("readxl")
#install.packages("dplyr")
#install.packages("tidyr")
#install.packages("lme4")
#install.packages("ggplot2")

# Load necessary libraries
library(readxl)
library(dplyr)
library(tidyr)
library(lme4)
library(ggplot2)

getwd() # to get working directory
setwd("~/Desktop/Leoni/Lectures/NANP/2024") # to set wd

# Read the data
data <- read_excel("2024NANP_EntericGasData.xlsx", sheet = "Study 1")

# Standardize column names. 
# All models were created with Response, trt, time, cow, period.
# Change according to variable name described in your data set.

data <- data %>%
  rename(
    Response = CH4,
    trt = Treatment,
    time = Time,
    cow = Cow,
    period = Period
  )

# Convert necessary columns to numeric
data$Response <- as.numeric(data$Response)

# Convert necessary columns to factors
data$trt <- as.factor(data$trt)
data$cow <- as.factor(data$cow)
data$period <- as.factor(data$period)

# Remove missing values in response variables
data <- na.omit(data)

# Calculate cosine and sine components for 24-hour and 4-hour cycles
data <- data %>%
  mutate(
    cos_24 = cos(2 * pi * time / 24),
    sin_24 = sin(2 * pi * time / 24),
    cos_12 = cos(2 * pi * time / 12),
    sin_12 = sin(2 * pi * time / 12)
  )

# Fit the initial mixed-effects model. (1|cow) sets cow as a random effect
initial_model <- lmer(Response ~ trt * cos_24 + trt * sin_24 + cos_12 + sin_12 + (1|cow), data = data)

# Print the summary of the initial model
summary(initial_model)

# Calculate residuals
residuals <- residuals(initial_model)

# Identify outliers based on studentized residuals >= 3.5
studentized_residuals <- residuals / sd(residuals)
outliers <- abs(studentized_residuals) >= 3.5
num_outliers <- sum(outliers)
cat("Number of outliers removed:", num_outliers, "\n")

# Remove outliers
data_cleaned <- data[!outliers, ]

# Fit the final mixed-effects model with cow as a random effect
final_model <- lmer(Response ~ trt * cos_24 + trt * sin_24 + cos_12 + sin_12 + (1|cow), data = data_cleaned)

# Print the summary of the final model
summary(final_model)

# Perform the zero-amplitude likelihood ratio test. 
# This determines if there is a rhythm in the Response variable
linear_model <- lmer(Response ~ 1 + trt + (1|cow), data = data_cleaned)
final_model <- lmer(Response ~ trt*cos_24 + trt*sin_24 + cos_12 + sin_12 + (1|cow), data = data_cleaned)

lr_stat <- 2 * (logLik(final_model) - logLik(linear_model))
p_value <- pchisq(lr_stat, df = 4, lower.tail = FALSE)

cat("Zero-Amplitude LRT: LR stat =", lr_stat, ", p-value =", p_value, "\n")

# Extract coefficients for each treatment
results <- data.frame()

# Example treatments - automatically fetch from data_cleaned
treatments <- unique(data_cleaned$trt)

# Extracting coefficient names to determine the reference treatment
coef_names <- names(fixef(final_model))
interaction_terms <- coef_names[grepl("trt", coef_names)]
trt_ref <- setdiff(treatments, gsub("trt", "", unique(sub("[:].*", "", interaction_terms))))


for (trt in treatments) {
  # Extract coefficients
  beta_cos_24 <- fixef(final_model)["cos_24"]
  beta_sin_24 <- fixef(final_model)["sin_24"]
  
  # Adjust for treatment-specific coefficients
  if (trt != trt_ref) {
    beta_cos_24 <- beta_cos_24 + ifelse(paste0("trt", trt, ":cos_24") %in% names(fixef(final_model)), fixef(final_model)[paste0("trt", trt, ":cos_24")], 0)
    beta_sin_24 <- beta_sin_24 + ifelse(paste0("trt", trt, ":sin_24") %in% names(fixef(final_model)), fixef(final_model)[paste0("trt", trt, ":sin_24")], 0)
  }
  
  amplitude_24 <- sqrt(beta_cos_24^2 + beta_sin_24^2)
  
  # Calculate SE for 24-hour amplitude
  var_cos_24 <- vcov(final_model)["cos_24", "cos_24"]
  var_sin_24 <- vcov(final_model)["sin_24", "sin_24"]
  if (trt != trt_ref) {
    if (paste0("trt", trt, ":cos_24") %in% colnames(vcov(final_model))) {
      var_cos_24 <- var_cos_24 + vcov(final_model)[paste0("trt", trt, ":cos_24"), paste0("trt", trt, ":cos_24")]
      var_cos_24 <- var_cos_24 + 2 * vcov(final_model)["cos_24", paste0("trt", trt, ":cos_24")]
    }
    if (paste0("trt", trt, ":sin_24") %in% colnames(vcov(final_model))) {
      var_sin_24 <- var_sin_24 + vcov(final_model)[paste0("trt", trt, ":sin_24"), paste0("trt", trt, ":sin_24")]
      var_sin_24 <- var_sin_24 + 2 * vcov(final_model)["sin_24", paste0("trt", trt, ":sin_24")]
    }
  }
  SE_amplitude_24 <- sqrt(var_cos_24 + var_sin_24)
  
  # Calculate acrophase_24
  if (!is.na(beta_cos_24) && !is.na(beta_sin_24) && !is.nan(beta_cos_24) && !is.nan(beta_sin_24) && beta_cos_24 != 0) {
    acrophase_24 <- atan2(beta_sin_24, beta_cos_24) * (12 / pi)
    acrophase_24 <- (acrophase_24 + 24) %% 24  # Ensure acrophase is within 0-24 hours
  } else {
    acrophase_24 <- NA  # Handle cases where acrophase cannot be calculated
  }
  SE_acrophase_24 <- sqrt(var_cos_24 + var_sin_24) * (12 / pi) # standard error for 24-hour acrophase
  
  # For the 12-hour cycle
  beta_cos_12 <- fixef(final_model)["cos_12"]
  beta_sin_12 <- fixef(final_model)["sin_12"]
  
  # Adjust for treatment-specific coefficients
  if (trt != trt_ref) {
    beta_cos_12 <- beta_cos_12 + ifelse(paste0("trt", trt, ":cos_12") %in% names(fixef(final_model)), fixef(final_model)[paste0("trt", trt, ":cos_12")], 0)
    beta_sin_12 <- beta_sin_12 + ifelse(paste0("trt", trt, ":sin_12") %in% names(fixef(final_model)), fixef(final_model)[paste0("trt", trt, ":sin_12")], 0)
  }
  
  amplitude_12 <- sqrt(beta_cos_12^2 + beta_sin_12^2)
  
  # Calculate SE for 12-hour amplitude
  var_cos_12 <- vcov(final_model)["cos_12", "cos_12"]
  var_sin_12 <- vcov(final_model)["sin_12", "sin_12"]
  if (trt != trt_ref) {
    if (paste0("trt", trt, ":cos_12") %in% colnames(vcov(final_model))) {
      var_cos_12 <- var_cos_12 + vcov(final_model)[paste0("trt", trt, ":cos_12"), paste0("trt", trt, ":cos_12")]
      var_cos_12 <- var_cos_12 + 2 * vcov(final_model)["cos_12", paste0("trt", trt, ":cos_12")]
    }
    if (paste0("trt", trt, ":sin_12") %in% colnames(vcov(final_model))) {
      var_sin_12 <- var_sin_12 + vcov(final_model)[paste0("trt", trt, ":sin_12"), paste0("trt", trt, ":sin_12")]
      var_sin_12 <- var_sin_12 + 2 * vcov(final_model)["sin_12", paste0("trt", trt, ":sin_12")]
    }
  }
  SE_amplitude_12 <- sqrt(var_cos_12 + var_sin_12)
  
  # Calculate acrophase_12
  if (!is.na(beta_cos_12) && !is.na(beta_sin_12) && !is.nan(beta_cos_12) && !is.nan(beta_sin_12) && beta_cos_12 != 0) {
    acrophase_12 <- atan2(beta_sin_12, beta_cos_12) * (6 / pi)
    acrophase_12 <- (acrophase_12 + 12) %% 12  # Ensure acrophase is within 0-12 hours
  } else {
    acrophase_12 <- NA  # Handle cases where acrophase cannot be calculated
  }
  SE_acrophase_12 <- sqrt(var_cos_12 + var_sin_12) * (6 / pi) # standard error for 4-hour acrophase
  
  # Calculate mesor
  mesor <- fixef(final_model)["(Intercept)"]
  if (trt != trt_ref) {
    mesor <- mesor + ifelse(paste0("trt", trt) %in% names(fixef(final_model)), fixef(final_model)[paste0("trt", trt)], 0)
  }
  
  # Confidence intervals
  CI_amplitude_24 <- c(amplitude_24 - 1.96 * SE_amplitude_24, amplitude_24 + 1.96 * SE_amplitude_24)
  CI_acrophase_24 <- c(acrophase_24 - 1.96 * SE_acrophase_24, acrophase_24 + 1.96 * SE_acrophase_24)
  CI_amplitude_12 <- c(amplitude_12 - 1.96 * SE_amplitude_12, amplitude_12 + 1.96 * SE_amplitude_12)
  CI_acrophase_12 <- c(acrophase_12 - 1.96 * SE_acrophase_12, acrophase_12 + 1.96 * SE_acrophase_12)
  
  results <- rbind(results, data.frame(
    Treatment = trt,
    Mesor = mesor,
    Amplitude_24 = amplitude_24,
    SE_Amplitude_24 = SE_amplitude_24,
    CI_Amplitude_24_Lower = CI_amplitude_24[1],
    CI_Amplitude_24_Upper = CI_amplitude_24[2],
    Acrophase_24 = acrophase_24,
    SE_Acrophase_24 = SE_acrophase_24,
    CI_Acrophase_24_Lower = CI_acrophase_24[1],
    CI_Acrophase_24_Upper = CI_acrophase_24[2],
    Amplitude_12 = amplitude_12,
    SE_Amplitude_12 = SE_amplitude_12,
    CI_Amplitude_12_Lower = CI_amplitude_12[1],
    CI_Amplitude_12_Upper = CI_amplitude_12[2],
    Acrophase_12 = acrophase_12,
    SE_Acrophase_12 = SE_acrophase_12,
    CI_Acrophase_12_Lower = CI_acrophase_12[1],
    CI_Acrophase_12_Upper = CI_acrophase_12[2]
  ))
}
# Print results
print(results)


# Hypothesis testing for overlapping CI based on Knezevic 2008
for (i in 1:(nrow(results) - 1)) {
  for (j in (i + 1):nrow(results)) {
    trt1 <- results$Treatment[i]
    trt2 <- results$Treatment[j]
    
    # For 24-hour amplitude
    t_amplitude_24 <- (results$Amplitude_24[i] - results$Amplitude_24[j]) / 
      sqrt(results$SE_Amplitude_24[i]^2 + results$SE_Amplitude_24[j]^2)
    p_value_amplitude_24 <- 2 * pt(-abs(t_amplitude_24), df = nrow(data_cleaned) - length(fixef(final_model)))
    
    # For 24-hour acrophase
    t_acrophase_24 <- (results$Acrophase_24[i] - results$Acrophase_24[j]) / 
      sqrt(results$SE_Acrophase_24[i]^2 + results$SE_Acrophase_24[j]^2)
    p_value_acrophase_24 <- 2 * pt(-abs(t_acrophase_24), df = nrow(data_cleaned) - length(fixef(final_model)))
    
    # For 12-hour amplitude
    t_amplitude_12 <- (results$Amplitude_12[i] - results$Amplitude_12[j]) / 
      sqrt(results$SE_Amplitude_12[i]^2 + results$SE_Amplitude_12[j]^2)
    p_value_amplitude_12 <- 2 * pt(-abs(t_amplitude_12), df = nrow(data_cleaned) - length(fixef(final_model)))
    
    # For 12-hour acrophase
    t_acrophase_12 <- (results$Acrophase_12[i] - results$Acrophase_12[j]) / 
      sqrt(results$SE_Acrophase_12[i]^2 + results$SE_Acrophase_12[j]^2)
    p_value_acrophase_12 <- 2 * pt(-abs(t_acrophase_12), df = nrow(data_cleaned) - length(fixef(final_model)))
    
    cat("Comparison of", trt1, "and", trt2, ":\n")
    cat("24-hour amplitude: t =", t_amplitude_24, ", p-value =", p_value_amplitude_24, "\n")
    cat("24-hour acrophase: t =", t_acrophase_24, ", p-value =", p_value_acrophase_24, "\n")
    cat("4-hour amplitude: t =", t_amplitude_12, ", p-value =", p_value_amplitude_12, "\n")
    cat("4-hour acrophase: t =", t_acrophase_12, ", p-value =", p_value_acrophase_12, "\n\n")
  }
}


# Plotting fitted curves and observed data points by treatment

# Generate time points
time_points <- seq(0, 24, length.out = 100)

# Function to calculate fitted curve for a given treatment
calculate_fitted_curve <- function(trt) {
  cos_24_points <- cos(2 * pi * time_points / 24)
  sin_24_points <- sin(2 * pi * time_points / 24)
  cos_12_points <- cos(2 * pi * time_points / 12)
  sin_12_points <- sin(2 * pi * time_points / 12)
  
  if (trt == trt_ref) {
    fitted_curve <- (fixef(final_model)["(Intercept)"] +
                       fixef(final_model)["cos_24"] * cos_24_points +
                       fixef(final_model)["sin_24"] * sin_24_points +
                       fixef(final_model)["cos_12"] * cos_12_points +
                       fixef(final_model)["sin_12"] * sin_12_points)
  } else {
    fitted_curve <- (fixef(final_model)["(Intercept)"] +
                       fixef(final_model)["cos_24"] * cos_24_points +
                       fixef(final_model)["sin_24"] * sin_24_points +
                       fixef(final_model)[paste0("trt", trt, ":cos_24")] * cos_24_points +
                       fixef(final_model)[paste0("trt", trt, ":sin_24")] * sin_24_points +
                       fixef(final_model)["cos_12"] * cos_12_points +
                       fixef(final_model)["sin_12"] * sin_12_points)
  }
  
  return(data.frame(time = time_points, fitted_curve = fitted_curve, trt = trt))
}

# Generate data frames for each treatment
fitted_curves <- lapply(treatments, calculate_fitted_curve)

# Combine data frames
fitted_curves_df <- do.call(rbind, fitted_curves)

# Plotting using ggplot
ggplot() +
  geom_line(data = fitted_curves_df, aes(x = time, y = fitted_curve, color = trt), linewidth = 1, alpha = 0.8) +
  #geom_point(data = data_cleaned, aes(x = time, y = Response, color = trt), size = 3, alpha = 0.6) +
  labs(title = "Fitted Curves and Observed Data Points by Treatment", x = "Time (hours)", y = "Response") +
  theme_minimal()


###########################################################


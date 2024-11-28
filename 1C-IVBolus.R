# Concentration-time data frame for 1 g IV bolus dose of drug X
time <- c(1, 2, 4, 6, 8, 12)
concentration <- c(93.5, 72.8, 44.1, 26.8, 16.2, 6.0)
ct_data <- data.frame(Time = time, Concentration = concentration)

# Visualizing the concentration-time data
ggplot(ct_data, aes(x = Time, y = Concentration)) +
  geom_point(color = "#62bec6") +
  geom_line(color = "gray") +
  labs(title = "Concentration-Time Data of Drug X", x = "Time (hours)", y = "Concentration (mg/L)") +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),  # Bold the title
    axis.title = element_text(face = "bold")               # Bold the axis labels
  )

# Visualizing the log-transformed concentration-time data
ggplot(ct_data, aes(x = Time, y = Concentration)) +
  geom_point(color = "#62bec6") +
  geom_line(color = "gray") +
  scale_y_log10() +
  labs(title = "Concentration-Time Data of Drug X", x = "Time (hours)", y = "logConcentration (mg/L)") +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),  # Bold the title
    axis.title = element_text(face = "bold") # Bold the axis labels
  )

# Adding time point 0 for AUC analysis
ct_data_0 <- rbind(data.frame(Time = 0, Concentration = 93.5), ct_data)

# Performing NCA using PK package
nca_result <- pk.calc.auc(
  time = ct_data_0$Time,
  conc = ct_data_0$Concentration,
  method = "lin up/log down"
)

# Calculating the elimination rate constant and half-life
k <- -coef(lm(log(ct_data$Concentration) ~ ct_data$Time))[2]  # Slope of log-linear portion
half_life <- log(2) / k

# Calculating Cmax and Tmax
Cmax <- max(ct_data_0$Concentration)  # Maximum concentration
Tmax <- ct_data$Time[which.max(ct_data$Concentration)]  # Time of maximum concentration

# Calculating AUC to infinity
last_time <- tail(ct_data_0$Time, 1)  # Last time point
last_conc <- tail(ct_data_0$Concentration, 1)  # Last concentration
AUC_inf <- nca_result + (last_conc / k)  # AUC(0-t) + extrapolated part

# Displaying NCA results
cat("Non-Compartmental Analysis Results:\n")
cat("AUC (0-t):", nca_result, "mg·h/L\n")
cat("AUC (0-∞):", AUC_inf, "mg·h/L\n")
cat("Cmax:", Cmax, "mg/L\n")
cat("Tmax:", Tmax, "hours\n")
cat("Elimination Rate Constant (k):", k, "1/h\n")
cat("Estimated Terminal Half-Life (t1/2):", half_life, "hours\n")

# Defining the one-compartment model
one_compartment_model <- function(C0, k, Time) {
  C0 * exp(-k * Time)
}

# Initial parameter guesses
# Assumed C0 = C1 because of IV route and used k from NCA
initial_params <- c(C0 = 93.5, k = 0.25)

# Fitting the model using nls
fit <- nls(
  Concentration ~ one_compartment_model(C0, k, Time),
  data = ct_data,
  start = initial_params
)

# Displaying parameter estimates
summary(fit)

# Extracting the fitted parameters
C0 <- coef(fit)["C0"]
k <- coef(fit)["k"]
V <- 1000 / C0 # Dose = 1 g = 1000 mg
CL <- V * k

# Displaying compartmental analysis results
cat("C0 (Initial Concentration):", C0, "mg/L\n")
cat("k (Elimination Rate Constant):", k, "1/h\n")
cat("V (Volume of Distribution):", V, "L\n")
cat("CL (Clearance):", CL, "L/h\n")

# Predicting concentrations using the fitted model
ct_data$Predicted <- predict(fit, newdata = ct_data)

# Plotting observed + predicted concentrations vs time
ggplot(ct_data, aes(x = Time)) +
  geom_point(aes(y = Concentration), color = "#62bec6") +
  geom_line(aes(y = Predicted), linetype = "dashed") +
  labs(title = "Observed & Predicted Concentrations vs Time",
       x = "Time (hours)",
       y = "Concentration (mg/L)") +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),  # Bold the title
    axis.title = element_text(face = "bold") # Bold the axis labels
  )

# Plotting observed vs predicted concentrations
ggplot(ct_data, aes(x = Predicted, y = Concentration)) +
  geom_point(color = "#62bec6") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  labs(
    title = "Observed vs Predicted Concentrations",
    x = "Predicted Concentration (mg/L)",
    y = "Observed Concentration (mg/L)"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),  # Bold the title
    axis.title = element_text(face = "bold") # Bold the axis labels
  )

# Residual analysis
ct_data$Residuals <- ct_data$Concentration - ct_data$Predicted
rss <- sum(ct_data$Residuals^2) # Residual Sum of Squares (RSS)
mse <- mean(ct_data$Residuals^2) # Mean Squared Error (MSE)

# Displaying residual analysis results
cat("Residual Sum of Squares (RSS):", rss, "\n")
cat("Mean Squared Error (MSE):", mse, "\n")

# Calculating and displaying R-squared result
ss_total <- sum((ct_data$Concentration - mean(ct_data$Concentration))^2)
ss_residual <- sum(ct_data$Residuals^2)
r_squared <- 1 - (ss_residual / ss_total)
cat("R-squared:", r_squared, "\n")

# Calculating and displaying AIC result
aic <- AIC(fit)
cat("Akaike Information Criterion (AIC):", aic, "\n")

# Visualing the residual analyses
# Residuals vs time
ggplot(ct_data, aes(x = Time, y = Residuals)) +
  geom_point(color = "#62bec6") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(
    title = "Residuals vs Time",
    x = "Time (hours)",
    y = "Residuals (mg/L)"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),  # Bold the title
    axis.title = element_text(face = "bold") # Bold the axis labels
  )

# Residuals vs predicted concentrations
ggplot(ct_data, aes(x = Predicted, y = Residuals)) +
  geom_point(color = "#62bec6") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(
    title = "Residuals vs Predicted Concentrations",
    x = "Predicted Concentration (mg/L)",
    y = "Residuals (mg/L)"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),  # Bold the title
    axis.title = element_text(face = "bold") # Bold the axis labels
  )
library(ggplot2)
library(forecast)
library(tseries)
library(gridExtra)
library(splines)
library(mgcv)
library(nlme)
library(dplyr)
library(DIMORA)
library(pdp)
library(ggcorrplot)
library(gbm)
library(caret)
library(tinytex)

# Load the dataset
df <- read.csv("Unemployment Rate.csv", header = TRUE)
df1 <- read.csv("Production in Industry.csv", header = TRUE)
df2 <- read.csv("EUR-USD Exchange Rate.csv", header = TRUE)
df3 <- read.csv("Retail Trade Volume.csv", header = TRUE)
df4 <- read.csv("Harmonized Index of Consumer Prices.csv", header = TRUE)

# Preprocess the data
date_column_in_csv <- "TIME_PERIOD"
unemployment_rate_column_in_csv <- "OBS_VALUE"

df$Date <- df[[date_column_in_csv]]

if (inherits(df$Date, "character")) {
  temp_date <- as.Date(df$Date, format="%Y-%m-%d")
  if (any(is.na(temp_date)) && !all(is.na(as.Date(paste0(df$Date, "-01"), format="%Y-%m-%d")))) {
    df$Date <- as.Date(paste0(df$Date, "-01"), format="%Y-%m-%d")
  } else {
    df$Date <- temp_date
  }
  if (all(is.na(df$Date))) {
    warning(paste0("Warning: Could not convert '", date_column_in_csv, "' to a date format. Using row index as time for plotting."))
    df$Date <- 1:nrow(df)
  }
} else if (!inherits(df$Date, c("Date", "POSIXct", "POSIXlt"))) {
  warning(paste0("Warning: Date column '", date_column_in_csv, "' is not in a recognized date format. Using row index as time for plotting."))
  df$Date <- 1:nrow(df)
}

df$Unemployment_Rate <- as.numeric(df[[unemployment_rate_column_in_csv]])

if (any(is.na(df$Unemployment_Rate))) {
  warning("Warning: NA values introduced in 'Unemployment_Rate' column during conversion. Please check your data for non-numeric entries.")
  df <- na.omit(df, cols=c("Unemployment_Rate"))
  message("NA values in Unemployment_Rate column removed for analysis.")
}

df$time_index <- 1:nrow(df)

if (inherits(df$Date, "Date")) {
  df$Month <- factor(format(df$Date, "%m"), levels = sprintf("%02d", 1:12), labels = month.abb)
  df$Month_Number <- as.numeric(format(df$Date, "%m")) # Add this line to get 1-12 numeric month
} else {
  warning("Date column not a proper date object. Cannot extract month for seasonality.")
  df$Month <- NULL
  df$Month_Number <- NULL # And this
}

start_year <- as.numeric(format(df$Date[1], "%Y"))
start_month <- as.numeric(format(df$Date[1], "%m"))

unemployment_ts <- ts(df$Unemployment_Rate, start = c(start_year, start_month), frequency = 12)

# additional datasets
preprocess_additional <- function(data, date_col = "TIME_PERIOD", value_col = "OBS_VALUE", value_name) {
  
  data$Date <- data[[date_col]]
  
  if (inherits(data$Date, "character")) {
    temp_date <- as.Date(data$Date, format = "%Y-%m-%d")
    if (any(is.na(temp_date)) && !all(is.na(as.Date(paste0(data$Date, "-01"), format = "%Y-%m-%d")))) {
      data$Date <- as.Date(paste0(data$Date, "-01"), format = "%Y-%m-%d")
    } else {
      data$Date <- temp_date
    }
    if (all(is.na(data$Date))) {
      warning(paste0("Could not convert '", date_col, "' to date format. Using row index."))
      data$Date <- 1:nrow(data)
    }
  } else if (!inherits(data$Date, c("Date", "POSIXct", "POSIXlt"))) {
    warning(paste0("Date column '", date_col, "' not in recognized format. Using row index."))
    data$Date <- 1:nrow(data)
  }
  
  data[[value_name]] <- as.numeric(data[[value_col]])
  
  if (any(is.na(data[[value_name]]))) {
    warning(paste("NA values introduced in", value_name, "during conversion"))
    data <- na.omit(data, cols = value_name)
    message(paste("NA values in", value_name, "removed"))
  }
  
  return(data[, c("Date", value_name)])
}


df_processed <- preprocess_additional(data = df,
                                      date_col = "TIME_PERIOD",
                                      value_col = "OBS_VALUE",
                                      value_name = "Unemployment_Rate")


df1_processed <- preprocess_additional(data = df1,
                                       date_col = "TIME_PERIOD",
                                       value_col = "OBS_VALUE",
                                       value_name = "Industrial_Production")


df2_processed <- preprocess_additional(data = df2,
                                       date_col = "TIME_PERIOD",
                                       value_col = "OBS_VALUE",
                                       value_name = "Exchange_Rate")


df3_processed <- preprocess_additional(data = df3,
                                       date_col = "TIME_PERIOD",
                                       value_col = "OBS_VALUE",
                                       value_name = "Retail_Trade")


df4_processed <- preprocess_additional(data = df4,
                                       date_col = "TIME_PERIOD",
                                       value_col = "OBS_VALUE",
                                       value_name = "Consumer_Price_Index")


# Merge datasets
merged_df <- df_processed %>%
  select(Date, Unemployment_Rate) %>%
  left_join(df1_processed, by = "Date") %>%
  left_join(df2_processed, by = "Date") %>%
  left_join(df3_processed, by = "Date") %>%
  left_join(df4_processed, by = "Date") %>%
  na.omit()

# Add monthly and yearly lag
merged_df <- merged_df %>%
  arrange(Date) %>%
  mutate(
    Unemployment_Rate_lag1 = lag(Unemployment_Rate, 1),
    Industrial_Production_lag1 = lag(Industrial_Production, 1),
    Exchange_Rate_lag1 = lag(Exchange_Rate, 1),
    Retail_Trade_lag1 = lag(Retail_Trade, 1),
    Consumer_Price_Index_lag1 = lag(Consumer_Price_Index, 1),
    Unemployment_Rate_lag12 = lag(Unemployment_Rate, 12),
    Industrial_Production_lag12 = lag(Industrial_Production, 12),
    Exchange_Rate_lag12 = lag(Exchange_Rate, 12),
    Retail_Trade_lag12 = lag(Retail_Trade, 12),
    Consumer_Price_Index_lag12 = lag(Consumer_Price_Index, 12)
  ) %>%
  na.omit() 

head(merged_df)
names(merged_df)



# Unemployment Rate Time Series, ACF, PACF

p1 <- ggplot(df, aes(x = Date, y = Unemployment_Rate)) +
  geom_line(color = "darkblue", linewidth = 0.6) +
  labs(title = "Unemployment Rate Time Series", x = "Date", y = "Unemployment Rate (%)") +
  theme_minimal() + theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 16))
print(p1)

p2 <- ggAcf(df$Unemployment_Rate, lag.max = 72) +
  geom_segment(aes(xend = lag, yend = 0), size = 20) +
  labs(title = "ACF of Unemployment Rate", x = "Lag", y = "ACF") +
  theme_minimal() + theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 16))
print(p2)

p3 <- ggPacf(df$Unemployment_Rate, lag.max = 72) +
  labs(title = "PACF of Unemployment Rate", x = "Lag", y = "PACF") +
  theme_minimal() + theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 16))
print(p3)

# Model 1: Linear Regression with Quadratic Trend + Seasonality

split_date <- as.Date("2020-01-01")

# Prepare data
train_data <- df %>%
  filter(Date < split_date, Date > as.Date("2008-01-01")) %>%
  mutate(time_index_temp = 1:n())

last_train_index <- max(train_data$time_index_temp)
forecast_data <- df %>%
  filter(Date >= split_date) %>%
  mutate(time_index_temp = last_train_index + 1:n())

# Fit the model
lm_model_quad_seasonal <- lm(Unemployment_Rate ~ time_index_temp + I(time_index_temp^2) + Month, data = train_data)

print(summary(lm_model_quad_seasonal))

# Predict train data
train_data$fitted_values <- predict(lm_model_quad_seasonal, newdata = train_data)

# Forecast
forecast_results <- predict(lm_model_quad_seasonal, newdata = forecast_data, interval = "confidence", level = 0.95)

forecast_data$forecasted_unemployment <- forecast_results[, "fit"]
forecast_data$lower_ci <- forecast_results[, "lwr"]
forecast_data$upper_ci <- forecast_results[, "upr"]

# Plotting
df_after_2008 <- df %>% filter(Date > as.Date("2008-01-01"))

p_qreg_forecast <- ggplot() +
  geom_line(data = df_after_2008, aes(x = Date, y = Unemployment_Rate, color = "Original Data"), linewidth = 0.8) +
  geom_line(data = train_data, aes(x = Date, y = fitted_values, color = "Predicted Values"), linewidth = 0.8) +
  geom_line(data = forecast_data, aes(x = Date, y = forecasted_unemployment, color = "Predicted Values"), linewidth = 0.8, linetype = "solid") +
  geom_ribbon(data = forecast_data, aes(x = Date, ymin = lower_ci, ymax = upper_ci, fill = "95% Confidence Interval"), alpha = 0.4) +
  
  geom_vline(xintercept = as.numeric(split_date), linetype = "dotted", color = "black", linewidth = 0.8) +
  
  labs(title = "Linear Regression \n (Quadratic Trend + Seasonality)",
       x = "Date",
       y = "Unemployment Rate (%)") +
  scale_color_manual(name = "Legend",
                     values = c("Original Data" = "red",
                                "Predicted Values" = "darkblue")) +
  scale_fill_manual(name = "Legend",
                    values = c("95% Confidence Interval" = "lightblue")) +
  guides(fill = guide_legend(order = 2),
         color = guide_legend(order = 1)) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 12),
        legend.position = "right",
        legend.box = "vertical")

print(p_qreg_forecast)


# Residual Analysis
train_data$residuals_quad_seasonal <- residuals(lm_model_quad_seasonal)

p4_quad_seasonal <- ggplot(train_data, aes(x = time_index_temp, y = residuals_quad_seasonal)) +
  geom_line(color = "darkblue", linewidth = 0.6) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  labs(title = "Linear Regression (Quadratic Trend + Seasonality)", x = "Time", y = "Residuals") +
  theme_minimal() + theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 16))
print(p4_quad_seasonal)

p5_quad_seasonal <- ggAcf(train_data$residuals_quad_seasonal, lag.max = 36) +
  labs(title = "Linear Regression (Quadratic Trend + Seasonality)", x = "Lag", y = "ACF") +
  theme_minimal() + theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 16))
print(p5_quad_seasonal)

p6_quad_seasonal <- ggplot(train_data, aes(x = residuals_quad_seasonal)) +
  geom_histogram(aes(y = after_stat(density)), binwidth = 0.2, fill = "lightgreen", color = "darkgreen", alpha = 0.7) +
  geom_density(aes(color = "Empirical Density"), linewidth = 1) +
  stat_function(aes(color = "Normal Distribution"), fun = dnorm,
                args = list(mean = mean(train_data$residuals_quad_seasonal, na.rm = TRUE), sd = sd(train_data$residuals_quad_seasonal, na.rm = TRUE)),
                linetype = "dashed", linewidth = 1) +
  scale_color_manual(name = "Distribution Type", values = c("Empirical Density" = "darkblue", "Normal Distribution" = "purple")) +
  labs(title = "Linear Regression (Quadratic Trend + Seasonality)", x = "Residuals", y = "Density") +
  theme_minimal() + theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 16), legend.position = "bottom", legend.title = element_text(face = "bold"))
print(p6_quad_seasonal)

# RMSE, MAPE, AIC, BIC

# rmse function
calculate_rmse <- function(actual, predicted) {
  sqrt(mean((actual - predicted)^2, na.rm = TRUE))
}

# mape function
calculate_mape <- function(actual, predicted) {
  # Avoid division by zero for actual values close to zero
  # Replace 0 with a small epsilon to prevent Inf or NaN
  actual_safe <- ifelse(actual == 0, .Machine$double.eps, actual)
  mean(abs((actual - predicted) / actual_safe), na.rm = TRUE) * 100
}

# Calculate RMSE for training set
rmse_train_lm <- calculate_rmse(train_data$Unemployment_Rate, train_data$fitted_values)
cat(paste0("Training RMSE: ", round(rmse_train_lm, 4), "\n"))

# Calculate MAPE for training set
mape_train_lm <- calculate_mape(train_data$Unemployment_Rate, train_data$fitted_values)
cat(paste0("Training MAPE: ", round(mape_train_lm, 2), "%\n"))

# Calculate RMSE for test/forecast set
rmse_test_lm <- calculate_rmse(forecast_data$Unemployment_Rate, forecast_data$forecasted_unemployment)
cat(paste0("Test/Forecast RMSE: ", round(rmse_test_lm, 4), "\n"))

# Calculate MAPE for test/forecast set
mape_test_lm <- calculate_mape(forecast_data$Unemployment_Rate, forecast_data$forecasted_unemployment)
cat(paste0("Test/Forecast MAPE: ", round(mape_test_lm, 2), "%\n"))

aic_lm <- AIC(lm_model_quad_seasonal)
bic_lm <- BIC(lm_model_quad_seasonal)

cat(paste0("AIC (on Training Data): ", round(aic_lm, 2), "\n"))
cat(paste0("BIC (on Training Data): ", round(bic_lm, 2), "\n"))


# Model 2: Local Regression (LOESS)

loess_model <- loess(Unemployment_Rate ~ time_index, data = df, span = 0.75)
df$fitted_loess <- predict(loess_model, newdata = df)
df$residuals_loess <- residuals(loess_model)


p_loess_fit <- ggplot(df, aes(x = Date, y = Unemployment_Rate)) +
  geom_point(aes(color = "Original Data")) +
  geom_line(aes(y = fitted_loess, color = "LOESS Fit"), linewidth = 0.8) +
  labs(title = "LOESS", x = "Date", y = "Unemployment Rate (%)") +
  scale_color_manual(name = "Data Type", values = c("Original Data" = "red", "LOESS Fit" = "darkblue")) +
  theme_minimal() + theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 16), legend.position = "bottom")
print(p_loess_fit)

p4_loess <- ggplot(df, aes(x = time_index, y = residuals_loess)) +
  geom_line(color = "darkblue", linewidth = 0.6) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  labs(title = "LOESS", x = "Time", y = "Residuals") +
  theme_minimal() + theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 16))
print(p4_loess)

p5_loess <- ggAcf(df$residuals_loess, lag.max = 36) +
  labs(title = "LOESS", x = "Lag", y = "ACF") +
  theme_minimal() + theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 16))
print(p5_loess)

p6_loess <- ggplot(df, aes(x = residuals_loess)) +
  geom_histogram(aes(y = after_stat(density)), binwidth = 0.2, fill = "lightgreen", color = "darkgreen", alpha = 0.7) +
  geom_density(aes(color = "Empirical Density"), linewidth = 1) +
  stat_function(aes(color = "Normal Distribution"), fun = dnorm,
                args = list(mean = mean(df$residuals_loess, na.rm = TRUE), sd = sd(df$residuals_loess, na.rm = TRUE)),
                linetype = "dashed", linewidth = 1) +
  scale_color_manual(name = "Distribution Type", values = c("Empirical Density" = "darkblue", "Normal Distribution" = "purple")) +
  labs(title = "LOESS", x = "Residuals", y = "Density") +
  theme_minimal() + theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 16), legend.position = "bottom", legend.title = element_text(face = "bold"))
print(p6_loess)



# Model 3: Regression Splines d=1

# Define Knot Points
knot_date_1 <- as.Date("2008-01-01")
knot_time_index_1 <- df$time_index[which.min(abs(df$Date - knot_date_1))]

knot_date_2 <- as.Date("2016-01-01")
knot_time_index_2 <- df$time_index[which.min(abs(df$Date - knot_date_2))]

knots_val <- c(knot_time_index_1, knot_time_index_2)


lm_model_spline1 <- lm(Unemployment_Rate ~ bs(time_index, degree = 1, knots = knots_val), data = df)
df$fitted_spline1 <- predict(lm_model_spline1, newdata = df)
df$residuals_spline1 <- residuals(lm_model_spline1)

print(summary(lm_model_spline1))

# Plot the Spline fit over the original data
p_spline1_fit <- ggplot(df, aes(x = Date, y = Unemployment_Rate)) +
  geom_point(aes(color = "Original Data")) +
  geom_line(aes(y = fitted_spline1, color = "Spline Fit (d = 1)"), linewidth = 0.8) +
  labs(title = "Regression Spline (d = 1)", x = "Date", y = "Unemployment Rate (%)") +
  scale_color_manual(name = "Data Type", values = c("Original Data" = "red", "Spline Fit (d = 1)" = "darkblue")) +
  theme_minimal() + theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 16), legend.position = "bottom")
print(p_spline1_fit)

p4_spline1 <- ggplot(df, aes(x = time_index, y = residuals_spline1)) +
  geom_line(color = "darkblue", linewidth = 0.6) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  labs(title = "Regression Spline (d = 1)", x = "Time", y = "Residuals") +
  theme_minimal() + theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 16))
print(p4_spline1)

p5_spline1 <- ggAcf(df$residuals_spline1, lag.max = 36) +
  labs(title = "Regression Spline (d = 1)", x = "Lag", y = "ACF") +
  theme_minimal() + theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 16))
print(p5_spline1)

p6_spline1 <- ggplot(df, aes(x = residuals_spline1)) +
  geom_histogram(aes(y = after_stat(density)), binwidth = 0.2, fill = "lightgreen", color = "darkgreen", alpha = 0.7) +
  geom_density(aes(color = "Empirical Density"), linewidth = 1) +
  stat_function(aes(color = "Normal Distribution"), fun = dnorm,
                args = list(mean = mean(df$residuals_spline1, na.rm = TRUE), sd = sd(df$residuals_spline1, na.rm = TRUE)),
                linetype = "dashed", linewidth = 1) +
  scale_color_manual(name = "Distribution Type", values = c("Empirical Density" = "darkblue", "Normal Distribution" = "purple")) +
  labs(title = "Regression Spline (d = 1)", x = "Residuals", y = "Density") +
  theme_minimal() + theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 16), legend.position = "bottom", legend.title = element_text(face = "bold"))
print(p6_spline1)


# Model 4: Regression Splines d=2

lm_model_spline2 <- lm(Unemployment_Rate ~ bs(time_index, degree = 2, knots = knots_val), data = df)
df$fitted_spline2 <- predict(lm_model_spline2, newdata = df)
df$residuals_spline2 <- residuals(lm_model_spline2)

print(summary(lm_model_spline2))

# Plot the Spline fit over the original data
p_spline2_fit <- ggplot(df, aes(x = Date, y = Unemployment_Rate)) +
  geom_point(aes(color = "Original Data")) +
  geom_line(aes(y = fitted_spline2, color = "Spline Fit (d = 2)"), linewidth = 0.8) +
  labs(title = "Regression Spline (d = 2)", x = "Date", y = "Unemployment Rate (%)") +
  scale_color_manual(name = "Data Type", values = c("Original Data" = "red", "Spline Fit (d = 2)" = "darkblue")) +
  theme_minimal() + theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 16), legend.position = "bottom")
print(p_spline2_fit)

p4_spline2 <- ggplot(df, aes(x = time_index, y = residuals_spline2)) +
  geom_line(color = "darkblue", linewidth = 0.6) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  labs(title = "Regression Spline (d = 2)", x = "Time", y = "Residuals") +
  theme_minimal() + theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 16))
print(p4_spline2)

p5_spline2 <- ggAcf(df$residuals_spline2, lag.max = 36) +
  labs(title = "Regression Spline (d = 2)", x = "Lag", y = "ACF") +
  theme_minimal() + theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 16))
print(p5_spline2)

p6_spline2 <- ggplot(df, aes(x = residuals_spline2)) +
  geom_histogram(aes(y = after_stat(density)), binwidth = 0.2, fill = "lightgreen", color = "darkgreen", alpha = 0.7) +
  geom_density(aes(color = "Empirical Density"), linewidth = 1) +
  stat_function(aes(color = "Normal Distribution"), fun = dnorm,
                args = list(mean = mean(df$residuals_spline2, na.rm = TRUE), sd = sd(df$residuals_spline2, na.rm = TRUE)),
                linetype = "dashed", linewidth = 1) +
  scale_color_manual(name = "Distribution Type", values = c("Empirical Density" = "darkblue", "Normal Distribution" = "purple")) +
  labs(title = "Regression Spline (d = 2)", x = "Residuals", y = "Density") +
  theme_minimal() + theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 16), legend.position = "bottom", legend.title = element_text(face = "bold"))
print(p6_spline2)

# Model 5: Regression Splines d=3

lm_model_spline3 <- lm(Unemployment_Rate ~ bs(time_index, degree = 3, knots = knots_val), data = df)
df$fitted_spline3 <- predict(lm_model_spline3, newdata = df)
df$residuals_spline3 <- residuals(lm_model_spline3)

print(summary(lm_model_spline3))

# Plot the Spline fit over the original data
p_spline3_fit <- ggplot(df, aes(x = Date, y = Unemployment_Rate)) +
  geom_point(aes(color = "Original Data")) +
  geom_line(aes(y = fitted_spline3, color = "Spline Fit (d = 3)"), linewidth = 0.8) +
  labs(title = "Regression Spline (d = 3)", x = "Date", y = "Unemployment Rate (%)") +
  scale_color_manual(name = "Data Type", values = c("Original Data" = "red", "Spline Fit (d = 3)" = "darkblue")) +
  theme_minimal() + theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 16), legend.position = "bottom")
print(p_spline3_fit)

p4_spline3 <- ggplot(df, aes(x = time_index, y = residuals_spline3)) +
  geom_line(color = "darkblue", linewidth = 0.6) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  labs(title = "Regression Spline (d = 3)", x = "Time", y = "Residuals") +
  theme_minimal() + theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 16))
print(p4_spline3)

p5_spline3 <- ggAcf(df$residuals_spline3, lag.max = 36) +
  labs(title = "Regression Spline (d = 3)", x = "Lag", y = "ACF") +
  theme_minimal() + theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 16))
print(p5_spline3)

p6_spline3 <- ggplot(df, aes(x = residuals_spline3)) +
  geom_histogram(aes(y = after_stat(density)), binwidth = 0.2, fill = "lightgreen", color = "darkgreen", alpha = 0.7) +
  geom_density(aes(color = "Empirical Density"), linewidth = 1) +
  stat_function(aes(color = "Normal Distribution"), fun = dnorm,
                args = list(mean = mean(df$residuals_spline3, na.rm = TRUE), sd = sd(df$residuals_spline2, na.rm = TRUE)),
                linetype = "dashed", linewidth = 1) +
  scale_color_manual(name = "Distribution Type", values = c("Empirical Density" = "darkblue", "Normal Distribution" = "purple")) +
  labs(title = "Regression Spline (d = 3)", x = "Residuals", y = "Density") +
  theme_minimal() + theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 16), legend.position = "bottom", legend.title = element_text(face = "bold"))
print(p6_spline3)

# Add all degrees
p_splines_fit <- ggplot(df, aes(x = Date, y = Unemployment_Rate)) +
  geom_point(aes(color = "Original Data")) +
  geom_line(aes(y = fitted_spline1, color = "d = 1"), linewidth = 0.8) +
  geom_line(aes(y = fitted_spline2, color = "d = 2"), linewidth = 0.8) +
  geom_line(aes(y = fitted_spline3, color = "d = 3"), linewidth = 0.8) +
  labs(title = "Regression Splines (d = 1,2,3)", x = "Date", y = "Unemployment Rate (%)") +
  scale_color_manual(name = "Data Type", values = c("Original Data" = "grey", "d = 1" = "black", "d = 2" = "blue", "d = 3" = "red")) +
  theme_minimal() + theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 16), legend.position = "bottom")
print(p_splines_fit)

# Model 6: Smoothing Splines

smoothspline_model <- smooth.spline(x = df$time_index, y = df$Unemployment_Rate, cv = TRUE)
df$fitted_smoothspline <- predict(smoothspline_model, x = df$time_index)$y
df$residuals_smoothspline <- df$Unemployment_Rate - df$fitted_smoothspline


# Calculating optimal lambda
optimal_spar <- smoothspline_model$spar
optimal_gcv <- smoothspline_model$crit

spar_range <- seq(0.1, 1.5, length.out = 100) # Adjust range as needed
gcv_scores <- numeric(length(spar_range))

for (i in seq_along(spar_range)) {
  temp_model <- smooth.spline(x = df$time_index, y = df$Unemployment_Rate, spar = spar_range[i])
  gcv_scores[i] <- temp_model$crit
}

plot_data <- data.frame(
  spar = spar_range,
  gcv = gcv_scores
)

ggplot(plot_data, aes(x = spar, y = gcv)) +
  geom_line(color = "darkblue", size = 1) +
  geom_vline(xintercept = optimal_spar, linetype = "dashed", color = "red", size = 1) +
  geom_point(x = optimal_spar, y = optimal_gcv, color = "red", size = 3, shape = 16) +
  labs(title = "Bias-Variance Trade-off: GCV vs. Spar",
       x = "Spar Value",
       y = "Generalized Cross-Validation Score") +
  theme_minimal() +
  annotate("text", x = optimal_spar + 0.1, y = max(gcv_scores) * 0.9,
           label = paste0("Optimal Spar: ", round(optimal_spar, 3)),
           color = "red", hjust = 0)

print(smoothspline_model)

# Plot the Smoothing Spline fit over the original data
p_smoothspline_fit <- ggplot(df, aes(x = Date, y = Unemployment_Rate)) +
  geom_point(aes(color = "Original Data")) +
  geom_line(aes(y = fitted_smoothspline, color = "Smoothing Spline Fit"), linewidth = 0.8) +
  labs(title = expression("Smoothing Spline (" * lambda * " = 0.176)"), x = "Date", y = "Unemployment Rate (%)") +
  scale_color_manual(name = "Data Type", values = c("Original Data" = "red", "Smoothing Spline Fit" = "darkblue")) +
  theme_minimal() + theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 16), legend.position = "bottom")
print(p_smoothspline_fit)

p4_smoothspline <- ggplot(df, aes(x = time_index, y = residuals_smoothspline)) +
  geom_line(color = "darkblue", linewidth = 0.6) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  labs(title = expression("Smoothing Spline (" * lambda * " = 0.176)"), x = "Time", y = "Residuals") +
  theme_minimal() + theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 16))
print(p4_smoothspline)

p5_smoothspline <- ggAcf(df$residuals_smoothspline, lag.max = 36) +
  labs(title = expression("Smoothing Spline (" * lambda * " = 0.176)"), x = "Lag", y = "ACF") +
  theme_minimal() + theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 16))
print(p5_smoothspline)

p6_smoothspline <- ggplot(df, aes(x = residuals_smoothspline)) +
  geom_histogram(aes(y = after_stat(density)), binwidth = 0.2, fill = "lightgreen", color = "darkgreen", alpha = 0.7) +
  geom_density(aes(color = "Empirical Density"), linewidth = 1) +
  stat_function(aes(color = "Normal Distribution"), fun = dnorm,
                args = list(mean = mean(df$residuals_smoothspline, na.rm = TRUE), sd = sd(df$residuals_smoothspline, na.rm = TRUE)),
                linetype = "dashed", linewidth = 1) +
  scale_color_manual(name = "Distribution Type", values = c("Empirical Density" = "darkblue", "Normal Distribution" = "purple")) +
  labs(title = expression("Smoothing Spline (" * lambda * " = 0.176)"), x = "Residuals", y = "Density") +
  theme_minimal() + theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 16), legend.position = "bottom", legend.title = element_text(face = "bold"))
print(p6_smoothspline)

# Data Splitting for Time Series Models (Holt-Winters, ARIMA)

cutoff_date <- as.Date("2020-01-01")

# Create the training dataset
train_df <- df[df$Date < cutoff_date, ]

# Create the forecast/test dataset
forecast_data <- df[df$Date >= cutoff_date, ]

start_year_train <- as.numeric(format(train_df$Date[1], "%Y"))
start_month_train <- as.numeric(format(train_df$Date[1], "%m"))

train_ts <- ts(train_df$Unemployment_Rate, 
               start = c(start_year_train, start_month_train), 
               frequency = 12)


# Model 7: Holt-Winters

# Fit HoltWinters model on training data
holt_winters_model <- hw(train_ts, seasonal = "multiplicative")

fitted_hw_model <- holt_winters_model$model
print(summary(fitted_hw_model))

# Residuals for Holt-Winters
df_residuals_hw <- data.frame(
  Date = df$Date[1:length(residuals(holt_winters_model))],
  time_index = df$time_index[1:length(residuals(holt_winters_model))],
  residuals_holt_winters = as.numeric(residuals(holt_winters_model))
)

p4_hw <- ggplot(df_residuals_hw, aes(x = time_index, y = residuals_holt_winters)) +
  geom_line(color = "darkblue", linewidth = 0.6) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  labs(title = "Holt-Winters", x = "Time", y = "Residuals") +
  theme_minimal() + theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 16))
print(p4_hw)

p5_hw <- ggAcf(df_residuals_hw$residuals_holt_winters[!is.na(df_residuals_hw$residuals_holt_winters)], lag.max = 36) +
  labs(title = "Holt-Winters", x = "Lag", y = "ACF") +
  theme_minimal() + theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 16))
print(p5_hw)

p6_hw <- ggplot(df_residuals_hw, aes(x = residuals_holt_winters)) +
  geom_histogram(aes(y = after_stat(density)), binwidth = 0.02, fill = "lightgreen", color = "darkgreen", alpha = 0.7) +
  geom_density(aes(color = "Empirical Density"), linewidth = 1) +
  stat_function(fun = dnorm, args = list(mean = mean(df_residuals_hw$residuals_holt_winters, na.rm = TRUE), sd = sd(df_residuals_hw$residuals_holt_winters, na.rm = TRUE)),
                aes(color = "Normal Distribution"), linetype = "dashed", linewidth = 1) +
  scale_color_manual(name = "Distribution Type", values = c("Empirical Density" = "darkblue", "Normal Distribution" = "purple")) +
  labs(title = "Holt-Winters", x = "Residuals", y = "Density") +
  theme_minimal() + theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 16), legend.position = "bottom", legend.title = element_text(face = "bold"))
print(p6_hw)


# Forecast

forecast_hw <- forecast(fitted_hw_model, h = nrow(forecast_data))

fitted_values_hw_df <- data.frame(
  Date = df$Date[1:length(fitted(holt_winters_model))],
  Value = as.numeric(fitted(holt_winters_model))
)

forecast_df <- data.frame(
  Date = forecast_data$Date,
  Value = as.numeric(forecast_hw$mean)
)

ci_hw_df <- data.frame(
  Date = forecast_data$Date,
  Lower_CI = as.numeric(forecast_hw$lower[, 2]),
  Upper_CI = as.numeric(forecast_hw$upper[, 2])
)

p_forecast_hw <- ggplot(mapping = aes(x = Date)) +
  geom_line(data = df, aes(y = Unemployment_Rate, color = "Original Data"), linewidth = 0.8) +
  geom_line(data = forecast_df, aes(y = Value, color = "Forecast Values"), linewidth = 0.8) +
  geom_line(data = fitted_values_hw_df, aes(y = Value, color = "Predicted Values"), linewidth = 0.8) +
  geom_ribbon(data = ci_hw_df, aes(ymin = Lower_CI, ymax = Upper_CI, fill = "95% Confidence Interval"), alpha = 0.2) +
  geom_vline(xintercept = as.numeric(as.Date("2020-01-01")), linetype = "dashed", color = "grey50", linewidth = 1) +
  labs(title = "Holt-Winters", x = "Date", y = "Unemployment Rate (%)") +
  scale_color_manual(name = "Data Type", values = c("Original Data" = "red", "Predicted Values" = "darkblue", "Forecast Values" = "darkcyan")) +
  scale_fill_manual(name = "", values = c("95% Confidence Interval" = "darkcyan")) +
  theme_minimal() + theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 16), legend.position = "bottom")
print(p_forecast_hw)


# RMSE, AIC, BIC, MAPE

# Calculate RMSE for training set
rmse_train_hw <- calculate_rmse(train_ts, fitted(holt_winters_model))
cat(paste0("Training RMSE: ", round(rmse_train_hw, 4), "\n"))

# Calculate MAPE for training set
mape_train_hw <- calculate_mape(train_ts, fitted(holt_winters_model))
cat(paste0("Training MAPE: ", round(mape_train_hw, 2), "%\n"))

# Calculate RMSE for test/forecast set
rmse_test_hw <- calculate_rmse(forecast_data$Unemployment_Rate, forecast_df$Value)
cat(paste0("Test/Forecast RMSE: ", round(rmse_test_hw, 4), "\n"))

# Calculate MAPE for test/forecast set
mape_test_hw <- calculate_mape(forecast_data$Unemployment_Rate, forecast_df$Value)
cat(paste0("Test/Forecast MAPE: ", round(mape_test_hw, 2), "%\n"))

# Calculate AIC and BIC for the fitted model (on training data)
aic_hw <- AIC(fitted_hw_model)
bic_hw <- BIC(fitted_hw_model)

cat(paste0("AIC (on Training Data): ", round(aic_hw, 2), "\n"))
cat(paste0("BIC (on Training Data): ", round(bic_hw, 2), "\n"))


# Model 8: auto ARIMA

arima_model <- auto.arima(train_ts, seasonal = TRUE)

print(summary(arima_model))

# Residuals for ARIMA
df_residuals_arima <- data.frame(
  Date = df$Date[1:length(residuals(arima_model))],
  time_index = df$time_index[1:length(residuals(arima_model))],
  residuals_arima = as.numeric(residuals(arima_model))
)

p4_arima <- ggplot(df_residuals_arima, aes(x = time_index, y = residuals_arima)) +
  geom_line(color = "darkblue", linewidth = 0.6) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  labs(title = "ARIMA (0,1,1)(0,1,1)[12]", x = "Time", y = "Residuals") +
  theme_minimal() + theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 16))
print(p4_arima)

p5_arima <- ggAcf(df_residuals_arima$residuals_arima[!is.na(df_residuals_arima$residuals_arima)], lag.max = 36) +
  labs(title = "ARIMA (0,1,1)(0,1,1)[12]", x = "Lag", y = "ACF") +
  theme_minimal() + theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 16))
print(p5_arima)

p6_arima <- ggplot(df_residuals_arima, aes(x = residuals_arima)) +
  geom_histogram(aes(y = after_stat(density)), binwidth = 0.2, fill = "lightgreen", color = "darkgreen", alpha = 0.7) +
  geom_density(aes(color = "Empirical Density"), linewidth = 1) +
  stat_function(fun = dnorm, args = list(mean = mean(df_residuals_arima$residuals_arima, na.rm = TRUE), sd = sd(df_residuals_arima$residuals_arima, na.rm = TRUE)),
                aes(color = "Normal Distribution"), linetype = "dashed", linewidth = 1) +
  scale_color_manual(name = "Distribution Type", values = c("Empirical Density" = "darkblue", "Normal Distribution" = "purple")) +
  labs(title = "ARIMA (0,1,1)(0,1,1)[12]", x = "Residuals", y = "Density") +
  theme_minimal() + theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 16), legend.position = "bottom", legend.title = element_text(face = "bold"))
print(p6_arima)


# Forecast
forecast_arima <- forecast(arima_model, h = nrow(forecast_data))

fitted_values_arima_df <- data.frame(
  Date = df$Date[1:length(fitted(arima_model))],
  Value = as.numeric(fitted(arima_model))
)

forecast_arima_df <- data.frame(
  Date = forecast_data$Date,
  Value = as.numeric(forecast_arima$mean)
)

ci_arima_df <- data.frame(
  Date = forecast_data$Date,
  Lower_CI = as.numeric(forecast_arima$lower[, 2]),
  Upper_CI = as.numeric(forecast_arima$upper[, 2])
)

p_forecast_arima <- ggplot(mapping = aes(x = Date)) +
  geom_line(data = df, aes(y = Unemployment_Rate, color = "Original Data"), linewidth = 0.8) +
  geom_line(data = forecast_arima_df, aes(y = Value, color = "Forecast Values"), linewidth = 0.8) +
  geom_line(data = fitted_values_arima_df, aes(y = Value, color = "Predicted Values"), linewidth = 0.8) +
  geom_ribbon(data = ci_arima_df, aes(ymin = Lower_CI, ymax = Upper_CI, fill = "95% Confidence Interval"), alpha = 0.2) +
  geom_vline(xintercept = as.numeric(as.Date("2020-01-01")), linetype = "dashed", color = "grey50", linewidth = 1) +
  labs(title = "ARIMA (0,1,1)(0,1,1)[12]", x = "Date", y = "Unemployment Rate (%)") +
  scale_color_manual(name = "Data Type", values = c("Original Data" = "red", "Predicted Values" = "darkblue", "Forecast Values" = "darkcyan")) +
  scale_fill_manual(name = "", values = c("95% Confidence Interval" = "darkcyan")) +
  theme_minimal() + theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 16), legend.position = "bottom")
print(p_forecast_arima)


# RMSE, AIC, BIC, MAPE

# Calculate RMSE for training set
rmse_train_arima <- calculate_rmse(train_ts, fitted(arima_model))
cat(paste0("Training RMSE: ", round(rmse_train_arima, 4), "\n"))

# Calculate MAPE for training set
mape_train_arima <- calculate_mape(train_ts, fitted(arima_model))
cat(paste0("Training MAPE: ", round(mape_train_arima, 2), "%\n"))

# Calculate RMSE for test/forecast set
rmse_test_arima <- calculate_rmse(forecast_data$Unemployment_Rate, forecast_arima_df$Value)
cat(paste0("Test/Forecast RMSE: ", round(rmse_test_arima, 4), "\n"))

# Calculate MAPE for test/forecast set
mape_test_arima <- calculate_mape(forecast_data$Unemployment_Rate, forecast_arima_df$Value)
cat(paste0("Test/Forecast MAPE: ", round(mape_test_arima, 2), "%\n"))

# Calculate AIC and BIC for the fitted model (on training data)
aic_arima <- AIC(arima_model)
bic_arima <- BIC(arima_model)

cat(paste0("AIC (on Training Data): ", round(aic_arima, 2), "\n"))
cat(paste0("BIC (on Training Data): ", round(bic_arima, 2), "\n"))

# Model 8_2: ARIMA

arima_model <- arima(train_ts, order = c(1, 1, 0), seasonal = c(1, 1, 0))

print(summary(arima_model))

# Residuals for ARIMA
df_residuals_arima <- data.frame(
  Date = df$Date[1:length(residuals(arima_model))],
  time_index = df$time_index[1:length(residuals(arima_model))],
  residuals_arima = as.numeric(residuals(arima_model))
)

p4_arima <- ggplot(df_residuals_arima, aes(x = time_index, y = residuals_arima)) +
  geom_line(color = "darkblue", linewidth = 0.6) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  labs(title = "ARIMA (1,1,0)(1,1,0)[12]", x = "Time", y = "Residuals") +
  theme_minimal() + theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 16))
print(p4_arima)

p5_arima <- ggAcf(df_residuals_arima$residuals_arima[!is.na(df_residuals_arima$residuals_arima)], lag.max = 36) +
  labs(title = "ARIMA (1,1,0)(1,1,0)[12]", x = "Lag", y = "ACF") +
  theme_minimal() + theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 16))
print(p5_arima)

p6_arima <- ggplot(df_residuals_arima, aes(x = residuals_arima)) +
  geom_histogram(aes(y = after_stat(density)), binwidth = 0.2, fill = "lightgreen", color = "darkgreen", alpha = 0.7) +
  geom_density(aes(color = "Empirical Density"), linewidth = 1) +
  stat_function(fun = dnorm, args = list(mean = mean(df_residuals_arima$residuals_arima, na.rm = TRUE), sd = sd(df_residuals_arima$residuals_arima, na.rm = TRUE)),
                aes(color = "Normal Distribution"), linetype = "dashed", linewidth = 1) +
  scale_color_manual(name = "Distribution Type", values = c("Empirical Density" = "darkblue", "Normal Distribution" = "purple")) +
  labs(title = "ARIMA (1,1,0)(1,1,0)[12]", x = "Residuals", y = "Density") +
  theme_minimal() + theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 16), legend.position = "bottom", legend.title = element_text(face = "bold"))
print(p6_arima)


# Forecast
forecast_arima <- forecast(arima_model, h = nrow(forecast_data))

fitted_values_arima_df <- data.frame(
  Date = df$Date[1:length(fitted(arima_model))],
  Value = as.numeric(fitted(arima_model))
)

forecast_arima_df <- data.frame(
  Date = forecast_data$Date,
  Value = as.numeric(forecast_arima$mean)
)

ci_arima_df <- data.frame(
  Date = forecast_data$Date,
  Lower_CI = as.numeric(forecast_arima$lower[, 2]),
  Upper_CI = as.numeric(forecast_arima$upper[, 2])
)

p_forecast_arima <- ggplot(mapping = aes(x = Date)) +
  geom_line(data = df, aes(y = Unemployment_Rate, color = "Original Data"), linewidth = 0.8) +
  geom_line(data = forecast_arima_df, aes(y = Value, color = "Forecast Values"), linewidth = 0.8) +
  geom_line(data = fitted_values_arima_df, aes(y = Value, color = "Predicted Values"), linewidth = 0.8) +
  geom_ribbon(data = ci_arima_df, aes(ymin = Lower_CI, ymax = Upper_CI, fill = "95% Confidence Interval"), alpha = 0.2) +
  geom_vline(xintercept = as.numeric(as.Date("2020-01-01")), linetype = "dashed", color = "grey50", linewidth = 1) +
  labs(title = "ARIMA (1,1,0)(1,1,0)[12]", x = "Date", y = "Unemployment Rate (%)") +
  scale_color_manual(name = "Data Type", values = c("Original Data" = "red", "Predicted Values" = "darkblue", "Forecast Values" = "darkcyan")) +
  scale_fill_manual(name = "", values = c("95% Confidence Interval" = "darkcyan")) +
  theme_minimal() + theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 16), legend.position = "bottom")
print(p_forecast_arima)


# RMSE, AIC, BIC, MAPE

# Calculate RMSE for training set
rmse_train_arima <- calculate_rmse(train_ts, fitted(arima_model))
cat(paste0("Training RMSE: ", round(rmse_train_arima, 4), "\n"))

# Calculate MAPE for training set
mape_train_arima <- calculate_mape(train_ts, fitted(arima_model))
cat(paste0("Training MAPE: ", round(mape_train_arima, 2), "%\n"))

# Calculate RMSE for test/forecast set
rmse_test_arima <- calculate_rmse(forecast_data$Unemployment_Rate, forecast_arima_df$Value)
cat(paste0("Test/Forecast RMSE: ", round(rmse_test_arima, 4), "\n"))

# Calculate MAPE for test/forecast set
mape_test_arima <- calculate_mape(forecast_data$Unemployment_Rate, forecast_arima_df$Value)
cat(paste0("Test/Forecast MAPE: ", round(mape_test_arima, 2), "%\n"))

# Calculate AIC and BIC for the fitted model (on training data)
aic_arima <- AIC(arima_model)
bic_arima <- BIC(arima_model)

cat(paste0("AIC (on Training Data): ", round(aic_arima, 2), "\n"))
cat(paste0("BIC (on Training Data): ", round(bic_arima, 2), "\n"))









# Model 9: Generalized Additive Model (GAM)

gam_model <- gam(Unemployment_Rate ~ s(time_index), data = df)
df$fitted_gam <- predict(gam_model, newdata = df)
df$residuals_gam <- residuals(gam_model)

print(summary(gam_model))

# Plotting the original data with the GAM fit
p_gam_fit <- ggplot(df, aes(x = Date, y = Unemployment_Rate)) +
  geom_point(aes(color = "Original Data")) +
  geom_line(aes(y = fitted_gam, color = "GAM Fit"), linewidth = 0.8) +
  labs(title = "GAM", x = "Date", y = "Unemployment Rate (%)") +
  scale_color_manual(name = "Data Type", values = c("Original Data" = "red", "GAM Fit" = "darkblue")) +
  theme_minimal() + theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 16), legend.position = "bottom")
print(p_gam_fit)

p4_gam <- ggplot(df, aes(x = time_index, y = residuals_gam)) +
  geom_line(color = "darkblue", linewidth = 0.6) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  labs(title = "GAM", x = "Time", y = "Residuals") +
  theme_minimal() + theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 16))
print(p4_gam)

p5_gam <- ggAcf(df$residuals_gam[!is.na(df$residuals_gam)], lag.max = 36) +
  labs(title = "GAM", x = "Lag", y = "ACF") +
  theme_minimal() + theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 16))
print(p5_gam)

p6_gam <- ggplot(df, aes(x = residuals_gam)) +
  geom_histogram(aes(y = after_stat(density)), binwidth = 0.2, fill = "lightgreen", color = "darkgreen", alpha = 0.7) +
  geom_density(aes(color = "Empirical Density"), linewidth = 1) +
  stat_function(fun = dnorm, args = list(mean = mean(df$residuals_gam, na.rm = TRUE), sd = sd(df$residuals_gam, na.rm = TRUE)),
                aes(color = "Normal Distribution"), linetype = "dashed", linewidth = 1) +
  scale_color_manual(name = "Distribution Type", values = c("Empirical Density" = "darkblue", "Normal Distribution" = "purple")) +
  labs(title = "GAM", x = "Residuals", y = "Density") +
  theme_minimal() + theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 16), legend.position = "bottom", legend.title = element_text(face = "bold"))
print(p6_gam)










# Model 9_2: Generalized Additive Model (GAM) + Complementary Data
train_data <- merged_df[train_idx, ]

train_data$time_index <- 1:nrow(train_data)

gam_model <- gam(Unemployment_Rate ~ s(time_index) + s(Industrial_Production) +
                   s(Exchange_Rate) + s(Retail_Trade) + 
                   s(Consumer_Price_Index) , data = train_data)
train_data$fitted_gam <- predict(gam_model, newdata = train_data)
train_data$residuals_gam <- residuals(gam_model)

print(summary(gam_model))

# Plotting the original data with the GAM fit
p_gam_fit <- ggplot(train_data, aes(x = Date, y = Unemployment_Rate)) +
  geom_point(aes(color = "Original Data")) +
  geom_line(aes(y = fitted_gam, color = "GAM Fit"), linewidth = 0.8) +
  labs(title = "GAM with Complementary Data", x = "Date", y = "Unemployment Rate (%)") +
  scale_color_manual(name = "Data Type", values = c("Original Data" = "red", "GAM Fit" = "darkblue")) +
  theme_minimal() + theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 16), legend.position = "bottom")
print(p_gam_fit)

p4_gam <- ggplot(train_data, aes(x = time_index, y = residuals_gam)) +
  geom_line(color = "darkblue", linewidth = 0.6) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  labs(title = "GAM with Complementary Data", x = "Time", y = "Residuals") +
  theme_minimal() + theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 16))
print(p4_gam)

p5_gam <- ggAcf(train_data$residuals_gam[!is.na(train_data$residuals_gam)], lag.max = 36) +
  labs(title = "GAM with Complementary Data", x = "Lag", y = "ACF") +
  theme_minimal() + theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 16))
print(p5_gam)

p6_gam <- ggplot(train_data, aes(x = residuals_gam)) +
  geom_histogram(aes(y = after_stat(density)), binwidth = 0.2, fill = "lightgreen", color = "darkgreen", alpha = 0.7) +
  geom_density(aes(color = "Empirical Density"), linewidth = 1) +
  stat_function(fun = dnorm, args = list(mean = mean(train_data$residuals_gam, na.rm = TRUE), sd = sd(train_data$residuals_gam, na.rm = TRUE)),
                aes(color = "Normal Distribution"), linetype = "dashed", linewidth = 1) +
  scale_color_manual(name = "Distribution Type", values = c("Empirical Density" = "darkblue", "Normal Distribution" = "purple")) +
  labs(title = "GAM with Complementary Data", x = "Residuals", y = "Density") +
  theme_minimal() + theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 16), legend.position = "bottom", legend.title = element_text(face = "bold"))
print(p6_gam)


# Model 10: Gradient Boosting

cor_data_combined <- df_processed %>%
  select(Date, Unemployment_Rate) %>%
  left_join(df1_processed, by = "Date") %>%
  left_join(df2_processed, by = "Date") %>%
  left_join(df3_processed, by = "Date") %>%
  left_join(df4_processed, by = "Date") %>%
  arrange(Date)

cor_data <- cor_data_combined[, c("Unemployment_Rate", "Industrial_Production",
                                  "Exchange_Rate", "Retail_Trade", "Consumer_Price_Index")]

cor_data <- na.omit(cor_data)

constant_cols <- sapply(cor_data, function(x) length(unique(x)) == 1)
if(any(constant_cols)) {
  warning(paste("Constant columns removed:", names(constant_cols)[constant_cols]))
  cor_data <- cor_data[, !constant_cols]
}

cor_matrix <- cor(cor_data)

cor_plot <- ggcorrplot(cor_matrix,
                       hc.order = FALSE,
                       type = "lower",
                       lab = TRUE,
                       lab_size = 3,
                       colors = c("#6D9EC1", "white", "#E46726"),
                       title = "Correlation Matrix (Response Variable: Unemployment Rate)",
                       ggtheme = theme_minimal()) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
print(cor_plot)


set.seed(123)

# Split data into training and testing sets (80/20)
train_idx <- createDataPartition(merged_df$Unemployment_Rate, p = 0.8, list = FALSE)
train_data <- merged_df[train_idx, ]
test_data <- merged_df[-train_idx, ]

gbm_model <- gbm(
  Unemployment_Rate ~ Industrial_Production +
    Exchange_Rate + Retail_Trade + 
    Consumer_Price_Index,
  data = train_data,
  distribution = "gaussian",
  n.trees = 5000,
  interaction.depth = 4,
  shrinkage = 0.01,
  cv.folds = 2,
  n.cores = 7,
  verbose = FALSE
)


# Find optimal number of trees
best_iter <- gbm.perf(gbm_model, method = "cv")
cat(paste0("Optimal number of trees: ", best_iter, "\n"))

predictions <- predict(gbm_model, test_data, n.trees = best_iter)

# Calculate RMSE
rmse <- sqrt(mean((test_data$Unemployment_Rate - predictions)^2))
cat(paste0("GBM Model RMSE: ", round(rmse, 4), "\n"))


n_trees_seq <- 1:gbm_model$n.trees
train_error <- gbm_model$train.error
cv_error <- gbm_model$cv.error

plot_data <- data.frame(
  Iterations = n_trees_seq,
  Training_Error = train_error,
  CV_Error = cv_error
)

plot_data_long <- melt(plot_data, id.vars = "Iterations",
                       variable.name = "Error_Type",
                       value.name = "Error_Value")

performance_plot <- ggplot(plot_data_long, aes(x = Iterations, y = Error_Value, color = Error_Type)) +
  geom_line(size = 0.8) +
  geom_vline(xintercept = best_iter, linetype = "dashed", color = "blue", size = 0.8) +
  annotate("text", x = best_iter + 500, y = max(plot_data_long$Error_Value, na.rm = TRUE) * 0.95,
           label = paste("Optimal Iterations:", best_iter), color = "blue", hjust = 0) +
  labs(
    title = "Gradient Boosting Model Performance",
    x = "Number of Trees (Iterations)",
    y = "MSE",
    color = "Error Type"
  ) +
  scale_color_manual(values = c("Training_Error" = "black", "CV_Error" = "red"),
                     labels = c("Training Error", "Cross-Validation Error")) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10),
    legend.position = "bottom"
  )

print(performance_plot)

vip <- summary(gbm_model, plotit = FALSE)

vip_plot <- ggplot(vip, aes(x = reorder(var, rel.inf), y = rel.inf)) +
  geom_col(fill = "steelblue", width = 0.7) +
  coord_flip() +
  labs(title = "Variable Importance",
       x = "Predictors",
       y = "Relative Influence (%)") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
        panel.grid.major.y = element_blank())

print(vip_plot)

results_df <- data.frame(
  Date = test_data$Date,
  Actual = test_data$Unemployment_Rate,
  Predicted = predictions
)

avp_plot <- ggplot(results_df, aes(x = Date)) +
  geom_point(aes(y = Actual, color = "Actual")) +
  geom_line(aes(y = Predicted, color = "Predicted"), linewidth = 0.8, linetype = "dashed") +
  labs(title = "Gradient Boosting",
       y = "Unemployment Rate (%)",
       color = "Legend") +
  scale_color_manual(values = c("Actual" = "red", "Predicted" = "darkblue")) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
        legend.position = "bottom")

print(avp_plot)


# Calculate residuals
residuals <- results_df$Actual - results_df$Predicted

results_df$Residuals <- residuals

residuals_plot <- ggplot(results_df, aes(x = Date, y = Residuals)) +
  geom_line(alpha = 0.6, linewidth = 0.6, color = "darkblue") +
  geom_hline(yintercept = 0, color = "red", linetype = "dashed") +
  labs(title = "Residuals Plot (GBM Model)",
       x = "Date",
       y = "Residuals") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 16))

print(residuals_plot)

acf_residuals <- acf(residuals, plot = FALSE, lag.max = 20)

acf_plot_data <- with(acf_residuals, data.frame(lag, acf))

acf_plot <- ggplot(acf_plot_data, aes(x = lag, y = acf)) +
  geom_bar(stat = "identity", position = "identity", fill = "black", width = 0.1) +
  geom_hline(yintercept = 0, color = "black", linetype = "solid") +
  geom_hline(yintercept = qnorm(0.975) / sqrt(length(residuals)), color = "blue", linetype = "dashed") +
  geom_hline(yintercept = -qnorm(0.975) / sqrt(length(residuals)), color = "blue", linetype = "dashed") +
  labs(title = "Gradient Boosting",
       x = "Lag",
       y = "Autocorrelation") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 16))

print(acf_plot)

pacf_residuals <- pacf(residuals, plot = FALSE, lag.max = 20)
pacf_plot_data <- with(pacf_residuals, data.frame(lag, acf))

pacf_plot <- ggplot(pacf_plot_data, aes(x = lag, y = acf)) +
  geom_bar(stat = "identity", position = "identity", fill = "black", width = 0.1) +
  geom_hline(yintercept = 0, color = "black", linetype = "solid") +
  # Significance bounds for PACF are also approx 2/sqrt(N)
  geom_hline(yintercept = qnorm(0.975) / sqrt(length(residuals)), color = "blue", linetype = "dashed") +
  geom_hline(yintercept = -qnorm(0.975) / sqrt(length(residuals)), color = "blue", linetype = "dashed") +
  labs(title = "Gradient Boosting",
       x = "Lag",
       y = "Partial Autocorrelation") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 16))

print(pacf_plot)

p_gbm_residuals_dist <- ggplot(results_df, aes(x = Residuals)) + # Use results_df and Residuals column
  geom_histogram(aes(y = after_stat(density)), binwidth = 0.2, fill = "lightgreen", color = "darkgreen", alpha = 0.7) +
  geom_density(aes(color = "Empirical Density"), linewidth = 1) +
  stat_function(fun = dnorm, args = list(mean = mean(results_df$Residuals, na.rm = TRUE), sd = sd(results_df$Residuals, na.rm = TRUE)), 
                aes(color = "Normal Distribution"), linetype = "dashed", linewidth = 1) +
  scale_color_manual(name = "Distribution Type", values = c("Empirical Density" = "darkblue", "Normal Distribution" = "purple")) +
  labs(title = "Gradient Boosting (test set)", # Updated title
       x = "Residuals",
       y = "Density") +
  theme_minimal() + theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 16), legend.position = "bottom", legend.title = element_text(face = "bold"))
print(p_gbm_residuals_dist)

top_vars <- vip$var[1:4]

pd_plots <- list()
for (var in top_vars) {
  pd <- partial(gbm_model, pred.var = var, n.trees = best_iter, 
                train = train_data, rug = TRUE)
  pd_plots[[var]] <- autoplot(pd) + 
    labs(y = "Effect") +
    theme_minimal()
}

grid.arrange(grobs = pd_plots, ncol = 2)

# RMSE, AIC, BIC, MAPE

train_predictions <- predict(gbm_model, train_data, n.trees = best_iter)

# Calculate RMSE for training set
rmse_train_gbm <- calculate_rmse(train_data$Unemployment_Rate, train_predictions)
cat(paste0("Training RMSE: ", round(rmse_train_gbm, 4), "\n"))

# Calculate MAPE for training set
mape_train_gbm <- calculate_mape(train_data$Unemployment_Rate, train_predictions)
cat(paste0("Training MAPE: ", round(mape_train_gbm, 2), "%\n"))

# Calculate RMSE for test/forecast set
rmse_test_gbm <- calculate_rmse(test_data$Unemployment_Rate, predictions)
cat(paste0("Test/Forecast RMSE: ", round(rmse_test_gbm, 4), "\n"))

# Calculate MAPE for test/forecast set
mape_test_gbm <- calculate_mape(test_data$Unemployment_Rate, predictions)
cat(paste0("Test/Forecast MAPE: ", round(mape_test_gbm, 2), "%\n"))

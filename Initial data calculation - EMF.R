# Load R packages
install.packages("tidyverse")
library(tidyverse)

# Read data
data <- read.csv("restoration_metrics_clean_wide_v2.csv", header = TRUE)

# Select functional columns
func_vars <- c("agb", "forage_rfv", "plant_c", "plant_n",
               "root_c_stock", "soil_cn_ratio", "soil_moisture", "soil_total_c_stock")

# Extract functional data and calculate Z-scores
func_data <- data[func_vars]
func_z <- scale(func_data)   # Standardization (Z-score)

# Convert to data frame
func_z <- as.data.frame(func_z)

# Calculate EMF (mean of all functional Z-scores)
data$EMF <- apply(func_z, 1, mean, na.rm = TRUE)

# Save results
write.csv(data, "EMF计算结果.csv", row.names = FALSE)

# View the first few rows
head(data)
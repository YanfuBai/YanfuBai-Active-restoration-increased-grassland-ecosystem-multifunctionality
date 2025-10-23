# Load required packages
library(ppcor)
library(stats)
library(dplyr)
library(data.table)

# Read data
data <- fread("restoration_metrics_clean_wide_v2.csv", data.table = FALSE)

# Define diversity metrics and functional metrics
diversity_metrics <- c("bacteria_hill_q0", "fungi_hill_q0", "plant_hill_q0")
function_metrics <- c("agb", "forage_rfv", "plant_c", "plant_n", "root_c_stock", "soil_cn_ratio", "soil_moisture", "soil_total_c_stock", "emf")
control_variables <- "treatment_group"

# Ensure all variables are numeric
data[diversity_metrics] <- lapply(data[diversity_metrics], as.numeric)
data[function_metrics] <- lapply(data[function_metrics], as.numeric)
data[control_variables] <- lapply(data[control_variables], as.numeric)


# Initialize an empty data.frame to store partial correlation results
results <- data.frame(
  from = character(),
  to = character(),
  correlation = numeric(),
  p_value = numeric(),
  stringsAsFactors = FALSE
)

# Perform partial correlation analysis for each pair of diversity and function metrics
for (div_var in diversity_metrics) {
  for (func_var in function_metrics) {
    # Remove rows with NA for calculation
    temp_data <- na.omit(data[, c(div_var, func_var, control_variables)])
    
    # Check if there are enough complete data points
    if (nrow(temp_data) > 2) {
      # Perform partial correlation using pcor.test()
      pcor_result <- pcor.test(temp_data[[div_var]], temp_data[[func_var]], temp_data[[control_variables]])
      
      # Extract results and add to data.frame
      results <- rbind(results, data.frame(
        from = div_var,
        to = func_var,
        correlation = pcor_result$estimate,
        p_value = pcor_result$p.value
      ))
    }
  }
}

# Apply FDR adjustment to p-values
results$p_adj <- p.adjust(results$p_value, method = "fdr")

# --- Create node and edge data (including all results) ---
# Unique names for nodes, ensuring all diversity and functional metrics are included
nodes_names <- unique(c(diversity_metrics, function_metrics))

# Create node data
nodes <- data.frame(
  id = nodes_names,
  label = nodes_names,
  group = ifelse(nodes_names %in% diversity_metrics, "Diversity", "Function"),
  stringsAsFactors = FALSE
)

# Create edge data (including all results)
edges <- data.frame(
  from = results$from,
  to = results$to,
  correlation = results$correlation,
  p_value = results$p_value,
  p_adj = results$p_adj,
  weight = abs(results$correlation), # Use absolute value of correlation coefficient as weight
  type = ifelse(results$correlation > 0, "positive", "negative"),
  stringsAsFactors = FALSE
)

# Save node and edge data to CSV files
write.csv(nodes, "nodes.csv", row.names = FALSE)
write.csv(edges, "edges.csv", row.names = FALSE)

# Print results overview
cat("Partial correlation analysis completed. All results saved, a preview is below:\n")
print(head(results))

cat("\nNode data saved to 'nodes.csv'.\n")

cat("Edge data saved to 'edges.csv'.\n")

# ============================================================
# Convert bacteria.txt to CSV file as is
# ============================================================

# Read txt file (tab-separated, first column is OTU ID)
bacteria_raw <- read.table("bacterial raw data.txt",
                           header = TRUE,
                           sep = "\t",
                           check.names = FALSE,
                           quote = "",
                           comment.char = "")

# Write to csv file
write.csv(bacteria_raw, "bacteria.csv", row.names = FALSE, na = "")

# Notification
cat("Conversion complete, output file is bacteria.csv\n")

# ============================================================
# Calculate bacterial community Hill diversity indices (q=0,1,2)
# Input: bacteria.csv (OTU × plot abundance matrix)
# Output: bacteria_diversity.csv (includes plot_name, bacteria_q0, bacteria_q1, bacteria_q2)
# ============================================================

# ---- Install and load hillR package ----
if (!requireNamespace("hillR", quietly = TRUE)) {
  install.packages("hillR")
}
library(hillR)

# ---- 1. Read data ----
# The first column of bacteria.csv is OTU ID, subsequent columns are plot abundance
bacteria_raw <- read.csv("bacteria.csv",
                         header = TRUE,
                         row.names = 1,
                         check.names = FALSE)

# ---- 2. Transform matrix structure ----
# hill_taxa requires input: plot × OTU (rows = plots, columns = species abundance)
comm <- t(bacteria_raw)

# ---- 3. Calculate Hill numbers ----
hill_q0 <- hill_taxa(comm, q = 0)  # Species richness
hill_q1 <- hill_taxa(comm, q = 1)  # Shannon true diversity
hill_q2 <- hill_taxa(comm, q = 2)  # Simpson true diversity

# ---- 4. Organize result table ----
hill_df <- data.frame(
  plot_name = rownames(comm),
  bacteria_q0 = hill_q0,
  bacteria_q1 = hill_q1,
  bacteria_q2 = hill_q2,
  row.names = NULL
)

# ---- 5. Export results ----
write.csv(hill_df, "bacteria_diversity.csv", row.names = FALSE)

# Print first few rows for check
print(head(hill_df))
cat("Results saved as bacteria_diversity.csv\n")
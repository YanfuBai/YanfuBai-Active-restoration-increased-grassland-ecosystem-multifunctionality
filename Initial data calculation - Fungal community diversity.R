# ============================================================
# Convert fungi.txt to CSV file as is
# ============================================================

# Read txt file (tab-separated, first column is OTU ID)
fungi_raw <- read.table("fungal raw data.txt",
                        header = TRUE,
                        sep = "\t",
                        check.names = FALSE,
                        quote = "",
                        comment.char = "")

# Write to csv file
write.csv(fungi_raw, "fungi.csv", row.names = FALSE, na = "")

# Notification
cat("Conversion complete, output file is fungi.csv\n")

# ============================================================
# Calculate fungal community Hill diversity indices (q=0,1,2)
# Input: fungi.csv (OTU × plot abundance matrix)
# Output: fungi_diversity.csv (includes plot_name, fungi_q0, fungi_q1, fungi_q2)
# ============================================================

# ---- Install and load hillR package ----
if (!requireNamespace("hillR", quietly = TRUE)) {
  install.packages("hillR")
}
library(hillR)

# ---- 1. Read data ----
# The first column of fungi.csv is OTU ID, subsequent columns are plot abundance
fungi_raw <- read.csv("fungi.csv",
                      header = TRUE,
                      row.names = 1,
                      check.names = FALSE)

# ---- 2. Transform matrix structure ----
# hill_taxa requires input: plot × OTU (rows = plots, columns = species abundance)
comm <- t(fungi_raw)

# ---- 3. Calculate Hill numbers ----
hill_q0 <- hill_taxa(comm, q = 0)  # Species richness
hill_q1 <- hill_taxa(comm, q = 1)  # Shannon true diversity
hill_q2 <- hill_taxa(comm, q = 2)  # Simpson true diversity

# ---- 4. Organize result table ----
hill_df <- data.frame(
  plot_name = rownames(comm),
  fungi_q0 = hill_q0,
  fungi_q1 = hill_q1,
  fungi_q2 = hill_q2,
  row.names = NULL
)

# ---- 5. Export results ----
write.csv(hill_df, "fungi_diversity.csv", row.names = FALSE)

# Print first few rows for check
print(head(hill_df))
cat("Results saved as fungi_diversity.csv\n")
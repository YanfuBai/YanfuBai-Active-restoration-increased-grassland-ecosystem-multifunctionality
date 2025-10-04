# ===============================
# Complete code: Multiple columns → Multiple plots → Export to PPT
# ===============================

# Load packages
library(ggplot2)
library(officer)
library(rvg)

# 1) Import data (CSV file, first column as row names)
# Assume the file format is:
# name,DG,ABC,XYZ
# agb,0.35,0.22,-0.11
# forage_rfv,-0.12,0.05,0.33
# plant_c,0.48,-0.20,0.10
# ...
df <- read.csv("ΔMFRI data.csv", row.names = 1, check.names = FALSE)

# 2) Create a new PPT document
doc <- read_pptx()

# 3) Generate plots for each column and add them to the PPT
for (colname in colnames(df)) {
  
  # Extract data
  data <- data.frame(
    cell_type = rownames(df),
    geneexp   = df[[colname]]
  )
  
  # Plot 1
  p1 <- ggplot(data, aes(geneexp, reorder(cell_type, geneexp), fill = geneexp)) +
    geom_bar(stat = "identity", color = "black", linewidth = 0.8) +
    geom_text(
      aes(x = ifelse(geneexp > 0, -0.06, 0.06),
          label = round(geneexp, 3)),
      color = "#89ABE3FF"
    ) +
    scale_fill_gradient2(low = "#e3596d", mid = "white", high = "#47d4df",
                         midpoint = 0) +
    labs(title = paste("Barplot of", colname),
         x = "Correlation coefficient", y = "Cell Type") +
    theme_classic() +
    theme(legend.position = "none",
          axis.text = element_text(size = 15),
          axis.title.x = element_text(vjust = 0, size = 15, face = "plain"),
          axis.title.y = element_text(vjust = 0, size = 15, face = "plain"),
          plot.title  = element_text(size = 15, face = "bold"),
          axis.line   = element_line(color = "black", linewidth = 0.8))
  
  # Plot 2
  p2 <- ggplot(data, aes(geneexp, reorder(cell_type, geneexp), fill = geneexp)) +
    geom_bar(stat = "identity", color = "black", linewidth = 0.8) +
    geom_text(
      aes(x = ifelse(geneexp > 0, geneexp + 0.06, geneexp - 0.06),
          label = round(geneexp, 3)),
      color = "#89ABE3FF"
    ) +
    geom_text(aes(x = ifelse(geneexp > 0, -0.02, 0.02),
                  label = cell_type),
              color = "black", size = 5,
              hjust = ifelse(data$geneexp > 0, 1, 0)) +
    scale_fill_gradient2(low = "#e3596d", mid = "white", high = "#47d4df",
                         midpoint = 0) +
    labs(title = paste("Labeled Barplot of", colname),
         x = "Correlation coefficient", y = "") +
    theme(legend.position = "none",
          panel.grid.major.y = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks  = element_blank(),
          panel.border = element_blank(),
          axis.text   = element_text(size = 15),
          axis.title.x = element_text(vjust = 0, size = 15, face = "bold"))
  
  # Add both plots to PPT
  doc <- add_slide(doc, layout = "Title and Content", master = "Office Theme")
  doc <- ph_with(doc, dml(ggobj = p1), location = ph_location_fullsize())
  
  doc <- add_slide(doc, layout = "Title and Content", master = "Office Theme")
  doc <- ph_with(doc, dml(ggobj = p2), location = ph_location_fullsize())
}

# 4) Save PPT
print(doc, target = "All_Plots.pptx")

# ===============================
# After running, "All_Plots.pptx" will be generated
# Each numeric column will produce two plots (Plot 1 and Plot 2)
# and will be added to the PPT in order
# Gradient color: negative = red (#e3596d), 0 = white, positive = blue-green (#47d4df)
# ===============================

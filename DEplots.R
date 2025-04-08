library(Seurat)
library(ggplot2)
# anxa2a
# 0dpa: 9.889958e-06
# 1dpa: 0.001490988
# 2dpa: 2.590675e-05
# 3dpa: 0.2268929
# 7dpa: 0.5970043
bec_no_untreated <- subset(zf_filtered, idents = c("Endothelial Cell"), 
                     subset = timepoint %in% c("untreated") == FALSE)
# Extract data from Seurat object
bec_data <- FetchData(bec_no_untreated, vars = c("olfml3b", "timepoint"))

# Convert timepoint to factor for correct ordering
bec_data$timepoint <- factor(bec_data$timepoint, levels = c("mock", "0dpa", "1dpa", "2dpa", "3dpa", "7dpa"))

# Define manually provided p-values
p_values <- data.frame(
  timepoint = c("0dpa", "1dpa", "2dpa", "3dpa", "7dpa"),
  p_value = c(9.889958e-06, 0.001490988, 2.590675e-05, 0.2268929, 0.5970043)
)

# Define y-position for significance bars
y_max <- max(bec_data$anxa2a) * 1.6  # Slightly above the highest violin plot
p_values$y_position <- seq(y_max, y_max * 0.8, length.out = nrow(p_values))  # Staggered heights

# Define starting & ending positions for significance bars
p_values$x_start <- "mock"
p_values$x_end <- p_values$timepoint

# Plot violin plot
ggplot(bec_data, aes(x = timepoint, y = olfml3b, fill = timepoint)) +
  geom_violin(scale = "width", trim = FALSE) +
  geom_boxplot(width = 0.1, outlier.shape = NA) +  # Add boxplot for clarity
  scale_fill_manual(values = c("mock" = "gray", "0dpa" = "#ccd5ae", "1dpa" = "#e9c46a", "2dpa" = "#2a9d8f", "3dpa" = "#e63946", "7dpa" = "#eab69f")) +
  labs(title = "olfml3b Expression Across Timepoints in Endothelial Cells", x = "Timepoint", y = "Expression") +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) #+
  # Add horizontal lines connecting "mock" to each timepoint
  #geom_segment(data = p_values, aes(x = x_start, xend = x_end, y = y_position, yend = y_position), size = 0.5) +
  # Add p-value labels above the lines
  #geom_text(data = p_values, aes(x = x_end, y = y_position + 0.05 * y_max, label = paste0("p = ", signif(p_value, 3))), size = 4)


# Plot violin plot
ggplot(bec_data, aes(x = timepoint, y = anxa2a, fill = timepoint)) +
  geom_violin(scale = "width", trim = FALSE) +
  geom_boxplot(width = 0.1, outlier.shape = NA) +  # Add boxplot for clarity
  scale_fill_manual(values = c("mock" = "gray", "0dpa" = "#ccd5ae", "1dpa" = "#e9c46a", "2dpa" = "#2a9d8f", "3dpa" = "#e63946", "7dpa" = "#eab69f")) +
  labs(title = "anxa2a Expression Across Timepoints in Endothelial Cells", x = "Timepoint", y = "Expression") +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) 

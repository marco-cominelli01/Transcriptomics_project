library(reshape2)
library(ggplot2)
library(gridExtra)


# Converti i dati in formato long, ma mantieni il dataset originale intatto
logcpm_before_melted <- data.frame(value = as.vector(as.matrix(logcpm_before)),
                                   group = rep(colnames(logcpm_before), each = nrow(logcpm_before)))
logcpm_after_melted <- data.frame(value = as.vector(as.matrix(logcpm_after)),
                                  group = rep(colnames(logcpm_after), each = nrow(logcpm_after)))

# Creare il primo boxplot per logcpm_before
p1 <- ggplot(logcpm_before_melted, aes(x = group, y = value, fill = group)) +
  geom_boxplot(notch = TRUE) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 20),
    axis.title.x = element_text(size = 15),
    axis.title.y = element_text(size = 20),
    legend.position = "none"
  ) +
  labs(x = "Sample", y = expression(log[2]*CPM))

# Creare il secondo boxplot per logcpm_after
p2 <- ggplot(logcpm_after_melted, aes(x = group, y = value, fill = group)) +
  geom_boxplot(notch = TRUE) +
  theme_minimal() +
  theme(
    axis.title.x = element_text(size = 15),
    legend.title = element_text(size = 15),
    legend.text = element_text(size = 12)
  ) +
  labs(x = "Sample") +
  scale_fill_discrete(name = "Group")

# Visualizzare i due grafici uno accanto all'altro
grid.arrange(p1, p2, ncol = 2)



#####################################


# Caricare le librerie necessarie
library(ggplot2)
library(SummarizedExperiment)

# Supponiamo che '31616' sia l'ID del gene di interesse
gene_id <- 31616

# Estrarre i dati per il gene specifico
brain_data <- assays(rse_brain)$TPM[gene_id,]
colon_data <- assays(rse_colon)$TPM[gene_id,]
muscle_data <- assays(rse_muscle)$TPM[gene_id,]

# Creare un data frame lungo per ggplot2
data <- data.frame(
  Expression = c(brain_data, colon_data, muscle_data),
  Tissue = factor(rep(c("Brain", "Colon", "Muscle"), c(length(brain_data), length(colon_data), length(muscle_data))))
)

# Creare il boxplot colorato con ggplot2
ggplot(data, aes(x = Tissue, y = Expression, fill = Tissue)) +
  geom_boxplot(outlier.shape = NA) +  # Per nascondere gli outlier
  labs(title = "Gene Expression across Tissues", x = "Tissue", y = "Expression (TPM)") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 14),
    axis.text.y = element_text(size = 14),
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 16),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12),
    plot.title = element_text(size = 18, hjust = 0.5)
  ) +
  scale_fill_brewer(palette = "Set3")  # Usare una palette di colori predefinita


### Barplot con barra sopra per geni upregolati e negativa per downregolati

library(ggplot2)
library(gridExtra)
library(grid)

# Create data frame
data <- data.frame(
  Comparison = c("Colon vs Brain", "Muscle vs Brain", "Muscle vs Colon"),
  Down = c(2406, 4376, 3180),
  Up = c(2026, 3888, 2889)
)

# Function to extract legend
get_legend <- function(myplot) {
  tmp <- ggplotGrob(myplot)
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}

# Create a list to store the plots
plots <- list()

# Set the same y-axis limits for all plots
y_limits <- c(-4500, 4500)

for (i in 1:nrow(data)) {
  comparison <- data$Comparison[i]
  down <- data$Down[i]
  up <- data$Up[i]
  
  # Prepare data for ggplot
  plot_data <- data.frame(
    Category = factor(c("Down", "Up"), levels = c("Down", "Up")),
    Genes = c(-down, up)
  )
  
  # Create the bar plot
  p <- ggplot(plot_data, aes(x = Category, y = Genes, fill = Category)) +
    geom_bar(stat = "identity", width = 0.5) +
    geom_text(aes(label = abs(Genes)), vjust = ifelse(plot_data$Genes < 0, 1.5, -0.5), size = 5) +  # Increase text size to 5
    geom_hline(yintercept = 0, linetype = "dashed") +
    coord_cartesian(ylim = y_limits) +  # Set the same y-axis limits
    scale_y_continuous(breaks = seq(-4500, 4500, by = 500),
                       labels = function(x) abs(x)) +  # Show the numbers as positive
    labs(title = comparison, x = "", y = if (i == 1) "# of genes" else NULL) +  # Change y-axis label
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5, size = 16),
          axis.text.y = if (i == 1) element_text(size = 10) else element_blank(),  # Hide y-axis text for central and right plots
          axis.ticks.y = if (i == 1) element_line(size = 0.5) else element_blank(),  # Hide y-axis ticks for central and right plots
          axis.title.y = if (i == 1) element_text(size = 18) else element_blank(),  # Increase y-axis label size for the first plot
          legend.position = "none",
          panel.grid = element_blank(),  # Remove grid background
          panel.background = element_blank()) +  # Ensure background is white
    scale_fill_manual(values = c("Down" = "blue", "Up" = "red"), 
                      labels = c("Down" = "Down", "Up" = "Up"))
  
  # Add plot to the list
  plots[[i]] <- p
}

# Extract the legend from the last plot
legend <- get_legend(
  ggplot(data.frame(Category = factor(c("Down", "Up"), levels = c("Down", "Up")),
                    Genes = c(-data$Down[3], data$Up[3])),
         aes(x = Category, y = Genes, fill = Category)) +
    geom_bar(stat = "identity", width = 0.5) +
    scale_fill_manual(values = c("Down" = "blue", "Up" = "red"), 
                      labels = c("Down" = "Down", "Up" = "Up")) +
    theme_minimal() +
    theme(legend.position = "right",
          legend.text = element_text(size = 12),
          legend.title = element_blank(),
          panel.grid = element_blank(),  # Remove grid background
          panel.background = element_blank())  # Ensure background is white
)

# Combine the plots with the legend
combined_plot <- grid.arrange(
  grobs = c(plots, list(legend)),
  ncol = 4,
  widths = c(1, 1, 1, 0.3)
)

# Plot the combined plot
grid.newpage()
grid.draw(combined_plot)

########
#######

library(ggplot2)

# Create data frame
data <- data.frame(
  Comparison = c("Colon vs others", "Muscle vs others", "Brain vs others"),
  Height = c(890, 2067, 1823)
)

# Create the bar plot
p <- ggplot(data, aes(x = Comparison, y = Height, fill = Comparison)) +
  geom_bar(stat = "identity", width = 0.5) +
  geom_text(aes(label = Height), vjust = -0.5, size = 5) +
  geom_hline(yintercept = 0, linetype = "dashed") +  # Add dashed line at y = 0
  scale_y_continuous(breaks = seq(0, max(data$Height) + 300, by = 300), limits = c(0, max(data$Height) + 300)) +
  labs(x = "", y = "# of upregulated genes") +  # Change y-axis label
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
    axis.title.y = element_text(size = 18),  # Increase y-axis label size
    legend.position = "none",
    panel.grid = element_blank(),  # Remove grid background
    panel.background = element_blank()  # Ensure background is white
  ) +
  scale_fill_manual(values = c("Colon vs others" = "skyblue", "Muscle vs others" = "tomato", "Brain vs others" = "gold"))

# Plot the bar plot
print(p)

##############
#single-cell #
##############









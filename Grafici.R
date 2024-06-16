library(reshape2)

# Converti i dati in formato long, ma mantieni il dataset originale intatto
logcpm_before_melted <- data.frame(value = as.vector(as.matrix(logcpm_before)),
                                   group = rep(colnames(logcpm_before), each = nrow(logcpm_before)))
logcpm_after_melted <- data.frame(value = as.vector(as.matrix(logcpm_after)),
                                  group = rep(colnames(logcpm_after), each = nrow(logcpm_after)))

par(mfrow = c(1,2))
x11()

# Creare boxplot per logcpm_before
ggplot(logcpm_before_melted, aes(x = group, y = value, fill = group)) +
  geom_boxplot(notch = TRUE) +
  theme_minimal() +
  labs(title = "Boxplot of Log CPM before normalization", x = "Sample", y = "Log CPM")

# Creare boxplot per logcpm_after
ggplot(logcpm_after_melted, aes(x = group, y = value, fill = group)) +
  geom_boxplot(notch = TRUE) +
  theme_minimal() +
  labs(title = "Boxplot of logcpm_after", x = "Sample", y = "Log CPM")


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
  scale_fill_brewer(palette = "Set3")  # Usare una palette di colori predefinita



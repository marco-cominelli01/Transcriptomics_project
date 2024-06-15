library(reshape2)

# Converti i dati in formato long, ma mantieni il dataset originale intatto
logcpm_before_melted <- data.frame(value = as.vector(as.matrix(logcpm_before)),
                                   group = rep(colnames(logcpm_before), each = nrow(logcpm_before)))
logcpm_after_melted <- data.frame(value = as.vector(as.matrix(logcpm_after)),
                                  group = rep(colnames(logcpm_after), each = nrow(logcpm_after)))

# Creare boxplot per logcpm_before
ggplot(logcpm_before_melted, aes(x = group, y = value, fill = group)) +
  geom_boxplot(notch = TRUE) +
  theme_minimal() +
  labs(title = "Boxplot of logcpm_before", x = "Group", y = "Value")

# Creare boxplot per logcpm_after
ggplot(logcpm_after_melted, aes(x = group, y = value, fill = group)) +
  geom_boxplot(notch = TRUE) +
  theme_minimal() +
  labs(title = "Boxplot of logcpm_after", x = "Group", y = "Value")

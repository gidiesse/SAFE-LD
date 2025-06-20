library(dplyr)

setwd("/scratch/giulia.desanctis/susie_results/result_spreadsheets/more_phenos_than_indivs/")

file_path <- "/scratch/giulia.desanctis/susie_results/5000_individuals/model_%d/susie_metrics_%d_more_phenos_than_indiv.rds"
df <- data.frame()
for ( i in 1:10) {
  new_data <- readRDS(sprintf(file_path, i, i))
  new_data$model_number <- i
  new_data$effect_size_1 <- c(0.005)
  new_data$effect_size_2 <- c(0.05)
  new_data$effect_size_3 <- c(0.2)
  new_data <- new_data[, c(1, 8, 9, 10, 11, 2, 3, 4, 5, 6, 7)]
  df <- rbind(df, new_data)
}
df_10k <- df
View(df)
write.csv2(df, file = "metrics_5k_individuals_5k_6k_7k_8k_9k_10k_phenos.csv")

## Do graphs

library(ggplot2)
library(patchwork)

# If your data frame is called df (adapt the name if needed):
# Example:
df <- read.csv2("metrics_10k_individuals_1k_2k_3k_4k_5k_6k_7k_8k_9k_10k_phenos.csv")

# Plot 1: Power vs model_number
p1 <- ggplot(df, aes(x = model_number, y = power, color = model_name, group = model_name)) +
  geom_line(size = 1) +
  geom_point(size = 2) +
  labs(title = "Power across model_number 10k individuals",
       x = "Model Number",
       y = "Power") +
  theme_minimal() +
  theme(text = element_text(size = 14),
        legend.title = element_blank(),
        legend.text = element_text(size = 8),       
        legend.key.size = unit(0.4, "cm")) +
  scale_x_continuous(breaks = unique(df$model_number)) +
  ylim(0.65,1)  # Assuming power is between 0 and 1

# Plot 2: Coverage vs model_number
p2 <- ggplot(df, aes(x = model_number, y = coverage, color = model_name, group = model_name)) +
  geom_line(size = 1) +
  geom_point(size = 2) +
  labs(title = "Coverage across model_number 10k individuals",
       x = "Model Number",
       y = "Coverage") +
  theme_minimal() +
  theme(text = element_text(size = 14),
        legend.title = element_blank(), 
        legend.text = element_text(size = 8),       
        legend.key.size = unit(0.4, "cm"))  +
  scale_x_continuous(breaks = unique(df$model_number)) +
  ylim(0.3,1)  # Assuming coverage is between 0 and 1

p1 + p2

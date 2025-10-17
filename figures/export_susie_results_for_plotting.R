library(dplyr)

setwd("/scratch/giulia.desanctis/susie_results_method_1/result_spreadsheets/")

## METHOD 1 ie SAFE-LD
file_path <- "/scratch/giulia.desanctis/susie_results_method_1/10000_individuals/model_%d/susie_metrics_%d.rds"
df <- data.frame()
for ( i in 1:10) {
  file_existance <- file.exists(sprintf(file_path, i, i))
  if (file_existance) {
    new_data <- readRDS(sprintf(file_path, i, i))
    new_data$model_number <- i
    new_data <- new_data[, c(1, 8, 2, 3, 4, 5, 6, 7)]
    df <- rbind(df, new_data) 
  } else {
    warning(sprintf("File not found: model %d", i))
  }
}
#View(df)
write.csv2(df, file = sprintf("metrics_10k_individuals.csv"))


## METHOD 2 ie SAFE-LDss
file_path <- "/scratch/giulia.desanctis/susie_results_method_2/10000_individuals/model_%d/%s_perc_phenos_with_effect/effect_%.2f/susie_metrics_%d_%.2f_effect.rds"
effect_size <- 0.01
perc_effect <- "10"
df <- data.frame()
for ( i in 1:10) {
  file_existance <- file.exists(sprintf(file_path, i, perc_effect, effect_size, i, effect_size))
  if (file_existance) {
    new_data <- readRDS(sprintf(file_path, i, perc_effect, effect_size, i, effect_size))
    new_data$model_number <- i
    new_data$effect_size <- effect_size
    new_data$perc_phenos_with_effect <- perc_effect
    new_data$effect_size_1 <- c(0.005)
    new_data$effect_size_2 <- c(0.05)
    new_data$effect_size_3 <- c(0.2)
    new_data <- new_data[, c(1, 8, 9, 10, 11, 12, 13, 2, 3, 4, 5, 6, 7)]
    df <- rbind(df, new_data) 
  } else {
    warning(sprintf("File not found: model %d", i))
  }
}
#View(df)
write.csv2(df, file = sprintf("metrics_10k_individuals_%s_perc_%.2f_effect.csv", perc_effect, effect_size))

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

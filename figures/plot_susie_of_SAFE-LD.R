# R script to load, tidy, summarize, and plot SAFE-LD simulation grid

setwd("/scratch/giulia.desanctis/susie_results_method_1/result_spreadsheets/greater_n_than_p/")

# 0. Install/load required packages
# install.packages(c("dplyr","stringr","ggplot2","tidyr"))
library(dplyr)
library(stringr)
library(tidyr)
library(ggplot2)

# 1. List & read all simulation CSV files in the working directory
files <- list.files(path = ".", pattern = "^metrics_.*\\.csv$", full.names = TRUE)
df_list <- lapply(files, read.table, sep=";", header=T, stringsAsFactors = FALSE)
df <- bind_rows(df_list)

# 2. Convert coverage & power from strings with commas to numeric
df <- df %>%
  mutate(
    coverage = as.numeric(str_replace(coverage, ",", ".")),
    power    = as.numeric(str_replace(power,    ",", "."))
  )

# 3. Flag & parse the simulated_LD rows: extract n_phenos & n_individuals
df <- df %>%
  mutate(
    is_sim = str_detect(model_name, "^simulated_LD_"),
    n_phenos = if_else(
      is_sim,
      as.integer(str_extract(model_name, "(?<=simulated_LD_)\\d+(?=_phenos)")),
      NA_integer_
    ),
    n_individuals = if_else(
      is_sim,
      as.integer(str_extract(model_name, "(?<=_phenos_)\\d+(?=_individuals)")),
      NA_integer_
    )
  )

# 4. Forward-fill n_individuals within each model_number so internal/external inherit it
df <- df %>%
  group_by(model_number) %>%
  fill(n_individuals, .direction = "up") %>%
  ungroup()

# 5. Build a lookup of all (model_number, n_individuals, n_phenos) from the simulated rows
sim_lookup <- df %>%
  filter(is_sim) %>%
  select(model_number, n_individuals, n_phenos) %>%
  distinct()

# 6. Extract internal_LD & external_LD, drop their NA n_phenos, then replicate them
df_int_ext <- df %>%
  filter(model_name %in% c("internal_LD", "external_LD")) %>%
  select(-n_phenos)

df_int_ext_rep <- df_int_ext %>%
  inner_join(sim_lookup,
             by = c("model_number", "n_individuals"),
             relationship = "many-to-many")

# 7. Rename simulated_LD → SAFE-LD and select columns
df_sim <- df %>%
  filter(is_sim) %>%
  mutate(model_name = "SAFE-LD") %>%
  select(model_name, model_number, n_phenos, n_individuals, coverage, power)

# 8. Combine SAFE-LD with replicated internal/external rows
df_int_ext_final <- df_int_ext_rep %>%
  select(model_name, model_number, n_phenos, n_individuals, coverage, power)

df_combined <- bind_rows(df_sim, df_int_ext_final)

# 9. Report missing counts by model_name
missing_summary <- df_combined %>%
  group_by(model_name) %>%
  summarize(
    missing_coverage = sum(is.na(coverage)),
    missing_power    = sum(is.na(power)),
    total_rows       = length(model_name),
    .groups = "drop"
  )
print(missing_summary)

# 10. Compute mean coverage & power over the 10 replicates for each (model, N, P)
summary_df <- df_combined %>%
  dplyr::group_by(model_name, n_individuals, n_phenos) %>%
  dplyr::summarize(
    mean_coverage = mean(coverage, na.rm = TRUE),
    mean_power    = mean(power,    na.rm = TRUE),
    .groups = "drop"
  )

# 11. Save the summary table
write.csv(summary_df, "summary_df.csv", row.names = FALSE)

# 12. Plot Coverage curves
p_cov <- ggplot(summary_df, aes(x = n_phenos, y = mean_coverage, color = model_name)) +
  geom_line(size = 1) +
  geom_point(size = 2) +
  scale_x_continuous(
    breaks = seq(1000, 10000, by = 1000),
    limits = c(1000, 10000),
    labels = function(x) paste0(x / 1000, "k")
  ) +
  scale_y_continuous(
    limits = c(0.6, 1),
    expand = expansion(mult = c(0.01, 0))  # optional: remove top padding
  ) +
  facet_wrap(~ n_individuals, scales = "free_y", labeller = label_both) +
  labs(
    x = "Number of Phenotypes",
    y = "Mean Coverage",
    color = "Method",
    title = "Coverage vs. Number of Phenotypes 0.20 effect"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    strip.text = element_text(face = "bold", size = 14),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 18),
    legend.title = element_text(face = "bold"),
    legend.position = "bottom",
    panel.spacing = unit(1.2, "lines"),
    axis.title = element_text(face = "bold"),
    axis.text = element_text(size = 11),
    axis.ticks.length = unit(0.25, "cm"),
    axis.title.x = element_text(face = "bold", margin = margin(t = 10)),
    axis.title.y = element_text(face = "bold", margin = margin(r = 10)),
  ) +
  scale_color_manual(values = c("external_LD" = "#F8766D", "internal_LD" = "#00BA38", "SAFE-LD" = "#619CFF"))


# 13. Plot Power curves
p_pow <- ggplot(summary_df, aes(x = n_phenos, y = mean_power, color = model_name)) +
  geom_line(size = 1) +
  geom_point(size = 2) +
  scale_x_continuous(
    breaks = seq(1000, 10000, by = 1000),
    limits = c(1000, 10000),
    labels = function(x) paste0(x / 1000, "k")
  ) +
  scale_y_continuous(
    limits = c(0.6, 1),
    expand = expansion(mult = c(0.01, 0))  # optional: remove top padding
  ) +
  facet_wrap(~ n_individuals, scales = "free_y", labeller = label_both) +
  labs(
    x = "Number of Phenotypes",
    y = "Mean Power",
    color = "Method",
    title = "Power vs. Number of Phenotypes 0.20 effect"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    strip.text = element_text(face = "bold", size = 14),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 18),
    legend.title = element_text(face = "bold"),
    legend.position = "bottom",
    panel.spacing = unit(1.2, "lines"),
    axis.title = element_text(face = "bold"),
    axis.text = element_text(size = 11),
    axis.ticks.length = unit(0.25, "cm"),
    axis.title.x = element_text(face = "bold", margin = margin(t = 10)),
    axis.title.y = element_text(face = "bold", margin = margin(r = 10)),
    
  ) +
  scale_color_manual(values = c("external_LD" = "#F8766D", "internal_LD" = "#00BA38", "SAFE-LD" = "#619CFF"))



# 14. Print the plots to the active graphics device
print(p_cov)
print(p_pow)

# 15. Put power and coverage on the same graph
df_long <- summary_df %>%
  pivot_longer(cols = c(mean_power, mean_coverage),
               names_to = "Metric",
               values_to = "Value") %>%
  mutate(Metric = recode(Metric,
                         mean_power = "Power",
                         mean_coverage = "Coverage"))

p_together <- ggplot(df_long, aes(x = n_phenos, y = Value, color = model_name, linetype = Metric)) +
  geom_line(size = 1.5) +
  geom_point(size = 4) +
  scale_x_continuous(
    breaks = seq(1000, 10000, by = 1000),
    limits = c(1000, 10000),
    labels = function(x) paste0(x / 1000, "k")
  ) +
  scale_y_continuous(limits = c(0.3, 1)) +
  facet_wrap(~ n_individuals, labeller = label_both) +
  #facet_wrap(~ n_individuals, labeller = label_both, scales = "free_x") + 
  labs(
    title = "Fine-mapping performance when N >= P",
    x = "Number of Phenotypes",
    y = "Value",
    color = "Method",
    linetype = "Metric"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 28),
    legend.title = element_text(face = "bold", size = 22),
    legend.text = element_text(face = "bold", size = 22),
    legend.position = "bottom",
    axis.title = element_text(face = "bold", size = 22),
    axis.text = element_text(face = "bold", size = 18),
    strip.text = element_text(face = "bold", size = 20),
    theme(axis.title.x = element_text(margin = margin(t = 20))),
  ) +
  scale_color_manual(values = c("external_LD" = "#F8766D",
                                "internal_LD" = "#00BA38",
                                "SAFE-LD" = "#619CFF"))

print(p_together)



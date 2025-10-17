## Script to create figure 1b.
## The true correlation matrix (LD) and simulated corrrelation matrix (SAFE-LD)
## should be provided as data frames in input.
## Number of figures, titles, labels etc. can all be personalised.

cor_mat_fp <- "/project/statgen/figures_for_poster/fig_1_nulls/%s_indiv/3000_phenos/100_iterations_3000_phenos_ukb_3k_phenos_basic_script/it_1/cor_mat.txt"
sim_mat_fp <- "/project/statgen/figures_for_poster/fig_1_nulls/%s_indiv/3000_phenos/100_iterations_3000_phenos_ukb_3k_phenos_basic_script/it_1/cor_zscore.txt"

cor_mat_200 <- read.table(sprintf(cor_mat_fp, "200"))
sim_mat_200 <- read.table(sprintf(sim_mat_fp, "200"))

# Sample 50k off-diagonal elements
off_diag_idx <- which(row(sim_mat_200) != col(sim_mat_200), arr.ind = TRUE)
sampled_idx <- off_diag_idx[sample(nrow(off_diag_idx), 50000), ]
sim_vals_200 <- sim_mat_200[sampled_idx]
cor_vals_200 <- cor_mat_200[sampled_idx]

# Plot data
df_200 <- data.frame(sim_vals_200, cor_vals_200)
p1 <- ggplot(df_200, aes(x = sim_vals_200, y = cor_vals_200)) +
  geom_point(size = 0.5, alpha = 0.8, color = "black") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
  labs(
    x = "Correlation of Z-scores",
    y = "Correlation of SNPs (True LD)",
    title = "Z-score correlation vs LD - 200 individuals"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10)
  ) +
  coord_fixed(xlim = c(-1, 1), ylim = c(-1, 1))

cor_mat_500 <- read.table(sprintf(cor_mat_fp, "500"))
sim_mat_500 <- read.table(sprintf(sim_mat_fp, "500"))

# Sample 50k off-diagonal elements
off_diag_idx <- which(row(sim_mat_500) != col(sim_mat_500), arr.ind = TRUE)
sampled_idx <- off_diag_idx[sample(nrow(off_diag_idx), 50000), ]
sim_vals_500 <- sim_mat_500[sampled_idx]
cor_vals_500 <- cor_mat_500[sampled_idx]

# Plot data
df_500 <- data.frame(sim_vals_500, cor_vals_500)
p2 <- ggplot(df_500, aes(x = sim_vals_500, y = cor_vals_500)) +
  geom_point(size = 0.5, alpha = 0.8, color = "black") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
  labs(
    x = "Correlation of Z-scores",
    y = "Correlation of SNPs (True LD)",
    title = "Z-score correlation vs LD - 500 individuals"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10)
  ) +
  coord_fixed(xlim = c(-1, 1), ylim = c(-1, 1))

cor_mat_1000 <- read.table(sprintf(cor_mat_fp, "1k"))
sim_mat_1000 <- read.table(sprintf(sim_mat_fp, "1k"))

# Sample 50k off-diagonal elements
off_diag_idx <- which(row(sim_mat_1000) != col(sim_mat_1000), arr.ind = TRUE)
sampled_idx <- off_diag_idx[sample(nrow(off_diag_idx), 50000), ]
sim_vals_1000 <- sim_mat_1000[sampled_idx]
cor_vals_1000 <- cor_mat_1000[sampled_idx]

# Plot data
df_1000 <- data.frame(sim_vals_1000, cor_vals_1000)
p3 <- ggplot(df_1000, aes(x = sim_vals_1000, y = cor_vals_1000)) +
  geom_point(size = 0.5, alpha = 0.8, color = "black") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
  labs(
    x = "Correlation of Z-scores",
    y = "Correlation of SNPs (True LD)",
    title = "Z-score correlation vs LD - 1000 individuals"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10)
  ) +
  coord_fixed(xlim = c(-1, 1), ylim = c(-1, 1))

cor_mat_5000 <- read.table(sprintf(cor_mat_fp, "5k"))
sim_mat_5000 <- read.table(sprintf(sim_mat_fp, "5k"))

# Sample 50k off-diagonal elements
off_diag_idx <- which(row(sim_mat_5000) != col(sim_mat_5000), arr.ind = TRUE)
sampled_idx <- off_diag_idx[sample(nrow(off_diag_idx), 50000), ]
sim_vals_5000 <- sim_mat_5000[sampled_idx]
cor_vals_5000 <- cor_mat_5000[sampled_idx]

# Plot data
df_5000 <- data.frame(sim_vals_5000, cor_vals_5000)
p4 <- ggplot(df_5000, aes(x = sim_vals_5000, y = cor_vals_5000)) +
  geom_point(size = 0.5, alpha = 0.8, color = "black") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
  labs(
    x = "Correlation of Z-scores",
    y = "Correlation of SNPs (True LD)",
    title = "Z-score correlation vs LD - 5000 individuals"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10)
  ) +
  coord_fixed(xlim = c(-1, 1), ylim = c(-1, 1))

cor_mat_30000 <- read.table(sprintf(cor_mat_fp, "30k"))
sim_mat_30000 <- read.table(sprintf(sim_mat_fp, "30k"))

# Sample 50k off-diagonal elements
off_diag_idx <- which(row(sim_mat_30000) != col(sim_mat_30000), arr.ind = TRUE)
sampled_idx <- off_diag_idx[sample(nrow(off_diag_idx), 50000), ]
sim_vals_30000 <- sim_mat_30000[sampled_idx]
cor_vals_30000 <- cor_mat_30000[sampled_idx]

# Plot data
df_30000 <- data.frame(sim_vals_30000, cor_vals_30000)
p5 <- ggplot(df_30000, aes(x = sim_vals_30000, y = cor_vals_30000)) +
  geom_point(size = 0.5, alpha = 0.8, color = "black") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
  labs(
    x = "Correlation of Z-scores",
    y = "Correlation of SNPs (True LD)",
    title = "Z-score correlation vs LD - 30000 individuals"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10)
  ) +
  coord_fixed(xlim = c(-1, 1), ylim = c(-1, 1))

library(patchwork)

# Create layout
common_theme <- theme(
  plot.margin = margin(6, 6, 6, 6),
  axis.title.x = element_text(margin = margin(t = 6)),
  axis.title.y = element_text(margin = margin(r = 6))
)

wrap_plots(
  p1 + common_theme, p2 + common_theme, p3 + common_theme,
  p4 + common_theme, p5 + common_theme, plot_spacer(),
  ncol = 3, byrow = TRUE
) +
  plot_layout(
    widths  = c(1, 1, 1),
    heights = c(1, 1),
    axes    = "collect",   # align axes = same label space
    guides  = "collect"    # (if any legends) keep them uniform
  ) +
  plot_annotation(
    title = "Z-score Correlation vs True LD across number of individuals",
    theme = theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold"))
  )

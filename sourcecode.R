library(dplyr)
library(ggplot2)
library(summarytools)
library(MASS)
library(fitdistrplus)
library(patchwork)
library(viridis)

df <- read.csv("deliveries.csv", stringsAsFactors = FALSE)

death_overs <- df %>%
  filter(!is.na(batsman),
         !is.na(batsman_runs),
         !is.na(over),
         over >= 16, over <= 20)

x <- death_overs$batsman_runs

mean_x <- mean(x, na.rm = TRUE)
var_x <- var(x, na.rm = TRUE)
sd_x <- sd(x, na.rm = TRUE)

skewness_x <- sum((x - mean_x)^3, na.rm = TRUE) / length(x) / (sd_x^3)

descriptive_stats <- data.frame(
  Statistic = c("Mean", "Variance", "SD", "Skewness"),
  Value = c(mean_x, var_x, sd_x, skewness_x)
)
print(descriptive_stats)

cat("\n=== Checking for Overdispersion ===\n")
cat("Mean:", round(mean_x, 2), "\n")
cat("Variance:", round(var_x, 2), "\n")
cat("Variance/Mean ratio:", round(var_x/mean_x, 2), "\n")
cat("If ratio >> 1, Negative Binomial is appropriate\n\n")

Q1 <- quantile(x, 0.25, na.rm = TRUE)
Q3 <- quantile(x, 0.75, na.rm = TRUE)
IQR_x <- Q3 - Q1
lower_bound <- Q1 - 1.5 * IQR_x
upper_bound <- Q3 + 1.5 * IQR_x

x_clean <- x[x >= lower_bound & x <= upper_bound]

cat("Outliers removed:", length(x) - length(x_clean), "out of", length(x), "\n")
cat("Clean data size:", length(x_clean), "\n\n")

p_raw <- ggplot(data.frame(x = x_clean), aes(x = x)) +
  geom_histogram(aes(y = after_stat(density)),
                 bins = 20,
                 fill = "skyblue",
                 color = "black",
                 alpha = 0.8) +
  geom_density(color = "red", linewidth = 1) +
  geom_vline(xintercept = mean(x_clean), color = "blue", linetype = "dashed", linewidth = 1) +
  annotate("text",
           x = mean(x_clean),
           y = Inf,
           label = paste0("Mean = ", round(mean(x_clean), 2)),
           vjust = 2,
           color = "blue",
           fontface = "bold") +
  annotate("text",
           x = min(x_clean),
           y = Inf,
           label = paste0("Var = ", round(var(x_clean), 2),
                          "\nSD = ", round(sd(x_clean), 2),
                          "\nSkew = ", round(skewness_x, 2)),
           hjust = 0,
           vjust = 2,
           color = "darkgreen",
           fontface = "bold") +
  labs(title = "Distribution of Batsman Runs in Death Overs",
       subtitle = "Empirical histogram with density curve (outliers removed)",
       x = "batsman_runs", y = "Density") +
  theme_minimal()

print(p_raw)
ggsave("raw_distribution.png", p_raw, width = 10, height = 6, dpi = 300)

cat("=== Fitting Negative Binomial Distribution ===\n")

fit_nbinom <- fitdist(x_clean, "nbinom")
print(fit_nbinom)

size_param <- fit_nbinom$estimate["size"]
prob_param <- fit_nbinom$estimate["prob"]

cat("\nNegative Binomial Parameters:\n")
cat("size (dispersion):", round(size_param, 4), "\n")
cat("prob (success probability):", round(prob_param, 4), "\n")

cat("\n=== Goodness-of-Fit Tests ===\n")
print(summary(fit_nbinom))

p_nbinom_dense <- denscomp(fit_nbinom, 
                           legendtext = "Negative Binomial", 
                           main = "Negative Binomial Fit for Death Overs Runs") + 
  theme_minimal()

p_nbinom_qq <- qqcomp(fit_nbinom, 
                      legendtext = "Negative Binomial",
                      main = "QQ Plot - Negative Binomial Fit") + 
  theme_minimal()

p_nbinom_pp <- ppcomp(fit_nbinom,
                      legendtext = "Negative Binomial",
                      main = "PP Plot - Negative Binomial Fit") + 
  theme_minimal()

print(p_nbinom_dense)
ggsave("nbinom_density_fit.png", p_nbinom_dense, width = 8, height = 6, dpi = 300)

print(p_nbinom_qq)
ggsave("nbinom_qq_plot.png", p_nbinom_qq, width = 8, height = 6, dpi = 300)

print(p_nbinom_pp)
ggsave("nbinom_pp_plot.png", p_nbinom_pp, width = 8, height = 6, dpi = 300)

clutch_sr <- death_overs %>%
  group_by(batsman) %>%
  summarise(
    balls_faced = n(),
    runs_scored = sum(batsman_runs, na.rm = TRUE),
    clutch_sr = 100 * runs_scored / balls_faced,
    .groups = "drop"
  ) %>%
  filter(balls_faced >= 30, !is.infinite(clutch_sr), !is.na(clutch_sr))

league_avg <- mean(clutch_sr$clutch_sr, na.rm = TRUE)
clutch_sr_sd <- sd(clutch_sr$clutch_sr, na.rm = TRUE)
shrink_prior <- 50

clutch_sr <- clutch_sr %>%
  mutate(
    weight = balls_faced / (balls_faced + shrink_prior),
    shrunk_sr = weight * clutch_sr + (1 - weight) * league_avg,
    cpi = (shrunk_sr - league_avg) / clutch_sr_sd
  )

cat("\n=== Top 10 Players by Shrunk Strike Rate ===\n")
print(head(clutch_sr %>% arrange(desc(shrunk_sr)), 10))

cat("\n=== Summary of Shrunk Strike Rate ===\n")
print(descr(clutch_sr$shrunk_sr,
            stats = c("mean", "sd", "min", "max", "n"),
            transpose = TRUE))

mean_shrunk <- mean(clutch_sr$shrunk_sr, na.rm = TRUE)
var_shrunk <- var(clutch_sr$shrunk_sr, na.rm = TRUE)
sd_shrunk <- sd(clutch_sr$shrunk_sr, na.rm = TRUE)

p_shrunk <- ggplot(clutch_sr, aes(x = shrunk_sr)) +
  geom_histogram(aes(y = after_stat(density)),
                 bins = 25,
                 fill = "steelblue",
                 alpha = 0.7,
                 color = "black") +
  geom_density(color = "red", linewidth = 1) +
  geom_vline(xintercept = mean_shrunk, color = "blue", linetype = "dashed", linewidth = 1) +
  annotate("text",
           x = mean_shrunk,
           y = Inf,
           label = paste0("Mean = ", round(mean_shrunk, 2)),
           vjust = 2,
           color = "blue",
           fontface = "bold") +
  annotate("text",
           x = min(clutch_sr$shrunk_sr),
           y = Inf,
           label = paste0("Var = ", round(var_shrunk, 2),
                          "\nSD = ", round(sd_shrunk, 2)),
           hjust = 0,
           vjust = 2,
           color = "darkgreen",
           fontface = "bold") +
  labs(title = "Distribution of Shrunk Clutch SR",
       subtitle = "Bayesian-adjusted performance in death overs",
       x = "shrunk_sr", y = "Density") +
  theme_minimal()

print(p_shrunk)
ggsave("shrunk_sr_distribution.png", p_shrunk, width = 10, height = 6, dpi = 300)

top_players <- clutch_sr %>%
  filter(balls_faced >= 100) %>%
  arrange(desc(cpi)) %>%
  slice_head(n = 100)

cat("\n=== Top Clutch Performers ===\n")
cat("TOP", nrow(top_players), "players generated!\n")

top20 <- top_players %>% slice_head(n = 20)

p_top20 <- ggplot(top20, aes(x = reorder(batsman, cpi), y = cpi)) +
  geom_col(fill = "darkgreen", alpha = 0.8) +
  coord_flip() +
  labs(title = "Top 20 Clutch Finishers",
       subtitle = "IPL Death Overs 16-20 | Min 100 balls",
       x = "Batsman", y = "CPI") +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 9, fontface = "bold"))

print(p_top20)
ggsave("slide1_top20.png", p_top20, width = 10, height = 8, dpi = 300)

top10 <- top_players %>% slice_head(n = 10)

p_top10 <- ggplot(top10, aes(x = reorder(batsman, cpi), y = cpi)) +
  geom_col(aes(fill = batsman), alpha = 0.9, width = 0.7, color = "black") +
  scale_fill_viridis_d(option = "plasma", direction = -1) +
  geom_text(aes(label = round(cpi, 2)), hjust = -0.05, fontface = "bold", size = 4, color = "white") +
  labs(title = "Top 10 Clutch Finishers (Bayesian CPI)",
       subtitle = "IPL Death Overs 16-20 | Min 100 balls",
       x = "Batsman", y = "CPI Score") +
  theme_minimal() +
  theme(legend.position = "none",
        axis.text.y = element_text(size = 12, fontface = "bold"),
        plot.title = element_text(size = 18, face = "bold"),
        plot.subtitle = element_text(size = 12)) +
  coord_flip()

print(p_top10)
ggsave("slide2_top10.png", p_top10, width = 12, height = 9, dpi = 300)

cat("\n=== Creating Final Multi-Panel Plot ===\n")

p1 <- p_raw
p2 <- p_shrunk
p3 <- p_nbinom_dense
p4 <- p_nbinom_qq
p5 <- p_top20

final_plot <- (p1 | p2) / (p3 | p4) / p5 +
  plot_annotation(
    title = "ClutchVision: Statistical Analysis of Death Overs Performance",
    subtitle = "IPL Death Overs (16-20): Negative Binomial Distribution, Bayesian Shrinkage & Top Clutch Performers",
    theme = theme(plot.title = element_text(size = 16, face = "bold"),
                  plot.subtitle = element_text(size = 12))
  )

print(final_plot)
ggsave("final_multi_panel_plot.png", final_plot, width = 14, height = 20, dpi = 300)

cat("\n=== Files Generated ===\n")
print(list.files(pattern = "*.png"))

cat("\n=== Top 5 Players Summary ===\n")
print(head(top_players, 5))

cat("\n=== Analysis Complete ===\n")
# Publication-style external validation overlay for the pMHC shortlist.

suppressPackageStartupMessages(library(ggplot2))

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- sub("^--file=", "", script_arg[[1]])
root <- normalizePath(file.path(dirname(script_path), ".."), mustWork = TRUE)

base_dir <- file.path(root, "processed", "external_validation_benchmark")
ann_path <- file.path(base_dir, "external_validation_pair_annotations.csv")
summary_path <- file.path(base_dir, "external_validation_rank_recovery_summary.csv")
out_dir <- file.path(root, "processed", "publication_figures")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

dat <- read.csv(ann_path, stringsAsFactors = FALSE)
summary <- read.csv(summary_path, stringsAsFactors = FALSE)
dat$rank <- as.numeric(dat$rank)
dat$review_priority_heuristic <- as.numeric(dat$review_priority_heuristic)

dat$plot_group <- ifelse(
  dat$pair_validation == "classic_BALF5_MBP_pair",
  "Classic BALF5--MBP pair",
  ifelse(
    dat$pair_validation == "drosu_2024_EBV_glycoprotein",
    "2024 EBV gB/gH positive",
    ifelse(dat$pair_validation == "classic_component_only", "Single classic component", "Background")
  )
)
dat$plot_group <- factor(
  dat$plot_group,
  levels = c("Classic BALF5--MBP pair", "2024 EBV gB/gH positive", "Single classic component", "Background")
)

label_row <- summary[summary$test == "classic_BALF5_MBP_pair_recovery", ][1, ]
new_row <- summary[summary$test == "strict_new_literature_overlay", ][1, ]
label_text <- sprintf(
  "Classic: %d/10 top-10, p < 0.001\n2024 gB/gH: %d/5 top-10, p = %.3f",
  label_row$observed_top10_positive_count,
  new_row$observed_top10_positive_count,
  new_row$empirical_p_top10_count_ge_observed
)

theme_pub <- theme_classic(base_size = 10, base_family = "Arial") +
  theme(
    axis.line = element_line(colour = "#333333", linewidth = 0.45),
    axis.ticks = element_line(colour = "#333333", linewidth = 0.45),
    axis.ticks.length = grid::unit(2.2, "pt"),
    axis.text = element_text(colour = "#333333"),
    axis.title = element_text(colour = "#333333"),
    plot.title = element_text(colour = "#333333", face = "bold", size = 11, hjust = 0),
    plot.subtitle = element_text(colour = "#333333", size = 8.2, hjust = 0),
    plot.caption = element_text(colour = "#555555", size = 7, hjust = 0),
    legend.title = element_blank(),
    legend.position = "bottom",
    legend.text = element_text(size = 7.2, colour = "#333333"),
    legend.key.size = grid::unit(0.15, "in"),
    legend.spacing.x = grid::unit(0.08, "in"),
    legend.margin = margin(0, 0, 0, 0),
    legend.box.margin = margin(-4, 0, -4, 0),
    panel.grid = element_blank(),
    plot.margin = margin(7, 7, 7, 7)
  )

p <- ggplot(dat, aes(x = rank, y = review_priority_heuristic, fill = plot_group)) +
  annotate("rect", xmin = 0.5, xmax = 10.5, ymin = -Inf, ymax = Inf,
           fill = "#F2F2F2", colour = NA) +
  annotate("text", x = 5.5, y = 0.70, label = "top 10",
           family = "Arial", size = 2.8, colour = "#555555") +
  geom_point(shape = 21, size = 3.2, stroke = 0.45, colour = "#333333", alpha = 0.95) +
  annotate("label", x = 17.8, y = 0.59, label = label_text,
           family = "Arial", size = 2.55, hjust = 0,
           fill = "white", colour = "#333333", linewidth = 0.25) +
  scale_fill_manual(values = c(
    "Classic BALF5--MBP pair" = "#009E73",
    "2024 EBV gB/gH positive" = "#0072B2",
    "Single classic component" = "#E69F00",
    "Background" = "#BDBDBD"
  )) +
  scale_x_continuous(breaks = c(1, 5, 10, 15, 20, 25, 30), limits = c(0.5, 34.5),
                     expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, 0.73), expand = c(0, 0)) +
  labs(
    title = "External validation overlay",
    subtitle = "Literature-positive candidates projected onto the pMHC priority shortlist",
    x = "Priority rank",
    y = "Candidate-priority score",
    caption = "Rank-recovery benchmark only; not TCR binding, activation, affinity, or patient mechanism evidence."
  ) +
  guides(fill = guide_legend(nrow = 2, byrow = TRUE)) +
  theme_pub

png_path <- file.path(out_dir, "figure_5_external_validation_overlay_300dpi.png")
pdf_path <- file.path(out_dir, "figure_5_external_validation_overlay.pdf")
ggsave(png_path, p, width = 5, height = 4, units = "in", dpi = 300,
       device = ragg::agg_png, bg = "white")
img <- png::readPNG(png_path)
pdf(pdf_path, width = 5, height = 4, useDingbats = FALSE)
par(mar = c(0, 0, 0, 0)); plot.new(); rasterImage(img, 0, 0, 1, 1)
dev.off()

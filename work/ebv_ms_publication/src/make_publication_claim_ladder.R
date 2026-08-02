# Publication claim/evidence ladder figure.
# This figure is generated from the claim matrix to keep manuscript claims tied
# to explicit evidence classes and caveats.

suppressPackageStartupMessages(library(ggplot2))

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- sub("^--file=", "", script_arg[[1]])
root <- normalizePath(file.path(dirname(script_path), ".."), mustWork = TRUE)

input <- file.path(root, "processed", "publication_claim_matrix.csv")
out_dir <- file.path(root, "processed", "publication_figures")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

dat <- read.csv(input, stringsAsFactors = FALSE)

plot_dat <- data.frame(
  level = c(
    "Context",
    "Experimental pMHC anchor",
    "Modeled pMHC screen",
    "Ternary TCR inference",
    "Patient-level mechanism"
  ),
  support = c(
    "External literature",
    "Supported",
    "Prioritization only",
    "Rejected as evidence",
    "Not tested"
  ),
  y = 5:1,
  fill = c("#56B4E9", "#009E73", "#0072B2", "#D55E00", "#999999")
)

theme_pub <- theme_classic(base_size = 10, base_family = "Arial") +
  theme(
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    axis.text.y = element_text(colour = "#333333", size = 9),
    axis.text.x = element_blank(),
    axis.title = element_blank(),
    plot.title = element_text(colour = "#333333", face = "bold", size = 11, hjust = 0),
    plot.subtitle = element_text(colour = "#333333", size = 8.5, hjust = 0),
    plot.caption = element_text(colour = "#555555", size = 7, hjust = 0),
    panel.grid = element_blank(),
    legend.position = "none",
    plot.margin = margin(7, 8, 7, 7)
  )

p <- ggplot(plot_dat, aes(x = 1, y = y)) +
  annotate("segment", x = 1, xend = 1, y = 1, yend = 5,
           linewidth = 0.7, colour = "#A8A8A8") +
  geom_point(aes(fill = fill), shape = 21, size = 6.2, stroke = 0.75,
             colour = "#333333") +
  geom_label(aes(label = support), x = 1.09, hjust = 0,
             family = "Arial", size = 3.0,
             fill = "white", colour = "#333333", linewidth = 0.25) +
  scale_fill_identity() +
  scale_y_continuous(breaks = plot_dat$y, labels = plot_dat$level,
                     limits = c(0.55, 5.45), expand = c(0, 0)) +
  scale_x_continuous(limits = c(0.7, 1.75), expand = c(0, 0)) +
  labs(
    title = "Evidence ladder",
    subtitle = "Each claim is capped by calibrated evidence",
    caption = sprintf("Source: claim matrix (%d claims).\nTCR inference rejected after calibration failure.", nrow(dat))
  ) +
  theme_pub

png_path <- file.path(out_dir, "figure_3_claim_ladder_300dpi.png")
pdf_path <- file.path(out_dir, "figure_3_claim_ladder.pdf")
ggsave(png_path, p, width = 5, height = 4, units = "in", dpi = 300,
       device = ragg::agg_png, bg = "white")
img <- png::readPNG(png_path)
pdf(pdf_path, width = 5, height = 4, useDingbats = FALSE)
par(mar = c(0, 0, 0, 0)); plot.new(); rasterImage(img, 0, 0, 1, 1)
dev.off()

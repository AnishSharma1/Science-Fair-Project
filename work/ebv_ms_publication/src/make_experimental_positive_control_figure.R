# Publication-ready experimental positive-control figure.
# Inputs are residue-level C-alpha distances calculated from PDB 1H15 and 1BX2
# after superposition of their HLA peptide-binding platforms.

suppressPackageStartupMessages(library(ggplot2))

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- sub("^--file=", "", script_arg[[1]])
root <- normalizePath(file.path(dirname(script_path), ".."), mustWork = TRUE)
input <- file.path(root, "processed", "experimental_positive_control", "experimental_core_position_distances.csv")
out_dir <- file.path(root, "processed", "experimental_positive_control")
dat <- read.csv(input, stringsAsFactors = FALSE)
dat$core_position <- as.numeric(dat$core_position)
dat$ca_distance_after_hla_fit_a <- as.numeric(dat$ca_distance_after_hla_fit_a)
rmsd <- sqrt(mean(dat$ca_distance_after_hla_fit_a^2))

theme_pub <- theme_classic(base_size = 10, base_family = "Arial") +
  theme(
    axis.line = element_line(colour = "#333333", linewidth = 0.45),
    axis.ticks = element_line(colour = "#333333", linewidth = 0.45),
    axis.ticks.length = grid::unit(2.2, "pt"),
    axis.text = element_text(colour = "#333333"),
    axis.title = element_text(colour = "#333333", face = "plain"),
    plot.title = element_text(colour = "#333333", face = "bold", size = 11, hjust = 0),
    plot.subtitle = element_text(colour = "#333333", size = 8.6, hjust = 0),
    plot.caption = element_text(colour = "#555555", size = 7, hjust = 0),
    panel.grid = element_blank(),
    plot.margin = margin(7, 8, 7, 7)
  )

p <- ggplot(dat, aes(x = core_position, y = ca_distance_after_hla_fit_a)) +
  geom_hline(yintercept = rmsd, linewidth = 0.55, linetype = "dashed", colour = "#666666") +
  geom_segment(aes(xend = core_position, y = 0, yend = ca_distance_after_hla_fit_a),
               linewidth = 0.7, colour = "#B8C4D8") +
  geom_point(shape = 21, size = 4.2, stroke = 0.75, fill = "#0072B2", colour = "#0072B2") +
  geom_text(aes(label = residue_pair), vjust = -1.15, family = "Arial", size = 3.0, colour = "#333333") +
  annotate("label", x = 5.35, y = 1.08,
           label = sprintf("core RMSD = %.3f Å", rmsd),
           family = "Arial", size = 3.0,
           fill = "white", colour = "#333333") +
  scale_x_continuous(breaks = 1:7, expand = expansion(mult = c(0.04, 0.08))) +
  scale_y_continuous(limits = c(0, max(dat$ca_distance_after_hla_fit_a) + 0.48),
                     expand = expansion(mult = c(0, 0))) +
  labs(
    title = "pMHC surface equivalence",
    subtitle = "DR2a–BALF5 versus DR2b–MBP: Cα distance after HLA-groove superposition",
    x = "Aligned peptide-core position",
    y = expression("C"[alpha]*" distance after HLA fit (Å)"),
    caption = "PDB 1H15 (DRB5*01:01–BALF5) and 1BX2 (DRB1*15:01–MBP). Residue labels: BALF5/MBP."
  ) + theme_pub

png_path <- file.path(out_dir, "figure_1_experimental_positive_control_300dpi.png")
pdf_path <- file.path(out_dir, "figure_1_experimental_positive_control.pdf")
ggsave(png_path, p, width = 5, height = 4, units = "in", dpi = 300,
       device = ragg::agg_png, bg = "white")
# The PDF preserves the exact Arial-rendered 300-DPI panel rather than allowing
# a platform-specific PDF font substitution.
img <- png::readPNG(png_path)
pdf(pdf_path, width = 5, height = 4, useDingbats = FALSE)
par(mar = c(0, 0, 0, 0)); plot.new(); rasterImage(img, 0, 0, 1, 1)
dev.off()

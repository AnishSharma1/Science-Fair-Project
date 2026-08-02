# Publication-ready pMHC screen shortlist figure.
# This is a candidate-prioritization plot, not a receptor-binding result.

suppressPackageStartupMessages(library(ggplot2))

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- sub("^--file=", "", script_arg[[1]])
root <- normalizePath(file.path(dirname(script_path), ".."), mustWork = TRUE)

input <- file.path(root, "processed", "fullscreen_tier1_ebv_myelin_shortlist.csv")
out_dir <- file.path(root, "processed", "publication_figures")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

dat <- read.csv(input, stringsAsFactors = FALSE)
dat$review_priority_heuristic <- as.numeric(dat$review_priority_heuristic)
dat$local_peptide_ca_rmsd_after_hla_fit <- as.numeric(dat$local_peptide_ca_rmsd_after_hla_fit)
dat$property_similarity <- as.numeric(dat$property_similarity)
dat$pair_label <- paste(
  sub("^EBV_TCELL_", "EBV ", dat$ebv_candidate_id),
  sub("^HUMAN_MYELIN_", "M ", dat$human_candidate_id),
  sep = " / "
)
top <- dat[order(-dat$review_priority_heuristic), ][seq_len(min(12, nrow(dat))), ]
top$pair_label <- factor(top$pair_label, levels = rev(top$pair_label))
top$highlight <- ifelse(top$ebv_candidate_id == "EBV_TCELL_63843", "BALF5-family", "Other")

theme_pub <- theme_classic(base_size = 9.5, base_family = "Arial") +
  theme(
    axis.line = element_line(colour = "#333333", linewidth = 0.45),
    axis.ticks = element_line(colour = "#333333", linewidth = 0.45),
    axis.ticks.length = grid::unit(2.2, "pt"),
    axis.text = element_text(colour = "#333333"),
    axis.text.y = element_text(size = 7.3),
    axis.title = element_text(colour = "#333333"),
    plot.title = element_text(colour = "#333333", face = "bold", size = 11, hjust = 0),
    plot.subtitle = element_text(colour = "#333333", size = 8.0, hjust = 0),
    plot.caption = element_text(colour = "#555555", size = 6.8, hjust = 0),
    legend.title = element_blank(),
    legend.position = c(0.76, 0.18),
    legend.background = element_rect(fill = "white", colour = "#CCCCCC", linewidth = 0.25),
    panel.grid = element_blank(),
    plot.margin = margin(7, 8, 7, 7)
  )

p <- ggplot(top, aes(x = review_priority_heuristic, y = pair_label, fill = highlight)) +
  geom_col(width = 0.62, colour = "#333333", linewidth = 0.18) +
  geom_text(aes(label = sprintf("%.2f A", local_peptide_ca_rmsd_after_hla_fit)),
            hjust = -0.08, family = "Arial", size = 2.25, colour = "#333333") +
  scale_fill_manual(values = c("BALF5-family" = "#0072B2", "Other" = "#E69F00")) +
  scale_x_continuous(limits = c(0, max(top$review_priority_heuristic) * 1.30),
                     expand = expansion(mult = c(0, 0))) +
  labs(
    title = "Modeled pMHC shortlist",
    subtitle = "Top EBV--myelin pairs by predeclared priority score",
    x = "Candidate-priority score",
    y = NULL,
    caption = "Priority only; not TCR evidence.\nActivation, affinity, and patient mechanism were not tested."
  ) +
  theme_pub

png_path <- file.path(out_dir, "figure_2_pmhc_screen_shortlist_300dpi.png")
pdf_path <- file.path(out_dir, "figure_2_pmhc_screen_shortlist.pdf")
ggsave(png_path, p, width = 5, height = 4, units = "in", dpi = 300,
       device = ragg::agg_png, bg = "white")
img <- png::readPNG(png_path)
pdf(pdf_path, width = 5, height = 4, useDingbats = FALSE)
par(mar = c(0, 0, 0, 0)); plot.new(); rasterImage(img, 0, 0, 1, 1)
dev.off()

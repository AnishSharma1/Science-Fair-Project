# Publication-style plot for the GSE190847 PyTorch expression classifier.

suppressPackageStartupMessages(library(ggplot2))

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- sub("^--file=", "", script_arg[[1]])
root <- normalizePath(file.path(dirname(script_path), ".."), mustWork = TRUE)

base_dir <- file.path(root, "processed", "geo_gse190847", "pytorch_expression_classifier")
pred_path <- file.path(base_dir, "gse190847_pytorch_loocv_predictions.csv")
metric_path <- file.path(base_dir, "gse190847_pytorch_classifier_metrics.csv")
out_dir <- file.path(root, "processed", "publication_figures")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

pred <- read.csv(pred_path, stringsAsFactors = FALSE)
metrics <- read.csv(metric_path, stringsAsFactors = FALSE)
pred$group_label <- factor(
  pred$group,
  levels = c("healthy_control", "PPMS"),
  labels = c("Healthy control", "Untreated PPMS")
)

theme_pub <- theme_classic(base_size = 10, base_family = "Arial") +
  theme(
    axis.line = element_line(colour = "#333333", linewidth = 0.45),
    axis.ticks = element_line(colour = "#333333", linewidth = 0.45),
    axis.ticks.length = grid::unit(2.2, "pt"),
    axis.text = element_text(colour = "#333333"),
    axis.title = element_text(colour = "#333333"),
    plot.title = element_text(colour = "#333333", face = "bold", size = 11, hjust = 0),
    plot.subtitle = element_text(colour = "#333333", size = 8.3, hjust = 0),
    plot.caption = element_text(colour = "#555555", size = 7, hjust = 0),
    legend.position = "none",
    panel.grid = element_blank(),
    plot.margin = margin(7, 8, 7, 7)
  )

label <- sprintf(
  "LOOCV AUC = %.3f\nPermutation p = %.3f",
  metrics$auc[1],
  metrics$empirical_auc_p_ge_observed[1]
)

p <- ggplot(pred, aes(x = group_label, y = loocv_ppms_probability, fill = group_label)) +
  geom_hline(yintercept = 0.5, linetype = "dashed", linewidth = 0.45, colour = "#666666") +
  geom_boxplot(width = 0.48, outlier.shape = NA, linewidth = 0.45, colour = "#333333", alpha = 0.75) +
  geom_jitter(width = 0.08, height = 0, shape = 21, size = 2.5,
              colour = "#333333", stroke = 0.35, alpha = 0.9) +
  annotate("label", x = 1.48, y = 0.93, label = label,
           family = "Arial", size = 3.0, fill = "white",
           colour = "#333333", linewidth = 0.25) +
  scale_fill_manual(values = c("Healthy control" = "#56B4E9", "Untreated PPMS" = "#D55E00")) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25), expand = c(0, 0)) +
  labs(
    title = "Transcriptomic context classifier",
    subtitle = "PyTorch linear probe on seven HLA-II/APC genes in GSE190847 B cells",
    x = NULL,
    y = "Leave-one-out PPMS probability",
    caption = "Expression context only; not EBV infection, pMHC presentation, TCR binding, or activation evidence."
  ) +
  theme_pub

png_path <- file.path(out_dir, "figure_4_gse190847_pytorch_expression_300dpi.png")
pdf_path <- file.path(out_dir, "figure_4_gse190847_pytorch_expression.pdf")
ggsave(png_path, p, width = 5, height = 4, units = "in", dpi = 300,
       device = ragg::agg_png, bg = "white")
img <- png::readPNG(png_path)
pdf(pdf_path, width = 5, height = 4, useDingbats = FALSE)
par(mar = c(0, 0, 0, 0)); plot.new(); rasterImage(img, 0, 0, 1, 1)
dev.off()

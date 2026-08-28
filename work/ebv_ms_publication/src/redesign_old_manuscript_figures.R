#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(stringr)
  library(forcats)
  library(scales)
})

`%||%` <- function(x, y) if (is.null(x) || length(x) == 0 || is.na(x)) y else x

script_arg <- commandArgs(trailingOnly = FALSE)
script_path <- sub("^--file=", "", script_arg[grepl("^--file=", script_arg)][1] %||% "")
root_dir <- if (nzchar(script_path)) normalizePath(file.path(dirname(script_path), ".."), mustWork = FALSE) else getwd()
if (!dir.exists(file.path(root_dir, "processed"))) {
  root_dir <- "/Users/anishsharma/Documents/Codex/2026-08-01/i-al/work/ebv_ms_publication"
}

out_dir <- file.path(root_dir, "processed", "manuscript_redesign_figures")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

theme_manuscript <- function(base_size = 10) {
  theme_classic(base_size = base_size, base_family = "sans") +
    theme(
      text = element_text(color = "#333333"),
      plot.title = element_text(face = "bold", size = base_size + 2, margin = margin(b = 4)),
      plot.subtitle = element_text(size = base_size - 0.5, color = "#555555", margin = margin(b = 8)),
      plot.caption = element_text(size = base_size - 2, color = "#666666", hjust = 0, margin = margin(t = 8), lineheight = 1.05),
      axis.title = element_text(face = "bold", size = base_size - 0.5),
      axis.text = element_text(color = "#333333", size = base_size - 1),
      axis.line = element_line(color = "#333333", linewidth = 0.45),
      axis.ticks = element_line(color = "#333333", linewidth = 0.4),
      legend.title = element_text(face = "bold", size = base_size - 1),
      legend.text = element_text(size = base_size - 1),
      legend.key.size = unit(0.35, "cm"),
      plot.margin = margin(10, 16, 10, 10)
    )
}

save_pub <- function(plot, filename, width = 5, height = 4) {
  png_path <- file.path(out_dir, paste0(filename, ".png"))
  pdf_path <- file.path(out_dir, paste0(filename, ".pdf"))
  ggsave(png_path, plot = plot, width = width, height = height, dpi = 300, bg = "white", device = ragg::agg_png)
  ggsave(pdf_path, plot = plot, width = width, height = height, bg = "white", device = function(file, width, height, ...) grDevices::pdf(file = file, width = width, height = height, useDingbats = FALSE, ...))
  invisible(c(png_path, pdf_path))
}

cb <- c(
  blue = "#0072B2",
  sky = "#56B4E9",
  green = "#009E73",
  orange = "#E69F00",
  vermillion = "#D55E00",
  purple = "#CC79A7",
  gray = "#7A7A7A",
  dark = "#333333"
)

label_pair <- function(ebv, human) {
  paste0(str_remove(ebv, "^EBV_(TCELL|MHC)_"), "  ->  ", str_remove(human, "^HUMAN_MYELIN_"))
}

# Figure 1 -----------------------------------------------------------------
core <- read_csv(file.path(root_dir, "processed/experimental_positive_control/experimental_core_position_distances.csv"), show_col_types = FALSE)
metrics <- read_csv(file.path(root_dir, "processed/experimental_positive_control/experimental_drb2_positive_control_metrics.csv"), show_col_types = FALSE)
mean_core <- mean(core$ca_distance_after_hla_fit_a)
groove_rmsd <- as.numeric(metrics$value[metrics$metric == "HLA_groove_CA_RMSD_A"][1])

fig1 <- ggplot(core, aes(core_position, ca_distance_after_hla_fit_a)) +
  geom_hline(yintercept = mean_core, linetype = "22", color = cb[["orange"]], linewidth = 0.55) +
  geom_segment(aes(xend = core_position, y = 0, yend = ca_distance_after_hla_fit_a), color = "#B9CFE3", linewidth = 1.1) +
  geom_point(aes(fill = residue_pair), shape = 21, size = 3.8, color = cb[["dark"]], stroke = 0.35, show.legend = FALSE) +
  geom_text(aes(label = residue_pair), nudge_y = 0.075, size = 3.1, family = "sans", color = cb[["dark"]]) +
  annotate("label", x = 4.35, y = max(core$ca_distance_after_hla_fit_a) * 0.92,
           label = sprintf("Mean core C-alpha RMSD = %.3f A\nHLA-groove C-alpha RMSD = %.3f A", mean_core, groove_rmsd),
           hjust = 0, size = 2.7, family = "sans",
           fill = "white", color = cb[["dark"]]) +
  scale_x_continuous(breaks = core$core_position) +
  scale_y_continuous(limits = c(0, max(core$ca_distance_after_hla_fit_a) + 0.22), expand = expansion(mult = c(0, 0.02))) +
  scale_fill_manual(values = rep(cb[["blue"]], nrow(core))) +
  labs(
    title = "Experimental pMHC geometry positive control",
    subtitle = "BALF5-DRB1*15:01 vs MBP-DRB1*15:01 after HLA-groove fitting",
    x = "Aligned peptide core position",
    y = "C-alpha distance after HLA fit (A)",
    caption = "Source: experimental positive-control distance table. Supports local geometric mimicry only, not TCR recognition."
  ) +
  theme_manuscript(9)

save_pub(fig1, "figure_1_experimental_positive_control_redesign", 6.2, 4.3)

# Figure 2 -----------------------------------------------------------------
shortlist <- read_csv(file.path(root_dir, "processed/fullscreen_tier1_ebv_myelin_shortlist.csv"), show_col_types = FALSE) %>%
  mutate(
    rank = row_number(desc(review_priority_heuristic)),
    pair_label = label_pair(ebv_candidate_id, human_candidate_id),
    evidence_family = if_else(ebv_candidate_id == "EBV_TCELL_63843", "BALF5-family anchor", "Other modeled pair")
  ) %>%
  arrange(desc(review_priority_heuristic)) %>%
  slice_head(n = 12) %>%
  mutate(pair_label = fct_reorder(pair_label, review_priority_heuristic))

fig2 <- ggplot(shortlist, aes(review_priority_heuristic, pair_label)) +
  geom_col(aes(fill = evidence_family), width = 0.68, color = "white", linewidth = 0.2) +
  geom_text(aes(label = sprintf("RMSD %.2f A", local_peptide_ca_rmsd_after_hla_fit)),
            hjust = -0.08, size = 2.75, family = "sans", color = cb[["dark"]]) +
  scale_fill_manual(values = c("BALF5-family anchor" = cb[["blue"]], "Other modeled pair" = cb[["gray"]])) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.22))) +
  labs(
    title = "pMHC screen prioritizes a compact shortlist",
    subtitle = "Top 12 of 32 EBV-CNS pairs ranked by property similarity and modeled local pMHC geometry",
    x = "Review-priority heuristic score",
    y = NULL,
    fill = "Pair class",
    caption = "Source: Tier 1 shortlist. Ranking is hypothesis-generating; it does not establish binding or pathogenicity."
  ) +
  theme_manuscript(9) +
  theme(legend.position = "top")

save_pub(fig2, "figure_2_pmhc_screen_shortlist_redesign", 7.7, 4.9)

# Figure 3 -----------------------------------------------------------------
claims <- read_csv(file.path(root_dir, "processed/publication_claim_matrix.csv"), show_col_types = FALSE)

ladder <- tibble::tribble(
  ~level, ~claim_boundary, ~status_label, ~status_group,
  1, "EBV-MS / HLA risk", "Background", "Supported",
  2, "BALF5/MBP pMHC", "Direct anchor", "Supported",
  3, "Modeled shortlist", "Hypothesis-generating", "Modeled",
  4, "RNA signal", "Independent support", "Supportive",
  5, "TCR / patient mechanism", "Not claimed", "Boundary"
)

claim_counts <- claims %>%
  count(evidence_class, status, name = "n") %>%
  mutate(summary = paste0(n, " claim", if_else(n == 1, "", "s"), ": ", evidence_class, " / ", status)) %>%
  summarise(note = paste(summary, collapse = "\n")) %>%
  pull(note)

fig3 <- ggplot(ladder, aes(level, 1)) +
  annotate("segment", x = 1, xend = 5, y = 1, yend = 1, color = "#C9C9C9", linewidth = 1) +
  geom_point(aes(fill = status_group), shape = 21, size = 7, color = cb[["dark"]], stroke = 0.45) +
  geom_text(aes(label = level), color = "white", fontface = "bold", size = 3.4, family = "sans") +
  geom_text(aes(label = claim_boundary), y = 1.18, angle = 35, hjust = 0, size = 3.0, family = "sans", color = cb[["dark"]]) +
  geom_text(aes(label = status_label), y = 0.78, size = 2.75, family = "sans", color = "#555555") +
  scale_fill_manual(values = c("Supported" = cb[["green"]], "Modeled" = cb[["blue"]], "Supportive" = cb[["orange"]], "Boundary" = cb[["gray"]])) +
  scale_x_continuous(breaks = 1:5, limits = c(0.7, 5.55)) +
  coord_cartesian(ylim = c(0.62, 1.72), clip = "off") +
  labs(
    title = "Manuscript claim ladder",
    subtitle = "Strongest defensible claim: pMHC geometric mimicry plus supportive transcriptomic signal",
    x = NULL,
    y = NULL,
    fill = "Claim status",
    caption = "Source: publication claim matrix. TCR binding, patient-level causality, and diagnosis remain outside the claim boundary."
  ) +
  theme_void(base_family = "sans") +
  theme(
    plot.title = element_text(face = "bold", size = 12, color = cb[["dark"]], margin = margin(b = 4)),
    plot.subtitle = element_text(size = 9.5, color = "#555555", margin = margin(b = 8)),
    plot.caption = element_text(size = 7.5, color = "#666666", hjust = 0, margin = margin(t = 10)),
    legend.position = "bottom",
    legend.title = element_text(face = "bold", size = 9),
    legend.text = element_text(size = 9),
    plot.margin = margin(12, 16, 18, 10)
  )

save_pub(fig3, "figure_3_claim_ladder_redesign", 7.4, 4.4)

# Figure 4 -----------------------------------------------------------------
pred <- read_csv(file.path(root_dir, "processed/geo_gse190847/pytorch_expression_classifier/gse190847_pytorch_loocv_predictions.csv"), show_col_types = FALSE) %>%
  mutate(group_label = recode(group, healthy_control = "Healthy control", ppms = "PPMS"))
met <- read_csv(file.path(root_dir, "processed/geo_gse190847/pytorch_expression_classifier/gse190847_pytorch_classifier_metrics.csv"), show_col_types = FALSE)

fig4 <- ggplot(pred, aes(group_label, loocv_ppms_probability, fill = group_label)) +
  geom_hline(yintercept = 0.5, color = "#888888", linetype = "22", linewidth = 0.45) +
  geom_boxplot(width = 0.52, alpha = 0.75, outlier.shape = NA, color = cb[["dark"]], linewidth = 0.45) +
  geom_jitter(width = 0.12, height = 0, size = 2.0, alpha = 0.78, shape = 21, color = "white", stroke = 0.2) +
  stat_summary(fun = mean, geom = "point", shape = 23, size = 3.2, fill = "white", color = cb[["dark"]], stroke = 0.5) +
  annotate("label", x = 1.5, y = 0.91,
           label = sprintf("LOOCV AUC = %.3f\nEmpirical p = %.3f\nn = %d healthy / %d PPMS",
                           met$auc[1], met$empirical_auc_p_ge_observed[1], met$healthy_n[1], met$ppms_n[1]),
           family = "sans", size = 3.0, fill = "white", color = cb[["dark"]]) +
  scale_fill_manual(values = c("Healthy control" = cb[["sky"]], "PPMS" = cb[["vermillion"]])) +
  scale_y_continuous(labels = percent_format(accuracy = 1), limits = c(0, 1), expand = expansion(mult = c(0.01, 0.04))) +
  labs(
    title = "RNA classifier shows supportive PPMS signal",
    subtitle = "Seven-gene PyTorch leave-one-out classifier on GSE190847",
    x = NULL,
    y = "Predicted PPMS probability",
    caption = "Source: GSE190847 LOOCV predictions and classifier metrics. Supportive, not diagnostic."
  ) +
  theme_manuscript(9) +
  theme(legend.position = "none")

save_pub(fig4, "figure_4_gse190847_pytorch_expression_redesign", 6.2, 4.3)

# Figure 5 -----------------------------------------------------------------
annot <- read_csv(file.path(root_dir, "processed/external_validation_benchmark/external_validation_pair_annotations.csv"), show_col_types = FALSE) %>%
  mutate(
    validation_group = case_when(
      pair_validation == "classic_BALF5_MBP_pair" ~ "Classic BALF5-MBP",
      is_strict_external_new_literature ~ "Strict external overlay",
      is_external_overlay ~ "Other external overlay",
      TRUE ~ "Unannotated screen pair"
    ),
    validation_group = factor(validation_group, levels = c("Classic BALF5-MBP", "Strict external overlay", "Other external overlay", "Unannotated screen pair")),
    label = if_else(rank <= 5 | validation_group == "Strict external overlay",
                    paste0(str_remove(ebv_candidate_id, "^EBV_(TCELL|MHC)_"), " -> ", str_remove(human_candidate_id, "^HUMAN_MYELIN_")),
                    NA_character_)
  )
summary <- read_csv(file.path(root_dir, "processed/external_validation_benchmark/external_validation_rank_recovery_summary.csv"), show_col_types = FALSE)
classic <- summary %>% filter(test == "classic_BALF5_MBP_pair_recovery")
strict <- summary %>% filter(test == "strict_new_literature_overlay")

fig5 <- ggplot(annot, aes(rank, review_priority_heuristic)) +
  annotate("rect", xmin = 0.5, xmax = 10.5, ymin = -Inf, ymax = Inf, fill = "#F2F2F2", alpha = 0.9) +
  geom_point(aes(fill = validation_group, size = validation_group), shape = 21, color = cb[["dark"]], stroke = 0.35, alpha = 0.95) +
  geom_text(data = filter(annot, !is.na(label)), aes(label = label), nudge_y = 0.012, size = 2.35, family = "sans", color = cb[["dark"]], check_overlap = TRUE) +
  annotate("text", x = 5.5, y = max(annot$review_priority_heuristic, na.rm = TRUE) * 0.98,
           label = "Top-10 priority zone", size = 3.0, family = "sans", color = "#555555") +
  annotate("label", x = 23, y = max(annot$review_priority_heuristic, na.rm = TRUE) * 0.94,
           label = sprintf("Classic recovery: %d/%d in top 10, p < %.4f\nStrict overlay: %d/%d in top 10, p = %.3f",
                           classic$observed_top10_positive_count[1], classic$positive_pair_count[1],
                           classic$empirical_p_top10_count_ge_observed[1] + 0.000001,
                           strict$observed_top10_positive_count[1], strict$positive_pair_count[1],
                           strict$empirical_p_top10_count_ge_observed[1]),
           hjust = 0.5, family = "sans", size = 2.85, fill = "white", color = cb[["dark"]]) +
  scale_fill_manual(values = c(
    "Classic BALF5-MBP" = cb[["blue"]],
    "Strict external overlay" = cb[["vermillion"]],
    "Other external overlay" = cb[["orange"]],
    "Unannotated screen pair" = "#BDBDBD"
  )) +
  scale_size_manual(values = c(
    "Classic BALF5-MBP" = 3.4,
    "Strict external overlay" = 3.2,
    "Other external overlay" = 2.8,
    "Unannotated screen pair" = 2.2
  ), guide = "none") +
  scale_x_continuous(breaks = c(1, 5, 10, 16, 24, 32), limits = c(0.5, 32.8)) +
  scale_y_continuous(expand = expansion(mult = c(0.06, 0.12))) +
  labs(
    title = "External annotations align with pMHC screen priority",
    subtitle = "Ranked EBV-CNS pairs, with classic BALF5-MBP and newer literature overlays marked",
    x = "Priority rank (1 = highest)",
    y = "Review-priority heuristic score",
    fill = "External annotation",
    caption = "Source: external-validation benchmark tables. Recovery validates prioritization behavior, not causal disease mechanism."
  ) +
  guides(fill = guide_legend(nrow = 2, byrow = TRUE)) +
  theme_manuscript(9) +
  theme(legend.position = "bottom")

save_pub(fig5, "figure_5_external_validation_overlay_redesign", 8.2, 4.8)

readme <- c(
  "# Redesigned older five-figure manuscript set",
  "",
  "Generated by `src/redesign_old_manuscript_figures.R` from existing processed project tables.",
  "",
  "## Outputs",
  "",
  "- `figure_1_experimental_positive_control_redesign`: experimental BALF5/MBP pMHC positive-control geometry.",
  "- `figure_2_pmhc_screen_shortlist_redesign`: top pMHC screen shortlist from the 32-pair Tier 1 table.",
  "- `figure_3_claim_ladder_redesign`: manuscript evidence-boundary ladder from the claim matrix.",
  "- `figure_4_gse190847_pytorch_expression_redesign`: GSE190847 seven-gene PyTorch LOOCV probability separation.",
  "- `figure_5_external_validation_overlay_redesign`: external-validation overlay against the ranked pMHC screen.",
  "",
  "Each figure is exported as vector PDF and 300 DPI PNG.",
  "",
  "## Manuscript guardrail",
  "",
  "These redesigned panels intentionally preserve the conservative interpretation: the project supports experimental/model-based pMHC geometric mimicry and independent transcriptomic support, but does not claim direct TCR binding, patient-level causality, or diagnostic performance."
)
writeLines(readme, file.path(out_dir, "README.md"))

message("Wrote redesigned figure set to: ", out_dir)

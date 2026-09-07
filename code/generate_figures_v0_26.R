#!/usr/bin/env Rscript

# Publication figures for the locked v0.26 schema-fixed aggregate outputs.
# Only base R is used so the script remains runnable with the documented R install.

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- if (length(script_arg)) sub("^--file=", "", script_arg[[1]]) else "code/generate_figures_v0_26.R"
repo_root <- normalizePath(Sys.getenv("PV_REPO_ROOT", unset = file.path(dirname(script_path), "..")), winslash = "/", mustWork = TRUE)
data_dir <- file.path(repo_root, "aggregate_outputs", "schema_fixed")
figure_dir <- file.path(repo_root, "figures", "v0_26")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

faers_blue <- "#2C7FB8"
canada_orange <- "#F28E2B"
jader_green <- "#59A14F"
ink <- "#222222"
muted <- "#777777"
light_grid <- "#EAEAEA"
sparse_grey <- "#B8B8B8"

read_output <- function(name) {
  read.csv(file.path(data_dir, name), stringsAsFactors = FALSE, check.names = FALSE)
}

save_figure <- function(stem, width, height, draw) {
  pdf(file.path(figure_dir, paste0(stem, ".pdf")), width = width, height = height, family = "sans", pointsize = 8, useDingbats = FALSE)
  draw()
  dev.off()
  svg(file.path(figure_dir, paste0(stem, ".svg")), width = width, height = height, pointsize = 8, onefile = FALSE)
  draw()
  dev.off()
  png(file.path(figure_dir, paste0(stem, ".png")), width = width * 300, height = height * 300, res = 300, pointsize = 8, type = "cairo")
  draw()
  dev.off()
}

panel_label <- function(label) {
  mtext(label, side = 3, line = 0.35, adj = 0, font = 2, cex = 1.05)
}

draw_box <- function(x0, y0, x1, y1, text, border = ink, fill = "#FFFFFF", cex = 0.78) {
  rect(x0, y0, x1, y1, border = border, col = fill, lwd = 1.1)
  text((x0 + x1) / 2, (y0 + y1) / 2, text, cex = cex, col = ink)
}

draw_arrow <- function(x0, y0, x1, y1) {
  arrows(x0, y0, x1, y1, length = 0.07, lwd = 1, col = muted)
}

agent_context <- read_output("agent_within_nonculprit_context_v0_26.csv")
agent_primary <- read_output("agent_primary_suspect_suspect_v0_26.csv")
volume_proxy <- read_output("agent_same_class_report_volume_proxy_v0_26.csv")
coexposure_faers <- read_output("faers_coexposure_counts.csv")
coexposure_other <- read_output("coexposure_counts_canada_jader.csv")
jader_audit <- read_output("jader_structure_and_vaccine_audit_v0_26.csv")
country <- read_output("canonical_agent_country_v0_26.csv")

# Figure 1: raw sources, retained units, and the analytical hierarchy.
save_figure("Figure_1_source_to_retained_unit_workflow", 10.2, 5.8, function() {
  par(mar = c(0.7, 0.7, 0.7, 0.7), xpd = NA)
  plot.new()
  plot.window(xlim = c(0, 10), ylim = c(0, 7))

  draw_box(0.35, 5.45, 2.55, 6.45, "FAERS official raw archives\n2004 Q1-2025 Q4", border = faers_blue, fill = "#EAF3F8")
  draw_box(3.90, 5.45, 6.10, 6.45, "Canada Vigilance\nfrozen local extract\n31 Jul 2025", border = canada_orange, fill = "#FDF1E5")
  draw_box(7.45, 5.45, 9.65, 6.45, "JADER\nfrozen local extract\nNov 2025", border = jader_green, fill = "#EDF6EC")

  draw_arrow(1.45, 5.42, 1.45, 4.72)
  draw_arrow(5.00, 5.42, 5.00, 4.72)
  draw_arrow(8.55, 5.42, 8.55, 4.72)
  draw_box(0.35, 3.55, 2.55, 4.55, "19,961,236 retained\ncase-version reports", border = faers_blue)
  draw_box(3.90, 3.55, 6.10, 4.55, "1,204,167 retained\nreport identifiers", border = canada_orange)
  draw_box(7.45, 3.55, 9.65, 4.55, "998,397 retained\nreport identifiers", border = jader_green)

  draw_arrow(1.45, 3.52, 4.15, 2.75)
  draw_arrow(5.00, 3.52, 5.00, 2.75)
  draw_arrow(8.55, 3.52, 5.85, 2.75)
  draw_box(3.10, 1.80, 6.90, 2.75, "Three-PT core-anaphylaxis outcome\nAnaphylaxis | Anaphylactic reaction\nAnaphylactic shock", fill = "#F8F8F8", cex = 0.70)
  draw_arrow(5.00, 1.77, 5.00, 1.35)
  draw_box(1.55, 0.20, 8.45, 1.35, "Medication-marker context -> same-class suspect-agent contrasts\nSupporting diagnostics: report-volume proxy; co-suspected classes; role fields\ncountry mix; JADER structure; terminology sensitivity", fill = "#FFFFFF", cex = 0.62)
})

# Figure 2: the primary reporting contrasts and the FAERS reference-setting movement.
save_figure("Figure_2_agent_contrasts_reference_setting", 10.3, 7.2, function() {
  par(mfrow = c(1, 2), mar = c(4.2, 6.9, 1.3, 1.1), oma = c(0, 0, 0, 0), xpd = NA)
  agents <- c("rocuronium", "succinylcholine", "cefazolin", "teicoplanin")
  agent_labels <- c("Rocuronium", "Succinylcholine", "Cefazolin", "Teicoplanin")
  faers_base <- agent_primary[agent_primary$database == "FAERS" & agent_primary$agent %in% agents, ]
  faers_context <- agent_context[agent_context$database == "FAERS" & agent_context$agent %in% agents, ]
  faers_base <- faers_base[match(agents, faers_base$agent), ]
  faers_context <- faers_context[match(agents, faers_context$agent), ]
  y <- rev(seq_along(agents))
  plot(NA, xlim = c(log10(0.1), log10(10)), ylim = c(0.5, 4.5), yaxt = "n", xaxt = "n", xlab = "Reporting odds ratio (log scale)", ylab = "", bty = "n")
  abline(v = log10(c(0.1, 0.2, 0.5, 1, 2, 5, 10)), col = light_grid, lwd = 0.8)
  axis(1, at = log10(c(0.1, 0.2, 0.5, 1, 2, 5, 10)), labels = c("0.1", "0.2", "0.5", "1", "2", "5", "10"))
  axis(2, at = y, labels = agent_labels, las = 1, tick = FALSE)
  abline(v = 0, col = ink, lty = 2, lwd = 1)
  for (i in seq_along(agents)) {
    segments(log10(faers_base$ror[i]), y[i], log10(faers_context$ror[i]), y[i], col = muted, lwd = 1.5)
    points(log10(faers_base$ror[i]), y[i], pch = 21, bg = "#FFFFFF", col = muted, cex = 1.25)
    points(log10(faers_context$ror[i]), y[i], pch = 21, bg = faers_blue, col = faers_blue, cex = 1.25)
  }
  legend("bottomleft", legend = c("Role-symmetric baseline", "Medication-marker context"), pch = 21, pt.bg = c("#FFFFFF", faers_blue), col = c(muted, faers_blue), bty = "n", cex = 0.77)
  panel_label("A")

  context <- agent_context[agent_context$agent %in% agents, ]
  plot(NA, xlim = c(log10(0.08), log10(15)), ylim = c(0.5, 4.5), yaxt = "n", xaxt = "n", xlab = "Medication-context ROR (log scale)", ylab = "", bty = "n")
  abline(v = log10(c(0.1, 0.2, 0.5, 1, 2, 5, 10)), col = light_grid, lwd = 0.8)
  axis(1, at = log10(c(0.1, 0.2, 0.5, 1, 2, 5, 10)), labels = c("0.1", "0.2", "0.5", "1", "2", "5", "10"))
  axis(2, at = y, labels = agent_labels, las = 1, tick = FALSE)
  abline(v = 0, col = ink, lty = 2, lwd = 1)
  db_col <- c(FAERS = faers_blue, Canada = canada_orange, JADER = jader_green)
  db_offset <- c(FAERS = -0.16, Canada = 0, JADER = 0.16)
  for (db in names(db_col)) {
    rows <- context[context$database == db, ]
    rows <- rows[match(agents, rows$agent), ]
    for (i in seq_along(agents)) {
      yy <- y[i] + db_offset[[db]]
      if (rows$sparse_agent_contrast[i] == 0 && is.finite(rows$ror[i])) {
        segments(log10(rows$ror_ci_low[i]), yy, log10(rows$ror_ci_high[i]), yy, col = db_col[[db]], lwd = 1.2)
        points(log10(rows$ror[i]), yy, pch = 21, bg = db_col[[db]], col = db_col[[db]], cex = 1.05)
      } else {
        points(log10(0.12), yy, pch = 21, bg = sparse_grey, col = sparse_grey, cex = 0.95)
      }
    }
  }
  text(log10(0.13), y[4] + 0.30, "grey = sparse", pos = 4, cex = 0.67, col = muted)
  legend("bottomleft", legend = c("FAERS", "Canada Vigilance", "JADER"), pch = 21, pt.bg = unname(db_col), col = unname(db_col), bty = "n", cex = 0.77)
  panel_label("B")
})

# Figure 3: report-volume proxy and co-suspected class overlap.
save_figure("Figure_3_report_volume_proxy_and_co_suspected_classes", 10.3, 6.5, function() {
  par(mfrow = c(1, 2), mar = c(4.4, 4.6, 1.3, 1.2), xpd = NA)
  agents <- c("rocuronium", "succinylcholine", "cefazolin", "teicoplanin")
  labels <- c("Rocuronium", "Succinylcholine", "Cefazolin", "Teicoplanin")
  dat <- volume_proxy[volume_proxy$database == "FAERS" & volume_proxy$agent %in% agents, ]
  dat <- dat[match(agents, dat$agent), ]
  x <- seq_along(agents)
  plot(NA, xlim = c(0.6, 4.4), ylim = c(0, 70), xaxt = "n", xlab = "", ylab = "Agent share among same-class suspect reports (%)", bty = "n")
  abline(h = seq(0, 70, 10), col = light_grid, lwd = 0.8)
  axis(1, at = x, labels = labels, las = 2, tick = FALSE)
  all_share <- dat$all_same_class_agent_report_share_proxy * 100
  core_share <- dat$core_agent_share * 100
  for (i in seq_along(x)) {
    segments(x[i], all_share[i], x[i], core_share[i], col = muted, lwd = 1.5)
  }
  points(x, all_share, pch = 21, bg = "#FFFFFF", col = muted, cex = 1.25)
  points(x, core_share, pch = 21, bg = faers_blue, col = faers_blue, cex = 1.25)
  legend("topleft", legend = c("All same-class suspect reports", "Core-anaphylaxis reports"), pch = 21, pt.bg = c("#FFFFFF", faers_blue), col = c(muted, faers_blue), bty = "n", cex = 0.75)
  panel_label("A")

  target_pair <- "NMBA__antibiotic_anchor"
  get_pair <- function(database, total) {
    source <- if (database == "FAERS") coexposure_faers else coexposure_other
    pair <- source[source$database == database & source$pair == target_pair, ]
    as.numeric(pair$suspect_pair_n[[1]]) / total * 100
  }
  totals <- c(FAERS = 3206, Canada = 79, JADER = 1183)
  pct <- c(FAERS = get_pair("FAERS", totals[["FAERS"]]), Canada = get_pair("Canada", totals[["Canada"]]), JADER = get_pair("JADER", totals[["JADER"]]))
  bars <- barplot(pct, names.arg = c("FAERS", "Canada\nVigilance", "JADER"), ylim = c(0, 40), col = c(faers_blue, canada_orange, jader_green), border = NA, ylab = "NMBA-suspect core reports also naming\nan antibiotic anchor (%)", las = 1)
  abline(h = seq(0, 40, 10), col = light_grid, lwd = 0.8)
  text(bars, pct + 1.5, labels = c("997/3,206", "22/79", "264/1,183"), cex = 0.77, col = ink)
  panel_label("B")
})

# Figure 4: database-structure and country-mixture diagnostics.
save_figure("Figure_4_structural_diagnostics", 10.3, 6.5, function() {
  par(mfrow = c(1, 2), mar = c(4.4, 7.2, 1.3, 1.1), xpd = NA)
  scenarios <- c("baseline", "exclude_vaccine_like", "exclude_2021_q1_q2", "exclude_vaccine_like_and_2021_q1_q2")
  audit <- jader_audit[jader_audit$scope == "nonculprit_medication_marker_context" & jader_audit$scenario %in% scenarios, ]
  audit <- audit[match(scenarios, audit$scenario), ]
  y <- rev(seq_along(scenarios))
  labels <- c("Baseline", "No vaccine-like", "No 2021 Q1-Q2", "Both exclusions")
  plot(NA, xlim = c(log10(5), log10(60)), ylim = c(0.5, 4.5), xaxt = "n", yaxt = "n", xlab = "Chlorhexidine ROR (log scale)", ylab = "", bty = "n")
  axis(1, at = log10(c(5, 10, 20, 40, 60)), labels = c("5", "10", "20", "40", "60"))
  axis(2, at = y, labels = labels, las = 1, tick = FALSE)
  abline(v = log10(c(5, 10, 20, 40, 60)), col = light_grid, lwd = 0.8)
  for (i in seq_along(scenarios)) {
    segments(log10(audit$ror_ci_low[i]), y[i], log10(audit$ror_ci_high[i]), y[i], col = jader_green, lwd = 1.4)
    points(log10(audit$ror[i]), y[i], pch = 21, bg = jader_green, col = jader_green, cex = 1.15)
  }
  panel_label("A")

  par(mar = c(4.4, 5.3, 1.3, 1.1), xpd = NA)
  tei <- country[country$agent == "teicoplanin", ]
  order_levels <- c("GB", "Non-GB observed", "Missing OCCR_COUNTRY")
  tei <- tei[match(order_levels, tei$country_stratum), ]
  core <- tei$a / (tei$a + tei$b) * 100
  noncore <- tei$c / (tei$c + tei$d) * 100
  y <- rev(seq_along(order_levels))
  plot(NA, xlim = c(0, 28), ylim = c(0.5, 3.5), yaxt = "n", xlab = "Teicoplanin share among antibiotic-anchor reports (%)", ylab = "", bty = "n")
  axis(2, at = y, labels = c("GB", "Observed non-GB", "Missing country"), las = 1, tick = FALSE)
  abline(v = seq(0, 25, 5), col = light_grid, lwd = 0.8)
  for (i in seq_along(y)) {
    col <- if (tei$sparse_display[i] == 1) sparse_grey else faers_blue
    segments(noncore[i], y[i], core[i], y[i], col = muted, lwd = 1.3)
    points(noncore[i], y[i], pch = 21, bg = "#FFFFFF", col = muted, cex = 1.2)
    points(core[i], y[i], pch = if (tei$sparse_display[i] == 1) 1 else 21, bg = if (tei$sparse_display[i] == 1) "#FFFFFF" else col, col = col, cex = 1.2)
    if (tei$sparse_display[i] == 1) text(1.0, y[i] + 0.23, "sparse core cell", pos = 4, cex = 0.67, col = muted)
  }
  legend("topright", legend = c("Core reports", "Leave-event-out non-core reports"), pch = c(21, 21), pt.bg = c(faers_blue, "#FFFFFF"), col = c(faers_blue, muted), bty = "n", cex = 0.74)
  panel_label("B")
})

writeLines(c(
  "v0.26 figures generated from aggregate_outputs/schema_fixed.",
  "Sparse agent cells are de-emphasized rather than converted with continuity corrections.",
  "Output formats: PDF, SVG, PNG preview."
), file.path(figure_dir, "FIGURE_BUILD_NOTE.txt"))

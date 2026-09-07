#!/usr/bin/env Rscript

# Publication figures derived from v0.20 canonical numeric artifacts plus the
# v0.23 display-only sparse annotation. No numeric estimate is changed.

root <- Sys.getenv("PV_REPO_ROOT", unset = normalizePath(file.path(getwd()), winslash = "/", mustWork = FALSE))
out_dir <- file.path(root, "figures")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

read_project_csv <- function(path) {
  d <- read.csv(path, stringsAsFactors = FALSE)
  # UTF-8 BOMs can be converted by read.csv/check.names into X...<name>.
  names(d)[1] <- sub("^X\\.+", "", names(d)[1])
  d
}

context <- read_project_csv(file.path(root, "aggregate_outputs", "canonical_agent_context_v0_20.csv"))
burden <- read_project_csv(file.path(root, "aggregate_outputs", "canonical_chlorhexidine_burden_v0_20.csv"))
country <- read_project_csv(file.path(root, "aggregate_outputs", "agent_country_strata_display_v0_23.csv"))
composition <- read_project_csv(file.path(root, "aggregate_outputs", "canonical_antibiotic_country_composition_v0_20.csv"))
roles <- read_project_csv(file.path(root, "aggregate_outputs", "canonical_role_profiles_v0_20.csv"))
legacy <- read_project_csv(file.path(root, "aggregate_outputs", "agent_primary_suspect_suspect_v0_14.csv"))

pal <- list(
  ink = "#1F2933", muted = "#697782", grid = "#D7DEE2", light = "#F3F6F7",
  context = "#007C83", legacy = "#8A969E", nmba = "#126782",
  antibiotic = "#C75B2A", teico = "#C75B2A", cefa = "#126782",
  ps = "#126782", ss = "#C75B2A", sparse = "#9AA5AB"
)

agent_label <- function(x) {
  map <- c(rocuronium = "Rocuronium", succinylcholine = "Succinylcholine",
           cefazolin = "Cefazolin", teicoplanin = "Teicoplanin")
  unname(map[x])
}

open_devices <- function(stem, width, height, draw_fun) {
  pdf_path <- file.path(out_dir, paste0(stem, ".pdf"))
  png_path <- file.path(out_dir, paste0(stem, ".png"))
  grDevices::cairo_pdf(pdf_path, width = width, height = height, family = "sans")
  draw_fun()
  grDevices::dev.off()
  grDevices::png(png_path, width = width * 300, height = height * 300, res = 300, type = "cairo")
  draw_fun()
  grDevices::dev.off()
}

draw_fig1 <- function() {
  par(mar = c(1.2, 1.2, 1.7, 1.2), xpd = NA)
  plot.new(); plot.window(xlim = c(0, 10), ylim = c(0, 10))
  title("Five diagnostics for suspected-drug attribution contrasts", line = 0.15,
        cex.main = 1.2, font.main = 2, col.main = pal$ink)
  boxes <- list(
    c(0.7, 6.5, 2.4, 8.2, "1  Role-consistent\nexposure"),
    c(3.0, 6.5, 4.7, 8.2, "2  Leave-event-out\nreference"),
    c(5.3, 6.5, 7.0, 8.2, "3  Medication\nburden"),
    c(7.6, 6.5, 9.3, 8.2, "4  Report-country\ncomposition"),
    c(2.55, 3.5, 4.45, 5.2, "5  Role construct\nPS / SS audit"),
    c(5.55, 3.5, 7.45, 5.2, "Bounded reporting\ndiagnostic")
  )
  for (i in seq_along(boxes)) {
    b <- boxes[[i]]; xy <- as.numeric(b[1:4]); fill <- if (i == 6) pal$context else pal$light
    rect(xy[1], xy[2], xy[3], xy[4], col = fill, border = if (i == 6) pal$context else pal$grid, lwd = 1.4)
    text((xy[1] + xy[3]) / 2, (xy[2] + xy[4]) / 2, b[5], cex = 0.82, col = if (i == 6) "white" else pal$ink, font = if (i == 6) 2 else 1)
  }
  arrows(2.48, 7.35, 2.92, 7.35, length = 0.08, col = pal$muted, lwd = 1.4)
  arrows(4.78, 7.35, 5.22, 7.35, length = 0.08, col = pal$muted, lwd = 1.4)
  arrows(7.08, 7.35, 7.52, 7.35, length = 0.08, col = pal$muted, lwd = 1.4)
  arrows(8.45, 6.42, 6.85, 5.30, length = 0.08, col = pal$muted, lwd = 1.4)
  arrows(4.50, 4.35, 5.47, 4.35, length = 0.08, col = pal$muted, lwd = 1.4)
  text(5, 1.9, "Output: a reporting-attribution diagnostic, not a clinical culprit ranking or risk estimate.",
       cex = 0.83, col = pal$ink)
}

draw_fig2 <- function() {
  key <- c("database", "drug_class", "agent")
  d <- merge(legacy[, c(key, "percentage_point_difference", "ror")],
             context[, c(key, "difference_pp", "ror")], by = key, suffixes = c("_unrestricted", "_context"), all = FALSE)
  d$label <- paste(d$database, agent_label(d$agent), sep = ": ")
  d <- d[order(match(d$database, c("FAERS", "Canada", "JADER")), match(d$agent, c("rocuronium", "succinylcholine", "cefazolin", "teicoplanin"))), ]
  y <- rev(seq_len(nrow(d)))
  par(mfrow = c(1, 2), mar = c(4.1, 10.2, 2.8, 1.2), oma = c(0, 0, 1.0, 0))
  plot(NA, xlim = c(-40, 52), ylim = c(0.5, nrow(d) + 0.5), yaxt = "n", xlab = "Core minus non-core agent share (percentage points)", ylab = "", main = "A. Absolute scale")
  abline(v = 0, col = pal$grid, lwd = 1.2)
  axis(2, at = y, labels = d$label, las = 1, cex.axis = 0.74)
  for (i in seq_len(nrow(d))) {
    segments(d$percentage_point_difference[i], y[i], d$difference_pp[i], y[i], col = pal$grid, lwd = 2)
    points(d$percentage_point_difference[i], y[i], pch = 21, bg = pal$legacy, col = pal$legacy, cex = 1.05)
    points(d$difference_pp[i], y[i], pch = 21, bg = pal$context, col = pal$context, cex = 1.05)
  }
  legend("bottomleft", inset = c(0, 0.035), bty = "n", horiz = TRUE, cex = 0.75,
         legend = c("Unrestricted background", "Target-independent context"),
         pch = 21, pt.bg = c(pal$legacy, pal$context), col = c(pal$legacy, pal$context))
  valid <- is.finite(d$ror_context) & is.finite(d$ror_unrestricted)
  plot(NA, log = "x", xlim = c(0.03, 15), ylim = c(0.5, nrow(d) + 0.5), yaxt = "n", xaxt = "n",
       xlab = "Reporting odds ratio (log scale)", ylab = "", main = "B. Relative scale")
  abline(v = 1, col = pal$grid, lwd = 1.2)
  axis(1, at = c(0.05, 0.1, 0.25, 0.5, 1, 2, 5, 10), labels = c("0.05", "0.1", "0.25", "0.5", "1", "2", "5", "10"), cex.axis = 0.65)
  for (i in which(valid)) {
    segments(d$ror_unrestricted[i], y[i], d$ror_context[i], y[i], col = pal$grid, lwd = 2)
    points(d$ror_unrestricted[i], y[i], pch = 21, bg = pal$legacy, col = pal$legacy, cex = 1.05)
    points(d$ror_context[i], y[i], pch = 21, bg = pal$context, col = pal$context, cex = 1.05)
  }
  if (any(!valid)) text(0.04, y[which(!valid)[1]], "Sparse\nnot estimated", pos = 4, cex = 0.62, col = pal$sparse)
  mtext("Agent contrasts before and after the target-independent medication context", outer = TRUE, line = -0.2, cex = 1.05, font = 2, col = pal$ink)
}

draw_fig3 <- function() {
  b <- burden[!burden$sparse, ]
  b$label <- b$drug_sequence_stratum
  chl <- roles[roles$drug_class == "chlorhexidine" & roles$role_profile %in% c("PS present", "SS present without PS"), ]
  chl$arm <- ifelse(chl$outcome_group == "core", "Core", "Non-core")
  ps <- chl$report_n[chl$role_profile == "PS present"]
  ss <- chl$report_n[chl$role_profile == "SS present without PS"]
  ps <- setNames(ps, chl$arm[chl$role_profile == "PS present"])
  ss <- setNames(ss, chl$arm[chl$role_profile == "SS present without PS"])
  par(mfrow = c(1, 2), mar = c(4.1, 12.3, 2.8, 1.2), oma = c(0, 0, 1.0, 0))
  y <- rev(seq_len(nrow(b)))
  plot(NA, log = "x", xlim = c(5, 35), ylim = c(0.5, nrow(b) + 0.5), yaxt = "n",
       xlab = "Reporting odds ratio (log scale)", ylab = "", main = "A. Medication-burden strata")
  abline(v = 1, col = pal$grid)
  axis(2, at = y, labels = paste0(b$label, "  a/b/c/d: ", b$a_target_suspect_core, "/", b$b_target_suspect_noncore, "/", b$c_non_target_core, "/", b$d_non_target_noncore), las = 1, cex.axis = 0.62)
  for (i in seq_len(nrow(b))) {
    segments(b$ror_ci_low[i], y[i], b$ror_ci_high[i], y[i], col = pal$context, lwd = 2)
    points(b$ror[i], y[i], pch = 21, bg = pal$context, col = pal$context, cex = 1.1)
    text(33.5, y[i], sprintf("%.1f", b$ror[i]), pos = 2, cex = 0.75, col = pal$ink)
  }
  mids <- barplot(rbind(ps[c("Core", "Non-core")], ss[c("Core", "Non-core")]), beside = FALSE,
                 col = c(pal$ps, pal$ss), border = NA, names.arg = c("Core", "Non-core"),
                 ylab = "Chlorhexidine suspect reports", main = "B. Reconstructed FAERS role profile", ylim = c(0, 1500))
  text(mids, colSums(rbind(ps[c("Core", "Non-core")], ss[c("Core", "Non-core")])) + 55,
       labels = colSums(rbind(ps[c("Core", "Non-core")], ss[c("Core", "Non-core")])), cex = 0.78, col = pal$ink)
  legend("topleft", bty = "n", cex = 0.78, fill = c(pal$ps, pal$ss), legend = c("Primary suspect present", "Secondary suspect only"))
  mtext("Chlorhexidine burden and role-construct diagnostics", outer = TRUE, line = -0.2, cex = 1.05, font = 2, col = pal$ink)
}

draw_fig4 <- function() {
  agent_order <- c("rocuronium", "succinylcholine", "cefazolin", "teicoplanin")
  stratum_order <- c("GB", "Non-GB observed", "Missing OCCR_COUNTRY")
  country$agent <- factor(country$agent, levels = agent_order)
  country$country_stratum <- factor(country$country_stratum, levels = stratum_order)
  country <- country[order(country$agent, country$country_stratum), ]
  y_base <- rep(rev(seq(2.0, 9.5, by = 2.5)), each = 3)
  offsets <- c(0.36, 0, -0.36)
  y <- y_base + rep(offsets, 4)
  par(mfrow = c(2, 1), mar = c(3.4, 7.6, 2.3, 1.2), oma = c(1.2, 0, 1.0, 0))
  plot(NA, log = "x", xlim = c(0.04, 30), ylim = c(0.7, 10.8), yaxt = "n", xaxt = "n", xlab = "", ylab = "", main = "")
  mtext("A", side = 3, line = 0.4, adj = 0, cex = 1.1, font = 2, col = pal$ink)
  abline(v = 1, col = pal$grid, lwd = 1.2)
  axis(1, at = c(0.05, 0.1, 0.25, 0.5, 1, 2, 5, 10, 20), labels = c("0.05", "0.1", "0.25", "0.5", "1", "2", "5", "10", "20"), cex.axis = 0.67)
  axis(2, at = rev(seq(2.0, 9.5, by = 2.5)), labels = agent_label(agent_order), las = 1, cex.axis = 0.83)
  for (i in seq_len(nrow(country))) {
    clr <- if (country$drug_class[i] == "NMBA") pal$nmba else pal$antibiotic
    pch <- c(21, 22, 24)[as.integer(country$country_stratum[i])]
    if (country$sparse_display[i] == 1) {
      points(country$ror[i], y[i], pch = pch, bg = "white", col = pal$sparse, cex = 1.05)
      text(23, y[i], "sparse", pos = 2, cex = 0.65, col = pal$sparse)
    } else {
      segments(country$ror_ci_low[i], y[i], country$ror_ci_high[i], y[i], col = clr, lwd = 1.8)
      points(country$ror[i], y[i], pch = pch, bg = clr, col = clr, cex = 1.0)
    }
  }
  legend("bottomleft", bty = "n", horiz = TRUE, cex = 0.72, pch = 21, pt.bg = c(pal$nmba, pal$antibiotic, pal$sparse), col = c(pal$nmba, pal$antibiotic, pal$sparse), legend = c("NMBA", "Antibiotic", "Sparse contrast"))
  legend("bottomright", bty = "n", cex = 0.7, pch = c(21, 22, 24), pt.bg = "white", col = pal$ink, legend = c("GB", "Observed non-GB", "Missing country"))
  comp <- composition[composition$outcome_group == "Core antibiotic reports" & composition$country_stratum %in% c("GB", "Non-GB observed"), ]
  comp$agent <- factor(comp$agent, levels = c("cefazolin", "teicoplanin"))
  comp$country_stratum <- factor(comp$country_stratum, levels = c("GB", "Non-GB observed"))
  mat <- matrix(0, nrow = 2, ncol = 2, dimnames = list(c("Cefazolin", "Teicoplanin"), c("GB", "Observed non-GB")))
  comp$plot_country <- ifelse(as.character(comp$country_stratum) == "Non-GB observed", "Observed non-GB", as.character(comp$country_stratum))
  for (i in seq_len(nrow(comp))) mat[agent_label(as.character(comp$agent[i])), comp$plot_country[i]] <- comp$agent_share_percent[i]
  mids <- barplot(mat, beside = TRUE, col = c(pal$cefa, pal$teico), border = NA, ylim = c(0, 52), ylab = "Agent share within core antibiotic reports (%)", main = "")
  mtext("B", side = 3, line = 0.4, adj = 0, cex = 1.1, font = 2, col = pal$ink)
  text(mids, as.vector(mat) + 2.1, labels = sprintf("%.1f%%", as.vector(mat)), cex = 0.75, col = pal$ink)
  legend("topright", bty = "n", cex = 0.78, fill = c(pal$cefa, pal$teico), legend = c("Cefazolin", "Teicoplanin"))
}

open_devices("Figure_1_five_diagnostics_v0_23", 10.0, 5.7, draw_fig1)
open_devices("Figure_2_reference_setting_contrasts_v0_23", 13.5, 7.1, draw_fig2)
open_devices("Figure_3_chlorhexidine_diagnostics_v0_23", 12.8, 6.8, draw_fig3)
open_devices("Figure_4_country_composition_v0_23", 10.5, 9.0, draw_fig4)

writeLines(c(
  "Public-release figure package",
  "All numerical panels are generated from the aggregate CSV artifacts; sparse annotations are display-only when a or c is below 5.",
  "Figure 2 unrestricted-background points are retained only as a reference-setting comparison; target-independent-context values are the current estimates.",
  "PDF files are submission-quality vector graphics. PNG files are 300 dpi review previews."
), file.path(out_dir, "README.md"))

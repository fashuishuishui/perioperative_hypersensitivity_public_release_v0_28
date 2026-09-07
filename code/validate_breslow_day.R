# Independent cross-check of the Breslow-Day family of statistics.
# The test reads the locked v0.26 country strata instead of retaining a second
# hard-coded copy of the 2 x 2 tables.

root <- Sys.getenv("PV_REPO_ROOT", unset = normalizePath(file.path(getwd()), winslash = "/", mustWork = FALSE))
out_dir <- Sys.getenv("PV_ANALYSIS_OUTPUT_DIR", unset = file.path(root, "aggregate_outputs", "schema_fixed"))
country_path <- file.path(out_dir, "canonical_agent_country_v0_26.csv")
if (!file.exists(country_path)) {
  stop("Missing locked country-strata input: ", country_path)
}
country <- read.csv(country_path, check.names = FALSE, stringsAsFactors = FALSE)
required_columns <- c("drug_class", "agent", "country_stratum", "a", "b", "c", "d")
if (!all(required_columns %in% names(country))) {
  stop("Locked country-strata file is missing required columns")
}
stratum_order <- c("GB", "Non-GB observed", "Missing OCCR_COUNTRY")
agents <- unique(country[, c("drug_class", "agent")])
country_tables <- lapply(seq_len(nrow(agents)), function(index) {
  subset <- country[country$drug_class == agents$drug_class[index] & country$agent == agents$agent[index], ]
  subset <- subset[match(stratum_order, subset$country_stratum), ]
  if (nrow(subset) != length(stratum_order) || any(is.na(subset$country_stratum))) {
    stop("Incomplete country strata for ", agents$agent[index])
  }
  lapply(seq_len(nrow(subset)), function(row) as.numeric(subset[row, c("a", "b", "c", "d")]))
})
names(country_tables) <- agents$agent

moments <- function(tab, theta) {
  a <- tab[1]; b <- tab[2]; c <- tab[3]; d <- tab[4]
  row1 <- a + b; row2 <- c + d; col1 <- a + c
  x <- max(0, col1 - row2):min(row1, col1)
  log_weights <- lchoose(row1, x) + lchoose(row2, col1 - x) + x * log(theta)
  weights <- exp(log_weights - max(log_weights)); weights <- weights / sum(weights)
  mean <- sum(x * weights)
  c(mean = mean, variance = sum((x - mean)^2 * weights))
}

conditional_mle_bd <- function(tables) {
  observed <- sum(vapply(tables, function(x) x[1], numeric(1)))
  score <- function(log_theta) sum(vapply(tables, function(x) moments(x, exp(log_theta))[1], numeric(1))) - observed
  theta <- exp(uniroot(score, interval = c(-30, 30))$root)
  m <- t(vapply(tables, function(x) moments(x, theta), numeric(2)))
  a <- vapply(tables, function(x) x[1], numeric(1))
  statistic <- sum((a - m[, "mean"])^2 / m[, "variance"])
  data.frame(conditional_mle_common_or = theta, conditional_mle_bd_statistic = statistic, conditional_mle_bd_p = pchisq(statistic, length(tables) - 1, lower.tail = FALSE))
}

desc_env <- new.env(parent = globalenv())
lazy_path <- file.path(
  Sys.getenv("PV_R_LIBS_USER", unset = ""),
  "DescTools", "R", "DescTools"
)
if (file.exists(paste0(lazy_path, ".rdb"))) {
  invisible(lazyLoad(lazy_path, desc_env))
} else if (requireNamespace("DescTools", quietly = TRUE)) {
  desc_env <- asNamespace("DescTools")
} else {
  stop("DescTools is required; set PV_R_LIBS_USER or install the declared package.")
}

out <- do.call(rbind, lapply(names(country_tables), function(agent) {
  tables <- country_tables[[agent]]
  array_input <- array(unlist(tables), dim = c(2, 2, length(tables)))
  uncorrected <- desc_env$BreslowDayTest(array_input, correct = FALSE)
  tarone <- desc_env$BreslowDayTest(array_input, correct = TRUE)
  cbind(agent = agent, conditional_mle_bd(tables),
        desctools_mh_bd_statistic = unname(uncorrected$statistic),
        desctools_mh_bd_p = uncorrected$p.value,
        desctools_tarone_statistic = unname(tarone$statistic),
        desctools_tarone_p = tarone$p.value)
}))
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
write.csv(out, file.path(out_dir, "breslow_day_variant_validation_v0_26.csv"), row.names = FALSE)
print(out)

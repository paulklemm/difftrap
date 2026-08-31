# Shared test fixtures.

make_col_data <- function(treatment_levels, source_levels = c("Input", "IP")) {
  grid <- expand.grid(
    treatment = treatment_levels,
    source = source_levels,
    stringsAsFactors = FALSE
  )
  tibble::tibble(
    sample_name = paste0("s", seq_len(nrow(grid))),
    treatment = factor(grid$treatment, levels = treatment_levels),
    source = factor(grid$source, levels = source_levels)
  )
}

make_settings <- function(treatment_groups, source_groups,
                          treatment_comparison, source_comparison) {
  list(
    treatment = treatment_groups,
    source = source_groups,
    treatment_comparison = as.list(treatment_comparison),
    source_comparison = as.list(source_comparison)
  )
}

# Negative-binomial counts over enough features that DESeq2's parametric
# dispersion fit converges; pure Poisson counts make it fail.
make_txi <- function(cd, n_genes = 200, n_zero = 50) {
  set.seed(1)
  mu <- rep(2^runif(n_genes, 6, 11), each = nrow(cd))
  counts <- matrix(stats::rnbinom(n_genes * nrow(cd), mu = mu, size = 5),
                   nrow = n_genes, byrow = TRUE,
                   dimnames = list(paste0("g", seq_len(n_genes)), cd$sample_name))
  counts[seq_len(n_zero), ] <- 0L
  list(counts = counts, abundance = counts,
       length = matrix(1000, nrow = n_genes, ncol = nrow(cd),
                       dimnames = dimnames(counts)),
       countsFromAbundance = "no")
}

make_replicated_col_data <- function(treatment_levels = c("WT", "KO")) {
  cd <- make_col_data(treatment_levels)
  cd <- cd[rep(seq_len(nrow(cd)), 3), ]
  cd$sample_name <- paste0("s", seq_len(nrow(cd)))
  cd
}

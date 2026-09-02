# Regression tests for the review findings on commit b54d4c8: coefficient
# naming, covariate validation, subset designs, and the low-count filter.

make_config_json <- function() {
  cfg <- list(
    treatment = list(WT = c("s1", "s2"), KO = c("s3", "s4")),
    source = list(Input = c("s1", "s3"), IP = c("s2", "s4")),
    treatment_comparison = c("KO", "WT"),
    source_comparison = c("IP", "Input")
  )
  path <- tempfile(fileext = ".json")
  jsonlite::write_json(cfg, path)
  path
}

# --- coefficient names ------------------------------------------------------

test_that("contrast coefs match DESeq2 resultsNames for digit-leading levels", {
  cd <- make_replicated_col_data(c("ctrl", "6h"))
  dds <- create_deseq_dataset(make_txi(cd), cd, level = "gene",
                              design = ~ source + treatment + source:treatment,
                              verbose = FALSE)
  ci <- define_contrasts(extract_factor_levels(cd, print_table = FALSE),
                         verbose = FALSE)
  expect_true(all(ci$coef %in% DESeq2::resultsNames(dds)))
})

# --- covariate validation ---------------------------------------------------

test_that("run_deseq2_pipeline rejects duplicated covariate rows", {
  covariates <- tibble::tibble(sample_name = c("s1", "s1", "s2", "s3", "s4"),
                               weight = 1:5)
  expect_error(
    run_deseq2_pipeline(make_config_json(), "nope.gtf",
                        covariates = covariates, verbose = FALSE),
    "duplicated"
  )
})

test_that("run_deseq2_pipeline rejects covariate columns that collide", {
  covariates <- tibble::tibble(sample_name = paste0("s", 1:4),
                               treatment = "oops")
  expect_error(
    run_deseq2_pipeline(make_config_json(), "nope.gtf",
                        covariates = covariates, verbose = FALSE),
    "collide"
  )
})

# --- design handling --------------------------------------------------------

test_that("subset_design keeps transformed and interacting covariate terms", {
  full <- ~ log(weight) + source + treatment + source:treatment
  expect_equal(subset_design(full, "treatment"), ~ log(weight) + source,
               ignore_formula_env = TRUE)
  full2 <- ~ sex:batch + source + treatment + source:treatment
  expect_equal(subset_design(full2, "treatment"), ~ source + sex:batch,
               ignore_formula_env = TRUE)
})

test_that("subset_design never reintroduces an absent factor", {
  expect_equal(subset_design(~ batch + treatment, "treatment"), ~ batch,
               ignore_formula_env = TRUE)
})

test_that("run_deseq2_pipeline rejects designs without the factorial terms", {
  expect_error(
    run_deseq2_pipeline("nope.json", "nope.gtf",
                        design = ~ batch + source + treatment),
    "source:treatment"
  )
})

# --- low-count filter -------------------------------------------------------

test_that("min_samples without min_count is an error, not a silent no-op", {
  cd <- make_replicated_col_data()
  expect_error(
    create_deseq_dataset(make_txi(cd), cd, level = "gene",
                         design = ~ source + treatment,
                         min_samples = 3, verbose = FALSE),
    "min_count"
  )
})

test_that("explicit min_samples is capped at the dataset's sample count", {
  cd <- make_replicated_col_data()
  dds <- create_deseq_dataset(make_txi(cd), cd, level = "gene",
                              design = ~ source + treatment,
                              min_count = 10, min_samples = 50, verbose = FALSE)
  expect_gt(nrow(dds), 0)
})

test_that("default min_samples derives cells from the design's factors", {
  cd <- make_replicated_col_data()
  cd2 <- tibble::tibble(sample_name = cd$sample_name,
                        condition = factor(rep(c("a", "b"), each = 6)))
  dds <- create_deseq_dataset(make_txi(cd2), cd2, level = "gene",
                              design = ~ condition,
                              min_count = 10, verbose = FALSE)
  expect_equal(nrow(dds), 150)
})

test_that("default min_samples errors when the design has no factors", {
  cd3 <- tibble::tibble(sample_name = paste0("s", 1:12),
                        dose = rnorm(12))
  expect_error(
    suppressWarnings(
      create_deseq_dataset(make_txi(cd3), cd3, level = "gene",
                           design = ~ dose, min_count = 10, verbose = FALSE)
    ),
    "min_samples"
  )
})

# --- contrast input validation ----------------------------------------------

test_that("define_contrasts rejects factor_levels without treatment_levels", {
  old_style <- list(treatment_a = "WT", treatment_b = "KO",
                    source_a = "Input", source_b = "IP")
  expect_error(define_contrasts(old_style, verbose = FALSE), "treatment_levels")
})

# --- shrinkage --------------------------------------------------------------

test_that("run_deseq2_pipeline rejects an unknown shrinkage method early", {
  expect_error(
    run_deseq2_pipeline(make_config_json(), "nope.gtf",
                        shrinkage = "foo", verbose = FALSE),
    "apeglm"
  )
})

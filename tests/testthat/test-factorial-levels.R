# `source` stays a two-level factor (IP/Input), but `treatment` may carry more
# than two levels. The two-level case must keep producing byte-identical
# contrast definitions, because published analyses filter on these effect names.

# Fixtures live in helper-data.R.

# --- two levels: the existing contract ---------------------------------------

test_that("extract_factor_levels keeps the two-level keys", {
  fl <- extract_factor_levels(make_col_data(c("WT", "KO")), print_table = FALSE)

  expect_equal(fl$treatment_a, "WT")
  expect_equal(fl$treatment_b, "KO")
  expect_equal(fl$source_a, "Input")
  expect_equal(fl$source_b, "IP")
  expect_equal(fl$treatment_levels, c("WT", "KO"))
})

test_that("define_contrasts is unchanged for a 2x2 design", {
  fl <- extract_factor_levels(make_col_data(c("WT", "KO")), print_table = FALSE)
  ci <- define_contrasts(fl, verbose = FALSE)

  expect_equal(nrow(ci), 5)
  expect_equal(
    ci$effect_name,
    c("IP_vs_Input_in_WT", "IP_vs_Input_in_KO",
      "KO_vs_WT_in_Input", "KO_vs_WT_in_IP", "interaction_effect")
  )
  expect_equal(
    ci$coef,
    c("source_IP_vs_Input", "source_IP_vs_Input",
      "treatment_KO_vs_WT", "treatment_KO_vs_WT", "sourceIP.treatmentKO")
  )
  expect_equal(
    ci$model_id,
    c("subset_treatment_WT", "subset_treatment_KO",
      "subset_source_Input", "subset_source_IP", "full")
  )
  expect_equal(ci$subset_variable,
               c("treatment", "treatment", "source", "source", NA_character_))
  expect_true(all(ci$shrinkage_method == "apeglm"))
})

# --- three levels: the extension ---------------------------------------------

test_that("extract_factor_levels accepts more than two treatment levels", {
  cd <- make_col_data(c("ctrl", "camk2a", "rosa26"), c("noad", "ad"))
  fl <- extract_factor_levels(cd, print_table = FALSE)

  expect_equal(fl$treatment_levels, c("ctrl", "camk2a", "rosa26"))
  expect_equal(fl$source_a, "noad")
  expect_equal(fl$source_b, "ad")
  # Ambiguous with more than two levels, so not reported.
  expect_null(fl$treatment_a)
  expect_null(fl$treatment_b)
})

test_that("extract_factor_levels still requires exactly two source levels", {
  cd <- make_col_data(c("WT", "KO"), c("a", "b", "c"))
  expect_error(extract_factor_levels(cd, print_table = FALSE), "source")
})

test_that("extract_factor_levels rejects a single treatment level", {
  cd <- make_col_data("WT")
  expect_error(extract_factor_levels(cd, print_table = FALSE), "treatment")
})

test_that("define_contrasts emits one interaction per non-reference level", {
  cd <- make_col_data(c("ctrl", "camk2a", "rosa26"), c("noad", "ad"))
  ci <- define_contrasts(extract_factor_levels(cd, print_table = FALSE),
                         verbose = FALSE)

  # 3 source-effects (one per treatment level)
  # + 4 treatment-effects (2 non-reference levels x 2 source levels)
  # + 2 interactions
  expect_equal(nrow(ci), 9)

  expect_equal(
    ci$effect_name,
    c("ad_vs_noad_in_ctrl", "ad_vs_noad_in_camk2a", "ad_vs_noad_in_rosa26",
      "camk2a_vs_ctrl_in_noad", "rosa26_vs_ctrl_in_noad",
      "camk2a_vs_ctrl_in_ad", "rosa26_vs_ctrl_in_ad",
      "interaction_effect_camk2a", "interaction_effect_rosa26")
  )
  expect_equal(
    ci$coef[ci$model_type == "full"],
    c("sourcead.treatmentcamk2a", "sourcead.treatmentrosa26")
  )
  # Both treatment effects within one source level come from the same subset
  # model: with three levels that model carries both coefficients.
  expect_equal(
    ci$model_id[ci$effect_name %in% c("camk2a_vs_ctrl_in_noad",
                                      "rosa26_vs_ctrl_in_noad")],
    c("subset_source_noad", "subset_source_noad")
  )
  expect_true(all(ci$shrinkage_method == "apeglm"))
})

# --- config validation --------------------------------------------------------

test_that("get_coldata accepts a treatment_comparison listing three levels", {
  settings <- make_settings(
    treatment_groups = list(ctrl = c("s1", "s2"), camk2a = c("s3", "s4"),
                            rosa26 = c("s5", "s6")),
    source_groups = list(noad = c("s1", "s3", "s5"), ad = c("s2", "s4", "s6")),
    treatment_comparison = c("rosa26", "camk2a", "ctrl"),
    source_comparison = c("ad", "noad")
  )
  cd <- get_coldata(settings)

  # Reversed, so the last entry becomes the DESeq2 reference.
  expect_equal(levels(cd$treatment), c("ctrl", "camk2a", "rosa26"))
  expect_equal(levels(cd$source), c("noad", "ad"))
})

test_that("get_coldata still requires source_comparison to name two levels", {
  settings <- make_settings(
    treatment_groups = list(ctrl = c("s1", "s2"), camk2a = c("s3", "s4")),
    source_groups = list(a = "s1", b = "s2", c = c("s3", "s4")),
    treatment_comparison = c("camk2a", "ctrl"),
    source_comparison = c("c", "b", "a")
  )
  expect_error(get_coldata(settings), "source_comparison")
})

test_that("get_coldata still requires every observed treatment level", {
  settings <- make_settings(
    treatment_groups = list(ctrl = c("s1", "s2"), camk2a = c("s3", "s4"),
                            rosa26 = c("s5", "s6")),
    source_groups = list(noad = c("s1", "s3", "s5"), ad = c("s2", "s4", "s6")),
    treatment_comparison = c("camk2a", "ctrl"),
    source_comparison = c("ad", "noad")
  )
  expect_error(get_coldata(settings), "rosa26")
})

# --- covariates ---------------------------------------------------------------

test_that("subset_design drops the held factor and keeps the other", {
  full <- ~ source + treatment + source:treatment
  expect_equal(subset_design(full, "treatment"), ~ source, ignore_formula_env = TRUE)
  expect_equal(subset_design(full, "source"), ~ treatment, ignore_formula_env = TRUE)
})

test_that("subset_design carries covariates into the subset models", {
  full <- ~ cp + source + treatment + source:treatment
  expect_equal(subset_design(full, "treatment"), ~ cp + source,
               ignore_formula_env = TRUE)
  expect_equal(subset_design(full, "source"), ~ cp + treatment,
               ignore_formula_env = TRUE)
})

test_that("subset_design keeps several covariates in order", {
  full <- ~ cp + quality + source + treatment + source:treatment
  expect_equal(subset_design(full, "treatment"), ~ cp + quality + source,
               ignore_formula_env = TRUE)
})

# --- low-count filter ---------------------------------------------------------

test_that("create_deseq_dataset keeps every feature by default", {
  cd <- make_replicated_col_data()
  dds <- create_deseq_dataset(make_txi(cd), cd, level = "gene",
                              design = ~ source + treatment, verbose = FALSE)
  expect_equal(nrow(dds), 200)
})

test_that("create_deseq_dataset drops low-count features when asked", {
  cd <- make_replicated_col_data()
  dds <- create_deseq_dataset(make_txi(cd), cd, level = "gene",
                              design = ~ source + treatment,
                              min_count = 10, min_samples = 2, verbose = FALSE)
  expect_equal(nrow(dds), 150)
  expect_false(any(paste0("g", 1:50) %in% rownames(dds)))
})

test_that("min_samples defaults to the smallest design cell", {
  cd <- make_replicated_col_data()
  # 2x2 with 3 replicates each, so the smallest cell is 3.
  expect_equal(min(table(interaction(cd$source, cd$treatment, drop = TRUE))), 3L)
  dds <- create_deseq_dataset(make_txi(cd), cd, level = "gene",
                              design = ~ source + treatment,
                              min_count = 10, verbose = FALSE)
  expect_equal(nrow(dds), 150)
})

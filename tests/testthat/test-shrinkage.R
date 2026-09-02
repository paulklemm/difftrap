# The interaction coefficient's apeglm MAP is bimodal at small n (the estimate
# flips between the MLE and ~0 depending on which samples are in the fit), so
# callers must be able to choose ashr or no shrinkage, and always get the MLE.

fit_full <- function() {
  cd <- make_replicated_col_data()
  create_deseq_dataset(make_txi(cd), cd, level = "gene",
                       design = ~ source + treatment + source:treatment,
                       verbose = FALSE)
}

test_that("define_contrasts passes the shrinkage method through", {
  fl <- extract_factor_levels(make_col_data(c("WT", "KO")), print_table = FALSE)
  expect_true(all(define_contrasts(fl, verbose = FALSE)$shrinkage_method == "apeglm"))
  expect_true(all(define_contrasts(fl, shrinkage = "ashr", verbose = FALSE)$shrinkage_method == "ashr"))
  expect_true(all(define_contrasts(fl, shrinkage = "none", verbose = FALSE)$shrinkage_method == "none"))
  expect_error(define_contrasts(fl, shrinkage = "foo", verbose = FALSE))
})

test_that("extract_de_results returns the plain Wald results for 'none'", {
  dds <- fit_full()
  ci <- tibble::tibble(effect_name = "KO_vs_WT", coef = "treatment_KO_vs_WT",
                       shrinkage_method = "none", model_id = "full")
  res <- extract_de_results(list(full = dds), ci)$KO_vs_WT
  ref <- DESeq2::results(dds, name = "treatment_KO_vs_WT")
  expect_equal(res$log2FoldChange, ref$log2FoldChange)
  expect_equal(res$padj, ref$padj)
})

test_that("shrunk results carry the MLE and its SE alongside", {
  dds <- fit_full()
  ci <- tibble::tibble(effect_name = "KO_vs_WT", coef = "treatment_KO_vs_WT",
                       shrinkage_method = "ashr", model_id = "full")
  res <- extract_de_results(list(full = dds), ci)$KO_vs_WT
  ref <- DESeq2::results(dds, name = "treatment_KO_vs_WT")
  expect_equal(res$log2FoldChange_MLE, ref$log2FoldChange)
  expect_equal(res$lfcSE_MLE, ref$lfcSE)
  expect_equal(res$padj, ref$padj)
  expect_false(isTRUE(all.equal(res$log2FoldChange, ref$log2FoldChange)))

  tbl <- results_to_tibble(list(KO_vs_WT = res), ci, level = "gene")
  expect_true(all(c("log2FoldChange_MLE", "lfcSE_MLE") %in% colnames(tbl)))
})

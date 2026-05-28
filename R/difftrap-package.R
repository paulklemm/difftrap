#' @keywords internal
"_PACKAGE"

#' @importFrom magrittr %>%
NULL

# Silence R CMD check NOTEs about non-standard evaluation in dplyr / tibble
# verbs. These names are column references inside tidy-eval contexts, not
# globals.
utils::globalVariables(c(
  "description", "effect", "effect_name", "level", "log2FoldChange",
  "model_id", "model_type", "padj", "sample_name", "shrinkage_method",
  "subset_level", "subset_variable", "treatment", "tx"
))

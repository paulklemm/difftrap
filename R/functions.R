#' Get col_data for DESeq2 from JSON config
#'
#' Builds the sample-level metadata table from a parsed JSON config and
#' relevels the `treatment` and `source` factors so that DESeq2 uses the
#' intended reference (denominator) level.
#'
#' Expected `settings` structure:
#' * `treatment`: named list mapping each treatment group to a vector of
#'   sample names (e.g. `list(WT = c("s1","s2"), KO = c("s3","s4"))`).
#' * `source`: named list mapping each source group (e.g. IP, Input) to a
#'   vector of sample names. Every sample must appear in exactly one
#'   treatment group **and** exactly one source group.
#' * `treatment_comparison`: length-2 character vector
#'   `c(numerator, denominator)`, e.g. `c("KO", "WT")` to test KO vs WT.
#' * `source_comparison`: length-2 character vector `c(numerator, denominator)`.
#'
#' Both `*_comparison` vectors must list exactly the levels present in the
#' corresponding group. The function fails fast if any required key is
#' missing, references an unknown level, or if a sample is missing from
#' either grouping.
#'
#' @param settings Configurations object from JSON
#' @return col_data tibble with treatment and source columns
#' @export
get_coldata <- function(settings) {
  required <- c("treatment", "source", "treatment_comparison", "source_comparison")
  missing_keys <- setdiff(required, names(settings))
  if (length(missing_keys) > 0) {
    stop("Config is missing required key(s): ", paste(missing_keys, collapse = ", "))
  }

  col_data_treatment <- dplyr::bind_rows(
    purrr::map(names(settings$treatment), function(group_name) {
      tibble::tibble(
        treatment = group_name,
        sample_name = unlist(settings$treatment[[group_name]])
      )
    })
  )

  col_data_source <- dplyr::bind_rows(
    purrr::map(names(settings$source), function(group_name) {
      tibble::tibble(
        source = group_name,
        sample_name = unlist(settings$source[[group_name]])
      )
    })
  )

  # Use full_join so unmatched samples surface as NAs we can report.
  col_data <- dplyr::full_join(
    x = col_data_source,
    y = col_data_treatment,
    by = "sample_name"
  )

  missing_treatment <- col_data$sample_name[is.na(col_data$treatment)]
  missing_source <- col_data$sample_name[is.na(col_data$source)]
  if (length(missing_treatment) > 0 || length(missing_source) > 0) {
    stop(
      "Every sample must appear in both 'source' and 'treatment' groupings.",
      if (length(missing_treatment) > 0)
        paste0("\n  Missing from treatment: ", paste(missing_treatment, collapse = ", ")),
      if (length(missing_source) > 0)
        paste0("\n  Missing from source: ", paste(missing_source, collapse = ", "))
    )
  }

  # DESeq2 uses the first factor level as the reference (denominator), so we
  # reverse the user-provided c(numerator, denominator) to put the denominator
  # in position 1.
  treatment_fct_levels <- rev(unlist(settings$treatment_comparison))
  source_fct_levels <- rev(unlist(settings$source_comparison))

  validate_comparison <- function(values, observed, name, fixed_length = TRUE) {
    if (fixed_length && length(values) != 2) {
      stop(name, " must be a length-2 vector c(numerator, denominator); got length ",
           length(values))
    }
    if (!fixed_length && length(values) < 2) {
      stop(name, " must list at least two levels, ordered so that the reference",
           " comes last; got length ", length(values))
    }
    unknown <- setdiff(values, observed)
    if (length(unknown) > 0) {
      stop(name, " references level(s) not present in the data: ",
           paste(unknown, collapse = ", "),
           "\n  Observed levels: ", paste(observed, collapse = ", "))
    }
    extra <- setdiff(observed, values)
    if (length(extra) > 0) {
      stop(name, " must list every observed level. Missing: ",
           paste(extra, collapse = ", "))
    }
  }
  validate_comparison(treatment_fct_levels, unique(col_data$treatment),
                      "treatment_comparison", fixed_length = FALSE)
  validate_comparison(source_fct_levels, unique(col_data$source), "source_comparison")

  col_data %>%
    dplyr::mutate(
      treatment = forcats::fct_relevel(forcats::as_factor(treatment), treatment_fct_levels),
      source = forcats::fct_relevel(forcats::as_factor(source), source_fct_levels)
    )
}



#' Add file paths to col_data
#'
#' @param col_data Column data tibble from get_coldata()
#' @param results_dir Directory that directly contains sample subfolders
#'   (e.g. "/nfcore-rnaseq-pipeline/results_grcm38/star_salmon/")
#' @param quant_file Quantification file name (default: "quant.sf")
#' @return col_data with filepath column added
#' @export
add_filepaths <- function(col_data,
                          results_dir = "/nfcore-rnaseq-pipeline/results_grcm38/star_salmon/",
                          quant_file = "quant.sf") {
  col_data %>%
    dplyr::mutate(
      filepath = file.path(results_dir, sample_name, quant_file)
    )
}



#' Load GTF annotation file
#'
#' @param gtf_path Path to GTF file
#' @return GTF object from rtracklayer
#' @export
load_gtf <- function(gtf_path) {
  if (!file.exists(gtf_path)) {
    stop("GTF file not found: ", gtf_path)
  }
  rtracklayer::import(gtf_path)
}



#' Create transcript-to-gene mapping from GTF
#'
#' @param gtf GTF object from rtracklayer::import()
#' @param verbose Print number of mappings (default: TRUE)
#' @return data.frame with tx and gene columns
#' @export
create_tx2gene <- function(gtf, verbose = TRUE) {
  tx2gene <- data.frame(
    tx = gtf$transcript_id,
    gene = gtf$gene_id
  ) %>%
    dplyr::filter(!is.na(tx)) %>%
    dplyr::distinct()

  # Downstream joins rely on a one-to-one tx -> gene mapping. Surface
  # ambiguous transcripts here with a clear message rather than letting a
  # generic dplyr "many-to-one" error fire deep in results_to_tibble().
  multi_gene_tx <- tx2gene$tx[duplicated(tx2gene$tx)]
  if (length(multi_gene_tx) > 0) {
    stop(length(unique(multi_gene_tx)),
         " transcript ID(s) map to more than one gene_id in the GTF, e.g.: ",
         paste(utils::head(unique(multi_gene_tx), 5), collapse = ", "),
         ". Resolve the ambiguity in the GTF before continuing.")
  }

  if (verbose) {
    message("tx2gene mapping: ", nrow(tx2gene), " transcript-gene pairs")
  }

  tx2gene
}



#' Import Salmon quantification data using tximport
#'
#' @param col_data Column data with filepath and sample_name columns
#' @param tx2gene Transcript-to-gene mapping data.frame
#' @param level Quantification level: "transcript" or "gene" (default: "transcript")
#' @param counts_from_abundance Abundance scaling method (default: "no").
#'   DESeq2 accounts for transcript length via the offset matrix that tximport
#'   provides, so scaling counts by abundance is unnecessary.
#' @return tximport object
#' @export
import_salmon_data <- function(col_data, tx2gene, 
                               level = "transcript",
                               counts_from_abundance = "no") {
  # Validate level parameter
  if (!level %in% c("transcript", "gene")) {
    stop("level must be 'transcript' or 'gene'")
  }
  
  # Extract files vector with names (required for tximport)
  files_vector <- col_data$filepath %>%
    purrr::set_names(col_data$sample_name)
  
  # Check if all files exist
  missing_files <- files_vector[!file.exists(files_vector)]
  if (length(missing_files) > 0) {
    stop("Missing quantification files:\n", paste(missing_files, collapse = "\n"))
  }

  tximport::tximport(
    files_vector,
    type = "salmon",
    tx2gene = tx2gene,
    txOut = level == "transcript",
    countsFromAbundance = counts_from_abundance
  )
}



#' Design formula for a subset model
#'
#' A subset model holds one of the two factors fixed, so it drops that factor
#' and keeps the other. Any covariates in the full design are kept, since a
#' nuisance term that needs adjusting in the full model needs adjusting here
#' too.
#'
#' @param design Full design formula.
#' @param subset_variable Factor held fixed: "treatment" or "source".
#' @return A one-sided formula.
#' @keywords internal
subset_design <- function(design, subset_variable) {
  kept <- if (subset_variable == "treatment") "source" else "treatment"
  covariates <- setdiff(all.vars(design), c("source", "treatment"))
  stats::reformulate(c(covariates, kept))
}

#' Extract factor levels from col_data
#'
#' `source` must have exactly two levels; `treatment` may have two or more.
#' The first level of each is the reference (denominator in DESeq2).
#'
#' `treatment_levels` holds every treatment level, reference first.
#' `treatment_a` and `treatment_b` are also returned when `treatment` has
#' exactly two levels, so existing two-level callers keep working; with more
#' than two they are ambiguous and are omitted.
#'
#' @param col_data Column data tibble with treatment and source factors
#' @param print_table Print summary table (default: TRUE)
#' @return Named list with treatment_levels, source_a, source_b, and -- for a
#'   two-level treatment -- treatment_a and treatment_b
#' @export
extract_factor_levels <- function(col_data, print_table = TRUE) {
  treatment_levels <- levels(col_data$treatment)
  source_levels <- levels(col_data$source)

  if (length(source_levels) != 2) {
    stop("Need exactly 2 levels for source (got ", length(source_levels), ")")
  }
  if (length(treatment_levels) < 2) {
    stop("Need at least 2 levels for treatment (got ",
         length(treatment_levels), ")")
  }

  factor_levels <- list(
    treatment_levels = treatment_levels,
    source_a = source_levels[1],
    source_b = source_levels[2]
  )
  if (length(treatment_levels) == 2) {
    factor_levels$treatment_a <- treatment_levels[1]
    factor_levels$treatment_b <- treatment_levels[2]
  }

  if (print_table) {
    summary_table <- tibble::tibble(
      treatment = treatment_levels,
      role = c("reference", rep("compared to reference",
                                length(treatment_levels) - 1))
    )
    print(summary_table)
    print(tibble::tibble(source = source_levels,
                         role = c("reference", "compared to reference")))
  }

  factor_levels
}



#' Create DESeq2 dataset from tximport object
#'
#' @param txi tximport object
#' @param col_data Column data for samples
#' @param level Quantification level used ("transcript" or "gene")
#' @param design Design formula (default: ~ source + treatment + source:treatment)
#' @param min_count Drop features without at least this many counts in
#'   `min_samples` samples. Default 0, which keeps every feature -- apeglm
#'   shrinks a noisy near-zero feature's LFC to 0 but leaves its Wald p-value
#'   untouched, so unfiltered results can carry "significant" features with a
#'   base mean near zero and no effect size.
#' @param min_samples Number of samples that must meet `min_count`. Default 0,
#'   in which case the smallest group in `col_data` is used.
#' @param verbose Print dimensions and size factors (default: TRUE)
#' @return DESeqDataSet object after running DESeq()
#' @export
create_deseq_dataset <- function(txi, col_data, 
                                 level = "transcript",
                                 design = ~ source + treatment + source:treatment,
                                 min_count = 0,
                                 min_samples = 0,
                                 verbose = TRUE) {
  dds <- DESeq2::DESeqDataSetFromTximport(
    txi,
    colData = col_data,
    design = design
  )

  if (min_count > 0) {
    if (min_samples <= 0) {
      min_samples <- min(table(interaction(col_data$source, col_data$treatment,
                                           drop = TRUE)))
    }
    keep <- rowSums(DESeq2::counts(dds) >= min_count) >= min_samples
    if (verbose) {
      message("Filtering: keeping ", sum(keep), " of ", length(keep),
              " features with >= ", min_count, " counts in >= ", min_samples,
              " samples")
    }
    dds <- dds[keep, ]
  }

  if (verbose) {
    feature_type <- if (level == "transcript") "transcripts" else "genes"
    message("DESeq2 object dimensions: ", dim(dds)[1], " ", feature_type,
            " x ", dim(dds)[2], " samples")
  }

  dds <- DESeq2::DESeq(dds)

  if (verbose) {
    # Check size factors (may be NULL when normalizationFactors are used, e.g. from tximport)
    size_factors <- DESeq2::sizeFactors(dds)
    if (!is.null(size_factors)) {
      message("Size factors (should be ~1): ",
              paste(round(size_factors, 3), collapse = " "))
    } else {
      message("Using normalization factors from tximport (avgTxLength correction)")
    }
  }

  dds
}



#' Define contrasts for differential expression analysis
#'
#' Defines contrasts using separate models for non-interaction effects and
#' the full interaction model only for the interaction terms. Non-interaction
#' contrasts (source effects within a treatment, treatment effects within a
#' source) use subset models that only include the relevant samples.
#' This avoids distorted dispersion estimates when the IP/Input difference is
#' large, and allows apeglm shrinkage for all contrasts.
#'
#' `source` must have exactly two levels. `treatment` may have more, in which
#' case every non-reference level gets its own treatment-effect contrast within
#' each source level, and its own interaction term. For k treatment levels that
#' is k source-effects + 2(k-1) treatment-effects + (k-1) interactions; at
#' k = 2 this is the original five contrasts, unchanged.
#'
#' With a single interaction the effect is named `interaction_effect`, as
#' before. With more than one they are named `interaction_effect_<level>`,
#' since one name can no longer identify them.
#'
#' Note: because the effect-within-group contrasts come from subset models
#' (per-subset dispersions and per-subset tximport offsets) while the
#' interaction contrasts come from the full model (pooled dispersions), an
#' interaction LFC is **not** equal to the algebraic difference of the two
#' source-in-treatment (or treatment-in-source) LFCs. Treat each interaction
#' as its own estimate from the full model rather than as a derived
#' difference.
#'
#' @param factor_levels List from [extract_factor_levels()], with
#'   treatment_levels, source_a and source_b
#' @param verbose Print contrast info (default: TRUE)
#' @return Tibble with contrast definitions including model_type, model_id,
#'   subset_variable, and subset_level columns
#' @export
define_contrasts <- function(factor_levels, verbose = TRUE) {
  treatment_levels <- factor_levels$treatment_levels
  treatment_ref <- treatment_levels[1]
  treatment_rest <- treatment_levels[-1]
  source_a <- factor_levels$source_a
  source_b <- factor_levels$source_b

  # DESeq2 creates coefficients as "level_b_vs_a" where a is the reference. We
  # use make.names() because DESeq2/model.matrix replaces spaces/chars with dots.
  source_coef_name <- paste0("source_", make.names(source_b), "_vs_",
                             make.names(source_a))
  source_coef_print <- paste0(source_b, "_vs_", source_a)

  # Source effect within each treatment level: one subset model per level,
  # each fitted as ~ source.
  source_effects <- tibble::tibble(
    effect_name = paste0(source_coef_print, "_in_", treatment_levels),
    description = paste("Source effect in", treatment_levels),
    coef = source_coef_name,
    shrinkage_method = "apeglm",
    model_type = "subset",
    model_id = paste0("subset_treatment_", make.names(treatment_levels)),
    subset_variable = "treatment",
    subset_level = treatment_levels
  )

  # Treatment effects within each source level. Both source subset models are
  # fitted as ~ treatment, so with k > 2 one model carries several coefficients.
  treatment_effects <- purrr::map_dfr(c(source_a, source_b), function(sl) {
    tibble::tibble(
      effect_name = paste0(treatment_rest, "_vs_", treatment_ref, "_in_", sl),
      description = paste0("Treatment effect (", treatment_rest, ") in ", sl),
      coef = paste0("treatment_", make.names(treatment_rest), "_vs_",
                    make.names(treatment_ref)),
      shrinkage_method = "apeglm",
      model_type = "subset",
      model_id = paste0("subset_source_", make.names(sl)),
      subset_variable = "source",
      subset_level = sl
    )
  })

  interactions <- tibble::tibble(
    effect_name = if (length(treatment_rest) == 1) {
      "interaction_effect"
    } else {
      paste0("interaction_effect_", treatment_rest)
    },
    description = paste0("Interaction term (", treatment_rest, ")"),
    coef = paste0("source", make.names(source_b), ".treatment",
                  make.names(treatment_rest)),
    shrinkage_method = "apeglm",
    model_type = "full",
    model_id = "full",
    subset_variable = NA_character_,
    subset_level = NA_character_
  )

  contrasts_info <- dplyr::bind_rows(source_effects, treatment_effects,
                                     interactions)

  if (verbose) {
    print(contrasts_info %>% dplyr::select(effect_name, description, model_type, model_id, shrinkage_method))
  }

  contrasts_info
}




#' Extract differential expression results with shrinkage
#'
#' Extracts results from the appropriate DESeq2 model for each contrast.
#' Each contrast specifies a model_id that maps to the corresponding
#' DESeqDataSet in dds_list.
#'
#' @param dds_list Named list of DESeqDataSet objects (e.g. "full" for
#'   interaction model, "subset_treatment_X" / "subset_source_X" for subset models)
#' @param contrasts_info Contrasts definition tibble from define_contrasts(),
#'   must contain coef, shrinkage_method, and model_id columns
#' @return List of DESeq2 results objects
#' @export
extract_de_results <- function(dds_list, contrasts_info) {
  results_list <- purrr::pmap(
    list(
      coef = contrasts_info$coef,
      method = contrasts_info$shrinkage_method,
      model_id = contrasts_info$model_id
    ),
    function(coef, method, model_id) {
      dds <- dds_list[[model_id]]
      if (is.null(dds)) {
        stop("No DESeq2 model found for model_id: ", model_id)
      }
      DESeq2::lfcShrink(dds, coef = coef, type = method)
    }
  )

  names(results_list) <- contrasts_info$effect_name
  results_list
}



#' Convert DESeq2 results to master tibble
#'
#' @param results_list List of DESeq2 results from extract_de_results()
#' @param contrasts_info Contrasts definition tibble
#' @param tx2gene Transcript-to-gene mapping (NULL for gene-level)
#' @param level Quantification level ("transcript" or "gene")
#' @return Combined tibble with all results
#' @export
results_to_tibble <- function(results_list, contrasts_info, tx2gene = NULL, level = "transcript") {
  if (level == "transcript" && is.null(tx2gene)) {
    stop("tx2gene mapping required for transcript-level results")
  }
  
  shrinkage_lookup <- stats::setNames(contrasts_info$shrinkage_method,
                                      contrasts_info$effect_name)

  purrr::imap(results_list, function(res, effect_name) {
    result_tbl <- tibble::as_tibble(res, rownames = "feature_id")

    if (level == "transcript") {
      result_tbl <- dplyr::left_join(
        result_tbl, tx2gene,
        by = c("feature_id" = "tx"),
        relationship = "many-to-one"
      )
    } else {
      result_tbl <- dplyr::rename(result_tbl, gene = "feature_id")
    }

    result_tbl %>%
      dplyr::mutate(
        effect = effect_name,
        shrinkage_method = shrinkage_lookup[[effect_name]],
        level = level
      ) %>%
      dplyr::arrange(padj)
  }) %>%
    dplyr::bind_rows()
}



#' Calculate summary statistics for DE results
#'
#' @param master_results Combined results tibble
#' @param alpha Significance threshold (default: 0.05)
#' @param verbose Print summary stats (default: TRUE)
#' @return Summary statistics tibble
#' @export
calculate_summary_stats <- function(master_results, alpha = 0.05, verbose = TRUE) {
  summary_stats <- master_results %>%
    dplyr::group_by(effect, shrinkage_method, level) %>%
    dplyr::summarise(
      significant = sum(padj < alpha, na.rm = TRUE),
      up = sum(log2FoldChange > 0 & padj < alpha, na.rm = TRUE),
      down = sum(log2FoldChange < 0 & padj < alpha, na.rm = TRUE),
      .groups = "drop"
    )
  
  if (verbose) {
    print(summary_stats)
  }

  summary_stats
}



#' Helper function to run single-level analysis
#'
#' Creates separate DESeq2 models for interaction and non-interaction contrasts.
#' The full interaction model (all samples) is used only for the interaction
#' term. Non-interaction contrasts use subset models with only the relevant
#' samples and a simpler design formula, providing better dispersion estimates
#' when IP/Input differences are large.
#'
#' @keywords internal
.run_single_level_analysis <- function(col_data, tx2gene, level,
                                       contrasts_info, design, alpha, verbose,
                                       min_count = 0, min_samples = 0) {
  if (verbose) message("=== Importing Salmon data (", level, "-level) ===")
  txi <- import_salmon_data(col_data, tx2gene, level = level)

  dds_list <- list()

  # --- Full interaction model (all samples) for interaction term ---
  full_contrasts <- contrasts_info %>% dplyr::filter(model_type == "full")
  if (nrow(full_contrasts) > 0) {
    if (verbose) message("=== Creating full DESeq2 dataset (interaction model, all samples) ===")
    dds_list[["full"]] <- create_deseq_dataset(
      txi, col_data, level = level, design = design,
      min_count = min_count, min_samples = min_samples, verbose = verbose
    )
  }

  # --- Subset models for non-interaction contrasts ---
  # Each subset model uses only the relevant samples (e.g. one treatment group)
  # with a simpler design formula (~ source or ~ treatment).
  subset_contrasts <- contrasts_info %>% dplyr::filter(model_type == "subset")
  if (nrow(subset_contrasts) > 0) {
    unique_subsets <- subset_contrasts %>%
      dplyr::distinct(model_id, subset_variable, subset_level)

    for (i in seq_len(nrow(unique_subsets))) {
      sub_var <- unique_subsets$subset_variable[i]
      sub_level <- unique_subsets$subset_level[i]
      mid <- unique_subsets$model_id[i]

      if (verbose) message("=== Creating subset DESeq2 dataset: ", sub_var, " = ", sub_level, " ===")

      col_data_sub <- col_data[col_data[[sub_var]] == sub_level, ]
      col_data_sub <- col_data_sub %>%
        dplyr::mutate(dplyr::across(dplyr::where(is.factor), droplevels))

      txi_sub <- import_salmon_data(col_data_sub, tx2gene, level = level)

      sub_design <- subset_design(design, sub_var)

      dds_list[[mid]] <- create_deseq_dataset(
        txi_sub, col_data_sub, level = level, design = sub_design,
        min_count = min_count, min_samples = min_samples, verbose = verbose
      )
    }
  }

  if (verbose) message("=== Extracting differential expression results ===")
  results_list <- extract_de_results(dds_list, contrasts_info)

  if (verbose) message("=== Converting results to tibble ===")
  master_results <- results_to_tibble(results_list, contrasts_info, tx2gene = tx2gene, level = level)

  if (verbose) message("=== Calculating summary statistics ===")
  summary_stats <- calculate_summary_stats(master_results, alpha = alpha, verbose = verbose)

  if (verbose) message("=== ", toupper(level), "-level analysis complete ===")

  list(
    dds_list = dds_list,
    master_results = master_results,
    summary_stats = summary_stats,
    col_data = col_data,
    tx2gene = tx2gene,
    txi = txi,
    factor_levels = extract_factor_levels(col_data, print_table = FALSE),
    contrasts_info = contrasts_info,
    results_list = results_list,
    level = level
  )
}



#' Complete DESeq2 pipeline from config to results
#'
#' End-to-end driver: reads the JSON config, loads the GTF, builds the
#' transcript-to-gene mapping, imports Salmon quantifications via tximport,
#' fits one full interaction model and four subset models, and returns
#' shrunk LFC results for all five contrasts plus summary statistics.
#'
#' # Assumptions
#'
#' The pipeline is restricted to a **2x2 factorial design**: exactly two
#' `treatment` levels and exactly two `source` levels (e.g. IP / Input).
#' `extract_factor_levels()` errors out otherwise.
#'
#' # Config JSON schema
#'
#' The config file must parse to a list with these keys:
#' \preformatted{
#' {
#'   "treatment": {
#'     "WT": ["sample1", "sample2", "sample3", "sample4"],
#'     "KO": ["sample5", "sample6", "sample7", "sample8"]
#'   },
#'   "source": {
#'     "Input": ["sample1", "sample2", "sample5", "sample6"],
#'     "IP":    ["sample3", "sample4", "sample7", "sample8"]
#'   },
#'   "treatment_comparison": ["KO", "WT"],
#'   "source_comparison":    ["IP", "Input"]
#' }
#' }
#' Both `*_comparison` vectors are `c(numerator, denominator)`, so the
#' example above produces LFCs of log2(KO/WT) and log2(IP/Input). Every
#' sample must appear in exactly one treatment and one source group.
#'
#' # Interaction LFC caveat
#'
#' Because the source-in-treatment and treatment-in-source contrasts come
#' from subset models while the interaction comes from the full interaction
#' model, the interaction LFC is **not** equal to the algebraic difference
#' of the two corresponding subset LFCs. See [define_contrasts()] for the
#' rationale.
#'
#' @param config_json_path Path to JSON configuration file
#' @param gtf_path Path to GTF annotation file
#' @param results_dir Directory that directly contains sample subfolders
#'   (one per `sample_name`, each holding a `quant.sf` file)
#' @param level Analysis level: "transcript", "gene", or "both" (default: "transcript")
#' @param design Design formula for the full interaction model
#'   (default: `~ source + treatment + source:treatment`)
#' @param alpha Significance threshold for the summary statistics (default: 0.05)
#' @param verbose Print progress messages (default: TRUE)
#' @return If level="both": list with "transcript" and "gene" elements, each containing analysis results.
#'         Otherwise: list containing dds_list (named list of DESeqDataSet objects:
#'         "full" for interaction model, "subset_*" for per-group models),
#'         master_results, summary_stats, and other objects
#' @export
run_deseq2_pipeline <- function(config_json_path,
                                gtf_path,
                                results_dir = "/nfcore-rnaseq-pipeline/results_grcm38/star_salmon/",
                                level = "transcript",
                                design = ~ source + treatment + source:treatment,
                                covariates = NULL,
                                min_count = 0,
                                min_samples = 0,
                                alpha = 0.05,
                                verbose = TRUE) {
  
  # Validate level parameter
  if (!level %in% c("transcript", "gene", "both")) {
    stop("level must be 'transcript', 'gene', or 'both'")
  }
  
  if (verbose) message("=== Loading configuration ===")
  dat_settings <- jsonlite::read_json(config_json_path)

  if (verbose) message("=== Creating column data ===")
  col_data <- get_coldata(dat_settings) %>%
    add_filepaths(results_dir = results_dir)

  if (!is.null(covariates)) {
    if (!"sample_name" %in% names(covariates)) {
      stop("covariates must have a 'sample_name' column")
    }
    missing <- setdiff(col_data$sample_name, covariates$sample_name)
    if (length(missing) > 0) {
      stop("covariates is missing ", length(missing), " sample(s): ",
           paste(utils::head(missing, 5), collapse = ", "))
    }
    col_data <- dplyr::left_join(col_data, covariates, by = "sample_name")
  }
  needed <- setdiff(all.vars(design), names(col_data))
  if (length(needed) > 0) {
    stop("design references column(s) not in col_data: ",
         paste(needed, collapse = ", "),
         "\n  Pass them via the `covariates` argument.")
  }

  if (verbose) message("=== Loading GTF annotation ===")
  gtf <- load_gtf(gtf_path)

  if (verbose) message("=== Creating transcript-to-gene mapping ===")
  tx2gene <- create_tx2gene(gtf, verbose = verbose)

  if (verbose) message("=== Extracting factor levels ===")
  factor_levels <- extract_factor_levels(col_data, print_table = verbose)

  if (verbose) message("=== Defining contrasts ===")
  contrasts_info <- define_contrasts(factor_levels, verbose = verbose)

  if (level == "both") {
    if (verbose) message("\n=== TRANSCRIPT-LEVEL ANALYSIS ===")
    results_tx <- .run_single_level_analysis(
      col_data, tx2gene, "transcript",
      contrasts_info, design, alpha, verbose, min_count, min_samples
    )

    if (verbose) message("\n=== GENE-LEVEL ANALYSIS ===")
    results_gene <- .run_single_level_analysis(
      col_data, tx2gene, "gene",
      contrasts_info, design, alpha, verbose, min_count, min_samples
    )

    list(
      transcript = results_tx,
      gene = results_gene
    )
  } else {
    .run_single_level_analysis(
      col_data, tx2gene, level,
      contrasts_info, design, alpha, verbose, min_count, min_samples
    )
  }
}


#' Plot dispersion estimates for each treatment group
#' 
#' Generates dispersion plots for each treatment group individually to check
#' for quality issues (e.g. issues with gene dispersions) between Input and IP 
#' samples within specific treatments.
#' 
#' @param dds A DESeqDataSet object
#' @param output_dir Directory to save plots. If NULL (default), plots are displayed.
#' @return Invisible list of subsetted DESeqDataSet objects
#' @export
plot_treatment_dispersions <- function(dds, output_dir = NULL) {
  
  # Check for treatment column
  if (is.null(dds$treatment)) {
    stop("dds object must contain a 'treatment' column in colData")
  }
  
  treatments <- levels(dds$treatment)
  results_list <- list()
  
  if (!is.null(output_dir)) {
    if (!dir.exists(output_dir)) {
      dir.create(output_dir, recursive = TRUE)
    }
  }
  
  for (trt in treatments) {
    # Subset dds for the current treatment
    dds_sub <- dds[, dds$treatment == trt]
    dds_sub$treatment <- droplevels(dds_sub$treatment)
    dds_sub$source <- droplevels(dds_sub$source)

    # The ~source design needs at least two source levels and at least two
    # samples per level to estimate dispersion. Skip otherwise with a warning.
    src_counts <- table(dds_sub$source)
    if (length(src_counts) < 2 || any(src_counts < 2)) {
      warning("Skipping treatment '", trt,
              "': needs >=2 source levels with >=2 samples each (got ",
              paste0(names(src_counts), "=", src_counts, collapse = ", "), ")")
      next
    }

    message("Generating dispersion plot for treatment: ", trt)

    DESeq2::design(dds_sub) <- ~ source
    
    # Re-run DESeq2 (estimates dispersions for the subset)
    # Using quiet=TRUE to suppress standard DESeq2 progress
    dds_sub <- DESeq2::DESeq(dds_sub, quiet = TRUE)
    
    # Store result
    results_list[[trt]] <- dds_sub
    
    # Plot
    if (!is.null(output_dir)) {
      pdf_file <- file.path(output_dir, paste0("dispersion_plot_", trt, ".pdf"))
      grDevices::pdf(pdf_file)
      DESeq2::plotDispEsts(dds_sub, main = paste("Dispersion Estimates:", trt))
      grDevices::dev.off()
    } else {
      DESeq2::plotDispEsts(dds_sub, main = paste("Dispersion Estimates:", trt))
    }
  }
  
  invisible(results_list)
}


#' DESeq2 RNA-seq Analysis Pipeline
#' 
#' A complete pipeline for differential expression analysis using DESeq2
#' with tximport for Salmon quantification data. Supports both transcript-level
#' and gene-level analysis.
#' 
#' @author Bioinformatics Pipeline
#' @date 2026-01-27


# Required packages ------------------------------------------------------------
# jsonlite, dplyr, purrr, magrittr, forcats, tibble, rtracklayer, 
# tximport, DESeq2, knitr



#' Get col_data for DESeq2 from JSON config
#'
#' @param settings Configurations object from JSON
#' @return col_data tibble with treatment and source columns
#' @export
#' @import purrr magrittr jsonlite dplyr forcats tibble
get_coldata <- function(settings) {
  # Iterate over group names and create treatment column
  col_data_treatment <- dplyr::bind_rows(
    settings$treatment %>%
      # Iterate over all treatment names
      names() %>%
      purrr::map(function(group_name) {
        # Create tibble for treatment and its samples
        tibble::tibble(
          treatment = group_name,
          sample_name = settings$treatment[group_name] %>% unlist()
        ) %>% return()
      })
  )
  
  # Attach IP/Input to col_data
  col_data_source <- dplyr::bind_rows(
    settings$source %>%
      names() %>%
      purrr::map(function(group_name) {
        tibble::tibble(
          source = group_name,
          sample_name = settings$source[group_name] %>% unlist()
        ) %>% return()
      })
  )
  
  # Join col_data source and treatment
  col_data <- dplyr::left_join(
    x = col_data_source,
    y = col_data_treatment,
    by = "sample_name"
  )
  
  # Get the proper factor levels of treatment and source in col_data
  # The first element in the comparison list is the "treatment" (numerator),
  # the second is the "reference" (denominator).
  # DESeq2 uses the first factor level as reference.
  # So we reverse the input list to set the second element as the first level.
  treatment_fct_levels <- rev(unlist(settings$treatment_comparison))
  source_fct_levels <- rev(unlist(settings$source_comparison))
  
  col_data <-
    col_data %>%
    dplyr::mutate(
      treatment = forcats::as_factor(treatment) %>%
        forcats::fct_relevel(treatment_fct_levels),
      source = forcats::as_factor(source) %>%
        forcats::fct_relevel(source_fct_levels)
    )
  
  return(col_data)
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
    tx = gtf$transcript_id[!is.na(gtf$transcript_id)],
    gene = gtf$gene_id[!is.na(gtf$transcript_id)]
  )
  
  if (verbose) {
    cat("tx2gene mapping:", nrow(tx2gene), "transcript-gene pairs\n")
  }
  
  return(tx2gene)
}



#' Import Salmon quantification data using tximport
#'
#' @param col_data Column data with filepath and sample_name columns
#' @param tx2gene Transcript-to-gene mapping data.frame
#' @param level Quantification level: "transcript" or "gene" (default: "transcript")
#' @param counts_from_abundance Abundance scaling method (default: "scaledTPM")
#' @return tximport object
#' @export
import_salmon_data <- function(col_data, tx2gene, 
                               level = "transcript",
                               counts_from_abundance = "scaledTPM") {
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
  
  tx_out <- ifelse(level == "transcript", TRUE, FALSE)
  
  txi <- tximport::tximport(
    files_vector,
    type = "salmon",
    tx2gene = tx2gene,
    txOut = tx_out,
    countsFromAbundance = counts_from_abundance
  )
  
  return(txi)
}



#' Extract factor levels from col_data
#'
#' @param col_data Column data tibble with treatment and source factors
#' @param print_table Print summary table (default: TRUE)
#' @return Named list with treatment_a, treatment_b, source_a, source_b
#' @export
extract_factor_levels <- function(col_data, print_table = TRUE) {
  treatment_levels <- levels(col_data$treatment)
  source_levels <- levels(col_data$source)
  
  if (length(treatment_levels) < 2 || length(source_levels) < 2) {
    stop("Need at least 2 levels for both treatment and source")
  }
  
  factor_levels <- list(
    treatment_a = treatment_levels[1],
    treatment_b = treatment_levels[2],
    source_a = source_levels[1],
    source_b = source_levels[2]
  )
  
  if (print_table) {
    summary_table <- tibble::tibble(
      condition = c("a", "b"),
      treatment = c(factor_levels$treatment_a, factor_levels$treatment_b),
      source = c(factor_levels$source_a, factor_levels$source_b)
    )
    print(knitr::kable(summary_table))
  }
  
  return(factor_levels)
}



#' Create DESeq2 dataset from tximport object
#'
#' @param txi tximport object
#' @param col_data Column data for samples
#' @param level Quantification level used ("transcript" or "gene")
#' @param design Design formula (default: ~ source + treatment + source:treatment)
#' @param verbose Print dimensions and size factors (default: TRUE)
#' @return DESeqDataSet object after running DESeq()
#' @export
create_deseq_dataset <- function(txi, col_data, 
                                 level = "transcript",
                                 design = ~ source + treatment + source:treatment,
                                 verbose = TRUE) {
  dds <- DESeq2::DESeqDataSetFromTximport(
    txi,
    colData = col_data,
    design = design
  )
  
  if (verbose) {
    feature_type <- ifelse(level == "transcript", "transcripts", "genes")
    cat("DESeq2 object dimensions:", dim(dds)[1], feature_type, "x", dim(dds)[2], "samples\n")
  }
  
  # Run DESeq2 analysis
  dds <- DESeq2::DESeq(dds)
  
  if (verbose) {
    # Check size factors
    size_factors <- DESeq2::sizeFactors(dds)
    cat("Size factors (should be ~1):", round(size_factors, 3), "\n")
  }
  
  return(dds)
}



#' Define contrasts for differential expression analysis
#'
#' @param factor_levels List with treatment_a, treatment_b, source_a, source_b
#' @param verbose Print contrast info (default: TRUE)
#' @return Tibble with contrast definitions
#' @export
define_contrasts <- function(factor_levels, verbose = TRUE) {
  treatment_a <- factor_levels$treatment_a
  treatment_b <- factor_levels$treatment_b
  source_a <- factor_levels$source_a
  source_b <- factor_levels$source_b
  
  # Construct dynamic coefficient names
  # DESeq2 creates coefficients as "level_b_vs_a" where a is reference (first level)
  # So the actual comparison in the results is b vs a (numerator vs denominator)
  # We use make.names() because DESeq2/model.matrix replaces spaces/chars with dots
  source_coef_name <- paste0("source_", make.names(source_b), "_vs_", make.names(source_a))
  treatment_coef_name <- paste0("treatment_", make.names(treatment_b), "_vs_", make.names(treatment_a))
  interaction_coef_name <- paste0("source", make.names(source_b), ".treatment", make.names(treatment_b))
  
  # Construct comparison names - Match actual DESeq2 output direction
  # The results are b vs a (log2(b/a)), so we should label them b_vs_a
  source_coef_print <- paste0(source_b, "_vs_", source_a)
  treatment_coef_print <- paste0(treatment_b, "_vs_", treatment_a)
  interaction_coef_print <- paste0(source_a, "_", treatment_a, "_vs_", source_b, "_", treatment_b)
  
  # Define contrasts with appropriate shrinkage method
  contrasts_info <- tibble::tribble(
    ~effect_name, ~description, ~coefs_to_use, ~shrinkage_method,
    paste0(source_coef_print, "_in_", treatment_a), paste("Source effect in", treatment_a), source_coef_name, "apeglm",
    paste0(source_coef_print, "_in_", treatment_b), paste("Source effect in", treatment_b), c(source_coef_name, interaction_coef_name), "ashr",
    paste0(treatment_coef_print, "_in_", source_a), paste("Treatment effect in", source_a), treatment_coef_name, "apeglm",
    paste0(treatment_coef_print, "_in_", source_b), paste("Treatment effect in", source_b), c(treatment_coef_name, interaction_coef_name), "ashr",
    "interaction_effect", "Interaction term", interaction_coef_name, "apeglm"
  )
  
  if (verbose) {
    print(contrasts_info %>% dplyr::select(effect_name, description, shrinkage_method))
  }
  
  return(contrasts_info)
}




#' Extract differential expression results with shrinkage
#'
#' @param dds DESeqDataSet object
#' @param contrasts_info Contrasts definition tibble from define_contrasts()
#' @param tx2gene Transcript-to-gene mapping (NULL for gene-level)
#' @param level Quantification level ("transcript" or "gene")
#' @return List of DESeq2 results objects
#' @export
extract_de_results <- function(dds, contrasts_info, tx2gene = NULL, level = "transcript") {
  results_list <- purrr::map2(
    contrasts_info$coefs_to_use,
    contrasts_info$shrinkage_method,
    function(coefs_vector, method) {
      if (length(coefs_vector) == 1) {
        # Simple coefficient - use apeglm (more accurate)
        DESeq2::lfcShrink(dds, coef = coefs_vector[1], type = method)
      } else {
        # Combined coefficients - must use ashr
        # ashr supports contrast list syntax
        DESeq2::lfcShrink(dds, contrast = list(coefs_vector, character(0)), type = method)
      }
    }
  )
  
  names(results_list) <- contrasts_info$effect_name
  
  return(results_list)
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
  
  master_results <- purrr::imap(results_list, function(res, effect_name) {
    result_tbl <- res %>%
      tibble::as_tibble(rownames = "feature_id")
    
    # For transcript level, add gene info; for gene level, rename column
    if (level == "transcript") {
      result_tbl <- result_tbl %>%
        dplyr::left_join(tx2gene, by = c("feature_id" = "tx"))
    } else {
      result_tbl <- result_tbl %>%
        dplyr::rename(gene = feature_id)
    }
    
    result_tbl %>%
      dplyr::mutate(
        effect = effect_name,
        shrinkage_method = contrasts_info$shrinkage_method[contrasts_info$effect_name == effect_name],
        level = level
      ) %>%
      dplyr::arrange(padj)
  }) %>%
    dplyr::bind_rows()
  
  return(master_results)
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
  
  return(summary_stats)
}



#' Helper function to run single-level analysis
#'
#' @keywords internal
.run_single_level_analysis <- function(col_data, tx2gene, level, 
                                       contrasts_info, design, alpha, verbose) {
  if (verbose) cat("=== Importing Salmon data (", level, "-level) ===\n", sep = "")
  txi <- import_salmon_data(col_data, tx2gene, level = level)
  
  if (verbose) cat("=== Creating DESeq2 dataset ===\n")
  dds <- create_deseq_dataset(txi, col_data, level = level, design = design, verbose = verbose)
  
  if (verbose) cat("=== Extracting differential expression results ===\n")
  results_list <- extract_de_results(dds, contrasts_info, tx2gene = tx2gene, level = level)
  
  if (verbose) cat("=== Converting results to tibble ===\n")
  master_results <- results_to_tibble(results_list, contrasts_info, tx2gene = tx2gene, level = level)
  
  if (verbose) cat("=== Calculating summary statistics ===\n")
  summary_stats <- calculate_summary_stats(master_results, alpha = alpha, verbose = verbose)
  
  if (verbose) cat("=== ", toupper(level), "-level analysis complete ===\n", sep = "")
  
  return(list(
    dds = dds,
    master_results = master_results,
    summary_stats = summary_stats,
    col_data = col_data,
    tx2gene = tx2gene,
    txi = txi,
    factor_levels = extract_factor_levels(col_data, print_table = FALSE),
    contrasts_info = contrasts_info,
    results_list = results_list,
    level = level
  ))
}



#' Complete DESeq2 pipeline from config to results
#'
#' @param config_json_path Path to JSON configuration file
#' @param gtf_path Path to GTF annotation file
#' @param results_dir Directory that directly contains sample subfolders
#' @param level Analysis level: "transcript", "gene", or "both" (default: "transcript")
#' @param design Design formula for DESeq2
#' @param alpha Significance threshold
#' @param verbose Print progress messages
#' @return If level="both": list with "transcript" and "gene" elements, each containing analysis results
#'         Otherwise: list containing dds, master_results, summary_stats, and other objects
#' @export
run_deseq2_pipeline <- function(config_json_path,
                                gtf_path,
                                results_dir = "/nfcore-rnaseq-pipeline/results_grcm38/star_salmon/",
                                level = "transcript",
                                design = ~ source + treatment + source:treatment,
                                alpha = 0.05,
                                verbose = TRUE) {
  
  # Validate level parameter
  if (!level %in% c("transcript", "gene", "both")) {
    stop("level must be 'transcript', 'gene', or 'both'")
  }
  
  if (verbose) cat("=== Loading configuration ===\n")
  dat_settings <- jsonlite::read_json(config_json_path)
  
  if (verbose) cat("=== Creating column data ===\n")
  col_data <- get_coldata(dat_settings) %>%
    add_filepaths(results_dir = results_dir)
  
  if (verbose) cat("=== Loading GTF annotation ===\n")
  gtf <- load_gtf(gtf_path)
  
  if (verbose) cat("=== Creating transcript-to-gene mapping ===\n")
  tx2gene <- create_tx2gene(gtf, verbose = verbose)
  
  if (verbose) cat("=== Extracting factor levels ===\n")
  factor_levels <- extract_factor_levels(col_data, print_table = verbose)
  
  if (verbose) cat("=== Defining contrasts ===\n")
  contrasts_info <- define_contrasts(factor_levels, verbose = verbose)
  
  # Run analysis for specified level(s)
  if (level == "both") {
    # Run both transcript and gene level analyses
    if (verbose) cat("\n=== TRANSCRIPT-LEVEL ANALYSIS ===\n")
    results_tx <- .run_single_level_analysis(
      col_data, tx2gene, "transcript", 
      contrasts_info, design, alpha, verbose
    )
    
    if (verbose) cat("\n=== GENE-LEVEL ANALYSIS ===\n")
    results_gene <- .run_single_level_analysis(
      col_data, tx2gene, "gene", 
      contrasts_info, design, alpha, verbose
    )
    
    return(list(
      transcript = results_tx,
      gene = results_gene
    ))
  } else {
    # Run single level analysis
    return(.run_single_level_analysis(
      col_data, tx2gene, level,
      contrasts_info, design, alpha, verbose
    ))
  }
}

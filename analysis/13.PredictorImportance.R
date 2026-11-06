# ---
# title: National Models 5.0 - predictor importance
# author: Elly Knight
# created: November 5, 2026
# ---

#NOTES################################

# This script extracts bootstrap-level predictor importance directly from the
# fitted landbird models. It preserves the native model-region identifiers and
# does not apply release decisions or distribution-stage region renaming.
#
# The completed bootstrap-level object is the source of truth for both
# 14.Summarize.R and downstream distribution products.

#PREAMBLE############################

#1. Load packages----
library(tidyverse)
library(gbm)
library(parallel)

#2. Set paths----
root <- "G:/Shared drives/BAM_NationalModels5"
input <- file.path(root, "output", "06_bootstraps")
output <- file.path(root, "output", "13_predictor_importance")
output_file <- file.path(
  output,
  "BAMV5-predictor-importance-bootstrap.rds"
)

if (!dir.exists(input)) {
  stop("Bootstrap-model directory not found: ", input)
}

dir.create(output, recursive = TRUE, showWarnings = FALSE)

#EXTRACT IMPORTANCE##################

#1. Index model files----
model_files <- list.files(
  input,
  pattern = "\\.RData$",
  recursive = TRUE,
  full.names = TRUE,
  ignore.case = TRUE
)

if (length(model_files) == 0) {
  stop("No bootstrap-model files were found in: ", input)
}

#2. Extract and complete one species-region model----
extract_model_importance <- function(model_file) {
  model_env <- new.env(parent = emptyenv())
  loaded_objects <- load(model_file, envir = model_env)

  if (!"b.list" %in% loaded_objects) {
    stop("Expected object `b.list` was not found in: ", model_file)
  }

  b.list <- model_env[["b.list"]]

  if (!is.list(b.list) || length(b.list) != 32) {
    stop("Expected 32 bootstrap models in: ", model_file)
  }

  file_stem <- tools::file_path_sans_ext(basename(model_file))
  file_parts <- strsplit(file_stem, "_", fixed = TRUE)[[1]]

  if (length(file_parts) < 2) {
    stop("Could not parse species and region from: ", model_file)
  }

  spp <- file_parts[length(file_parts) - 1]
  region <- file_parts[length(file_parts)]

  importance <- lapply(seq_along(b.list), function(boot) {
    fit <- b.list[[boot]]

    if (is.null(fit)) {
      stop("Bootstrap ", boot, " is NULL in: ", model_file)
    }

    fit_importance <- summary(fit, plotit = FALSE)

    if (!all(c("var", "rel.inf") %in% names(fit_importance))) {
      stop("Unexpected GBM importance columns in: ", model_file)
    }

    total_importance <- sum(fit_importance$rel.inf, na.rm = TRUE)

    if (!is.finite(total_importance) || total_importance <= 0) {
      stop(
        "Non-positive total importance for bootstrap ",
        boot,
        " in: ",
        model_file
      )
    }

    data.frame(
      boot = as.integer(boot),
      predictor = as.character(fit_importance$var),
      rel_inf = fit_importance$rel.inf / total_importance * 100,
      stringsAsFactors = FALSE
    )
  }) |>
    bind_rows()

  # A predictor absent from a bootstrap has zero relative influence. Complete
  # every predictor observed in any bootstrap for this species-region model.
  importance <- expand_grid(
    boot = seq_along(b.list),
    predictor = unique(importance$predictor)
  ) |>
    left_join(importance, by = c("boot", "predictor")) |>
    mutate(
      spp = spp,
      region = region,
      boot = as.integer(boot),
      rel_inf = replace_na(rel_inf, 0)
    ) |>
    dplyr::select(spp, region, boot, predictor, rel_inf)

  importance
}

#3. Extract files in parallel----
n_workers <- suppressWarnings(
  as.integer(Sys.getenv("BAM_N_WORKERS", unset = "5"))
)

if (is.na(n_workers) || n_workers < 1) {
  stop("BAM_N_WORKERS must be a positive integer.")
}

n_workers <- min(n_workers, length(model_files))

if (n_workers > 1) {
  cl <- makePSOCKcluster(n_workers, type = "PSOCK")
  importance_list <- tryCatch(
    {
      invisible(clusterEvalQ(cl, {
        library(tidyverse)
        library(gbm)
      }))
      clusterExport(cl, "extract_model_importance")
      parLapply(cl, model_files, extract_model_importance)
    },
    finally = stopCluster(cl)
  )
} else {
  importance_list <- lapply(model_files, extract_model_importance)
}

bam_predictor_importance_bootstrap_v5 <- importance_list |>
  bind_rows() |>
  arrange(spp, region, boot, desc(rel_inf)) |>
  as.data.frame()

rm(importance_list)
gc()

#VALIDATE############################

#1. Confirm each model-bootstrap-predictor row is unique----
duplicate_rows <- bam_predictor_importance_bootstrap_v5 |>
  count(spp, region, boot, predictor) |>
  dplyr::filter(n > 1)

if (nrow(duplicate_rows) > 0) {
  stop("Duplicate predictor-importance rows were produced.")
}

#2. Confirm relative influence sums to 100 within every bootstrap----
invalid_totals <- bam_predictor_importance_bootstrap_v5 |>
  group_by(spp, region, boot) |>
  summarize(total = sum(rel_inf), .groups = "drop") |>
  dplyr::filter(abs(total - 100) > 1e-8)

if (nrow(invalid_totals) > 0) {
  stop("Relative influence does not sum to 100 within every model bootstrap.")
}

#SAVE################################

#1. Save the native bootstrap-level source object----
saveRDS(
  bam_predictor_importance_bootstrap_v5,
  file = output_file,
  compress = "xz"
)

print(paste("Predictor importance saved to", output_file))

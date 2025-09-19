
# R/load_sets.R
# Read manifest + datasets back. No L/tau reattachment.

#' Load datasets from a recipe-sets manifest
#'
#' Reads a \code{manifest.rds} created by \code{generate_recipe_sets()}, loads one
#' dataset per row (preferring \code{rds} → \code{rdata} → \code{csv} → \code{txt}),
#' restores attribute \code{"achieved_censoring"}, and returns a named list
#' of \code{list(data = <data.frame>, meta = <list>)}.
#' @param manifest_path Path to \code{manifest.rds}.
#' @return A named list where each element is \code{list(data=..., meta=...)}.
#' @examples
#' \dontrun{
#' sets <- load_recipe_sets("checks/manifest.rds")
#' names(sets)
#' str(sets[[1]]$meta)
#' }
#' @export
load_recipe_sets <- function(manifest_path) {
  man <- readRDS(manifest_path)
  if (nrow(man) == 0) stop("manifest has 0 rows; rebuild it with rebuild_manifest()")
  out <- vector("list", nrow(man))

  read_one <- function(row) {
    if (!is.na(row$file_rds))       df <- readRDS(row$file_rds)
    else if (!is.na(row$file_rdata)) { e <- new.env(parent = emptyenv()); load(row$file_rdata, envir = e); df <- get("dat_obj", envir = e) }
    else if (!is.na(row$file_csv))  df <- utils::read.csv(row$file_csv, check.names = FALSE)
    else if (!is.na(row$file_txt))  df <- utils::read.table(row$file_txt, header = TRUE, sep = "\t", check.names = FALSE)
    else stop("No data file path in manifest row.")

    ac <- if ("achieved_censoring" %in% names(row)) row$achieved_censoring else NA_real_
    if (is.na(ac)) ac <- mean(df$status == 0, na.rm = TRUE)
    attr(df, "achieved_censoring") <- ac

    meta <- if ("meta" %in% names(row)) row$meta[[1]] else as.list(row)
    if (is.null(meta) || length(meta) == 0L || is.null(meta$model)) {
      meta <- list(
        dataset_id = sprintf("sc%03d_r%02d", row$scenario_id, row$rep),
        scenario_id = row$scenario_id,
        rep = row$rep,
        seed_used = row$seed,
        n = nrow(df),
        event_rate = mean(df$status == 1L, na.rm = TRUE),
        achieved_censoring = ac
      )
    }
    list(data = df, meta = meta)
  }

  for (i in seq_len(nrow(man))) {
    out[[i]] <- read_one(man[i, , drop = FALSE])
  }
  names(out) <- sprintf("sc%03d_r%02d", man$scenario_id, man$rep)
  out
}

#' Rebuild manifest for an existing output directory (no re-simulation)
#'
#' Reads datasets already written in \code{out_dir} and reconstructs a
#' \code{manifest.rds} with rich metadata (model, baseline, effects, etc.).
#'
#' @param base_recipe A validated recipe list (use \code{validate_recipe()} if needed).
#' @param vary Named list used originally in \code{generate_recipe_sets()} for the grid.
#' @param out_dir Directory that already contains the datasets.
#' @param filename_template The same template you used when writing files
#'   (default \code{"sc{scenario_id}_r{rep}"}). It may also include tokens for dotted paths from \code{vary}.
#' @return The rebuilt manifest (also writes \code{manifest.rds} in \code{out_dir}).
#' @export
rebuild_manifest <- function(base_recipe, vary, out_dir,
                             filename_template = "sc{scenario_id}_r{rep}") {
  `%||%` <- function(x, y) if (is.null(x)) y else x
  sanitize <- function(x) {
    x <- gsub("[^A-Za-z0-9._-]+", "_", x)
    x <- gsub("_+", "_", x)
    x <- sub("^_", "", x)
    x <- sub("_$", "", x)
    x
  }
  build_name <- function(template, tokens) {
    out <- template
    for (nm in names(tokens)) {
      ph  <- paste0("{", nm, "}")
      val <- tokens[[nm]]
      if (is.numeric(val) && length(val) == 1)
        val <- format(val, digits = 6, trim = TRUE, scientific = FALSE)
      if (length(val) > 1) val <- paste(val, collapse = "-")
      out <- gsub(ph, as.character(val), out, fixed = TRUE)
    }
    sanitize(out)
  }
  get_by_path <- function(obj, path) {
    parts <- strsplit(path, ".", fixed = TRUE)[[1]]
    cur <- obj
    for (p in parts) cur <- cur[[p]]
    cur
  }
  collect_meta <- function(dat, recipe, files, seed, scenario_id, rep, sc_params) {
    has_arm <- "arm" %in% names(dat)
    list(
      dataset_id = sprintf("sc%03d_r%02d", scenario_id, rep),
      scenario_id = scenario_id,
      rep = rep,
      seed_used = seed %||% NA_integer_,
      n = nrow(dat),
      n_treat = if (has_arm) sum(dat$arm == 1L) else NA_integer_,
      n_control = if (has_arm) sum(dat$arm == 0L) else NA_integer_,
      event_rate = mean(dat$status == 1L),
      achieved_censoring = attr(dat, "achieved_censoring", exact = TRUE) %||% mean(dat$status == 0L),
      model = recipe$event_time$model,
      baseline = recipe$event_time$baseline,
      effects  = recipe$event_time$effects,
      treatment = recipe$treatment %||% NULL,
      censoring = recipe$censoring %||% NULL,
      covariates = if (!is.null(recipe$covariates$defs))
        lapply(recipe$covariates$defs, function(d) d[c("name","type","dist","params")]) else NULL,
      allocation = recipe$treatment$allocation %||% NA_character_,
      params = sc_params,
      files = files,
      created_at = as.character(Sys.time())
    )
  }

  if (!dir.exists(out_dir))
    stop("out_dir does not exist: ", out_dir)

  base_recipe <- validate_recipe(base_recipe)
  grids <- recipe_grid(base_recipe, vary)
  if (!length(grids)) grids <- list(base_recipe)

  rows <- list()
  sc_id <- 0L
  matched_any <- FALSE

  for (r in grids) {
    sc_id <- sc_id + 1L
    sc_params <- list()
    if (length(vary)) for (k in names(vary)) sc_params[[k]] <- get_by_path(r, k)

    consecutive_gaps <- 0L
    for (rep in 1:999) {
      base_name <- build_name(filename_template, c(list(scenario_id = sc_id, rep = rep), sc_params))
      fp <- list(
        rds   = file.path(out_dir, paste0(base_name, ".rds")),
        csv   = file.path(out_dir, paste0(base_name, ".csv")),
        txt   = file.path(out_dir, paste0(base_name, ".txt")),
        rdata = file.path(out_dir, paste0(base_name, ".RData"))
      )
      exists_any <- file.exists(fp$rds) | file.exists(fp$csv) | file.exists(fp$txt) | file.exists(fp$rdata)
      if (!exists_any) {
        consecutive_gaps <- consecutive_gaps + 1L
        if (rep == 1) break
        if (consecutive_gaps >= 3L) break
        next
      }
      matched_any <- TRUE
      consecutive_gaps <- 0L

      if (file.exists(fp$rds)) {
        dat <- readRDS(fp$rds)
      } else if (file.exists(fp$rdata)) {
        e <- new.env(parent = emptyenv()); load(fp$rdata, envir = e); dat <- get("dat_obj", envir = e)
      } else if (file.exists(fp$csv)) {
        dat <- utils::read.csv(fp$csv, check.names = FALSE)
      } else {
        dat <- utils::read.table(fp$txt, header = TRUE, sep = "\t", check.names = FALSE)
      }

      meta_full <- collect_meta(
        dat = dat, recipe = r, files = fp, seed = NA_integer_,
        scenario_id = sc_id, rep = rep, sc_params = sc_params
      )

      row <- data.frame(
        scenario_id = sc_id, rep = rep, seed = NA_integer_,
        achieved_censoring = meta_full$achieved_censoring, n = nrow(dat),
        file_txt = if (file.exists(fp$txt)) fp$txt else NA_character_,
        file_csv = if (file.exists(fp$csv)) fp$csv else NA_character_,
        file_rds = if (file.exists(fp$rds)) fp$rds else NA_character_,
        file_rdata = if (file.exists(fp$rdata)) fp$rdata else NA_character_,
        stringsAsFactors = FALSE
      )
      if (length(sc_params)) {
        add <- as.list(sc_params); names(add) <- paste0("p__", sanitize(names(add)))
        for (nm in names(add)) row[[nm]] <- as.character(add[[nm]])
      }
      row$meta <- I(list(meta_full))
      rows[[length(rows) + 1L]] <- row
    }
  }

  if (!matched_any) {
    stop(
      "rebuild_manifest(): no files matched the expected names under out_dir = '", out_dir, "'.\n",
      "Check filename_template and vary, and ensure files exist."
    )
  }

  man <- do.call(rbind, rows)
  saveRDS(man, file.path(out_dir, "manifest.rds"))
  man
}
